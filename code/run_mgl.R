########################################################
# Running the Tiny Overlap Finder (mgl) on our Hungarian data
########################################################

# -- head -- #

library(tidyverse)
library(magrittr)
library(glue)
library(stringi)
library(tictoc)

setwd('~/Github/Racz2024b')

source('code/mgl.R')

flag = F

# -- fun -- #

# get predictions for various outputs across a set of items in a test df
# the minimal generalisation learner [puts money in jar] tiny overlap finder outputs rules
# but in order to find out how these rules apply to test forms, you need to cross them first
# take test set, rule set, return test words with rules that apply to them, with rule metric
#
# this fn will absolutely break if you have a test input / output pair that no rule matches. ho boy will it ever. but in this project every test pair is from a forced choice task and there are always A -> B / every context rules so this won't happen. but if you're using this for something else, be careful.
#
getPredictions = function(test, rules){
  # rules was rowwise, bad
  rules = ungroup(rules)
  # test has two outputs. we want them to have their own rows so we can join
  matches = test %>% 
    pivot_longer(c(output1,output2), names_to = 'variant', values_to = 'output') %>% 
    mutate(
      input = glue('#{input}#'),
      output = glue('#{output}#')
    ) %>% 
    # we cross with rules and then keep rows where the rule actually applies to the test input and generates the relevant output
    crossing(rules) %>% 
    filter(
      str_detect(input, glue('{C}{A}{D}$')),
      str_detect(output, glue('{C}{B}{D}$')),
    ) %>% 
    group_by(
      base, variant, resp1, resp2, log_odds
    ) %>% 
    # for each output variant, we keep the BEST RULE and take its MOST SOPHISTICATED confidence value, the impugned lower confidence.
    arrange(
      -impugned_lower_confidence_limit
    ) %>%
    slice(1)
  return(matches)
}

# okay we have the best rule for each test input and its both existing output variants.
# we want to get a mgl score that expresses the confidence of these two rules, following Rácz Beckner Hay Pierrehumbert 2020.
# in principle, if there's no rule for one of the variants then the mgl confidence in the other variant is 1.
# if there's no rule for either that's bad
getScores = . %>% 
  select(base,resp1,resp2,log_odds,variant,impugned_lower_confidence_limit) %>% 
  pivot_wider(names_from = variant, values_from = impugned_lower_confidence_limit) %>% 
  mutate(
    mgl_score = case_when(
      output1 == 0 ~ 0,
      output2 == 0 ~ 1,
      output1 != 0 & output2 != 0 ~ output1 / (output1 + output2)
    )
  ) %>% 
  select(-output1,-output2)

# -- read -- #

# read in training and test
# a lot of finnicking is involved in hammering the training and test sets into the appropriate format, see training_preprocessor
# most notably, for the two variations where vowel harmony doesn't matter, all vowel harmony is pruned from the suffixes so the mgl doesn't try to learn it
# this is fine: for cselekszenek, we want to know if there's vowel deletion or not (so doháṉz-o-tok or doháṉoz-tok) and don't care about the vowel (for now). for lakok, we want to know if the suffix is k or m (so lak-ok or lak-om) and again don't care about the vowel (for now).
# we also split training data into suffixes because the mgl shouldn't have to learn morphology proper.
training1 = read_tsv('dat/training_sets/training_mgl.tsv')
test1 = read_tsv('dat/training_sets/test_mgl.tsv')

# -- setup -- #

# we try some lower and upper parameters
# alpha_lower tweaks the LOWER confidence limit of the rule's reliability. it sets the penalty on a rule that applies to few forms but applies to all of them versus a rule that applies to heaps of forms but has more exceptions.
# if a rule applies to 3 forms and 2 of them are correct, that is inuitively worse than a rule that applies to 3000 forms and 2000 are correct. you can tune this with this one.
# alpha_upper tweaks the UPPER confidence limit of the residue rule in rule impungment.
# basically it expresses how much you are bothered by a subset rule doing how much of the work of the big rule. 
parameters = crossing(
  alpha_upper = c(.25,.5,.75,.9),
  alpha_lower = c(.25,.5,.75,.9)
)

# all variations have different suffixes, lakok is always 3sg but I have to add a col or functions break down the line.
test1 %<>%
  mutate(
    suffix = ifelse(
      variation == 'lakok/lakom', '3sg', suffix
    )
  )

# take trainings, rename base so it works with the mgl() fun, nest by var and suffix (so a separate mgl runs on each) cross with parameter settings (so we have one mgl for each var, suff, parameter setting combo)
trainings = training1 %>% 
  rename(orth = base) %>% 
  filter(!is.na(suffix)) %>% 
  mutate(
    input = ifelse(variation == 'hotelban/hotelben', input, str_replace(input, 'ik$', 'iK')) # this is legacy
  ) %>% 
  group_by(variation,suffix) %>% 
  nest() %>% 
  crossing(parameters) 

# -- run mgl I mean mgl to get rules -- #

# this takes a while so I put it in a dummy if.
if (flag == T){
  # set up empty models list for the mgls
  models = as.list(NULL)
  for (i in 1:nrow(trainings)) {
    tic('running current iteration: ')
    # build each mgl using data and parameters from trainings
    models[[i]] = mgl(training = trainings$data[[i]], alpha_lower = trainings$alpha_lower[[i]], alpha_upper = trainings$alpha_upper[[i]])
    toc()
    print(i)
  }
  # save resulting models
  save(models, file = 'dat/mgl/trainings_mgl_final.rda')
} else {
  # load models
  load('dat/mgl/trainings_mgl_final.rda')
}

# -- build master file -- #

# add the following things to trainings:
trainings %<>%
  mutate(
    # the models for each var, suff, parameter setting combo
    mgl = map(models, ~ ungroup(.)),
    # the relevant test set
    test = map2(variation, suffix, ~
                  {
                    test1 %>%
                      filter(variation == .x, suffix == .y) # see? we filter for var and suff
                  }
                  ),
    prediction = map2(test, mgl, ~ getPredictions(.x, .y)), # best rule for each output for each input in test
    score = map(prediction, ~ getScores(.)) # aggr score from best rules for input (how output1 or output2 it is)
  )

# we keep the var, suff, settings, scores
master = trainings %>%
  select(variation,suffix,alpha_upper,alpha_lower,score)

# we unnest score and then nest again because we want to get total scores for var, not var + suff. suff fit separately to put our finger on the scale for the mgl so it has a fighting chance against the much simpler knn and gcm.
results_a = master %>%
  unnest(score) %>%
  group_by(variation,alpha_upper,alpha_lower) %>%
  nest() %>%
  # how well did the mgl do? we use the same method as for the gcm and knn, fit a glm to predict response counts for var1 / var2 using mgl_score (how much does this test form like var1 according to the mgl)
  mutate(
    # fit glm on each nested table
    glm = map(data, ~
                glm(cbind(resp1,resp2) ~ mgl_score, data = ., family = binomial)
                ),
    # get glm summary
    glm_summary = map(glm, ~ broom::tidy(.))
  ) %>%
  select(-data) %>%
  # keep summary and unnest it
  unnest(glm_summary) %>%
  filter(term == 'mgl_score')

results = results_a %>%
  group_by(variation) %>%
  # get highest abs t value for each var
  arrange(-abs(statistic)) %>%
  slice(1)

results
# https://media1.tenor.com/m/pWZZ9gzx_p0AAAAC/face-melting-indiana-jones.gif

# here are the best models
best_mgls = results %>%
  select(variation,alpha_upper,alpha_lower) %>%
  left_join(trainings)

best_mgls %<>%
  select(variation,suffix,alpha_upper,alpha_lower,mgl) %>%
  unnest(mgl)

predictions = master %>% 
  inner_join(results) %>% 
  select(variation,score) %>% 
  unnest(score)

predictions %>% 
  group_by(variation) %>% 
  ggplot(aes(mgl_score, log_odds)) +
  geom_point() +
  geom_smooth(method = mgcv::gam) +
  theme_bw() +
  facet_wrap( ~ variation, scales = 'free')
  
# -- illustrations for the readme -- #

i_lakok_pred = trainings %>% 
  filter(alpha_upper == .25, alpha_lower == .25, variation == 'lakok/lakom') %>% 
  select(prediction) %>% 
  unnest(prediction) %>% 
  arrange(base)

i_lakok_score = trainings %>% 
  filter(alpha_upper == .25, alpha_lower == .25, variation == 'lakok/lakom') %>% 
  select(score) %>% 
  unnest(score) %>% 
  arrange(base)

# -- write -- #

write_tsv(results_a, 'dat/mgl/mgl_stats.tsv')
write_tsv(predictions, 'dat/mgl/best_mgl_predictions.tsv')
write_tsv(best_mgls, 'dat/mgl/best_mgls.tsv')
write_tsv(i_lakok_score, 'dat/mgl/i_lakok_score.tsv')
write_tsv(i_lakok_pred, 'dat/mgl/i_lakok_pred.tsv')
save(results, file = 'dat/mgl/mgl_results.rda')
save(trainings, file = 'dat/mgl/all_trainings.rda')
