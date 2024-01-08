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

getPredictions = function(test, rules){
  rules = ungroup(rules)
  matches = test %>% 
    pivot_longer(c(output1,output2), names_to = 'variant', values_to = 'output') %>% 
    mutate(
      input = glue('#{input}#'),
      output = glue('#{output}#')
    ) %>% 
    crossing(rules) %>% 
    filter(
      str_detect(input, glue('{C}{A}{D}$')),
      str_detect(output, glue('{C}{B}{D}$')),
    ) %>% 
    group_by(
      base, variant, resp1, resp2, log_odds
    ) %>% 
    arrange(
      -impugned_lower_confidence_limit
    ) %>%
    slice(1)
  return(matches)
}

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

training1 = read_tsv('dat/training_sets/training_mgl.tsv')
test1 = read_tsv('dat/training_sets/test_mgl.tsv')

# -- setup -- #
parameters = crossing(
  alpha_upper = c(.25,.5,.75,.9),
  alpha_lower = c(.25,.5,.75,.9)
)

test1 %<>%
  mutate(
    suffix = ifelse(
      variation == 'lakok/lakom', '3sg', suffix
    )
  )

trainings = training1 %>% 
  rename(orth = base) %>% 
  filter(!is.na(suffix)) %>% 
  mutate(
    input = ifelse(variation == 'hotelban/hotelben', input, str_replace(input, 'ik$', 'iK'))
  ) %>% 
  group_by(variation,suffix) %>% 
  nest() %>% 
  crossing(parameters) 

# -- run mgl to get rules -- #

if (flag == T){
  models = as.list(NULL)
  for (i in 1:nrow(trainings)) {
    tic('running current iteration: ')
    models[[i]] = mgl(training = trainings$data[[i]], alpha_lower = trainings$alpha_lower[[i]], alpha_upper = trainings$alpha_upper[[i]])
    toc()
    print(i)
  }
  save(models, file = 'dat/training_sets/trainings_mgl_cs.rda')
} else {
  # load('dat/training_sets/trainings_mgl.rda')
  # models_other = models[49:112]
  # load('dat/training_sets/trainings_mgl_cs.rda')
  # models_cs = models
  # models = c(models_cs,models_other)
  # save(models, file = 'dat/training_sets/trainings_mgl_final.rda')
  load('dat/training_sets/trainings_mgl_final.rda')
}

# -- build master file -- #

trainings %<>%
  mutate(
    mgl = map(models, ~ ungroup(.)),
    test = map2(variation, suffix, ~
                  {
                    test1 %>%
                      filter(variation == .x, suffix == .y)
                  }
                  ),
    prediction = map2(test, mgl, ~ getPredictions(.x, .y)),
    score = map(prediction, ~ getScores(.))
  )

master = trainings %>% 
  select(variation,suffix,alpha_upper,alpha_lower,score)

results = master %>% 
  unnest(score) %>% 
  group_by(variation,alpha_upper,alpha_lower) %>% 
  nest() %>%
  mutate(
    glm = map(data, ~ 
                glm(cbind(resp1,resp2) ~ mgl_score, data = ., family = binomial)
                ),
    glm_summary = map(glm, ~ broom::tidy(.))
  ) %>% 
  select(-data,-glm) %>% 
  unnest(glm_summary) %>% 
  group_by(variation) %>% 
  arrange(-statistic) %>% 
  slice(1)

results
# https://media1.tenor.com/m/pWZZ9gzx_p0AAAAC/face-melting-indiana-jones.gif

best_mgls = results %>% 
  select(variation,alpha_upper,alpha_lower) %>% 
  left_join(trainings)
