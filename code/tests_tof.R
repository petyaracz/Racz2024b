setwd('~/Github/Racz2024b')

library(tidyverse)
library(viridis)

source('code/tof.R')

# let's make up some rules
# i. x -> y / c _
# ii. x -> y / bc _
# iii. x -> y / abc _
# iv. x -> z / bbc _

training = tibble(
  input = c(
    'vacx',
    'gucx',
    'pocx',
    'necx',
    'penicx',
    'glonabcx',
    'flanabcx',
    'binabcx',
    'frobbcx',
    'knebbcx',
    'kribbcx',
    'klubcx',
    'flobcx',
    'blibcx'
  ),
  output = c(
    'vacy',
    'gucy',
    'pocy',
    'necz',
    'penicz',
    'glonabcy',
    'flanabcy',
    'binabcz',
    'frobbcz',
    'knebbcz',
    'kribbcz',
    'klubcy',
    'flobcy',
    'blibcy'
  )
)

d = formatTraining(training)

flag1 = all(d$a == 'x')
flag2 = all(d$b %in% c('y','z'))

d2 = findLargerRules(d)

flag3 = all(str_detect(d2$rule, '(y|z) \\/'))
flag4 = all(str_detect(d2$rule, '^x'))

d3 = getRulesAndWords(d,d2)

flag5 = all(str_detect(d3$c, glue('{d3$C}$')))
flag6 = all(str_detect(d3$output, glue('{d3$b}#$')))

alphas = crossing(
  alpha_lower = seq(.1,1,.1),
  alpha_upper = seq(.1,1,.1)
) %>% 
  mutate(
    d4 = map2(alpha_upper, alpha_lower, ~ getRuleStats(d3, alpha_upper = .x, alpha_lower = .y))
    )

alpha1 = alphas[alphas$alpha_upper == .1,]$d4[[1]]
alpha9 = alphas[alphas$alpha_upper == .9,]$d4[[1]]

flag7 = all(alpha1$lower_confidence_limit <= alpha1$reliability)

alphas %>% 
  unnest(d4) %>% 
  mutate(
    rule_info = glue('{rule} (scope: {scope}, hits: {hits})') %>% 
      fct_reorder(-scope)
    ) %>%
  distinct(alpha_lower,rule_info,lower_confidence_limit) %>% 
  ggplot(aes(x = alpha_lower, y = lower_confidence_limit, color = rule_info)) +
  geom_line() +
  theme_bw() +
  scale_colour_viridis_d()

# this looks about right

alpha11 = impugnRules(alpha1, d3)
alpha91 = impugnRules(alpha9, d3)

flag8 = (setdiff(alpha1$rule,alpha11$rule) == 0)

alphas2 = alphas %>% 
  # filter(alpha_lower == .1) %>% 
  mutate(
    d5 = map(d4, ~ impugnRules(., d3))
  )

alphas2rows = alphas2 %>% 
  mutate(
    nd4 = map(d4, ~ nrow(.)),
    nd5 = map(d5, ~ nrow(.))
  ) %>% 
  select(-d4,-d5)

alphas2 %>% 
  unnest(d5) %>% 
  mutate(
    rule_info = glue('{rule} (rel: {round(reliability,2)})') %>% 
      fct_reorder(-reliability)
  ) %>%
  distinct(alpha_lower,alpha_upper,rule_info,impugned_lower_confidence_limit) %>% 
  ggplot(aes(x = alpha_upper, y = impugned_lower_confidence_limit, color = rule_info)) +
  geom_line() +
  theme_bw() +
  scale_colour_viridis_d() +
  facet_wrap( ~ alpha_lower)

# ho hum

alpha2 = alphas2$d5[[1]]
alpha3 = alphas2$d5[[2]]
alpha4 = alphas2$d5[[3]]

# only big rules remain. fixed
# one rule gets duplicated somehow.