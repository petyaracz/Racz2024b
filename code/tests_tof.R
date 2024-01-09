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
    'trumabcx',
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
    'trumabcy',
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

d

flag1 = all(d$a == 'x')
flag2 = all(d$b %in% c('y','z'))

d2 = findLargerRules(d)

d2

flag3 = all(str_detect(d2$rule, '(y|z) \\/'))
flag4 = all(str_detect(d2$rule, '^x'))

d3 = getRulesAndWords(d,d2)

d3

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

# this looks about right. the farther lower confidence is pushed out, the worse small rules are doing.

alpha11 = impugnRules(alpha1, d3)
alpha91 = impugnRules(alpha9, d3)

flag8 = all(sort(alpha1$rule)==sort(alpha11$rule))
flag9 = all(sort(alpha1$rule)==sort(alpha91$rule))

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
  select(-d4,-d5) %>% 
  unnest(c(nd4,nd5))

flag9 = all(alphas2rows$nd4 == alphas2rows$nd5)

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

# lgtm

flags = c(flag1,flag2,flag3,flag4,flag5,flag6,flag7,flag8,flag9)

if (all(flags)) {
  print('all tests passed')
} else {
  print('some tests failed')
}