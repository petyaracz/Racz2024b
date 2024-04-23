# -- head -- #

library(tidyverse)
library(magrittr)
library(glue)
library(stringi)
library(tictoc)

setwd('~/Github/Racz2024b')

source('code/mgl.R')

# -- read -- #

mglrules = read_tsv('dat/mgl/CELEXTop200_no_stress.rules')
training = read_tsv('dat/mgl/CELEXTop200_no_stress.in', skip = 1, col_names = F, n_max = 199)

# -- format -- #

# mgl setup: .75, .75, no features, doppelgaengers, impungment, phonology.

training %<>% 
  rename(input = X1, output = X2, orth = X4) %>% 
  select(orth, input, output)

mglrules = mgl(training, .75, .75)

# -- compare rules -- #

mglrules2 = mglrules %>% 
  mutate(
    A = ifelse(is.na(A), '', A),
    B = ifelse(is.na(B), '', B),
    C = ifelse(is.na(C), '', C),
    D = ifelse(is.na(D), '', D),
    rule = glue('{A} -> {B} / {C} _ {D}')
  ) %>% 
  select(A,B,C,D,rule,Scope,Hits,Reliability)

mglrules2 = mglrules %>% 
  select(A,B,C,D,rule,scope,hits,reliability)

# 1. the mgl makes "everything rules" for all alternations. the mgl doesn't, why would it?
# 2. the scopes are larger in the mgl than in the mgl
formatTraining(training) %>% 
  filter(str_detect(c, 'n$'))
# this is effectively the same issue as (1). the mgl counts scopes EVEN IF they don't belong to the given a->b alternation. I think that's the main thing?
# 3. the hits look about right
# 4. the mgl has more rules, even though it hasn't got the everything rules
mglrules3 = mglrules2 %>% 
  select(A,B,C,D) %>% 
  # replace # with '' in all cols
  mutate_all(~str_replace(., '#', '')) %>% 
  mutate(model = 'mgl')
mglrules3 = mglrules2 %>%
  select(A,B,C,D) %>% 
  mutate_all(~str_replace(., 'X ', '')) %>% 
  mutate_all(~str_replace(., 'X', '')) %>% 
  mutate(model2 = 'mgl')
models = full_join(mglrules3,mglrules3)
models %>% 
  filter(!is.na(model),!is.na(model2))
# 73 rules overlap!
models %>% 
  filter(is.na(model),!is.na(model2))
#  the mgl-specific rules look like the mgl also has them and this is only a formatting diff
models %>% 
  filter(!is.na(model),is.na(model2))
# the mgl-specific rules look like supersets of other rules (I to E)
training %>% 
  filter(str_detect(input, 'ɪ'), str_detect(output, 'æ'))
# also MGL seems to have "examples" and not all words listed for the rule

# taken together, one clear difference is that something is going on with the scopes. the mgl has larger scopes.
# but there are other differences that are unclear to me.
# welp