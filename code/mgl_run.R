library(tidyverse)
library(magrittr)
library(glue)
library(stringi)
library(tictoc)

setwd('~/Github/Racz2024b')

source('code/mgl.R')

training1 = read_tsv('dat/training_sets/training_mgl.tsv')
test1 = read_tsv('dat/training_sets/test_set.tsv')

parameters = crossing(
  alpha_upper = c(.25,.5,.75,.9),
  alpha_lower = c(.25,.5,.75,.9)
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

  # mutate(
  #   model =  pmap(across(c(data,alpha_lower,alpha_upper)), mgl)
  # )

tic('dplyr')
mgl(training = trainings$data[[1]], alpha_lower = trainings$alpha_lower[[1]], alpha_upper = trainings$alpha_upper[[1]])
toc()