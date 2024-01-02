setwd('~/Github/Racz2024b')

library(tidyverse)

source('code/helper.R')

# toy feature matrix

fm = tibble(
  segment = c('a','b','c'),
  feature1 = c(1,0,1),
  feature2 = c(0,1,1),
  feature3 = c(1,1,0)
)

# toy test set

test = tibble(
  string = c('abc','abb','aac','bbc','bba','bca','b')
)

# toy training set

training = tibble(
  string = c('a','b','c','abc'),
  category = c('high','low','low','high')
)

# a. generate natural classes

nc = fm |>
  generateNaturalClasses()

flag1 = all(
any(str_detect(nc$feature_bundle, 'feature1')),
any(str_detect(nc$feature_bundle, 'feature2')),
any(str_detect(nc$feature_bundle, 'feature3')),
any(str_detect(nc$segments, 'a')),
any(str_detect(nc$segments, 'b')),
any(str_detect(nc$segments, 'c'))
)

# b. build distance table for pairwise segment comparisons

lookup = buildDistTable(fm, nc) %>% 
  addLevenshtein()

n_a_and_b = nc %>% 
  filter(str_detect(segments,'a'),str_detect(segments,'b')) %>% 
  nrow()
n_a_or_b = nc %>% 
  filter(str_detect(segments,'a|b')) %>% 
  nrow()

a_b_dist = 1 - n_a_and_b / (n_a_or_b)

a_b_dist_l = lookup %>% 
  filter(segment1 == 'a', segment2 == 'b') %>% 
  pull(dist)

flag2 = all(
  a_b_dist == a_b_dist_l,
  nrow(fm) * nrow(fm) + nrow(fm) * 2 == nrow(lookup)  
)

# c. align test and target words to find best phon-based alignment, here, for the 'lakok/lakom' variation

alignments = runLookup(test,training,lookup)

alignments

# d. get distance based on best alignment, here, for lakok. this takes ages.

word_distance = alignments |>
  dplyr::distinct(test,training,phon_dist) |>
  dplyr::left_join(training, join_by('training' == 'string'))

bad_lines = word_distance %>% 
  mutate(lv = stringdist::stringdist(test,training,method = 'lv')) %>% 
  filter(phon_dist > lv) %>% 
  nrow()

flag3 = bad_lines == 0

# e. use the paired data to fit some sort of a phon distance based learning model, here, a KNN

knn_out1 = KNN(dat = word_distance, distance_type = 'phon', var_s = .1, var_k = 3, var_p = 1)

# you can use a wrapper function to do a-e:
knn_out2 = KNNwrapper(test = test, training = training, feature_matrix = fm, my_distance = 'phon', my_s = .1, my_k = 3, my_p = 1)

flag4 = all(knn_out1 == knn_out2)

# you can use the wrapper function to use some other distance and skip the whole alignment bit:
knn_out3 = KNNwrapper(test = test, training = training, feature_matrix = fm, my_distance = 'jaccard', my_s = .1, my_k = 3, my_p = 1)

flag5 = any(knn_out2 != knn_out3)

# make sure the same things work with the gcm
# GCM fun with phon distance
gcm_out1 = GCM(dat = word_distance, distance_type = 'phon', var_s = .1, var_p = 1)
# GCM fun with edit distance
gcm_out2 = GCM(dat = word_distance, distance_type = 'edit', var_s = .1, var_p = 1)
# GCM wrapper
gcm_out3 = GCMwrapper(test = test, training = training, feature_matrix = fm, my_distance = 'phon', my_s = .1, my_p = 1)

flag6 = all(gcm_out1 == gcm_out3)

flags = all(flag1,flag2,flag3,flag4,flag5,flag6)

if (flags){
  'checks completed successfully'
} else {
  'ruh roh'
}
