setwd('~/Github/Racz2024b')

library(tidyverse)

source('code/helper.R')

flag = F

if(flag == T){

# training and test data
test = read_tsv('dat/training_sets/test_set.tsv')
training = read_tsv('dat/training_sets/training_set.tsv')
# segment lookup
lookup_h = read_tsv('dat/segmental_distances/siptar_torkenczy_toth_racz_hungarian_dt.tsv')

# set up test and training data
testl = filter(test, variation == 'lakok/lakom')
trainingl = filter(training, variation == 'lakok/lakom')
testcs = filter(test, variation == 'cselekszenek/cselekednek')
trainingcs = filter(training, variation == 'cselekszenek/cselekednek')
testh = filter(test, variation == 'hotelban/hotelben')
trainingh = filter(training, variation == 'hotelban/hotelben')

# run the lookup for all variations
alignments_lakok = runLookup(testl,trainingl,lookup_h)
alignments_cselekszenek = runLookup(testcs,trainingcs,lookup_h)
alignments_hotelban = runLookup(testh,trainingh,lookup_h)

# save alignments
write_tsv(alignments_lakok, 'dat/alignments/alignments_lakok.tsv')
write_tsv(alignments_cselekszenek, 'dat/alignments/alignments_cselekszenek.tsv')
write_tsv(alignments_hotelban, 'dat/alignments/alignments_hotelban.tsv')

# save distances
alignments_lakok_2 = alignments_lakok %>% 
  mutate(variation = 'lakok/lakom') %>% 
  distinct(test,training,phon_dist,variation)

alignments_hotelban_2 = alignments_hotelban %>% 
  mutate(variation = 'hotelban/hotelben') %>% 
  distinct(test,training,phon_dist,variation)

alignments_cselekszenek_2 = alignments_cselekszenek %>% 
  mutate(variation = 'cselekszenek/cselekednek') %>%
  distinct(test,training,phon_dist,variation)

word_distance = bind_rows(
  alignments_lakok_2,
  alignments_hotelban_2,
  alignments_cselekszenek_2
)

write_tsv(word_distance, 'dat/alignments/word_distances.tsv')
}