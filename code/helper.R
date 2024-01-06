##
# build natural classes
##

# summarise w/o the groups message which would resurface, even with the global option set to F, if I use furrr.
summarise2 <- function(.data, ..., .groups = "drop") {
  dplyr::summarise(.data, ..., .groups = .groups)
}

# take feature matrix as tibble, return all possible feature bundles as list
makeCombinations = function(dat){
  # possible feature bundles
  my_features = names(dat)[names(dat) != 'segment']
  my_combinations = purrr::map(1:length(my_features), ~ combn(my_features, ., simplify = F)) |>
    unlist(recursive = F)
  return(my_combinations)  
}

# take feature matrix and feature combinations list (output of makeCombinations), return every possible feature combination as natural class in tibble
getSegments = function(dat, my_combination){
  
  # select features in combination out of wide dat 
  my_table = dat |>
    dplyr::select(segment,dplyr::all_of(my_combination))
  
  # keep only rows where there are no NA values
  my_table = my_table |>
    dplyr::filter(
      rowSums(sapply(my_table, is.na)) == 0
    )
  
  # list segments per unique feature combination (= natural class)
  my_classes = my_table |>
    # group by all cols except segment
    dplyr::group_by(across(-segment)) |>
    summarise2(
      segments = paste(segment, collapse = ', ')
    )
  
  # get names of features which is colnames except "segments"
  class_names = names(my_classes)[names(my_classes) != 'segments']
  
  # paste together colnames (feature name) and cells (feature value) and then turn it into neat feature bundle
  natural_classes = my_classes |>
    dplyr::group_by(segments) |>
    dplyr::summarise(
      # paste together all columns in a list
      feature_values = list(dplyr::across(dplyr::all_of(class_names))),
      feature_names = list(class_names),
      # paste together feature_values and feature_names
      feature = purrr::map2(feature_values, feature_names, ~ glue::glue('{.x}{.y}')),
      # collapse feature into one string
      feature_bundle = purrr::map(feature, ~ glue::glue_collapse(.x, sep = ', ')) |>
        unlist() |>
        stringr::str_replace_all(c('^' = '[', '$' = ']', '1' = '+', '0' = '-')) |>
        stringr::str_replace('\\+(?=[,\\]])', '1') |> # turn these back if they are in feature names
        stringr::str_replace('\\-(?=[,\\]])', '0')
    ) |>
    dplyr::ungroup() |>
    dplyr::select(feature_bundle,segments)
  
  return(natural_classes)
}

# helper function to map getSegments across dat and then bind things together into tibble
mapGetSegments = function(my_combinations,dat){
  my_combinations |>
    furrr::future_map(~ getSegments(dat, .)) |># this is the slow one
    dplyr::bind_rows(.id = 'class')
}

# helper function to take the natural classes and return the, uh, most natural ones
getNaturalClasses = function(dat){
  dat |>
    dplyr::mutate(
      nsegments = stringr::str_count(segments, ', ') + 1, # nice
      nfeatures = stringr::str_count(feature_bundle, ', ') + 1
    ) |>
    dplyr::group_by(segments) |>
    dplyr::arrange(+nfeatures) |>
    dplyr::slice(1) |>
    dplyr::arrange(nsegments) |>
    dplyr::ungroup() |>
    dplyr::select(feature_bundle,segments)
} 

# run nc building pipeline
generateNaturalClasses = function(dat){
  dat |>
    makeCombinations() |>
    mapGetSegments(dat) |>
    getNaturalClasses()
}

##
# calc segment distance
##

# take natural classes, return natural classes with s1 and s2, return shared and non-shared classes in a list
calcNC = function(nc, s1, s2){
  nc_shared = nc |> 
    dplyr::filter(stringr::str_detect(segments, s1) & stringr::str_detect(segments, s2))
  nc1 = nc |> 
    dplyr::filter(stringr::str_detect(segments, s1))
  nc2 = nc |>
    dplyr::filter(stringr::str_detect(segments, s2))
  nc_else = dplyr::bind_rows(nc1,nc2) |> 
    dplyr::distinct() |> 
    dplyr::anti_join(nc_shared, by = dplyr::join_by(feature_bundle, segments))
  list(nc_shared = nc_shared, nc_else = nc_else)
}

# take shared and non-shared classes, calc dist
calcDist = function(ncl){
  similarity = nrow(ncl$nc_shared) / (nrow(ncl$nc_shared) + nrow(ncl$nc_else))
  distance = 1 - similarity
  return(distance)
}

# find segment distance in lookup table
getDist = function(s1, s2, dist_lookup){
  dist_lookup |>
    dplyr::filter(segment1 == s1, segment2 == s2) |>
    dplyr::pull(dist)
}

## build long lookup table of segmental distances from feature matrix and natural class
buildDistTable = function(fm,nc){
  tidyr::crossing(
    segment1 = fm$segment,
    segment2 = fm$segment
  ) |> 
    dplyr::mutate(
      ncl = purrr::map2(segment1, segment2, ~ calcNC(nc, .x, .y)),
      dist = purrr::map(ncl, ~ calcDist(.x))
    ) |>
    dplyr::select(segment1,segment2,dist) |>
    tidyr::unnest(cols = dist)
}

# the levenshtein distance between any segment and nothing is 1. I wire this in the lookup table so you have to do this only once
addLevenshtein = function(lukap){
  segments = lukap |>
    dplyr::distinct(segment1) |>
    dplyr::pull()
  void = rep(' ', length(segments))
  dist = rep(1, length(segments))
  lev1 = tidyr::tibble(segment1 = segments, segment2 = void, dist = dist)
  lev2 = tidyr::tibble(segment1 = void, segment2 = segments, dist = dist)
  dplyr::bind_rows(lukap,lev1,lev2)
}

##
# calc word distance
##


# replace characters in pairs with their IPA equivalents
transcribeIPA = function(string, direction){
  if (direction == 'single'){
    stringr::str_replace_all(string, c(
      'ccs' = 'cscs', 'ssz' = 'szsz', 'zzs' = 'zszs', 'tty' = 'tyty', 'ggy' = 'gygy', 'nny' = 'nyny', 'lly' = 'jj', 'cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's'))
  } else if (direction == 'double'){
    stringr::str_replace_all(string, c('s' = 'ß', 'š' = 's', 'ṉ' = 'ny', 'ḏ' = 'gy', 'ṯ' = 'ty', 'ž' = 'zs', 'ß' = 'sz', 'č' = 'cs'))
  }
}

# take word broken up into letters (v) and length of window (l) and create all possible alignments, return as list
createAlignmentWindows = function(v, l) {
  combinations = combn(l, length(v))
  alignments = purrr::map(seq_len(ncol(combinations)), ~ {
    out = rep(' ', l)
    indices = combinations[, .x]
    out[indices] = v
    out
  })
  return(alignments)
}

# big aligner function, takes two strings, returns best alignment
alignWords = function(s1,s2,lookup_dist){
  
  # check if s1 and s2 are character vectors
  if (!is.character(s1) | !is.character(s2)){
    stop('s1 and s2 must be character vectors')
  }
  
  # make string into vector
  st1 = unlist(strsplit(s1, ''))
  st2 = unlist(strsplit(s2, ''))
  
  # check for lookup_dist format
  if (!all(c('segment1','segment2','dist') %in% colnames(lookup_dist))){
    stop('lookup_dist must have columns segment1, segment2, dist')
  }
  # check for lookup_dist NA
  if (any(is.na(lookup_dist))){
    glue::glue('lookup_dist cannot have NA values. did you change them to " "? Im gonna do it for you one last time.')
    lookup_dist = lookup_dist |> 
      dplyr::mutate(
        across(everything(), ~ ifelse(is.na(.x), ' ', .x))
      )
  }
  # check for segments not in lookup_dist
  if (any(!st1 %in% lookup_dist$segment1)){
    stop('s1 segment ', paste(st1[!st1 %in% lookup_dist$segment1], collapse = ', '), ' not in lookup table')
  }
  # check for segments not in lookup_dist for st2
  if (any(!st2 %in% lookup_dist$segment1)){
    stop('s2 segment ', paste(st2[!st2 %in% lookup_dist$segment1], collapse = ', '), ' not in lookup table')
  }
  
  # possible total alignment window length is l1 + l2 - 1
  window_length = length(st1) + length(st2) - 1
  # n empties: how many empty slots do we permit in the alignment? n empties is window length minus length of longer string. we add a toggle
  scale_n_empties = 0.5
  n_empties = window_length - max(length(st1), length(st2))
  window_length = window_length - round( n_empties * scale_n_empties, 0)
  
  # create all possible alignment windows for st1
  a1 = createAlignmentWindows(st1,window_length)
  # create all possible alignment windows for st2
  a2 = createAlignmentWindows(st2,window_length)
  
  # get all combinations of alignments1 and alignments2
  a_combinations = tidyr::crossing(a1, a2)
  
  # we now have all segments in s1 facing off all segments in s2 with holes in the middle as well. but it's in an ugly nested list that would take ages to churn through (I um happen to know that!)
  # we grab the possible arrangements from s1 and flatten them
  a1b = a_combinations |>
    dplyr::mutate(id = 1:n()) |>
    dplyr::select(id,a1) |>
    tidyr::unnest(a1) |>
    dplyr::mutate(pos = rep(1:window_length, nrow(a_combinations)))
  # we do the same for s2
  a2b = a_combinations |>
    dplyr::mutate(id = 1:n()) |>
    dplyr::select(id,a2) |>
    tidyr::unnest(a2) |>
    dplyr::mutate(pos = rep(1:window_length, nrow(a_combinations)))
  # now we can join the two into a big df and then join that df with our lookup table
  best_alignments = dplyr::left_join(a1b, a2b, by = c('id', 'pos')) |>
    # remove repetitions
    dplyr::filter(a1 != ' ' | a2 != ' ') |>
    dplyr::rename(segment1 = a1, segment2 = a2) |>
    dplyr::left_join(lookup_dist, by = c('segment1', 'segment2')) |>
    # we group alignments by id and sum the distances in each alignment
    dplyr::group_by(id) |>
    dplyr::mutate(
      phon_dist = sum(dist),
      length = n()
    ) |>
    dplyr::ungroup()
  best_alignment = best_alignments |># several ones are best. but these are equivalent except things are in different pos. so we take the shortest ones, pick the first one and remove pos 
    dplyr::filter(phon_dist == min(phon_dist)) |>
    dplyr::filter(length == min(length)) |>
    dplyr::filter(id == min(id)) |>
    dplyr::select(segment1, segment2, dist, phon_dist)
  
  # this is the best alignment. hooray!
  return(best_alignment)
}

# loop the aligner through a paired table of test and training forms
runLookup = function(test,training,lookup_dist){
  test_string = test |>
    dplyr::pull(string)
  training_string = training |>
    dplyr::pull(string)
  string_distance_lookup = tidyr::crossing(
    test = test_string,
    training = training_string
  )
  
  string_distance_res = string_distance_lookup |>
    dplyr::mutate(res = furrr::future_map2(test, training, ~ alignWords(.x, .y, lookup_dist = lookup_dist)))

  string_distance_res = string_distance_res |>
  dplyr::select(test,training,res) |>
  tidyr::unnest(res)
  
  return(string_distance_res)
}

##
# algo
##

# GCM (Nosofsky 1988) takes the word_distance data frame, distance type, variation type, and variation parameters and returns test words with predictions
GCM = function(dat, distance_type, var_s, var_p){

if(!distance_type %in% c('phon','jaccard','edit','hamming')){
  stop('distance_type must be one of phon, jaccard, edit, hamming')
}

  dists = dat %>% 
    mutate(
      dist = case_when(
        distance_type == 'phon' ~ phon_dist,
        distance_type == 'jaccard' ~ stringdist::stringdist(test,training,method = 'jaccard'),
        distance_type == 'edit' ~ stringdist::stringdist(test,training,method = 'lv'),
        distance_type == 'qgram' ~ stringdist::stringdist(test,training,method = 'qgram'),
        distance_type == 'hamming' ~ stringdist::stringdist(test,training,method = 'hamming')
      )
    ) %>% 
    mutate(
      pairwise_sim = exp ( - dist / var_s )^var_p,
      ) %>% 
    group_by(test) %>%
    mutate(
      total_sim = sum(pairwise_sim),
      ) %>% 
    group_by(test,category) %>% 
    mutate(
      category_sim = sum(pairwise_sim),
      ) %>% 
    ungroup() %>% 
    mutate(
      category_high = category_sim / total_sim, # cheating, because this is category_low for low. but we're skipping low as there's only two categories and they are complementary.
    ) %>% 
    filter(category == 'high') %>% 
    distinct(
      test,category_high
    )
  return(dists)
}

# K-nearest neighbours. takes the word_distance data frame, distance type, variation type, and variation parameters and returns test words with predictions
KNN = function(dat,distance_type,var_p,var_s,var_k){
  
if(!distance_type %in% c('phon','jaccard','edit','hamming')){
  stop('distance_type must be one of phon, jaccard, edit, hamming')
}
  
  knn = dat %>% 
    mutate(
      dist = case_when(
        distance_type == 'phon' ~ phon_dist,
        distance_type == 'jaccard' ~ stringdist::stringdist(test,training,method = 'jaccard'),
        distance_type == 'edit' ~ stringdist::stringdist(test,training,method = 'lv'),
        distance_type == 'qgram' ~ stringdist::stringdist(test,training,method = 'qgram'),
        distance_type == 'hamming' ~ stringdist::stringdist(test,training,method = 'hamming')
      ),
      category_high = category == 'high'
    ) %>% 
    mutate(
      pairwise_sim = exp ( - dist / var_s )^var_p,
    ) %>% 
    group_by(test) %>% 
    arrange(-pairwise_sim) %>% 
    slice(1:var_k) %>% 
    summarise(category_high = mean(category_high))
  
  return(knn)
  
}

##
# wrappers
##

# wrapper fun for knn
KNNwrapper = function(test,training,feature_matrix,my_distance,my_s,my_k,my_p){
if(my_distance == 'phon'){
  nc = feature_matrix |>
    generateNaturalClasses()
  
  lookup = buildDistTable(feature_matrix, nc) |>
    addLevenshtein()
  
  alignments = runLookup(test,training,lookup)
  
  distances = alignments |>
    dplyr::distinct(test,training,phon_dist) |>
    dplyr::left_join(training, join_by('training' == 'string'))
} else {
  test = test |>
    dplyr::rename(
      test = string
    )
  training = training |>
    dplyr::rename(
      training = string
    )
  distances = tidyr::crossing(test, training) |>
    mutate(
      phon_dist = NA
    )
}

KNN(dat = distances, distance_type = my_distance, var_s = my_s, var_k = my_k, var_p = my_p)
  
}

# wrapper fun for gcm
GCMwrapper = function(test,training,feature_matrix,my_distance,my_s,my_p){
  if(my_distance == 'phon'){
    nc = feature_matrix |>
      generateNaturalClasses()
    
    lookup = buildDistTable(feature_matrix, nc) |>
      addLevenshtein()
    
    alignments = runLookup(test,training,lookup)
    
    distances = alignments |>
      dplyr::distinct(test,training,phon_dist) |>
      dplyr::left_join(training, join_by('training' == 'string'))
  } else {
    test = test |>
    dplyr::rename(
      test = string
    )
  training = training |>
    dplyr::rename(
      training = string
    )
  distances = tidyr::crossing(test, training) |>
    mutate(
      phon_dist = NA
    )
  }
  
  GCM(dat = distances, distance_type = my_distance, var_s = my_s, var_p = my_p)
  
}
