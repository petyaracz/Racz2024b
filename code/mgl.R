# -- fun -- #

# replace characters in pairs with their IPA equivalents
transcribeIPA = function(string, direction){
  if (direction == 'single'){
    stringr::str_replace_all(string, c(
      'ccs' = 'cscs', 'ssz' = 'szsz', 'zzs' = 'zszs', 'tty' = 'tyty', 'ggy' = 'gygy', 'nny' = 'nyny', 'lly' = 'jj', 'cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's'))
  } else if (direction == 'double'){
    stringr::str_replace_all(string, c('s' = 'ß', 'š' = 's', 'ṉ' = 'ny', 'ḏ' = 'gy', 'ṯ' = 'ty', 'ž' = 'zs', 'ß' = 'sz', 'č' = 'cs'))
  }
}

# turn training set into mgl-compatible input
pairsIntoT = function(dat){
  dat |> 
  dplyr::mutate(
    category = ifelse(log_odds >= mean(log_odds), 'high', 'low'),
    input = transcribeIPA(base, 'single'),
    output = dplyr::case_when(
      category == 'high' ~ transcribeIPA(form_1, 'single'),
      category == 'low' ~ transcribeIPA(form_2, 'single')
    ),
    output = glue::glue('#{output}#'),
    input = glue::glue('#{input}#')
  ) |> 
  dplyr::select(input,output)
}

# find string overlap from left
str_overlap_left = function(str1, str2) {
  if (nchar(str1) == 0){
    return('')
  } else if (nchar(str2) == 0){
    return('')
  }
  else if (str1 == str2){
    return(str1)
  } else {
  v1 = str_split_1(str1, '')
  v2 = str_split_1(str2, '')
  max_length = min(length(v1), length(v2))
  v1 = v1[1:max_length]
  v2 = v2[1:max_length]
  overlap = as.list(NULL)
  for (i in 1:max_length){
    if (v1[i] == v2[i]){
      overlap[[i]] = v1[i]
    } else {
      break
    }
  }
  overlap = stringr::str_c(overlap, collapse = '')
  return(overlap)
  }
}

# find string overlap from right
str_overlap_right = function(a, b){
  a2 = stringi::stri_reverse(a)
  b2 = stringi::stri_reverse(b)
  rev_right = str_overlap_left(a2, b2)
  stringi::stri_reverse(rev_right)
}

# delete part b of string a
str_delete = function(a, b){
  if (nchar(b) == 0){
    return(a)
  } else {
    stringr::str_replace(a, b, '') 
  }
}

# format training input to build rules
formatTraining = function(dat){
  dat |>
  dplyr::rowwise() |>
  dplyr::mutate(
    input = glue::glue('#{input}#'),
    output = glue::glue('#{output}#'),
    c = str_overlap_left(input, output),
    rem1 = str_delete(input, c),
    rem2 = str_delete(output, c),
    d = str_overlap_right(rem1, rem2),
    a = str_delete(rem1, d),
    b = str_delete(rem2, d),
  ) |> 
  dplyr::select(-rem1,-rem2)
}

# cross rules on output of formatTraining

findLargerRules = function(t){
  r1 = t |> 
    dplyr::select(a,b,c,d) |> 
    dplyr::rename_all(~ glue::glue('{.}1'))
  
  r2 = t |> 
    dplyr::select(a,b,c,d) |> 
    dplyr::rename_all(~ glue::glue('{.}2'))
  
  all_rules = tidyr::crossing(r1,r2) |> 
    dplyr::filter(a1 == a2 & b1 == b2) |> 
    dplyr::filter(!(c1 == c2 & d1 == d2)) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      C = str_overlap_right(c1,c2),
      D = str_overlap_left(d1,d2)
    ) |>
    dplyr::rename(
      A = a1, 
      B = b1
    ) |> 
    dplyr::distinct(A,B,C,D) |> 
    dplyr::mutate(
      rule = glue::glue('{A} → {B} / {C}__{D}')
    )
  
  return(all_rules)
}

# take output of formatTraining and output of findLargerRules and make rules and words set
getRulesAndWords = function(t,all_rules){
  
  rules_and_words = tidyr::crossing(t,all_rules) |> 
    dplyr::mutate(
      scope = a == A & stringr::str_detect(c, glue::glue('{C}$')) & stringr::str_detect(d, glue::glue('^{D}')),
      hit = scope & b == B
    ) |> 
    dplyr::filter(scope)
}

# take rules and words set and dplyr::count matches and scope for rules, then calc confidence and so on
getRuleStats = function(rules_and_words,alpha_upper,alpha_lower){
  
  all_rules_stats = rules_and_words |> 
    dplyr::count(rule, A, B, C, D, hit) |> 
    tidyr::pivot_wider(names_from = hit, values_from = n, values_fill = 0) |>
    dplyr::rename(hits = `TRUE`, misses = `FALSE`) |> 
    dplyr::mutate(
      scope = hits + misses,
      reliability = hits / scope
    )
  
  all_rules_stats = all_rules_stats |>
    dplyr::rowwise() |>
    dplyr::mutate(
      adjusted_reliability = (hits + .5) / (scope + 1),
      t_lower = qt(p = (1 - alpha_lower)/2, df = scope - 1, lower.tail = F),
      t_upper = qt(p = (1 - alpha_upper)/2, df = scope - 1, lower.tail = F),
      sample_variance = sqrt((adjusted_reliability * (1 - adjusted_reliability)) / scope),
      lower_confidence_limit = adjusted_reliability - t_lower * sample_variance,
      upper_confidence_limit = adjusted_reliability + t_upper * sample_variance
    )
  
  return(all_rules_stats)
}

# find overlap between c1 and c2 but if c1 is empty they do in fact overlap
compareContextsC = function(c1,c2){
  if(nchar(c1) == 0){
    return(T)
  } else {
    return(stringr::str_detect(c1, glue::glue('{c2}$')))
  }
}

# for superset rule1 and its subset rule2 find residue reliability of rule1 for positions where rule2 doesn't apply
getResidueRelability = function(rule1,rule2,rules_and_words){
  
  r1_scope = rules_and_words |> 
    dplyr::filter(rule == rule1) |> 
    dplyr::select(input,output)
  
  r1_hits = rules_and_words |> 
    dplyr::filter(
      rule == rule1, 
      hit
    ) |> 
    dplyr::select(input,output)
  
  r2_scope = rules_and_words |> 
    dplyr::filter(rule == rule2) |> 
    dplyr::select(input,output)
  
  r2_hits = rules_and_words |>
    dplyr::filter(
      rule == rule2, 
      hit
    ) |> 
    dplyr::select(input,output)
  
  scope_residue = dplyr::anti_join(r1_scope,r2_scope, by = join_by(input, output))
  hits_residue = dplyr::anti_join(r1_hits,r2_hits, by = join_by(input, output))
  reliability_residue = (nrow(hits_residue) + .5) / (nrow(scope_residue) + 1)
  return(reliability_residue)
}

# impugn rules. 
impugnRules = function(all_rules_stats,rules_and_words){
  
  all_rules_stats = all_rules_stats |>
    dplyr::mutate(
      context = glue::glue('{C}__{D}')
    )
  
  all_rules_stats_2 = all_rules_stats |> 
    dplyr::rename_all(~ glue::glue('{.}_2'))
  
  rule_cross = tidyr::crossing(all_rules_stats,all_rules_stats_2) |>
    dplyr::filter(
      A == A_2, # same thing turns into...
      B == B_2, # same other thing
      nchar(context) == min(nchar(context)) | nchar(context) < nchar(context_2), # we want the broader rules on the left. so we keep smaller string contexts (smaller context = applies to more forms) or the "everything context"*
      nchar(context_2) != min(nchar(context_2))
    ) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      keep = compareContextsC(C_2,C)
    ) |> 
    dplyr::filter(keep)
  
  # *if like a rule 1 is a -> b / cd_e and rule 2 is a -> b / d_e then rule 2 will apply everywhere where rule 1 applies, since words ending in 'cd' is a subset of words ending in 'd'
  
  # rule 1 is always the bigger rule, since context 'a$' will always superset context 'ba$'
  # now I need to go back to rb and calc the rel(C1-C2) = (hits(C1) - hits(C2))/ (scope(C1) - scope(C2))
  # for that I need to grab hits and scope words for C1 and C2 from all rules w/ words
  
  rule_cross = rule_cross |>
    dplyr::mutate(
      adjusted_residue_reliability = purrr::map2_dbl(rule, rule_2, ~ getResidueRelability(.x,.y,rules_and_words))
    )
  
  rules = rule_cross |> 
    dplyr::mutate(
      upper_residue_confidence_limit = adjusted_residue_reliability + t_upper * sample_variance,
      impugned_lower_confidence_limit = ifelse(
        upper_residue_confidence_limit < lower_confidence_limit, upper_residue_confidence_limit, upper_confidence_limit)
    ) |> 
    dplyr::distinct(A,B,C,D,rule,scope,hits,reliability,adjusted_reliability,lower_confidence_limit,impugned_lower_confidence_limit)
  
  return(rules)
}

# add examples and exceptions
addWordsToRules = function(rules,rules_and_words){
  
  rules_w_e = rules_and_words |>
    dplyr::select(orth,hit,rule) |>
    dplyr::group_by(rule,hit) |>
    tidyr::nest() |>
    dplyr::mutate(
        cat = dplyr::case_when(
          hit ~ 'examples',
          !hit ~ 'exceptions'
          ),
        words = purrr::map_chr(data, ~ .x |> dplyr::pull(orth) |> paste(collapse = ', '))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(rule,cat,words) |>
    tidyr::pivot_wider(names_from = cat, values_from = words, values_fill = '') |>
    dplyr::select(rule,examples,exceptions)
  
  left_join(rules,rules_w_e, by = join_by(rule))
}

# merge rules w/ test and get predictions from best rules

getPredictions = function(test,rules){
  test |> 
    dplyr::mutate(input = glue::glue('#{input}#')) |> 
    tidyr::crossing(rules) |> 
    dplyr::filter(str_detect(input, glue::glue('{C}{A}{D}$'))) |> #!
    dplyr::group_by(input,B) |> 
    dplyr::arrange(-impugned_lower_confidence_limit) |>
    dplyr::slice(1) |> 
    dplyr::select(input,log_odds,B,impugned_lower_confidence_limit) |> 
    tidyr::pivot_wider(names_from = B, values_from = impugned_lower_confidence_limit, values_fill = 0) |> 
    magrittr::set_colnames(c('base','log_odds','cat1','cat2')) |> 
    dplyr::mutate(
      cat1_weight  = dplyr::case_when(
        cat1 == 0 ~ 0,
        cat2 == 0 ~ 1,
        cat1 != 0 & cat2 != 0 ~ cat1 / (cat1 + cat2)
      )
    )
  
}

# -- run the whole thing -- #

mgl = function(training,alpha_lower,alpha_upper){
  
  # check for input column in training
  if(!'input' %in% colnames(training)){
    stop('training data must have a column named "input"')
  }
  # check for output column in training
  if(!'output' %in% colnames(training)){
    stop('training data must have a column named "output"')
  }
  # check for alpha_lower
  if(!is.numeric(alpha_lower)){
    stop('you have to define alpha_lower as a number between 0-1')
  }
  # check for alpha_upper
  if(!is.numeric(alpha_upper)){
    stop('you have to define alpha_upper as a number between 0-1')
  }
  
  # Modeling English Past Tense Intuitions with Minimal Generalization Albright and Hayes 2002
  # also: Automatic Rule Induction for Unknown-Word Guessing Mikheev 1997
  
  # "Formally, the structural change can be represented in the format A → B, and the context in the format / C__D, to yield word-specific rules like those in (2). (The symbol ‘#’ stands for a word boundary.)"
  
  # "The exact procedure for finding a word-specific rule is as follows: given an input pair (X, Y), the model first finds the maximal left-side substring shared by the two forms (e.g., #ms), to create the C term (left side context). The model then exam- ines the remaining material and finds the maximal substring shared on the right side, to create the D term (right side context). The remaining material is the change; the non-shared string from the first form is the A term, and from the second form is the B term."
  
  d = formatTraining(training)
  
  # [Now you need to find overlaps between structural descriptions of rules by comparing every rule to every other rule]
  
  d2 = findLargerRules(d)
  
  # combine rules with words
  
  d3 = getRulesAndWords(d,d2)
  
  # [Now you find scope, and hits for each rule and then calculate lower conf following Mikhailev 1997]
  
  # First, we determine the number of forms in the training data that meet the structural description of the rule (for A → B / C__D, these are the forms that contain CAD). This number is the scope of the rule. The hits of the rule is the number of forms that it actually derives correctly. The reliability of a rule is simply the ratio of its hits to its scope.
  
  # Mikhaleev:
  # The estimate is a good indicator of the rule accuracy but it frequently suffers from large estimation error due to insufficient training data. For example, if a rule was found to apply just once and the total number of observations was also one, its estimate p has the maximal value (1) but clearly this is not a very reliable estimate. We tackle this problem by calculating the lower confidence limit 71"L for the rule estimate, which can be seen as the minimal expected value of/~ for the rule if we were to draw a large number of samples. 
  
  d4 = getRuleStats(d3, alpha_upper, alpha_lower)
  
  # [ we skip phonological rules. for this particular set of uh datasets I've done morpho-phonology beforehand (like, I've removed VH from -ik, split others into suffix classes)]
  
  # formalizing this intuition, we propose a re- finement of the way that confidence is calculated, in order to diagnose when a subpart of a generali- zation is doing most of the work of the larger gen- eralization. When we consider the confidence of a context C associated with a change A → B, we must consider every other context C′ associated with A → B, checking to see whether C′ covers a subset of the cases that C covers. In the present case, when we assess the confidence of adding [-t] after any consonant, we would check all of the other rules adding [-t], including the one that adds [-t] after obstruents. For each C′ that covers a sub- set of C,we must ask whether the rule A→B/C′ is actually “doing most of the work” of the larger rule A→B/C.
  
  # To find out if the smaller rule is doing most of the work, we calculate how well the larger rule (C) performs outside the area covered by the smaller rule (C′). The reliability of the residue area (C -C') is calculated as follows: rel(C-C') = (hits(C) - hits(C'))/ (scope(C) - scope(C')).
  
  # From the reliability of this residue area (C – C′), we can then calculate its confidence, using confidence limit statistics in a way similar to that described above in section 3.2. However, there is a crucial difference: when we are assessing whether a rule explains enough cases to be trustable, we are interested in the denseness of cases within the gen- eralization. But when we are assessing whether a rule offers an improvement over a subpart, we are interested in the sparseness of cases in the residue outside of the subpart. Therefore, when calculating the confidence of the residue, we must use the up- per confidence limit rather than the lower confi-dence limit.
  
  # If the upper confidence limit of the reliability of the residue (C – C′) is lower than the lower confi- dence limit of the reliability of the larger context (C), then we can infer that the smaller rule (A→B/C′) is doing most of the work of the larger rule (A → B / C). Therefore, we penalize the larger rule by replacing its confidence value (Lower confidence(C)) with the confidence value of the residue (Upper confidence(C – C′)). We call this penalty impugnment, because the validity of the larger rule is being called into question by the smaller rule. Impugnment is carried out for all contexts of all rules.
  
  d5 = impugnRules(d4, d3)
  
  # add examples and exceptions to each rule
  
  d6 = addWordsToRules(d5, d3)
  
  return(d6)
}
