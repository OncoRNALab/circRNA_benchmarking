
#' ---
#' title: "Calculate validation for 1500 circRNAs"
#' author: "Marieke Vromman"
#' output: 
#'    html_document:
#'       toc: TRUE
#'       toc_float: TRUE
#'       theme: paper
#'       df_print: paged
#'       highlight: tango
#' ---
#' 

#' # File set-up

#' ## Set working directory to current directory
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

#' ## Load standard libraries and resolves conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')


#' ## Read data
#' table with all measured Cq values
cq = read_tsv("../data/details/selected_circ_cq.txt")
cq

#' # Label each circ for three methods
#' ## qPCR validation

cq = cq %>%
  # qPCR_val is 'pass' for all circ that have at least one replicate > 10
  mutate(qPCR_val = ifelse((Cq_min_untreated < 10) & (Cq_max_untreated < 10), 'fail', 'pass'),
         qPCR_val = ifelse((is.na(Cq_min_untreated)) & (is.na(Cq_max_untreated)), 'fail', qPCR_val),
         # qPCR_val_detail describes the reason why a circ fails this validation step if it does
         qPCR_val_detail = ifelse((Cq_min_untreated < 10) & (Cq_max_untreated < 10), 'all_Cq_under_10', 'pass'),
         qPCR_val_detail = ifelse(is.na(Cq_min_untreated) & is.na(Cq_max_untreated), 'all_NA', qPCR_val_detail))

cq

#' statistics
cq %>% count(qPCR_val)
cq %>% count(qPCR_val_detail)


#' ## RR val

cq = cq %>%
  # calculate delta cq
  # first make temp 47 for NAs
  mutate(Cq_min_treated_tmp = Cq_min_treated,
         Cq_min_treated = ifelse(is.na(Cq_min_treated), 47, Cq_min_treated),
         Cq_diff = Cq_min_treated - Cq_max_untreated) %>% 
  mutate(# fail if diff > 3
         RR_val = ifelse(Cq_diff > 3, "fail", 'pass'),
         RR_val_detail = ifelse(RR_val == 'fail', "FP (Cq_diff > 3)", 'pass'),
         # when Cq untreated is high => put NA
         RR_val = ifelse(Cq_min_untreated > 32, NA, RR_val),
         RR_val_detail = ifelse(Cq_min_untreated > 32, 'out_of_range', RR_val_detail),
         # if qPCR already failed => RR also becomes 'fail'
         RR_val = ifelse(qPCR_val == "fail", 'fail', RR_val),
         RR_val_detail = ifelse(qPCR_val == "fail", "fail_qPCR", RR_val_detail)) %>%
  # clean up NAs for 47
  mutate(Cq_diff = ifelse(Cq_min_treated == 47, NA, Cq_diff),
         Cq_min_treated = Cq_min_treated_tmp) %>%
  select(-Cq_min_treated_tmp) %>%
  #' to clean up dataframe: also remove Cq_diff for 'all_Cq_under_10' group and out_of_range group
  mutate(Cq_diff = ifelse(qPCR_val_detail == 'all_Cq_under_10', NA, Cq_diff),
         Cq_diff = ifelse(RR_val_detail == 'out_of_range', NA, Cq_diff))

cq

#' statistics
cq %>% count(RR_val)
cq %>% count(RR_val_detail)


#' ## Amplicon sequencing
cq_amp = read_tsv('../data/details/amp_seq.txt')

#' statistics
cq_amp %>% count(amp_seq_val)
cq_amp %>% count(amp_seq_val_detail)


#' add to cq dataframe
cq = cq %>% 
  left_join(cq_amp %>% unique()) %>%
  mutate(amp_seq_val_detail = ifelse(is.na(amp_seq_val_detail), "not_included", amp_seq_val_detail))

#' statistics
cq %>% count(amp_seq_val)
cq %>% count(amp_seq_val_detail)


#' ## Add circRNA annotation and save as Supplementary Table 3
all_circ = read_tsv("../data/Supplementary_Table_2_all_circRNAs.txt")

all_circ

cq %>%
  left_join(all_circ) %>%
  write_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')

#' # Calculate val rates (precision) per method per tool
perc_val = cq %>% 
  select(tool, circ_id, cell_line, count_group, qPCR_val, RR_val, RR_val_detail, amp_seq_val, amp_seq_val_detail) %>%
  group_by(tool, count_group) %>%
  summarise(nr_qPCR_total = n(),
            nr_qPCR_fail = sum(qPCR_val == 'fail'),
            nr_qPCR_val = sum(qPCR_val == 'pass'),
            nr_RR_total = n() - sum(is.na(RR_val)),  # here NA are the ones that have are 'out_of_range'
            nr_RR_fail = sum(RR_val == "fail", na.rm = T),
            nr_RR_val = sum(RR_val == "pass", na.rm = T),
            nr_amp_total = n() - sum(is.na(amp_seq_val)),  # here NA are the ones 'not_included'
            nr_amp_fail = sum(amp_seq_val == "fail", na.rm = T),
            nr_amp_val = sum(amp_seq_val == "pass", na.rm = T)) %>%
  mutate(perc_qPCR_val = nr_qPCR_val/nr_qPCR_total,
         perc_RR_val = nr_RR_val/nr_RR_total,
         perc_amp_val = nr_amp_val/nr_amp_total) %>%
  ungroup()

perc_val

#' # Calculate compound val rate per tool
#' ignore NAs and multiply
perc_val_compound = cq %>% 
  select(tool, circ_id, cell_line, count_group, qPCR_val, RR_val, RR_val_detail, amp_seq_val, amp_seq_val_detail) %>%
  group_by(tool, count_group) %>%
  summarise(nr_qPCR_total = n(),
            nr_qPCR_fail = sum(qPCR_val == 'fail'),
            nr_qPCR_val = sum(qPCR_val == 'pass'),
            nr_RR_total = n() - sum(is.na(RR_val)) - sum(RR_val_detail == "fail_qPCR"), # NA = 'out_of_range' + 'qPCR_fail'
            nr_RR_fail = sum(RR_val == "fail", na.rm = T),
            nr_RR_val = sum(RR_val == "pass", na.rm = T),
            nr_amp_total = n() - sum(is.na(amp_seq_val) | RR_val_detail == "fail_qPCR"), # NA = 'not_included' + 'qPCR_fail'
            nr_amp_fail = sum(amp_seq_val == "fail", na.rm = T),
            nr_amp_val = sum(amp_seq_val == "pass", na.rm = T)) %>%
  mutate(perc_qPCR_val = nr_qPCR_val/nr_qPCR_total,
         perc_RR_val = nr_RR_val/nr_RR_total,
         perc_amp_val = nr_amp_val/nr_amp_total,
         perc_compound_val = perc_qPCR_val*perc_RR_val*perc_amp_val) %>%
  ungroup()

perc_val_compound %>% select(tool, count_group, perc_compound_val) %>% view()
 
#' add to other dataframe with separate validation rates
perc_val = perc_val %>%
  left_join(perc_val_compound %>% select(tool, count_group, perc_compound_val))

perc_val

#' ## Per tool val rates => stats
perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_RR_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_qPCR_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_amp_val) %>%
  quantile()

#' # Calculate sensitivity
sens = cq %>% 
  # get set of uniquely validated circRNAs
  filter(qPCR_val == 'pass',
         RR_val == 'pass',
         amp_seq_val == 'pass') %>%
  select(circ_id, cell_line) %>% unique() %>%
  # check witch tools have detected these
  left_join(all_circ %>%
              select(tool, circ_id, cell_line) %>% unique()) %>%
  group_by(tool) %>% 
  summarise(nr_detected = n(),
            sensitivity = nr_detected/957)

sens

#' add sensitivity to dataframe
perc_val = perc_val %>%
  left_join(sens %>% select(tool, sensitivity))


#' # Save dataframe (Sup Table 4)
perc_val %>% write_tsv('../data/Supplementary_Table_4_precision_values.txt')



#' # Overall validation rates (ignoring tools)
#' ## Calculate overall validation rates
#' RT-qPCR
cq %>%
  select(circ_id, cell_line, qPCR_val) %>% unique() %>%
  count(qPCR_val)

round(1479/(1479+37), digits = 3)

#' RNase R
cq %>%
  select(circ_id, cell_line, RR_val) %>% unique() %>%
  count(RR_val)

round(1319/(1319+85), digits = 3)
round(112/(1319+85+112), digits = 3)

#' amplicon sequencing
cq %>%
  select(circ_id, cell_line, amp_seq_val) %>% unique() %>%
  count(amp_seq_val)

cq %>%
  select(circ_id, cell_line, amp_seq_val_detail) %>% unique() %>%
  count(amp_seq_val_detail)
round(1014/(1014+165), digits = 3)
round(337/(1014+165+337), digits = 3)


#' ## Calculate perc of circRNAs upsetplot

#' nr of unique circ that pass all val methods
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  filter(qPCR_val == 'pass', RR_val == 'pass', amp_seq_val == 'pass')

round(957/1516, digits = 3)

#' nr of unique circ that fail all val methods
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  filter(qPCR_val == 'fail', RR_val == 'fail', amp_seq_val == 'fail')

round(18/1516, digits = 3)


#' nr of unique circ that fail 1 or 2 val methods

cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  filter(!is.na(RR_val), !is.na(amp_seq_val)) %>%
  unique() %>%
  mutate(val_nr = str_count(paste(qPCR_val, RR_val, amp_seq_val),
                              pattern = 'pass')) %>%
  #filter(val_nr == 0)
  #filter(val_nr == 3)
  filter(val_nr < 3, val_nr > 0)

round(128/1516, digits = 3)


#' nr of cicrc that have at least one NA
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  mutate(NA_nr = str_count(paste(qPCR_val, RR_val, amp_seq_val),
                            pattern = 'NA')) %>%
  #filter(NA_nr == 0)
  #filter(NA_nr == 3)
  filter(NA_nr > 0)

round(413/1516, digits = 3)


#' nr of unique circ that have no NAs
cq %>% 
  filter(!is.na(RR_val), !is.na(amp_seq_val)) %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line)  %>%
  unique()

#' => recalculate previous numbers based on no NAs
#' all pass
round(957/1103, digits = 3)
#' all fail
round(18/1103, digits = 3)
#' one or two fail
round(128/1103, digits = 3)


#' ## Calculate  val rate < 5
#' RT-qPCR
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, qPCR_val) %>% unique() %>%
  count(qPCR_val)

round(243/(243+17), digits = 3)

#' RNase R
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, RR_val) %>% unique() %>%
  count(RR_val)

round(165/(165+25), digits = 3)

#' amplicon sequencing
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, amp_seq_val) %>% unique() %>%
  count(amp_seq_val)

round((63+136)/(63+136+61), digits = 3)

round(136/(63+136), digits = 3)

