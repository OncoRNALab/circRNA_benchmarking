#' ---
#' title: "Link between circRNA annotation and validation rates"
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

#' ## Load standard libraries and resolve conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')


#' ## Read data
cq = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')
all_circ = read_tsv("../data/Supplementary_Table_2_all_circRNAs.txt")
val = read_tsv('../data/Supplementary_Table_4_precision_values.txt')

cq
all_circ


#' # Presence in db
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, count_group_median, n_db) %>%
  unique() %>%
  # only keep high abundant circRNAs
  filter(count_group_median == 'count ≥ 5') %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 'pass', 'fail')) %>%
  # change nr of databases to binary
  mutate(n_db = ifelse(is.na(n_db), 0, 1)) %>%
  group_by(n_db) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-n_db)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("novel", "in_db")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR


#' # Link with detected by multiple tools
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, count_group_median, n_detected) %>%
  unique() %>%
  filter(count_group_median == 'count ≥ 5') %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 'pass', 'fail')) %>%
  # change nr of databases to binary
  mutate(n_detected = ifelse(n_detected == 1, 0, 1)) %>%
  group_by(n_detected) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-n_detected)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("unique", "multiple_tools")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)


#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR



#' # Link with nr of exons
#' ## Single exon vs multi-exon
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, count_group_median, nr_exons) %>%
  unique() %>%
  filter(count_group_median == 'count ≥ 5') %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 'pass', 'fail')) %>%
  # also remove the ones that do not have a exon annotation
  filter(!is.na(nr_exons), !nr_exons == "ambiguous") %>%
  mutate(exon_bin = ifelse(nr_exons == 1, 0, 1)) %>%
  group_by(exon_bin) %>%
  count(all_val) %>% 
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  mutate(fail = ifelse(is.na(fail), 0, fail)) %>%
  select(-exon_bin)

cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("single_exon", "multi_exons")
cont_table


#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)


#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR

#' ## Is 'exon 1' included in the circRNA?
cq_start = cq %>% 
  mutate(start_exon_nr = substr(start_match, 19, 27),
         start_exon_nr = ifelse(substr(start_exon_nr, 1, 1) == '_',
                                substr(start_exon_nr, 2, 10),
                                start_exon_nr))
cq_start 

cont_table = cq_start %>%
  filter(count_group_median == 'count ≥ 5') %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, start_exon_nr) %>%
  unique() %>%
  filter(!is.na(start_exon_nr)) %>%
  mutate(exon_1 = ifelse(start_exon_nr == "exon_1", 1, 0)) %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 
                          'pass', 'fail')) %>%
  # change nr of databases to binary
  group_by(exon_1) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() 


cont_table

#' => not enough values
#'

#' # Canonical splicing

cont_table = cq %>%
  filter(count_group_median == 'count ≥ 5',
         !strand == 'unknown') %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, ss_motif) %>%
  unique() %>%
  mutate(ss_can = ifelse(ss_motif == "AGNGT", 1, 0)) %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 
                          'pass', 'fail')) %>%
  # change nr of databases to binary
  group_by(ss_can) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-ss_can)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("non_cannonical", "cannonical")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR



#' # Known annotation

cont_table = cq %>%
  filter(count_group_median == "count ≥ 5") %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, ENST) %>%
  unique() %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 
                          'pass', 'fail'),
         ENST_group = ifelse(is.na(ENST), 'no_match', "match")) %>% 
  # change nr of databases to binary
  group_by(ENST_group) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-ENST_group)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("match", "no_match")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR

#' # CircRNA detection tool approach

cont_table = cq %>%
  filter(count_group_median == "count ≥ 5") %>%
  left_join(read_tsv('../data/details/tool_annotation.txt')) %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, approach) %>%
  filter(!approach == 'integrative') %>%
  unique() %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 
                          'pass', 'fail')) %>% 
  # change nr of databases to binary
  group_by(approach) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-approach)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("candidate-based", "segmented read-based")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' OR
OR = (cont_table[1,2]/cont_table[1,1])/(cont_table[2,2]/cont_table[2,1])
OR

#' # BSJ count group
cont_table = cq %>%
  select(circ_id, cell_line, qPCR_val, RR_val, amp_seq_val, count_group_median) %>%
  unique() %>%
  # to use all val together
  filter(!is.na(amp_seq_val), !is.na(RR_val)) %>%
  mutate(all_val = ifelse(qPCR_val == RR_val & qPCR_val == amp_seq_val & qPCR_val == 'pass', 
                          'pass', 'fail')) %>%
  # change nr of databases to binary
  group_by(count_group_median) %>%
  count(all_val) %>%
  pivot_wider(values_from = n, names_from = all_val) %>%
  ungroup() %>%
  select(-count_group_median)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("count < 5", "count ≥ 5")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' OR
OR = (cont_table[2,2]/cont_table[2,1])/(cont_table[1,2]/cont_table[1,1])
OR

#' # Using sensitivity
#' 

tool_anno = read_tsv('../data/details/tool_annotation.txt')
sens = read_tsv('../data/Supplementary_Table_5_sensitivity_values.txt')

#' add annotation to sensitivity
sens_anno = sens %>% 
  left_join(tool_anno) %>%
  select(-tool_lt) %>% ungroup() %>%
  filter(count_group_median == 'count ≥ 5')

sens_anno

#' ## Detection approach: segmented read-based VS candidate-based
wilcox.test(sensitivity ~ approach, data=sens_anno %>% filter(!approach == 'integrative')) 

kruskal.test(sensitivity ~ approach, data=sens_anno)

#' ## Detection based on annotation: known splice sites VS entire genome
wilcox.test(sensitivity ~ lin_annotation, data=sens_anno)

#' ## Strand annotation method
kruskal.test(sensitivity ~ strand_anno, data=sens_anno %>% 
               filter(!strand_anno == "no strand reported"))


#' ## Splicing: canonical VS non-canonical
wilcox_ss = wilcox.test(sensitivity ~ splicing, data=sens_anno) 

wilcox_ss

sens_anno %>%
  group_by(splicing) %>%
  summarise(mean_sens = mean(sensitivity))

#' ### method 1
N = nrow(sens_anno)
N
Z = qnorm(wilcox_ss$p.value/2)
Z
r = abs(Z)/sqrt(N)
r

#' ### method 2
library(rstatix)
library(ggpubr)
sens_anno %>% rstatix::wilcox_test(sensitivity ~ splicing)

sens_anno %>% rstatix::wilcox_effsize(sensitivity ~ splicing)

#sens_anno %>% rstatix::wilcox_effsize(sensitivity ~ splicing, ci = T)



#' ### Median difference in sensitivity

median_diff = tibble()

for (tool_can in sens_anno %>% filter(splicing == "canonical") %>% pull(tool)) {
  sens_can = sens_anno %>% filter(tool == tool_can) %>% pull(sensitivity)
  
  for (tool_non_can in sens_anno %>% filter(splicing == "non-canonical") %>% pull(tool)) {
    sens_non_can = sens_anno %>% filter(tool == tool_non_can) %>% pull(sensitivity)
    
    median_diff = median_diff %>%
      bind_rows(tibble(tool_can, sens_can, tool_non_can, sens_non_can))
    
  }
}

median_diff = median_diff %>%
  mutate(sens_diff = sens_can - sens_non_can)

median_diff

median_diff %>% pull(sens_diff) %>% median()


#' ## Filtering
wilcox.test(sensitivity ~ BSJ_filter, data=sens_anno) 


#' ## Correlation estimated sensitivity and theoretical nr of TPs

val_cor = val %>%
  # use perc_compound_val
  select(tool, count_group, perc_compound_val) %>%
  # get the number of detected circRNAs for each cell line and tool, per count group
  left_join(all_circ %>%
              group_by(cell_line, tool, count_group) %>%
              summarise(total_n_ut = n())) %>%
  # calculate the theoretical nr of validated circ
  mutate(theo_TP_all = perc_compound_val * total_n_ut) %>%
  left_join(sens, by = c('count_group' = 'count_group_median', 'tool')) %>%
  group_by(tool, count_group, perc_compound_val, sensitivity) %>%
  summarize(theo_TP_all_cl = sum(theo_TP_all)) %>% ungroup()

val_cor


#' below 5
val_cor_tmp = val_cor %>% filter(count_group == "count < 5")
cor.test(val_cor_tmp$sensitivity, val_cor_tmp$theo_TP_all_cl, method = 'spearman')

#' above 5
val_cor_tmp = val_cor %>% filter(count_group == "count ≥ 5")
cor.test(val_cor_tmp$sensitivity, val_cor_tmp$theo_TP_all_cl, method = 'spearman')
