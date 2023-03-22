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
val = read_tsv('../data/Supplementary_Table_6A_precision_values.txt')

cq
all_circ


#' # Precision
#' ## Presence in db
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, compound_val, count_group_median, n_db) %>%
  # only keep high abundant circRNAs
  filter(count_group_median == 'count ≥ 5') %>%
  # filter out circ that are NAs for amplicon sequencing validation
  filter(!is.na(compound_val)) %>%
  # change nr of databases to binary
  mutate(n_db = ifelse(is.na(n_db), 0, 1)) %>%
  group_by(n_db) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val)  %>%
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


#' ## Link with detected by multiple tools
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, compound_val, count_group_median, n_detected) %>%
  unique() %>%
  filter(count_group_median == 'count ≥ 5') %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  # change nr of tools to binary
  mutate(n_detected = ifelse(n_detected == 1, 0, 1)) %>%
  group_by(n_detected) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
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



#' ## Link with nr of exons
#' ### Single exon vs multi-exon
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, compound_val, count_group_median, nr_exons) %>%
  unique() %>%
  filter(count_group_median == 'count ≥ 5') %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  # also remove the ones that do not have a exon annotation
  filter(!is.na(nr_exons), !nr_exons == "ambiguous") %>%
  mutate(exon_bin = ifelse(nr_exons == 1, 0, 1)) %>%
  group_by(exon_bin) %>%
  count(compound_val) %>% 
  pivot_wider(values_from = n, names_from = compound_val) %>%
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

#' ### Is 'exon 1' included in the circRNA?
cq_start = cq %>% 
  mutate(start_exon_nr = substr(start_match, 19, 27),
         start_exon_nr = ifelse(substr(start_exon_nr, 1, 1) == '_',
                                substr(start_exon_nr, 2, 10),
                                start_exon_nr))
cq_start 

cont_table = cq_start %>%
  filter(count_group_median == 'count ≥ 5') %>%
  select(circ_id, cell_line, compound_val, start_exon_nr) %>%
  unique() %>%
  filter(!is.na(start_exon_nr)) %>%
  mutate(exon_1 = ifelse(start_exon_nr == "exon_1", 1, 0)) %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  group_by(exon_1) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
  ungroup() 


cont_table

#' => not enough values
#'

#' ## Canonical splicing

cont_table = cq %>%
  filter(count_group_median == 'count ≥ 5',
         !strand == 'unknown') %>%
  select(circ_id, cell_line, compound_val, ss_motif) %>%
  unique() %>%
  mutate(ss_can = ifelse(ss_motif == "AGNGT", 1, 0)) %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  # change nr of databases to binary
  group_by(ss_can) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
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



#' ## Known annotation

cont_table = cq %>%
  filter(count_group_median == "count ≥ 5") %>%
  select(circ_id, cell_line, compound_val, ENST) %>%
  unique() %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  # make ENST match binary
  mutate(ENST_group = ifelse(is.na(ENST), 'no_match', "match")) %>% 
  # change nr of databases to binary
  group_by(ENST_group) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
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
OR = (cont_table[1,2]/cont_table[1,1])/(cont_table[2,2]/cont_table[2,1])
OR

#' ## CircRNA detection tool approach

cont_table = cq %>%
  filter(count_group_median == "count ≥ 5") %>%
  left_join(read_tsv('../data/details/tool_annotation.txt')) %>%
  select(circ_id, cell_line, compound_val, approach) %>%
  filter(!approach == 'integrative') %>%
  unique() %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  group_by(approach) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
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

#' ## BSJ count group
cont_table = cq %>%
  select(circ_id, cell_line, compound_val, count_group_median) %>%
  unique() %>%
  # filter out NAs
  filter(!is.na(compound_val)) %>%
  group_by(count_group_median) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
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

#' # Sensitivity

tool_anno = read_tsv('../data/details/tool_annotation.txt')
sens = read_tsv('../data/Supplementary_Table_6B_sensitivity_values.txt')

#' add annotation to sensitivity
sens_anno = sens %>% 
  left_join(tool_anno) %>%
  select(-tool_lt) %>% ungroup() %>%
  filter(count_group_median == 'count ≥ 5')

sens_anno

cq %>% select(circ_id, compound_val, cell_line) %>% unique() %>% count(compound_val)

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



#' # Cell lines - statistics
#' ## Chi-squared
#' make count table
cont_table = cq %>%
  select(circ_id, cell_line, compound_val, count_group_median) %>%
  unique() %>%
  # only keep high abundant circRNAs
  filter(count_group_median == 'count ≥ 5') %>%
  # to use all val together
  filter(!is.na(compound_val)) %>%
  # change nr of databases to binary
  group_by(cell_line) %>%
  count(compound_val) %>%
  pivot_wider(values_from = n, names_from = compound_val) %>%
  ungroup() %>%
  select(-cell_line)


cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("HLF", "NCI-H23", "SW480")

cont_table

#' plot data
mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

#' if chisquare expected < 5 => then Fisher should be used
chisq.test(cont_table)$expected

#' chisquare test
chisq.test(cont_table)

#' ## ANOVA - Sensitivity
#' first calculate the total nr of validated circ per count group

nr_val_cl = cq %>% 
  # get set of uniquely validated circRNAs
  filter(compound_val == 'pass') %>%
  select(circ_id, cell_line, count_group_median) %>% unique() %>%
  group_by(count_group_median, cell_line) %>%
  summarise(nr_expected = n())

nr_val_cl

#' then calculate sensitivity by dividing nr of circ found by total
sens_cl = cq %>% 
  # get set of uniquely validated circRNAs
  filter(compound_val == 'pass') %>%
  select(circ_id, cell_line, count_group_median) %>% unique() %>%
  # check witch tools have detected these
  left_join(all_circ %>%
              select(tool, circ_id, cell_line) %>% unique()) %>%
  group_by(tool, count_group_median, cell_line) %>% 
  summarise(nr_detected = n()) %>%
  left_join(nr_val_cl) %>%
  mutate(sensitivity = nr_detected/nr_expected) %>%
  ungroup()

#' then ANOVA
sens_cl_anova = sens_cl %>%
  filter(count_group_median == 'count ≥ 5') %>%
  select(tool, cell_line, sensitivity)

sens_cl_anova

sens_cl_anova$cell_line = factor(sens_cl_anova$cell_line, levels = c('HLF', 'NCI-H23', 'SW480'))
sens_cl_anova$tool = factor(sens_cl_anova$tool, levels = c("CIRCexplorer3", "CirComPara2", "circRNA_finder", "circseq_cup",  "CircSplice", "circtools","CIRI2", "CIRIquant", "ecircscreen","find_circ",  "KNIFE",  "NCLcomparator",  "NCLscan", "PFv2","Sailfish-cir", "segemehl"))

res.aov = aov(sensitivity ~ cell_line+tool, data = sens_cl_anova)

summary(res.aov)

0.0036/(0.0036+1.8085+0.0427)
1.8085/(0.0036+1.8085+0.0427)


#' ## ANOVA - precision

#' recalculate the precision per cell line
val_cl = cq %>% 
  select(tool, circ_id, cell_line, count_group, qPCR_val, RR_val, RR_val_detail, amp_seq_val, amp_seq_val_detail, compound_val) %>%
  group_by(tool, count_group, cell_line) %>%
  summarise(nr_qPCR_total = n(),
            nr_qPCR_fail = sum(qPCR_val == 'fail'),
            nr_qPCR_val = sum(qPCR_val == 'pass'),
            nr_RR_total = n() - sum(is.na(RR_val)),  # here NA are the ones that have are 'out_of_range'
            nr_RR_fail = sum(RR_val == "fail", na.rm = T),
            nr_RR_val = sum(RR_val == "pass", na.rm = T),
            nr_amp_total = n() - sum(is.na(amp_seq_val)),  # here NA are the ones 'not_included'
            nr_amp_fail = sum(amp_seq_val == "fail", na.rm = T),
            nr_amp_val = sum(amp_seq_val == "pass", na.rm = T),
            nr_compound_total = n() - sum(is.na(amp_seq_val)), # here NA are the ones 'not_included
            nr_compound_fail = sum(compound_val == "fail", na.rm = T),
            nr_compound_val = sum(compound_val == "pass", na.rm = T)) %>%
  mutate(perc_qPCR_val = nr_qPCR_val/nr_qPCR_total,
         perc_RR_val = nr_RR_val/nr_RR_total,
         perc_amp_val = nr_amp_val/nr_amp_total,
         perc_compound_val = nr_compound_val/nr_compound_total) %>%
  ungroup()

#' ANOVA
val_cl_anova = val_cl %>%
  filter(count_group == 'count ≥ 5') %>%
  select(tool, cell_line, perc_qPCR_val, perc_RR_val, perc_amp_val, perc_compound_val)

val_cl_anova

val_cl_anova$cell_line = factor(val_cl_anova$cell_line, levels = c('HLF', 'NCI-H23', 'SW480'))
val_cl_anova$tool = factor(val_cl_anova$tool, levels = c("CIRCexplorer3", "CirComPara2", "circRNA_finder", "circseq_cup",  "CircSplice", "circtools","CIRI2", "CIRIquant", "ecircscreen","find_circ",  "KNIFE",  "NCLcomparator",  "NCLscan", "PFv2","Sailfish-cir", "segemehl"))


#' ### RT-qPCR precision
res.aov = aov(perc_qPCR_val ~ cell_line+tool, data = val_cl_anova)
summary(res.aov)

0.008126/(0.008126+0.024932+0.017273)
0.024932/(0.008126+0.024932+0.017273)


#' ### RR precision
res.aov = aov(perc_RR_val ~ cell_line+tool, data = val_cl_anova)
summary(res.aov)

0.026300/(0.026300+0.011436+0.003043)
0.011436/(0.026300+0.011436+0.003043)


#' ### ampl seq precision
res.aov = aov(perc_amp_val ~ cell_line+tool, data = val_cl_anova)
summary(res.aov)

0.00912/(0.00912+0.09004+0.00279)
0.09004/(0.00912+0.09004+0.00279)


#' ### compound precision
res.aov = aov(perc_compound_val ~ cell_line+tool, data = val_cl_anova)
summary(res.aov)

0.06205/(0.06205+0.09406+0.00359)
0.09406/(0.06205+0.09406+0.00359)

#' # RNase R enrichment in sequencing data

treatment = read_tsv('../data/Supplementary_Table_5_RNase_R_enrichment_seq.txt')


cont_table = treatment %>%
  # remove Sailfish-cir as it does not report counts
  filter(!tool == 'Sailfish-cir') %>%
  # remove NAs
  filter(!is.na(count_UT), !is.na(count_T)) %>%
  # add median count group
  left_join(all_circ %>% 
              select(circ_id_strand, count_group_median, cell_line) %>% unique()) %>%
  select(circ_id_strand, cell_line, tool, enrichment_bin, count_group_median) %>% 
  group_by(count_group_median) %>%
  count(enrichment_bin) %>%
  pivot_wider(values_from = n, names_from = enrichment_bin) %>%
  ungroup() %>%
  select(-count_group_median)

cont_table = as.data.frame(cont_table)
rownames(cont_table) <- c("count < 5", "count ≥ 5")

cont_table

mosaicplot(cont_table,
           main = "Mosaic plot",
           color = TRUE)

chisq.test(cont_table)$expected

chisq.test(cont_table)
