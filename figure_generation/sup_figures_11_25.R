#' ---
#' title: "Generation of sup figures"
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

#' ## Load specific libraries
library(gplots)
library(ggpubr)
library("quantreg") 
library(ggrepel)

#' ## Set figure theme
mytheme = theme_bw(base_size = 10) + 
  theme(text = element_text(size=10, colour='black'),
        title = element_text(size=10, colour='black'),
        line = element_line(size=0.5),
        axis.title = element_text(size=10, colour='black'),
        axis.text = element_text(size=10, colour='black'),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size=0.5),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.8, 0.8), 
        legend.text=element_text(size=10)) 

mytheme_discrete_x = mytheme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# a colorblind-friendly palette 
# colorblind.palette.grey = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#' ## Read data and order
cq = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')
all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')
val = read_tsv('../data/Supplementary_Table_4_precision_values.txt')
sens = read_tsv('../data/Supplementary_Table_5_sensitivity_values.txt')

all_circ$tool = factor(all_circ$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))

cq$tool = factor(cq$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))

val$tool = factor(val$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))

sens$tool = factor(sens$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))


cq
all_circ


#' # Sup Figure 11: pilot exp for detectability of low abundance circRNAs
cq_test = read_tsv('../data/details/abundance_exp.txt') 

cq_test

cq_cum = cq_test %>%
  select(Cq_mean_PCR, BSJ_count, circ_id) %>%
  unique() %>%
  group_by(BSJ_count) %>%
  arrange(desc(BSJ_count)) %>%
  summarize(total_n = n(), val_n = sum(Cq_mean_PCR < 32)) %>%
  ungroup() %>%
  arrange(desc(BSJ_count)) %>%
  mutate(total_n_cum = cumsum(total_n), 
         val_n_cum = cumsum(val_n)) %>%
  mutate(perc_val = val_n_cum/total_n_cum)

cq_cum

cq_cum %>% 
  ggplot(aes(BSJ_count, perc_val)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 5, color = "#D55E00") +
  geom_smooth(color = '#5AB4E5') +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  mytheme

#ggsave('sup_figures/sup_figure_11.pdf', width = 10, height = 8.5, units = "cm")


#' # Sup Figure 12: primer design succes rate

primer_design = read_tsv('../data/details/primer_design.txt')

# showing how many primers failed

primer_design$tool = factor(primer_design$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

primer_design %>%
  mutate(fail_reason = ifelse(total_primer_pairs == failed_spec, 'failed specificity', NA),
         fail_reason = ifelse(total_primer_pairs == failed_SNP, 'failed SNPs', fail_reason),
         fail_reason = ifelse(total_primer_pairs == failed_sec_str_amp, 'failed sec str amp', fail_reason),
         fail_reason = ifelse(total_primer_pairs == failed_sec_str_temp, 'failed sec str temp', fail_reason),
         fail_reason = ifelse(is.na(fail_reason), 'failed sec str amp', fail_reason),
         fail_reason = ifelse(design == 0, 'no_design', fail_reason),
         fail_reason = ifelse(primer_found == 1, 'succesful_design', fail_reason)) %>%
  mutate(count_group = ifelse(BSJ_count > 4, 'count ≥ 5', 'count < 5'),
         count_group = ifelse(tool == 'Sailfish-cir', 'no_counts', count_group)) %>%
  ggplot(aes(tool, fill = fail_reason)) +
  scale_fill_manual(values = c("#E69F00", '#CC79A7', '#00B9F2', "#009E73")) +
  geom_bar(position = 'fill') +
  facet_grid(~count_group, scales = 'free', space = 'free') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x + 
  xlab('') + 
  ylab('% of circRNAs') +
  theme(legend.position = 'right')

#ggsave('sup_figures/sup_figure_12.pdf',  width = 20, height = 10, units = "cm")

# numbers
primer_design %>%
  filter(primer_found == 0) %>%
  mutate(fail_reason = ifelse(total_primer_pairs == failed_spec, 'failed specificity', NA),
         fail_reason = ifelse(total_primer_pairs == failed_SNP, 'failed SNPs', fail_reason),
         fail_reason = ifelse(total_primer_pairs == failed_sec_str_amp, 'failed sec str amp', fail_reason),
         fail_reason = ifelse(total_primer_pairs == failed_sec_str_temp, 'failed sec str temp', fail_reason),
         fail_reason = ifelse(is.na(fail_reason), 'failed sec str amp', fail_reason),
         fail_reason = ifelse(design == 0, 'no_design', fail_reason),
         fail_reason = ifelse(primer_found == 1, 'succesful_design', fail_reason)) %>%
  group_by(tool) %>%
  summarise(perc_ok = sum(fail_reason == 'failed specificity')/n()) %>%
  pull(perc_ok) %>% quantile()

#' # Sup Figure 13: distribution of BSJ counts of selected circRNAs

cq %>%
  filter(!tool == 'Sailfish-cir') %>%
  ggplot(aes(tool, BSJ_count)) +
  geom_boxplot() +
  #geom_point(size = 0.5) +
  #geom_jitter() +
  #facet_wrap(~count_group, scales = 'free') + 
  mytheme_discrete_x +
  scale_y_log10() +
  ylab('BSJ count') +
  xlab('')

#ggsave('sup_figures/sup_figure_13.pdf', width = 3, height = 10, units = "cm")

#' # Sup Figure 14: schematic overview of selected circRNAs
#' see illustator
#' 

#' # Sup Figure 15: summary of validation results (heatmap style)
#' make one validation label for each circ (that combines the 3 val methods)
#' this will be used later for ordering
#' and make the dataframe longer
hm = cq %>%
  select(tool, circ_id, count_group, qPCR_val, RR_val, amp_seq_val) %>% 
  mutate(all_val = paste(qPCR_val, RR_val, amp_seq_val, sep = '_')) %>%
  pivot_longer(cols = c(qPCR_val, RR_val, amp_seq_val),
               values_to = 'val', names_to = "val_type") 

hm

#' add 3 * 20 empty values for the 2 tools that don't report circRNA with BSJ count < 5
#' (this would probably have been more efficient to do before the pivot_longer)

hm = hm %>%
  bind_rows(tibble(tool = "circRNA_finder", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "qPCR_val")) %>%
  bind_rows(tibble(tool = "segemehl", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "qPCR_val")) %>%
  bind_rows(tibble(tool = "circRNA_finder", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "RR_val")) %>%
  bind_rows(tibble(tool = "segemehl", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "RR_val")) %>%
  bind_rows(tibble(tool = "circRNA_finder", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "amp_seq_val")) %>%
  bind_rows(tibble(tool = "segemehl", circ_id = paste('test', 1:20), count_group = 'count < 5', val = 'NANANA', 
                   val_type = "amp_seq_val"))

hm


#' add one extra empty line for the tools that do report BSJ count < 5 that can then separate both groups
for (tool_name in cq %>% #filter(count_group == 'count < 5') %>% 
     pull(tool) %>% unique()) {
  hm = hm %>% 
    bind_rows(tibble(tool = tool_name, circ_id = 'extra_line', count_group = 'extra_line', 
                     val = "NANANA", val_type = "qPCR_val")) %>%
    bind_rows(tibble(tool = tool_name, circ_id = 'extra_line', count_group = 'extra_line', 
                     val = "NANANA", val_type = "RR_val")) %>%
    bind_rows(tibble(tool = tool_name, circ_id = 'extra_line', count_group = 'extra_line', 
                     val = "NANANA", val_type = "amp_seq_val"))
}

hm

#' check if every tool has the same number of lines
hm %>% count(tool)

#' change label of specific line sailfish-cir (does not report counts)
hm = hm %>%
  mutate(count_group = ifelse(tool == 'sailfish-cir' & count_group == "extra_line",
                              'extra_line2', count_group))

#' set everything as factors in the right order

hm$val_type = factor(hm$val_type, levels = c("qPCR_val", 'RR_val', 'amp_seq_val')) #"pass_NA_NA",
hm$all_val = factor(hm$all_val, levels = c("fail_fail_fail", "fail_fail_NA", "fail_fail_pass", "pass_fail_fail", "pass_fail_NA", "pass_NA_fail", 
                                           "pass_NA_NA", "pass_fail_pass", "pass_pass_fail","pass_NA_pass", "pass_pass_NA", "pass_pass_pass"))
hm$count_group = factor(hm$count_group, levels = c('count ≥ 5', "extra_line",'count < 5', 'no_counts', 'extra_line2'))
hm$tool = factor(hm$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl", "Sailfish-cir")) 


#' order per tool and per count group all circ according to their validation (all_val) and add a number to each circRNA
hm = hm %>% 
  group_by(tool) %>%
  arrange(count_group, all_val) %>%
  mutate(nr = row_number()) %>%
  ungroup() %>%
  # to have only one number per circ (instead of 3), take the mean (this step can be skipped and gives same results, but I though it might solve the random points)
  group_by(circ_id, tool) %>%
  mutate(nr = mean(nr)) %>%
  ungroup()

hm$nr = factor(hm$nr, 
               levels = c("2", "5", "8", "11", "14", "17", "20", "23", "26", "29",
                          "32", "35", "38", "41", "44", "47", "50", "53", "56", 
                          "59", "62", "65", "68", "71", "74", "77", "80", "83", 
                          "86", "89", "92", "95","98", "101", "104", "107", "110", 
                          "113", "116", "119", "122", "125", "128", "131", "134", 
                          "137", "140", "143", "146", "149", "152", "155", "158", 
                          "161", "164", "167", "170", "173", "176", "179", "182", 
                          "185", "188", "191", "194", "197", "200", "203", "206", 
                          "209", "212", "215", "218", "221", "224", "227", "230", 
                          "233", "236", "239", "242", "245", "248", "251", "254", 
                          "257", "260", "263", "266", "269", "272", "275", "278", 
                          "281", "284", "287", "290", "293", "296", "299", "302",
                          "264", "265", "267", "268", "270", "271", "273", "274",
                          "276", "277", "279", "280", "282", "283"))


#' plot the whole thing
hm %>%
  ggplot(aes(val_type, nr)) + 
  geom_tile(aes(fill = val), colour = "white") + 
  scale_fill_manual(values=c('#CC79A7', "#ffffff", '#009E73', '#999999')) +
  facet_wrap(~tool, ncol = 16, scales = 'free_y') +
  #facet_grid(rows = vars(count_group), cols = vars(tool), scales = 'free', space = 'free') +
  mytheme_discrete_x +
  ylab('') +
  xlab('') +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right")

hm %>% count(tool)

#ggsave('sup_figures/sup_figure_15.pdf',  width = 23, height = 13, units = "cm")


#' # Sup Figure 16 & 17: same circRNAs in different cell lines

#' number of circRNAs measured in more than one cell line
circ_cl = cq %>% inner_join(cq %>% 
                              select(circ_id, cell_line) %>%
                              group_by(circ_id) %>%
                              unique() %>%
                              filter(n() > 1) %>%
                              select(circ_id) %>% unique()) %>%
  select(circ_id, qPCR_val, qPCR_val_detail, RR_val,
         RR_val_detail, amp_seq_val, amp_seq_val_detail) %>%
  unique() %>%
  #group_by(circ_id) %>% filter(n() == 1)
  select(circ_id) %>% unique()

circ_cl

# 58
round(54/58, 3)
round(4/58, 3)

# plot all Cq values
circ_cl = cq %>% inner_join(circ_cl) %>%
  select(circ_id, cell_line, Cq_min_untreated, Cq_max_untreated) %>% unique() %>%
  mutate(Cq_mean = (Cq_min_untreated + Cq_max_untreated)/2) %>%
  select(-Cq_min_untreated, -Cq_max_untreated)


circ_cl

circ_cl = circ_cl %>%
  # between HLF and NCI-H23
  filter(cell_line == "HLF") %>%
  rename(cell_line_1 = cell_line, Cq_mean_1 = Cq_mean) %>%
  inner_join(circ_cl %>%
               filter(cell_line == "NCI-H23") %>%
               rename(cell_line_2 = cell_line, Cq_mean_2 = Cq_mean)) %>%
  bind_rows(circ_cl %>%
              # between HLF and SW480
              filter(cell_line == "HLF") %>%
              rename(cell_line_1 = cell_line, Cq_mean_1 = Cq_mean) %>%
              inner_join(circ_cl %>%
                           filter(cell_line == "SW480") %>%
                           rename(cell_line_2 = cell_line, Cq_mean_2 = Cq_mean))) %>%
  bind_rows(circ_cl %>%
              # between NCI-H23 and SW480
              filter(cell_line == "NCI-H23") %>%
              rename(cell_line_1 = cell_line, Cq_mean_1 = Cq_mean) %>%
              inner_join(circ_cl %>%
                           filter(cell_line == "SW480") %>%
                           rename(cell_line_2 = cell_line, Cq_mean_2 = Cq_mean)))


circ_cl %>% 
  ggplot(aes(Cq_mean_1, Cq_mean_2)) +
  geom_point() +
  geom_smooth(method = "lm", color = '#00B9F2') +
  facet_wrap(~cell_line_1+cell_line_2, scales = 'free') +
  xlim(17, 32) +
  ylim(17, 32) +
  stat_regline_equation(label.y = 32, label.x=17, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 30, label.x=17, aes(label = ..rr.label..)) +
  xlab('Cq') +
  ylab('Cq') +
  #geom_abline(color = "#00B9F2") +
  mytheme +
  theme(aspect.ratio = 1) 

all_circ %>% select(circ_id) %>% unique()


#ggsave('sup_figure_16.pdf', width = 19, height = 9, units = "cm")

#' plot those 4 circ => sup figure 17
cq_4 = cq %>% filter(circ_id %in% c('chr10:96155325-96160402',
                                    'chr5:177209635-177212195',
                                    'chr5:36982164-36986301',
                                    'chr5:618989-655584')) %>%
  select(circ_id, cell_line, Cq_max_untreated, Cq_min_treated) %>%
  pivot_longer(cols = c(Cq_max_untreated, Cq_min_treated), values_to = 'Cq', names_to = 'sample') %>%
  mutate(Cq = as.numeric(Cq)) 

cq_4$circ_id = factor(cq_4$circ_id, levels = c('chr5:618989-655584', 'chr5:36982164-36986301','chr5:177209635-177212195', 'chr10:96155325-96160402'))
cq_4 %>%
  ggplot(aes(sample, Cq, color = cell_line)) +
  geom_point() +
  scale_color_manual(name = "",
                     values = c("#0072B2", "#E69F00", "#CC79A7")) +
  geom_line(aes(group = cell_line)) +
  facet_grid(~circ_id) +
  mytheme_discrete_x +
  xlab('')

#ggsave('sup_figures/sup_figure_17.pdf', width = 21, height = 10, units = "cm")


#' # Sup Figure 18: relation between BSJ counts and Cq values
count_Cq = all_circ %>%
  inner_join(cq %>%
               select(circ_id, Cq_min_untreated, Cq_max_untreated) %>%
               unique() %>%
               mutate(Cq_mean = rowMeans(select(.,c(Cq_min_untreated, Cq_max_untreated)), na.rm = T)) %>%
               filter(!is.na(Cq_mean)))

count_Cq %>% filter(!tool == "Sailfish-cir") %>% 
  pull(BSJ_count) %>% max()

count_Cq %>%
  filter(!tool == "Sailfish-cir") %>% 
  ggplot(aes(log2(BSJ_count), Cq_mean)) +
  #geom_density(adjust = 2, stat = 'identity') +
  geom_point(size = 0.5, alpha = 0.2, fill=alpha("black", 0.2)) +
  #geom_density(stat = 'identity') +
  geom_smooth(method = "lm", color = '#00B9F2') +
  facet_wrap(~tool, scales = 'free') +
  ylab('mean of Cq values (qPCR replicates)') +
  xlab('log2(BSJ count)') +
  #stat_regline_equation(label.y = 37, label.x=6, aes(label = ..eq.label..)) +
  #stat_regline_equation(label.y = 40, label.x=6, aes(label = ..rr.label..)) +
  #stat_regline_equation(label.y = 42, label.x=6, aes(label = ..p.label..)) +
  stat_regline_equation(label.x = 0, label.y = 10, size = 3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
           label.y= 5, label.x = 0, size = 3) +
  mytheme +
  theme(aspect.ratio = 1) +
  xlim(0,11) + 
  #xlim(-25,17) + 
  ylim(0,45)


#ggsave('sup_figures/sup_figure_18.pdf',  width = 21, height = 21, units = "cm")


#' # Sup Figure 19: cummulative percentage of on-target amplification
map_perc_cum_tool = cq %>%
  select(circ_id_strand, tool, perc_on_target, count_group) %>%
  unique() %>%
  filter(!is.na(perc_on_target)) %>%
  group_by(tool, count_group) %>%
  mutate(total_n = n()) %>% 
  arrange(desc(perc_on_target)) %>%
  mutate(n_cum = 1:n(),
         n_cum_perc = n_cum/total_n)

map_perc_cum_tool

#' plot
map_perc_cum_tool$tool = factor(map_perc_cum_tool$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl", "Sailfish-cir")) 


map_perc_cum_tool %>%
  ggplot(aes(perc_on_target, n_cum_perc, color = count_group)) +
  geom_line() +
  xlab('minimum % on-target amplification') +
  ylab('% of circRNAs') +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme + 
  facet_wrap(~tool) +
  scale_color_manual(values=c( '#56B4E9', '#E69F00', '#999999')) +
  theme(legend.position = 'right') +
  geom_vline(xintercept = 0.5)

#ggsave('sup_figures/sup_figure_19.pdf',  width = 30, height = 20, units = "cm")



#' # Sup Figure 20: compound precision values per tool

val %>% 
  ggplot(aes(tool, perc_compound_val, fill = count_group)) +
  geom_bar(stat = 'identity') +
  mytheme_discrete_x +
  facet_grid(~count_group, scales = 'free_x', space = 'free') +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c('#00B9F2', '#E69F00' , '#999999')) +
  xlab('') +
  ylab('compound precision metric') +
  theme(legend.position = "")

#ggsave('separate_figures/sup_figure_20.pdf', width = 21, height = 10, units = "cm")


#' # Sup Figure 21: theoratical nr of TP circRNAs
#' see script panel 3
#' 

#' # Sup Figure 22: sensitivity per tool

sens %>%
  mutate(margin = qnorm(0.975)*sqrt(sensitivity*(1-sensitivity)/957), #957 is total nr of val circ
         CI_low = sensitivity - margin,
         CI_up = sensitivity + margin) %>%
  ggplot(aes(tool, sensitivity, fill = count_group_median)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = scales::percent_format()) +
  mytheme_discrete_x +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), 
                width=.3, color = 'grey45') +
  facet_wrap(~count_group_median) +
  scale_fill_manual(values = c('#00B9F2', '#E69F00')) +
  xlab('') +
  theme(legend.position = '')

#ggsave('separate_figures/sup_figure_22.pdf', width =20, height = 10, units = "cm")

#' # Sup Figure 23: validation of circRNA compared to their precense in databases
val_db = cq %>% 
  select(circ_id, qPCR_val, RR_val, amp_seq_val, n_db, cell_line) %>% unique() %>%
  pivot_longer(cols = c(qPCR_val, RR_val, amp_seq_val), names_to = "val_type", values_to = "val") %>%
  mutate(n_db = ifelse(is.na(n_db), 0, n_db))

val_db

val_db$val_type = factor(val_db$val_type, levels = c('qPCR_val', 'RR_val', 'amp_seq_val'))

val_db %>%
  ggplot(aes(n_db, fill = val)) +
  geom_bar() +
  facet_grid(~val_type) +
  mytheme +
  scale_fill_manual(values = c('#CC79A7', '#00A875' ,'#999999')) +
  scale_x_continuous(breaks=c(0,5,10,15)) +
  ylab('number of circRNAs') +
  xlab('number of databases in which the circRNAs is reported')

#ggsave('sup_figure_23.pdf', width = 21, height = 8.5, units = "cm")



#' # Sup Figure 24: combo of two tools
#' see panel 4 script
#' 

#' # Sup Figure 25: combo of three tools

#' all possible combo's, only circRNAs ≥ 5

simple_union = read_tsv('../data/Supplementary_Table_6_combo_2tools.txt')
simple_union_3 = read_tsv('../data/Supplementary_Table_7_combo_3tools.txt')


#' ### simple version

#' add total nr of circ for that cell line
total_cl = all_circ %>%
  group_by(cell_line) %>%
  filter(count_group == 'count ≥ 5') %>%
  select(circ_id, cell_line) %>%
  unique() %>%
  count(cell_line) %>%
  rename(total_cell_line = n) %>% ungroup()

total_cl

#' make subset
union_sub_3 = simple_union_3 %>%
  #filter(nr_union > 4000) %>%
  filter(nr_union > 9000) %>% # is 70% of HLF : 0.7*13,087
  group_by(cell_line) %>%
  #filter(nr_union - pmax(total_n_1, total_n_2, total_n_3) > 2499) %>%
  ungroup() %>%
  filter(perc_compound_val_1 > 0.9,
         perc_compound_val_2  > 0.9,
         perc_compound_val_3 > 0.9)

union_sub_3


#' add individual tools of interest
# union_sub_3 = union_sub_3 %>%
#   bind_rows(simple_union_3 %>%
#               filter(tool_1 == tool_2 & tool_1 == tool_3,
#                      tool_1 %in% (union_sub_3 %>% pull(tool_1, tool_2) %>% unique())))

#' remove doubles
union_sub_3 = union_sub_3 %>%
  group_by(tool_combo) %>%
  sample_n(1) %>%
  ungroup() 

union_sub_3

#' give combo nr and add info

union_sub_3 = union_sub_3 %>%
  mutate(combo = paste("combo", 1:12, sep = "_"))

union_sub_3 = union_sub_3 %>%
  bind_rows(simple_union_3 %>%
              filter(nchar(tool_combo) == 1) %>% 
              unique() %>%
              right_join(union_sub_3 %>%
                           select(tool_lt_1, tool_lt_2, tool_lt_3, combo) %>%
                           pivot_longer(cols = -combo, names_to = 'tmp', values_to = 'tool_lt_1') %>%
                           select(-tmp)))


#' plot
union_sub_3 %>%
  filter(cell_line == 'HLF') %>%
  left_join(all_circ %>% select(circ_id, cell_line, count_group) %>%
              filter(count_group == "count ≥ 5") %>%
              unique() %>% count(cell_line) %>% rename(total_n = n)) %>%
  mutate(perc_union = nr_union/total_n) %>%
  ggplot(aes(w_val_rate, perc_union)) +
  geom_point(aes(color = (tool_1 == tool_2 & tool_2 == tool_3))) +
  geom_text_repel(aes(label=tool_combo, color = (tool_1 == tool_2 & tool_2 == tool_3)),
                  max.overlaps = 20) +
  facet_wrap(~combo, scales = 'free') +
  mytheme +
  theme(legend.position = 'NA') +
  scale_color_manual(values = c('#CC79A7', '#0072B2')) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('validation rate') +
  ylab('number of circRNAs')

#' ### also adding combo's?

union_sub_3 = union_sub_3 %>%
  bind_rows(simple_union %>%
              right_join(union_sub_3 %>% select(tool_combo, combo) %>%
                           filter(nchar(tool_combo) == 3) %>%
                           mutate(tool_combo1 = paste(substr(tool_combo, 1, 1), substr(tool_combo, 2, 2), sep = ''),
                                  tool_combo2 = paste(substr(tool_combo, 1, 1), substr(tool_combo, 3, 3), sep = ''),
                                  tool_combo3 = paste(substr(tool_combo, 2, 2), substr(tool_combo, 3, 3), sep = '')) %>%
                           select(combo, tool_combo1, tool_combo2, tool_combo3) %>%
                           pivot_longer(cols = c(tool_combo1, tool_combo2, tool_combo3), names_to = "tmp", values_to = 'tool_combo') %>%
                           select(-tmp) %>%
                           unique()))

union_sub_3 = union_sub_3 %>%
  group_by(combo, tool_combo, cell_line) %>%
  sample_n(1)


union_sub_3 %>%
  filter(cell_line == 'HLF') %>%
  left_join(all_circ %>% select(circ_id, cell_line, count_group) %>%
              filter(count_group == "count ≥ 5") %>%
              unique() %>% count(cell_line) %>% rename(total_n = n)) %>%
  mutate(perc_union = nr_union/total_n) %>%
  ggplot(aes(w_val_rate, perc_union)) +
  geom_point(aes(color = as.character(nchar(tool_combo)))) +
  geom_text_repel(aes(label=tool_combo, color = as.character(nchar(tool_combo))),
                  max.overlaps = 20) +
  facet_wrap(~combo, scales = 'free') +
  mytheme +
  theme(legend.position = 'NA') +
  scale_color_manual(values = c('#CC79A7', '#E69F00', '#0072B2')) + 
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('(weighted) compound precision value') +
  ylab('percentage of all predicted circRNAs')

#ggsave('separate_figures/sup_figure_25.pdf', width = 21, height = 18, units = "cm")


#' check mean increase in perc
simple_union_3 %>%
  filter(!(tool_1 == tool_2 & tool_1 == tool_3),
         !tool_1 == tool_2,
         !tool_2 == tool_3,
         !tool_1 == tool_3,
         perc_compound_val_1 >= 0.9,
         perc_compound_val_2 >= 0.9,
         perc_compound_val_3 >= 0.9) %>%
  mutate(increase = nr_union - total_n_1,
         increase_prec = increase/total_n_1) %>%
  pull(increase_prec) %>%
  quantile()
