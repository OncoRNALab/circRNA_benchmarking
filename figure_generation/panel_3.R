#' ---
#' title: "Generation of figures panel 3"
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

#' ## Loads standard libraries and resolve conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')

#' ## Load specific libraries
library(UpSetR)

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
all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')
cq = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')
val = read_tsv('../data/Supplementary_Table_6A_precision_values.txt')

val$tool = factor(val$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

cq
val
all_circ

#' # Figure 3A

val_longer = val %>% pivot_longer(cols = c(perc_qPCR_val, perc_RR_val, perc_amp_val),
                     names_to = 'val_type', values_to = 'perc')  %>%
  left_join(val %>% select(tool, count_group, nr_qPCR_total, nr_RR_total, nr_amp_total) %>%
              rename(perc_qPCR_val = nr_qPCR_total, perc_RR_val = nr_RR_total, perc_amp_val = nr_amp_total) %>%
              pivot_longer(cols = c(perc_qPCR_val, perc_RR_val, perc_amp_val), values_to = 'nr_in_group', 
                           names_to = "val_type")) %>%
  select(tool, count_group, val_type, perc, nr_in_group) %>%
  mutate(margin = qnorm(0.975)*sqrt(perc*(1-perc)/nr_in_group),
         CI_low = perc - margin,
         CI_up = perc + margin) %>%
  mutate(count_group == ifelse(count_group == "count < 5", 'BSJ count < 5', count_group),
         count_group == ifelse(count_group == "count ≥ 5", 'BSJ count ≥ 5', count_group),
         count_group == ifelse(count_group == "no_counts", 'no BSJ count', count_group)) 

val_longer$tool = factor(val_longer$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 
val_longer$val_type = factor(val_longer$val_type, levels = c('perc_qPCR_val', 'perc_RR_val', 'perc_amp_val'))

val_longer %>%
  filter(count_group == 'count < 5', val_type == "perc_amp_val") %>%
  pull(nr_in_group) %>% quantile()

val_longer %>%
  ggplot(aes(tool, perc, fill = count_group)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), 
                width=.3, color = 'grey45') +
  facet_grid(scales = 'free_x', space = "free", 
             rows = vars(val_type), cols = vars(count_group)) +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0,0.25, 0.5, 0.75, 1)) +
  mytheme_discrete_x +
  #theme(strip.text.y = element_text(size = 10), legend.position = "none") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#00B9F2', '#E69F00' , '#999999')) +
  xlab('') +
  ylab('')

#ggsave('../tmp_figures//figure_3A.pdf',  width = 21, height = 13, units = "cm")
#val_longer %>% write_tsv('../tmp_figures/source_data_fig_3A.txt')

#' # Figure 3B
#' upset
upset = list()
upset['qPCR_val'] = cq %>% filter(qPCR_val == "pass") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['RR_val'] = cq %>% filter(RR_val == "pass") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['amp_seq_val'] = cq %>% filter(amp_seq_val == "pass") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['qPCR_fail'] = cq %>% filter(qPCR_val == "fail") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['RR_fail'] = cq %>% filter(RR_val == "fail") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['amp_seq_fail'] = cq %>% filter(amp_seq_val == "fail") %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['RR_out_of_range'] = cq %>% filter(is.na(RR_val)) %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()
upset['amp_seq_not_included'] = cq %>% filter(is.na(amp_seq_val)) %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% pull(id_cl) %>% unique() %>% list()

# get number y-axis to add to plot
cq %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% select(id_cl, qPCR_val) %>% unique() %>% count(qPCR_val)
cq %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% select(id_cl, RR_val) %>% unique() %>% count(RR_val)
cq %>% mutate(id_cl = paste(circ_id, cell_line, sep = "_")) %>% select(id_cl, amp_seq_val) %>% unique() %>% count(amp_seq_val)


# pdf(file="../tmp_figures/figure_3B_upset_small.pdf", onefile=FALSE,
#     width = 8, height = 2.5)
upset(fromList(upset), 
      order.by = "freq",
      sets = c('qPCR_val', 'RR_val', 'amp_seq_val', 'qPCR_fail', 'RR_fail', 'amp_seq_fail', 'RR_out_of_range', 'amp_seq_not_included'),
      mainbar.y.label = "number of circRNAs", sets.x.label = "number of circRNAs",
      keep.order = TRUE,
      number.angles = 30,
      point.size = 2.5, line.size = 1)

# cq %>% select(qPCR_val, RR_val, amp_seq_val) %>% 
#    write_tsv('../tmp_figures/source_data_fig_3B.txt')

# dev.off()

#' # Figure 3C & Sup Figure 25
#' add total nr of circRNAs and calculate theoretical nr of TP
val_cl = val %>%
  # use perc_compound_val
  select(tool, count_group, perc_compound_val) %>%
  # get the number of detected circRNAs for each cell line and tool, per count group
  left_join(all_circ %>%
              group_by(cell_line, tool, count_group) %>%
              summarise(total_n_ut = n())) %>%
  # calculate the theoretical nr of validated circ
  mutate(theo_TP_all = perc_compound_val * total_n_ut)

val_cl

val_cl$tool = factor(val_cl$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

val_cl %>%
  filter(cell_line == 'HLF') %>%
  ggplot() +
  geom_bar(aes(tool, total_n_ut), stat = "identity", fill = 'grey') +
  geom_bar(aes(tool, theo_TP_all, fill = count_group), stat = "identity") +
  #facet_wrap(~cell_line+count_group, scales = 'free') +
  facet_wrap(~count_group, scales = 'free') +
  #facet_grid(rows = vars(cell_line), cols = vars(count_group), scales = 'free', space = 'free') +
  mytheme_discrete_x +
  scale_fill_manual(values = c('#00B9F2', '#E69F00' , '#999999')) +
  xlab('') +
  ylab('') +
  theme(legend.position = NULL) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(legend.position="none")

#ggsave('separate_figures/sup_figure_25.pdf',  width = 21, height = 24.5, units = "cm")
#ggsave('separate_figures/figure_3C_alt2.pdf',  width = 21, height = 8, units = "cm")

#val_cl %>% write_tsv('../tmp_figures/source_data_fig_3C.txt')

