#' ---
#' title: "Generation of figures pannel 4"
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

#' ## loads standard libraries and resolves conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')

#' ## load specific libraries
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
all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')
cq = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')
val = read_tsv('../data/Supplementary_Table_4_precision_values.txt')

val$tool = factor(val$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

all_circ
cq
val

#' # Figure 4A
#' first how many per group in selected circ?
nr_detected = all_circ %>% 
  group_by(chr, start, end, circ_id, cell_line) %>%
  summarise(n_detected = n()) %>%
  ungroup()

val_df = cq %>% left_join(nr_detected) %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, n_detected, cell_line) %>% unique() %>%
  pivot_longer(cols = c(qPCR_val, RR_val, amp_seq_val), names_to = "val_type", values_to = "val")

val_df

val_df$val_type = factor(val_df$val_type, levels = c('qPCR_val', 'RR_val', 'amp_seq_val'))

val_df %>%
  ggplot(aes(n_detected, fill = val)) +
  geom_bar() +
  facet_grid(~val_type) +
  mytheme +
  scale_fill_manual(values = c('#CC79A7', '#00A875' ,'#999999')) +
  ylab('number of circRNAs') +
  xlab('number of tools the circRNAs was detected by')


#ggsave('separate_figures/figure_4A.pdf',  width = 21, height = 8.5, units = "cm")



#' # Figure 4B & Sup Figure 24
simple_union = read_tsv('../data/Supplementary_Table_6_combo_2tools.txt')
simple_union

#' add total nr of circ for that cell line
total_cl = all_circ %>%
  group_by(cell_line) %>%
  filter(count_group == 'count ≥ 5') %>%
  select(circ_id, cell_line) %>%
  unique() %>%
  count(cell_line) %>%
  rename(total_cell_line = n) %>% ungroup()

total_cl

union_sub = simple_union %>%
  # filter based on percentage increase
  left_join(total_cl) %>%
  filter((nr_union - pmax(total_n_1, total_n_2))/total_cell_line > 0.075) %>%
  #filter(nr_union - pmax(total_n_1, total_n_2) > 999) %>%
  filter(perc_compound_val_1 >= 0.9,
         perc_compound_val_2 >= 0.9)

#' add individual tools of interest
union_sub = union_sub %>%
  bind_rows(simple_union %>%
              filter(tool_1 == tool_2,
                     tool_1 %in% (union_sub %>% pull(tool_1, tool_2) %>% unique())) %>%
              left_join(total_cl) )

#' remove doubles
union_sub = union_sub %>%
  group_by(tool_combo, cell_line) %>%
  sample_n(1) %>%
  ungroup()

union_sub

#' as percentage of all circ in that cell line

union_sub %>% 
  left_join(all_circ %>% select(circ_id, cell_line, count_group) %>%
              filter(count_group == "count ≥ 5") %>%
              unique() %>% count(cell_line) %>% rename(total_n = n)) %>%
  mutate(perc_union = nr_union/total_n) %>%
  #filter(cell_line == "HLF") %>%
  ggplot(aes(w_val_rate, perc_union)) +
  geom_point(aes(color = (tool_1 == tool_2))) +
  geom_text_repel(aes(label=str_remove(tool_combo, '-'), color = (tool_1 == tool_2)), max.overlaps = 20, size = 3) +
  mytheme +
  theme(legend.position = 'NA') +
  facet_wrap(~cell_line, nrow = 2, scales = 'free') + 
  scale_color_manual(values = c('#E69F00', '#0072B2')) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('(weighted) compound precision value') +
  ylab('% of all predicted circRNAs')

#ggsave('separate_figures/figure_4B.pdf',  width = 12, height = 10, units = "cm")
#ggsave('separate_figures/sup_figure_24.pdf',  width = 20, height = 20, units = "cm")

  
#' check mean increase in perc
simple_union %>%
  filter(count_group == "count ≥ 5",
         !tool_1 == tool_2,
         perc_compound_val_1 >= 0.9,
         perc_compound_val_2 >= 0.9) %>%
  mutate(increase = nr_union - total_n_1,
         increase_prec = increase/total_n_1) %>% #view()
  pull(increase_prec) %>%
  quantile()
