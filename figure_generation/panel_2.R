
#' ---
#' title: "Generation of figures pannel 2"
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
library(ggseqlogo)

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
all_circ$tool = factor(all_circ$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

all_circ

#' # Figure 2B & Sup Figure 1

all_circ %>% 
  group_by(tool, cell_line, count_group) %>% 
  summarise(n = n()) %>%
  group_by(tool, cell_line) %>% mutate(total_n = sum(n)) %>%
  filter(cell_line == "HLF") %>%
  ggplot(aes(reorder(tool, -total_n), n, fill = count_group)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #facet_wrap(~cell_line, nrow = 3) +
  mytheme_discrete_x +
  theme(axis.title.x=element_blank(), legend.title=element_blank()) +
  scale_fill_manual(values = c('#00B9F2', '#E69F00' , '#999999')) +
  ylab("number of circRNAs")
  

#ggsave('separate_figures/figure_2B.pdf',  width = 10, height = 8.5, units = "cm")
#ggsave('../supplemental/sup_figures/sup_figure_1.pdf',  width = 20, height = 12, units = "cm")

# all_circ %>% 
#   group_by(tool, cell_line, count_group) %>% 
#   summarise(n = n()) %>%
#   group_by(tool, cell_line) %>% mutate(total_n = sum(n)) %>%
#   write_tsv('source_data_fig_2B.txt')


#' # Figure 2C & Sup Figure 4

n_detected_per_tool = all_circ %>% 
  mutate(n_detected_group = NA,
         n_detected_group = ifelse(n_detected == 1, 'unique', n_detected_group),
         n_detected_group = ifelse(n_detected > 1, '2 tools', n_detected_group),
         n_detected_group = ifelse(n_detected > 2, '2-5 tools', n_detected_group),
         n_detected_group = ifelse(n_detected > 5, '6-9 tools', n_detected_group),
         n_detected_group = ifelse(n_detected > 9, '≥ 10 tools', n_detected_group),
         n_detected_group = ifelse(n_detected == 16, 'all tools', n_detected_group))

n_detected_per_tool$n_detected_group = factor(n_detected_per_tool$n_detected_group, 
                                              levels = c('all tools', '≥ 10 tools', '6-9 tools', "2-5 tools", '2 tools', 'unique'))

n_detected_per_tool$tool = factor(n_detected_per_tool$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 


n_detected_per_tool %>%
  #filter(cell_line == "HLF") %>%
  ggplot(aes(tool, fill = n_detected_group)) +
  geom_bar() +
  mytheme_discrete_x +
  ylab("number of circRNAs") +
  scale_fill_manual(name = "circRNAs detected by",
                    values = c('#D55E00',"#0072B2", '#00B9F2', '#00A875', "#E69F00", "#CC79A7")) +
  theme(axis.title.x=element_blank()) +
  #facet_wrap(~cell_line) +
  scale_y_continuous(labels = scales::comma_format())

#ggsave('separate_figures/figure_2C.pdf',  width = 10, height = 9, units = "cm")
#ggsave('../supplemental/sup_figures/sup_figure_4.pdf',  width = 21, height = 12, units = "cm")

# n_detected_per_tool %>% 
#   group_by(tool, cell_line, n_detected_group) %>%
#   tally() %>% ungroup() %>%
#   write_tsv('source_data_fig_2C.txt')

#' # Figure 2D

circ_ss_list = list()

for (tool_name in all_circ %>% pull(tool) %>% unique()) {
  ss = all_circ %>% filter(tool == tool_name) %>% pull(ss_motif) 
  if (length(ss) > 0 ){
    circ_ss_list[tool_name] = all_circ %>% filter(tool == tool_name) %>% pull(ss_motif) %>% list()
    
  }
}


# Create custom colour scheme
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                      cols=c('#00A875', '#F15A22', '#0072BC', '#F7941D'))

ggplot() + 
  geom_logo(circ_ss_list, method = 'prob', col_scheme = cs1) + 
  theme_logo() +
  guides(scale = "none") +
  facet_wrap(~seq_group, nrow = 2) +
  theme(axis.text.x=element_blank(),
        text = element_text(size=10, colour='gray30'),
        title = element_text(size=10, colour='gray30'))

#ggsave('separate_figures/figure_2D.pdf',  width = 20, height = 4.5, units = "cm")

# circ_ss_list %>% write_tsv('source_data_fig_2D.txt')


