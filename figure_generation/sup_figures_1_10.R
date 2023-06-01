#' ---
#' title: "Generation of sup figures and sup tables"
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
val = read_tsv('../data/Supplementary_Table_6A_precision_values.txt')

all_circ$tool = factor(all_circ$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))

cq$tool = factor(cq$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))


val$tool = factor(val$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl"))


cq
all_circ


#' # Sup Figure 1: number of circRNAs detected per tool
#' see panel 2 scripts
#'


#'# Sup Figure 2: circRNA BSJ count distribution per tool

all_circ %>% 
  filter(!tool == "Sailfish-cir") %>%
  ggplot(aes(tool, BSJ_count)) +
  geom_boxplot() +
  scale_y_log10(labels = scales::comma_format()) +
  ylab("BSJ count") +
  xlab('') +
  mytheme_discrete_x

#ggsave('separate_figures/sup_figure_2.pdf',  width = 20, height = 12, units = "cm")

#' numbers

all_circ %>% 
  count(tool)


#' # Sup Figure 3: correlation between BSJ counts from one circRNA detected by different tools

count_conc = tibble()
for (tool_1 in val %>% pull(tool) %>% unique()){
  for (tool_2 in val %>% pull(tool) %>% unique()){
    count_conc = count_conc %>%
      bind_rows(all_circ %>% filter(tool == tool_1) %>%
                  select(cell_line, circ_id, BSJ_count, tool) %>% unique() %>%
                  rename(BSJ_count_1 = BSJ_count, tool_1 = tool) %>%
                  inner_join(all_circ %>% filter(tool == tool_2) %>%
                               select(cell_line, circ_id, BSJ_count, tool) %>% unique() %>%
                               rename(BSJ_count_2 = BSJ_count, tool_2 = tool)))
  }
  
}

count_conc


# only keep unique combo's (every combo is calculated twice) => remove this step!
count_conc = count_conc %>%
  mutate(tool_comb = paste(tool_1, tool_2, sep = "/")) #,
  #       tool_comb = map_chr(str_split(tool_comb, "/"), 
  #                           ~str_c(str_sort(unique(.x)), collapse = "/"))) #%>%
  # group_by(tool_comb, circ_id, cell_line) %>% 
  # arrange(tool_1) %>%
  # filter(row_number()==1) %>%
  # ungroup()

count_conc %>%
  filter(!tool_1 == tool_2) %>%
  filter(!tool_1 == "Sailfish-cir", !tool_2 == 'Sailfish-cir') %>%
  #filter(tool_1 == "Sailfish-cir" | tool_2 == 'Sailfish-cir') %>%
  ggplot(aes(BSJ_count_1, BSJ_count_2)) +
  geom_smooth(method = "lm", color = '#00B9F2') +
  stat_regline_equation(label.x = -Inf, label.y = Inf, vjust = 1.5, hjust = -0.1, size = 3) +
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
  #          label.y= Inf, label.x = -Inf, vjust = 2, hjust = -0.1, size = 3) +
  stat_cor(aes(label = ..rr.label..),
           label.y= Inf, label.x = -Inf, vjust = 3, hjust = -0.1, size = 3) +
  stat_cor(aes(label = ..p.label..),
           label.y= Inf, label.x = -Inf, vjust = 5, hjust = -0.1, size = 3) +
  mytheme_discrete_x +
  theme(aspect.ratio = 1) +
  #coord_fixed(ratio = 1) +
  xlim(0,4000) + ylim(0,4000) +
  #xlim(0,130000) + ylim(0,130000) +
  facet_wrap(~tool_comb, nrow = 4) +
  xlab('') +
  ylab('')

#' ## slope and R2 plot
#' get same data in table
count_conc_lm = count_conc %>%
  filter(!tool_1 == tool_2) %>%
  filter(!tool_1 == "Sailfish-cir", !tool_2 == 'Sailfish-cir') %>% 
  group_by(tool_comb) %>% 
  do({
    mod = lm(log10(BSJ_count_2) ~ log10(BSJ_count_1), data = .)
    #mod = lm(BSJ_count_2 ~ BSJ_count_1, data = .)
    data.frame(intercept = coef(mod)[1],
               slope = coef(mod)[2],
               R_squared = summary(mod)$r.squared)
  })

count_conc_lm

count_conc_lm %>%
  #mutate(delta_slope = abs(1-slope)) %>%
  ggplot(aes(R_squared, slope)) +
  geom_point() +
  mytheme +
  geom_hline(yintercept = 1, color = '#999999', linetype = 'dashed') +
  geom_vline(xintercept = 1, color = '#999999', linetype = 'dashed') +
  coord_fixed(ratio = 1)

#ggsave('separate_figures/sup_figure_3.pdf', width = 10, height = 15, units = "cm")

# get median value
count_conc_lm %>% pull(slope) %>% quantile()
count_conc_lm %>% pull(R_squared) %>% quantile()
count_conc_lm %>% mutate(delta_slope = abs(1-slope)) %>%
  pull(delta_slope) %>% quantile()
  

#' # Sup Figure 4: number of tools by which the circRNA is detected
#' see panel 2 script
#' 

#' # Sup Figure 5: Jaccard distance between tools (heatmap)
#' Heatmap based on counts
#' here only for count ≥ 5 (in sup figure for all circ)
#' calculate Jaccard index 

jac = read_tsv("../data/Supplementary_Table_7_combo_2tools.txt") %>%
  filter(cell_line == 'SW480') %>%
  mutate(jac_index = nr_intersection/nr_union,
         jac_dist = 1 - jac_index) # calculate jac index and distance

jac

count_table = jac %>%
  select(tool_1, tool_2, jac_dist) %>%
  pivot_wider(values_from = jac_dist, names_from = tool_2)

count_table
rn = count_table$tool_1
count_table = count_table %>%
  select(-tool_1)
count_matrix = as.matrix(count_table)
rownames(count_matrix) = rn
 

#pdf("separate_figures/sup_figure_5_SW480.pdf", height=8, width=12)
heatmap.2(count_matrix,
          main = "SW480 - bin - Jaccard",
          trace = "none", density.info = "none")
#dev.off()


#' # Sup Figure 6: number of databases that also report that circRNA
all_circ_db = all_circ %>%
  mutate(n_db_group = ifelse(n_db == 1, '1 database', '0 databases'),
         n_db_group = ifelse(n_db > 1, '2-5 databases', n_db_group),
         n_db_group = ifelse(n_db > 5, '6-9 databases', n_db_group),
         n_db_group = ifelse(n_db > 9, '≥ 10 databases', n_db_group),
         n_db_group = ifelse(n_db == 16, 'all databases', n_db_group),
         n_db_group = ifelse(is.na(n_db), '0 databases', n_db_group))

all_circ_db$n_db_group = factor(all_circ_db$n_db_group, levels = c('all databases', '≥ 10 databases', '6-9 databases', "2-5 databases", '1 database', "0 databases"))
all_circ_db$tool = factor(all_circ_db$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 

all_circ_db %>% 
  ggplot(aes(tool, fill = n_db_group)) +
  geom_bar() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~cell_line) +
  ylab('number of circRNAs') +
  xlab('') +
  scale_fill_manual(name = "circRNAs reported by",
                    values = c("#0072B2", '#00B9F2', '#00A875', "#E69F00", "#CC79A7")) +
  mytheme_discrete_x

#ggsave('sup_figures/sup_figure_6.pdf',  width = 20, height = 9, units = "cm")


#' numbers for paper
all_circ_db %>% 
  select(circ_id, n_db) %>% unique() %>%
  filter(is.na(n_db))

round(215885/315312, digits = 3)

all_circ_db %>%
  group_by(tool) %>%
  summarise(perc_db = sum(is.na(n_db)) / n()) %>% arrange(perc_db)

all_circ_db %>%
  group_by(tool) %>%
  summarise(perc_db = sum(is.na(n_db)) / n()) %>% arrange(perc_db) %>%
  pull(perc_db) %>% quantile()

#' # Sup Figure 7: number of circRNAs with canonical linear annotation match
all_circ_enst = all_circ  %>%
  mutate(ENST_group = ifelse(is.na(ENST), NA, 'one match'),
         ENST_group = ifelse(ENST == 'ambiguous', 'ambiguous', ENST_group)) 

all_circ_enst$ENST_group = factor(all_circ_enst$ENST_group, levels = c('one match', 'ambiguous', NA))

all_circ_enst %>%
  ggplot(aes(tool, fill = ENST_group)) +
  geom_bar(position = 'fill') +
  mytheme_discrete_x +
  theme(legend.position = "right") +
  ylab('% of circRNAs') + 
  xlab('') +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(name = "",
                    values = c('#F0E442', '#CC79A7', '#99999'))

#ggsave('sup_figures/sup_figure_7.pdf',  width = 18, height = 13, units = "cm")

#' calculations
all_circ_enst %>%
  select(chr, start, end, strand, ENST_group) %>% unique() %>%
  count(ENST_group)

round(100 * 216483 / (216483 + 41683 + 144782), 1)
round(100 * 41683 / (216483 + 41683 + 144782), 1)
round(100 * 144782 / (216483 + 41683 + 144782), 1)

#' # Sup Figure 8: distribution of estimated circRNA length
all_circ %>% 
  mutate(estim_len_ex = end - start,
         len_type = 'including all exons and introns') %>%
  bind_rows(all_circ %>% 
              filter(!estim_len_ex == 'ambiguous',
                     !is.na(estim_len_ex)) %>%
              mutate(estim_len_ex = as.numeric(estim_len_ex),
                     len_type = 'including all exons')) %>%
  ggplot(aes(tool, estim_len_ex, color = len_type)) +
  geom_boxplot(outlier.shape=NA, lwd = 1/.pt) +
  mytheme_discrete_x +
  scale_color_manual(values = c( '#CC79A7', '#E69F00')) +
  theme(legend.title = element_blank(), legend.position = 'right') + 
  scale_y_log10(labels = scales::comma_format(), expand = expansion(mult = 0.5)) +
  coord_cartesian(ylim = c(50, 110000)) +
  ylab("estimated length of circRNAs")  +
  theme(axis.title.x=element_blank())


#ggsave('sup_figures/sup_figure_8.pdf', width = 20, height = 10, units = "cm")

#' numbers

all_circ %>% 
  mutate(estim_len_ex = end - start,
         len_type = 'including all exons and introns') %>%
  bind_rows(all_circ %>% 
              filter(!estim_len_ex == 'ambiguous',
                     !is.na(estim_len_ex)) %>%
              mutate(estim_len_ex = as.numeric(estim_len_ex),
                     len_type = 'including all exons')) %>% 
  #filter(len_type == 'including all exons and introns') %>%
  filter(len_type == 'including all exons') %>%
  count(tool)

#' # Sup Figure 9: number of exons per circRNA per tool
all_circ_exons = all_circ %>%
  filter(!nr_exons == 'ambiguous') %>%
  mutate(nr_exons = as.numeric(nr_exons),
         exon_group = ifelse(nr_exons == 1, 'single exon', NA),
         exon_group = ifelse(nr_exons > 1, '2-5 exons', exon_group),
         exon_group = ifelse(nr_exons > 5, '6-9 exons', exon_group),
         exon_group = ifelse(nr_exons > 9, '≥ 10 exons', exon_group),
         nr_exons = as.character(nr_exons)) %>%
  bind_rows(all_circ %>%
              filter(nr_exons == 'ambiguous') %>%
              mutate(exon_group = ifelse(nr_exons == "ambiguous", 'ambiguous', exon_group))) %>%
  bind_rows(all_circ %>%
              filter(is.na(nr_exons)) %>%
              mutate(exon_group = NA))

all_circ_exons %>% filter(exon_group == "ambiguous")

all_circ_exons$exon_group = factor(all_circ_exons$exon_group, 
                                   levels = c('≥ 10 exons', '6-9 exons', '2-5 exons', 'single exon', 'ambiguous', NA))
all_circ_exons$tool = factor(all_circ_exons$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl")) 


all_circ_exons %>%
  #filter(!exon_group == 'ambiguous', !is.na(exon_group)) %>%
  ggplot(aes(tool, fill = exon_group)) +
  geom_bar(position = 'fill') +
  mytheme_discrete_x +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(name = "circRNAs with ...",
                    values = c("#0072B2", '#00B9F2', '#00A875', "#E69F00", '#CC79A7', '#99999')) +
  theme(legend.position = NULL) +
  xlab('') +
  ylab('')

#ggsave('sup_figures/sup_figure_9.pdf',  width = 18, height = 13, units = "cm")


#' # Sup Figure 10: nr of circRNAs on both strands per tool
#' first calculate strand distribution linear genes
ens = read_tsv('~/Documents/PhD/indexes/Homo_sapiens.GRCh38.103.genes.gtf') %>%
  select(gene_id, strand) %>%
  rename(circ_id = gene_id) %>%
  mutate(tool = 'lineair genes',
         tool_group = 'lin')

ens

all_circ_strand = all_circ %>%
  mutate(tool_group = 'circ') %>%
  bind_rows(ens) 

all_circ_strand$tool = factor(all_circ_strand$tool, levels = c("circseq_cup", "CIRI2", "CIRIquant", "CircSplice", "find_circ", "CirComPara2",  "CIRCexplorer3", "circtools", "Sailfish-cir", "NCLscan", "NCLcomparator", "PFv2", "ecircscreen", "KNIFE",  "circRNA_finder", "segemehl", 'lineair genes')) 

all_circ_strand %>%
  ggplot(aes(tool, fill = strand)) +
  geom_bar(position = 'fill') +
  ylab('% of circRNAs') +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab('') +
  facet_grid(~tool_group, scales = 'free_x', space = 'free') +
  scale_fill_manual(values = c('#009E73', '#CC79A7' , '#5AB4E5')) +
  mytheme_discrete_x +
  theme(legend.position = 'right')

#ggsave('sup_figures/sup_figure_10.pdf',  width = 20, height = 12, units = "cm")