
#' ---
#' title: "Generation of figures panel 1"
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
library(europepmc)

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

#' # Figure 1A
#' source: https://cran.r-project.org/web/packages/europepmc/vignettes/evergreenreviewgraphs.html
prop_circ = europepmc::epmc_hits_trend(query = "circRNA", period = 2005:2022)

prop_circ

ggplot(prop_circ, aes(year, query_hits / all_hits)) + 
  geom_point() + 
  geom_line() +
  xlab("Year published") + 
  ylab("% of circRNA full-texts in Europe PubMed Central") +
  mytheme + 
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels = scales::percent_format())


# ggsave('separate_figures/figure_1A.pdf',  width = 10, height = 8.5, units = "cm")

#' calculate CAGR
#' CAGR = 100[(Ending Value + Beginning Value)1/n  1 ]
prop_circ %>% filter(query_hits > 0) %>%
  mutate(previous_year = year - 1) %>%
  left_join(prop_circ %>% rename(previous_year = year, query_hits_py = query_hits) %>%
              select(previous_year, query_hits_py)) %>%
  mutate(AGR = (query_hits - query_hits_py) / query_hits_py,
         nr_y = year - 2011,
         CAGR = 100 * ((query_hits + query_hits_py)^(1/ nr_y) -1)) %>%
  filter(!year == 2022) %>%
  top_n(5, nr_y) %>%
  pull(CAGR) %>% median()


100 * ((4222 + 215) ^ (1/5) - 1)

((4177-215)^(1/5))
