
library(visNetwork)
library(tidyverse)
library(rstatix)
library(ggsignif)
library(ggpubr)

network <- function(aa) {
node_list <- unique(c(unlist(aa$gene_CLSA), unlist(aa$metab_name_CLSA)))
nodes <- data.frame(id = unique(c(aa$gene_CLSA, aa$metab_name_CLSA )), 
                    shape = c(rep("square", length(unique(aa$gene_CLSA))), rep("circle", length(unique(aa$metab_name_CLSA ) ))))

edges <- data.frame(from = nodes$id[c(match(aa$gene_CLSA, nodes$id))],
                    to = nodes$id[c(match(aa$metab_name_CLSA, nodes$id))])

set.seed(42)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges,  directed = FALSE)
routes_tidy %>% 
  mutate(community=as.factor(group_louvain(weights = NULL))) -> routes_tidy

node_list <-
  routes_tidy %>%
  activate(nodes) %>%
  data.frame()

col_vector <- c(  "#e7b800", "black", "#FF2400",  "#2b6d65", "#5abce7" , "#e6009e", 
                  "#aa00c7", "#002fe7",  "#ff8f00", "green","#9ba39b" ,"#FFED6F", "darkred")

col_vector <- data.frame(color = col_vector, community = c(1:13))
additional <- data.frame(color = rep("lightgrey", length(c(14:25))), community = c(14:25))
col_vector <- rbind(col_vector, additional)

node_list <- merge(node_list, col_vector, by = "community")

visNetwork(node_list, edges) %>%
  visEdges(color = "black") %>%
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visOptions(selectedBy = list(variable = "community", multiple = TRUE)) %>%
  visLegend()

}