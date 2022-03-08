## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

library(graph4lg)
library(igraph)


## ---- echo = FALSE, eval = TRUE-----------------------------------------------
data("data_tuto")

mat_dps <- data_tuto[[1]]
mat_pg <- data_tuto[[2]]
graph_ci <- data_tuto[[3]]
dmc <- data_tuto[[4]]
land_graph <- data_tuto[[5]]
mat_ld <- data_tuto[[6]]

## -----------------------------------------------------------------------------
land_graph <- gen_graph_topo(mat_w = mat_ld,
                             mat_topo = mat_ld,
                             topo = "comp")

# Plot the histogram of its link weights
plot_w_hist(graph = land_graph)

## -----------------------------------------------------------------------------
miw_lg <- compute_node_metric(graph = land_graph, metrics = "miw")
head(miw_lg)

## -----------------------------------------------------------------------------
land_graph <- add_nodes_attr(graph = land_graph,
                             input = "df",
                             data = miw_lg,
                             index = "ID")

## ---- eval = FALSE------------------------------------------------------------
#  mat_dps <- mat_gen_dist(x = data_simul_genind, dist = "DPS")

## -----------------------------------------------------------------------------
gen_comp_graph <- gen_graph_topo(mat_w = mat_dps,
                                 mat_topo = mat_dps,
                                 topo = "comp")

## -----------------------------------------------------------------------------
plot_w_hist(graph = gen_comp_graph, 
            fill = "darkblue")

## -----------------------------------------------------------------------------
miw_comp <- compute_node_metric(graph = gen_comp_graph, metrics = "miw")
gen_comp_graph <- add_nodes_attr(graph = gen_comp_graph,
                                 input = "df",
                                 data = miw_comp,
                                 index = "ID")


## -----------------------------------------------------------------------------
graph_node_compar(x = land_graph, y = gen_comp_graph,
                  metrics = c("miw", "miw"), method = "spearman",
                  weight = TRUE, test = TRUE)

## -----------------------------------------------------------------------------

mat_geo <- mat_geo_dist(data = pts_pop_simul,
                        ID = "ID", x = "x", y = "y")
mat_geo <- reorder_mat(mat_geo, order = row.names(mat_dps))

gen_gab_graph <- gen_graph_topo(mat_w = mat_dps,
                                mat_topo = mat_geo,
                                topo = "gabriel")

# Associate the values of miw from the complete graph to this graph

gen_gab_graph <- add_nodes_attr(gen_gab_graph,
                                data = miw_comp,
                                index = "ID")

# Plot the graph with node sizes proportional to MIW

plot_graph_lg(graph = gen_gab_graph, 
              crds = pts_pop_simul,
              mode = "spatial",
              node_size = "miw",
              link_width = "inv_w")


## -----------------------------------------------------------------------------

land_graph_thr <- gen_graph_thr(mat_w = mat_ld, mat_thr = mat_ld,
                             thr = 2000, mode = "larger")

plot_graph_lg(land_graph_thr, 
              mode = "spatial", 
              crds = pts_pop_simul,
              link_width = "inv_w", 
              pts_col = "#80C342")


## -----------------------------------------------------------------------------

graph_topo_compar(obs_graph = land_graph_thr, 
                  pred_graph = gen_gab_graph,
                  mode = "mcc", 
                  directed = FALSE)

## -----------------------------------------------------------------------------
graph_plot_compar(x = land_graph_thr, y = gen_gab_graph,
                  crds = pts_pop_simul)


## -----------------------------------------------------------------------------
graph_modul_compar(x = land_graph_thr, 
                   y = gen_gab_graph)

## -----------------------------------------------------------------------------
module_land <- compute_graph_modul(graph = land_graph_thr, 
                                   algo = "fast_greedy",
                                   node_inter = "distance")

land_graph_thr <- add_nodes_attr(graph = land_graph_thr,
                                 data = module_land,
                                 index = "ID")

module_gen <- compute_graph_modul(graph = gen_gab_graph, 
                                   algo = "fast_greedy",
                                   node_inter = "distance")

gen_gab_graph <- add_nodes_attr(graph = gen_gab_graph,
                                 data = module_gen,
                                 index = "ID")

## -----------------------------------------------------------------------------
plot_graph_lg(graph = land_graph_thr,
              mode = "spatial",
              crds = pts_pop_simul,
              module = "module")

## -----------------------------------------------------------------------------
plot_graph_lg(graph = gen_gab_graph,
              mode = "spatial",
              crds = pts_pop_simul,
              module = "module")

