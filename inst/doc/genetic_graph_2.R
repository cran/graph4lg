## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

library(graph4lg)
library(igraph)


## ---- echo = TRUE, eval = TRUE------------------------------------------------
data("data_tuto")

mat_dps <- data_tuto[[1]]
mat_pg <- data_tuto[[2]]
graph_ci <- data_tuto[[3]]
dmc <- data_tuto[[4]]
land_graph <- data_tuto[[5]]
mat_ld <- data_tuto[[6]]

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  dmc <- dist_max_corr(mat_gd = mat_dps, mat_ld = mat_ld,
#                       interv = 500, pts_col = "black")

## -----------------------------------------------------------------------------
# DMC value
dmc[[1]]
# Correlation coefficients
dmc[[2]]
# Threshold distances tested
dmc[[3]]

## ---- eval = TRUE, echo = FALSE-----------------------------------------------

vec_t <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000,
           6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000, 10230.05)

cc_val <- c(NA, 0.2986565, 0.3154498, 0.5188747, 0.7059633, 0.7559539, 
            0.7850267, 0.7947691, 0.8038470, 0.7853646, 0.7760106, 
            0.7641339, 0.7530264, 0.7462445, 0.7386713, 0.7333936, 
            0.7305631, 0.7226695, 0.7137972, 0.7110962, 0.7041702)


dat <- data.frame(vec_t = vec_t, cc_val = cc_val)
if(any(is.na(dat$cc_val))){
  dat <- dat[-which(is.na(dat$cc_val)), ]
}

plot_dmc <- ggplot2::ggplot(data = dat, ggplot2::aes(x = vec_t, y = cc_val)) +
  ggplot2::geom_point(color = "#999999", size = 1, shape = 16) +
  ggplot2::geom_line(color = "black") +
  ggplot2::labs(x = "Distance threshold",
                y = "Correlation coefficient") +
  ggplot2::theme_bw()
plot_dmc


## -----------------------------------------------------------------------------
scatter_dist(mat_gd = mat_dps, mat_ld = mat_ld, 
             pts_col = "black")

## -----------------------------------------------------------------------------
# First compute the geographical distance between populations
mat_geo <- mat_geo_dist(data = pts_pop_simul,
                        ID = "ID", x = "x", y = "y",
                        crds_type = "proj")
# Reorder the matrix
mat_geo <- reorder_mat(mat_geo, order = row.names(mat_dps))

# Create the thresholded graph
graph_thr <- gen_graph_thr(mat_w = mat_dps, mat_thr = mat_geo,
                           thr = 12000, mode = "larger")
graph_thr

## -----------------------------------------------------------------------------
graph_gab_geo <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_geo,
                                topo = "gabriel")
graph_gab_geo
graph_gab_gen <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                                topo = "gabriel")

## -----------------------------------------------------------------------------
graph_mst <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                            topo = "mst")
graph_mst

## -----------------------------------------------------------------------------
graph_percol <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                               topo = "percol")

## -----------------------------------------------------------------------------
graph_k3 <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                               topo = "knn", k = 3)

## -----------------------------------------------------------------------------
graph_comp <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                             topo = "comp")

## ---- echo = FALSE------------------------------------------------------------
g_plan <- graph_plan(crds = pts_pop_simul,
                     ID = "ID", x = "x", y = "y",
                     weight = TRUE)
g_plan

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  graph_ci <- gen_graph_indep(x = data_genind,
#                              dist = "PCA",
#                              cov = "sq",
#                              adj = "holm")

## -----------------------------------------------------------------------------
graph_ci

## -----------------------------------------------------------------------------
df_metric <- compute_node_metric(graph = graph_percol)
head(df_metric)

## -----------------------------------------------------------------------------
graph_percol <- add_nodes_attr(graph = graph_percol,
                               data = df_metric,
                               index = "ID",
                               include = "all")
graph_percol

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  graph_percol <- add_nodes_attr(graph_percol,
#                               input = "shp",
#                               dir_path = system.file('extdata', package = 'graph4lg'),
#                               layer = "patches",
#                               index = "Id",
#                               include = "Area")

## -----------------------------------------------------------------------------
df_modul <- compute_graph_modul(graph = graph_percol, 
                    algo = "fast_greedy",
                    node_inter = "distance")

head(df_modul)
# Unique values of module ID
unique(df_modul$module)

## -----------------------------------------------------------------------------
graph_percol <- add_nodes_attr(graph = graph_percol, 
                               input = "df",
                               data = df_modul,
                               index = "ID")

## -----------------------------------------------------------------------------
p <- plot_graph_lg(graph = graph_mst, 
                   mode = "spatial",
                   crds = pts_pop_simul,
                   link_width = "inv_w")
p

## ---- eval = TRUE, echo = TRUE------------------------------------------------
# Compute the metrics
df_metric_mst <- compute_node_metric(graph = graph_mst)

# Associate them to the graph
graph_mst <- add_nodes_attr(graph = graph_mst,
                               data = df_metric_mst,
                               index = "ID",
                               include = "all")
# Compute the modules
df_module_mst <- compute_graph_modul(graph = graph_mst, 
                    algo = "fast_greedy",
                    node_inter = "distance")
# Associate them to the graph
graph_mst <- add_nodes_attr(graph = graph_mst,
                               data = df_module_mst,
                               index = "ID",
                               include = "all")


# Plot the graph
# Link width is inversely proportional to genetic distance
# Node size is proportional to MIW metric
# Node color depends on the node module

plot_graph_lg(graph = graph_mst,
              mode = "spatial",
              crds = pts_pop_simul,
              link_width = "inv_w",
              node_size = "miw",
              module = "module")


## -----------------------------------------------------------------------------
p <- plot_graph_lg(graph = graph_mst, 
                   mode = "aspatial", 
                   node_inter = "distance", 
                   link_width = "inv_w",
                   node_size = "miw",
                   module = "module")
p

## -----------------------------------------------------------------------------
scatter_dist_g(mat_y = mat_dps , 
               mat_x = mat_ld, 
               graph = graph_gab_geo)

## -----------------------------------------------------------------------------
p <- plot_w_hist(graph = graph_gab_gen)
p

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  graph_to_shp(graph = graph_mst,
#               crds = pts_pop_simul,
#               mode = "both",
#               layer = "test_shp_mst",
#               dir_path = "wd",
#               metrics = TRUE,
#               crds_crs = 2154)

