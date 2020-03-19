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


## -----------------------------------------------------------------------------
data_genind <- genepop_to_genind(path = paste0(system.file('extdata', 
                                                           package = 'graph4lg'), "/gpop_51_sim22_01_25.txt"),
                                 n.loci = 20, pop_names = as.character(1:50))
data_genind

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  genind_to_genepop(x = data_genind, output = "data_gpop_test.txt")

## -----------------------------------------------------------------------------
loci_names <- c("DkiD104", "DkiD124", "DkiD102", "CAM19",
                "DkiC118", "DkiD128", "DkiB12",  "Lswmu7",
                "DkiD109", "Lswmu5",  "TG12_15", "DkiD12" )
data(data_pc_genind)
ind_names <- row.names(data_pc_genind@tab)
pop_names <- c("BT-1", "BT-10", "BT-11", "BT-12", "BT-13", "BT-2",
               "BT-3",  "BT-4",  "BT-5",  "BT-6",  "BT-7",  "BT-8",
               "BT-9",  "GT-1",  "GT-2", "GT-3",  "GT-4",  "GT-5",
               "GT-6", "GT-7")
data_paru <- structure_to_genind(path = paste0(system.file('extdata', 
                                                           package = 'graph4lg'), 
                                               "/data_PC_str.txt"),
                                 loci_names = loci_names,
                                 pop_names = pop_names,
                                 ind_names = ind_names)
data_paru

## -----------------------------------------------------------------------------
head(data_pc_gstud)

## -----------------------------------------------------------------------------
gstud_to_genind(x = data_pc_gstud, pop_col = "Cluster",
                ind_col = "ID")

## ----eval=FALSE, echo =TRUE, message = FALSE, warning = FALSE-----------------
#  mat_dps <- mat_gen_dist(x = data_genind, dist = "DPS")

## ---- message = FALSE, warning = FALSE----------------------------------------
mat_dps[1:5, 1:5]

## ---- eval=FALSE, echo =TRUE, message = FALSE, warning = FALSE----------------
#  mat_pg <- mat_gen_dist(x = data_genind, dist = "PG")

## ---- message = FALSE, warning = FALSE----------------------------------------
mat_pg[1:5, 1:5]

## -----------------------------------------------------------------------------
land_graph <- graphab_to_igraph(dir_path = system.file('extdata', 
                                                       package = 'graph4lg'), 
                                nodes = "patches", 
                                links = "liens_rast2_1_11_01_19-links",
                                weight = "cost", fig = FALSE, crds = TRUE)

## -----------------------------------------------------------------------------
crds_patches <- land_graph[[2]]
land_graph <- land_graph[[1]]

## -----------------------------------------------------------------------------
mat_ld <- as_adjacency_matrix(land_graph, attr = "weight", 
                              type = "both", sparse = FALSE)
order <- row.names(mat_ld)[order(c(as.character(row.names(mat_ld))))]
mat_ld <- reorder_mat(mat = mat_ld, order = order)

## -----------------------------------------------------------------------------
head(crds_patches)

## -----------------------------------------------------------------------------
mat_geo <- mat_geo_dist(data = crds_patches, ID = "ID", x = "x", y = "y")
mat_geo <- reorder_mat(mat = mat_geo, order = order) 

## -----------------------------------------------------------------------------
mat_geo[1:5, 1:5]

## ---- out.width = '40%'-------------------------------------------------------
convert_res <- convert_cd(mat_euc = mat_geo, mat_ld = mat_ld, 
                          to_convert = 10000, fig = TRUE, 
                          method = "log-log", pts_col = "grey")
convert_res

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  dmc <- dist_max_corr(mat_gd = mat_dps, mat_ld = mat_ld,
#                       interv = 500, pts_col = "black")

## -----------------------------------------------------------------------------
# DMC value
dmc[[1]]
# Correlation coefficients
dmc[[2]]
# Threshold distances tested
dmc[[3]]

## ---- out.width='40%', eval = TRUE, echo = FALSE------------------------------

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


## ---- out.width='40%'---------------------------------------------------------
scatter_dist(mat_gd = mat_dps, mat_ld = mat_ld, 
             pts_col = "black")

## -----------------------------------------------------------------------------
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
graph_comp <- gen_graph_topo(mat_w = mat_dps, mat_topo = mat_dps,
                             topo = "comp")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  graph_ci <- gen_graph_indep(x = data_genind,
#                              dist = "PCA",
#                              cov = "sq",
#                              adj = "holm")

## -----------------------------------------------------------------------------
graph_ci

## ---- out.width = '60%'-------------------------------------------------------
p <- plot_graph_lg(graph = graph_mst, crds = crds_patches,
                   mode = "spatial", weight = TRUE, width = "inv")
p

## ---- out.width = '60%'-------------------------------------------------------
p <- plot_graph_lg(graph = graph_mst, crds = crds_patches,
                   mode = "aspatial", weight = TRUE, width = "inv")
p

## ---- out.width = '80%'-------------------------------------------------------
plot_graph_modul(graph = graph_gab_geo, crds = crds_patches)

## ---- out.width = '40%'-------------------------------------------------------
scatter_dist_g(mat_y = mat_dps , mat_x = mat_ld, graph = graph_gab_geo)

## ---- out.width = '40%'-------------------------------------------------------
p <- plot_w_hist(graph = graph_gab_gen)
p

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  graph_to_shp(graph = graph_mst, crds = crds_patches, mode = "both",
#               layer_name = "test_shp_mst",
#               dir_path = "wd",
#               metrics = TRUE,
#               crds_crs = "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3
#                        +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs")

## -----------------------------------------------------------------------------
igraph::degree(graph_percol)

## -----------------------------------------------------------------------------
igraph::cluster_fast_greedy(graph_thr, 
                            weights = 1/igraph::E(graph_thr)$weight)

## -----------------------------------------------------------------------------
df_met <- data.frame(ID = V(graph_percol)$name)
df_met$deg <- igraph::degree(graph_percol)
df_met$modul <- igraph::cluster_fast_greedy(graph_thr, 
                                            weights = 1/igraph::E(graph_thr)$weight)$membership

## -----------------------------------------------------------------------------
graph_percol <- add_nodes_attr(graph_percol,
                               input = "df",
                               data = df_met,
                               index = "ID")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  land_graph <- add_nodes_attr(land_graph,
#                               input = "shp",
#                               dir_path = system.file('extdata', package = 'graph4lg'),
#                               layer = "patches",
#                               index = "Id",
#                               include = "Area")

## -----------------------------------------------------------------------------

land_graph2 <- gen_graph_thr(mat_w = mat_ld, mat_thr = mat_ld,
                             thr = 2000, mode = "larger")

graph_topo_compar(obs_graph = land_graph2, 
                  pred_graph = graph_gab_geo,
                  mode = "mcc", 
                  directed = FALSE)

## ---- out.width = '70%'-------------------------------------------------------
graph_plot_compar(x = land_graph2, y = graph_gab_geo, 
                  crds = crds_patches)


## -----------------------------------------------------------------------------
graph_node_compar(x = graph_gab_geo, y = land_graph2,
                  metrics = c("btw", "btw"), method = "spearman",
                  weight = TRUE, test = TRUE)

## -----------------------------------------------------------------------------
graph_modul_compar(x = land_graph2, y = graph_gab_geo)

