## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

library(graph4lg)
library(igraph)
library(ggplot2)


## ----eval = FALSE-------------------------------------------------------------
# get_graphab()

## ----echo = FALSE, eval = TRUE------------------------------------------------
load(file = paste0(system.file('extdata', package = 'graph4lg'), 
                   "/", "res_g.RDa"))

## -----------------------------------------------------------------------------
rast <- terra::rast(paste0(system.file('extdata', package = 'graph4lg'), 
                   "/", "rast_simul50.tif"))
# convert the raster into a df while keeping coords
r_df <- as.data.frame(rast,
                      xy = TRUE)
# add a categorical field with the raster codes
r_df$code <- as.factor(r_df$rast_simul50)

# Plot it with ggplot2
g <- ggplot(r_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = code)) + coord_equal()+
  theme_bw()+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("#396D35", "#FB9013", "#EDC951", "#80C342", "black", "#396D35"),
                    labels = c("0 - Forest", "1 - Shrublands", "3 - Crops",
                               "4 - Grasslands","5 - Artificial areas", "6 - Forest"),
                    name = "Land use type")+
  labs(x="Longitude",y="Latitude")
g


## ----eval = FALSE-------------------------------------------------------------
# proj_name <- "graphab_example"
# 
# graphab_project(proj_name = proj_name,
#                 raster = "rast_simul50.tif")
# 

## ----echo = FALSE, eval = TRUE------------------------------------------------
proj_name <- "graphab_example"

## ----eval = FALSE-------------------------------------------------------------
# # Habitat creation
# graphab_habitat(proj_name = proj_name,
#                 name = "forest",
#                 type = "raster",
#                 rast_codes = c(0, 5),
#                 minarea = 200)
# 

## ----eval = FALSE-------------------------------------------------------------
# # Bushland
# graphab_habitat(proj_name = proj_name,
#                 name = "shrubland",
#                 type = "raster",
#                 rast_codes = 1,
#                 minarea = 20)

## ----eval = FALSE-------------------------------------------------------------
# graphab_show(proj_path = "graphab_example/graphab_example.xml")

## ----echo = FALSE, eval = FALSE-----------------------------------------------
# graphab_show(proj_path = paste0(system.file('extdata',
#                                             package = 'graph4lg'),
#                    "/", "graphab_example/graphab_example.xml"))
# 

## ----eval = FALSE-------------------------------------------------------------
# graphab_habitat(proj_name = proj_name,
#                 name = "grassland",
#                 type = "vector",
#                 vec_layer = "grassland_patches.gpkg",
#                 vec_capa_field = "capa")

## -----------------------------------------------------------------------------
cost <- data.frame(code = 0:5,
                   cost = c(1, 5, 60, 40, 1000, 1))

print(cost)

## ----eval=FALSE---------------------------------------------------------------
# 
# graphab_link(proj_name = proj_name,
#              distance = "cost",
#              name = "forest_link_planar",
#              habitat = "forest",
#              cost = cost,
#              topo = "planar")
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# graphab_link(proj_name = proj_name,
#              distance = "cost",
#              name = "shrub_link_planar",
#              habitat = "shrubland",
#              cost = cost,
#              topo = "planar")
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# graphab_link(proj_name = proj_name,
#              distance = "cost",
#              name = "inter_forest_shrub",
#              habitat = c("forest", "shrubland"),
#              inter = TRUE,
#              cost = cost,
#              topo = "planar")
# 

## ----eval = FALSE-------------------------------------------------------------
# link_forest <- get_graphab_linkset(proj_name = proj_name,
#                                    linkset = "forest_link_planar")
# print(link_forest[1:6, ])

## ----echo = FALSE, eval = FALSE-----------------------------------------------
# link_forest <- get_graphab_linkset(proj_name = proj_name,
#                                    linkset = "forest_link_planar",
#                                    proj_path = system.file('extdata',
#                                             package = 'graph4lg'))
# print(link_forest[1:6, ])

## ----eval = FALSE-------------------------------------------------------------
# link_param <- get_graphab_linkset_cost(proj_name = proj_name,
#                                       linkset = "forest_link_planar")
# 
# print(link_param)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
# link_param <- get_graphab_linkset_cost(proj_name = proj_name,
#                                        proj_path = system.file('extdata',
#                                             package = 'graph4lg'),
#                     linkset = "forest_link_planar")
# 
# print(link_param)

## ----eval = FALSE-------------------------------------------------------------
# graphab_project_desc(proj_name = proj_name)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
# graphab_project_desc(proj_name = proj_name,
#                      proj_path = system.file('extdata',
#                                             package = 'graph4lg'))

## ----eval=FALSE---------------------------------------------------------------
# graphab_graph(proj_name = proj_name,
#               linkset = "forest_link_planar",
#               name = "graph_forest")

## ----eval = FALSE-------------------------------------------------------------
# # Forest graph
# graphab_graph(proj_name = proj_name,
#               linkset = "forest_link_planar",
#               name = "graph_forest_1")
# # Shrub graph
# graphab_graph(proj_name = proj_name,
#               linkset = "shrub_link_planar",
#               name = "graph_shrub_2")
# # Inter graph
# graphab_graph(proj_name = proj_name,
#               linkset = "inter_forest_shrub",
#               name = "graph_inter_3")
# # Merge the 3 graphs
# graphab_merge_graph(proj_name = proj_name,
#               name = "graph_multi",
#               graphs = c("graph_forest_1",
#                          "graph_shrub_2",
#                          "graph_inter_3"))

## ----eval=FALSE---------------------------------------------------------------
# # Global metric: PC
# pc <- graphab_metric(proj_name = proj_name,
#                graph = "graph_forest",
#                metric = "PC",
#                dist = 10000,
#                prob = 0.05,
#                beta = 1,
#                cost_conv = TRUE)
# pc

## ----echo = FALSE-------------------------------------------------------------
res_g[["PC"]]

## ----eval=FALSE---------------------------------------------------------------
# f <- graphab_metric(proj_name = proj_name,
#                     graph = "graph_forest",
#                     metric = "F",
#                     dist = 10000,
#                     prob = 0.05,
#                     beta = 1,
#                     cost_conv = FALSE)

## ----echo = FALSE-------------------------------------------------------------
print(res_g[["F"]][1:6, ])

## ----eval = FALSE-------------------------------------------------------------
# metric_table <- get_graphab_metric(proj_name = proj_name,
#                    graph = "graph_forest")
# 

## ----eval = FALSE-------------------------------------------------------------
# ec_multi <- graphab_metric(proj_name = proj_name,
#                      graph = "graph_multi",
#                      multihab = "all",
#                      metric = "EC",
#                      dist = 10000,
#                      prob = 0.05,
#                      beta = 1,
#                      cost_conv = TRUE)
# ec_multi

## ----echo = FALSE-------------------------------------------------------------
res_g[["EC_multi"]]

## ----eval = FALSE-------------------------------------------------------------
# f_multi <- graphab_metric(proj_name = proj_name,
#                     graph = "graph_multi",
#                     multihab = "all",
#                     metric = "F",
#                     dist = 10000,
#                     prob = 0.05,
#                     beta = 1,
#                     cost_conv = FALSE)
# f_multi

## ----echo = FALSE-------------------------------------------------------------
print(res_g[["F_multi"]][1:6, ])

## ----eval=FALSE---------------------------------------------------------------
# graphab_modul(proj_name = proj_name,
#               graph = "graph_forest",
#               dist = 10000,
#               prob = 0.05,
#               beta = 1)

## ----eval=FALSE---------------------------------------------------------------
# get_graphab_linkset(proj_name = proj_name,
#                     linkset = "forest_link_planar")

## ----echo = FALSE-------------------------------------------------------------
print(res_g[["LK"]][1:6, ])

## ----eval=FALSE---------------------------------------------------------------
# get_graphab_metric(proj_name = proj_name,
#                    habitat = "forest")

## ----echo = FALSE-------------------------------------------------------------
print(res_g[["MET"]][1:6, ])

## ----eval=FALSE---------------------------------------------------------------
# land_graph <- graphab_to_igraph(proj_name = proj_name,
#                                 linkset = "forest_link_planar",
#                                 habitat = "forest",
#                                 weight = "cost",
#                                 fig = FALSE,
#                                 crds = TRUE)
# 
# crds_patches <- land_graph[[2]]
# land_graph <- land_graph[[1]]

## ----echo = FALSE-------------------------------------------------------------
crds_patches <- res_g[["CRDS"]]
land_graph <- res_g[["LGRAPH"]]

## -----------------------------------------------------------------------------
plot_graph_lg(land_graph,
              crds = crds_patches,
              mode = "spatial",
              node_size = "area")

