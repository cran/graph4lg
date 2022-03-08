## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

library(graph4lg)
library(igraph)
library(ggplot2)


## ---- eval = FALSE------------------------------------------------------------
#  get_graphab()

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
load(file = paste0(system.file('extdata', package = 'graph4lg'), 
                               "/", "res_g.RDa"))


## -----------------------------------------------------------------------------
load(file = paste0(system.file('extdata', package = 'graph4lg'), 
                               "/", "rast_simul50.RDa"))

r.spdf <- as(rast, "SpatialPixelsDataFrame")
r.df <- as.data.frame(r.spdf)

r.df$layer <- as.factor(r.df$rast_simul50)

g <- ggplot(r.df, aes(x=x, y=y)) + geom_tile(aes(fill = layer)) + coord_equal()+
  theme_bw()+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("#396D35", "#FB9013", "#EDC951", "#80C342", "black", "#396D35"),
                    labels = c("0 - Forest", "1 - Shrublands", "2 - Crops", 
                               "3 - Grasslands","4 - Artificial areas", "5 - Forest"),
                    name = "Land use type")+
  labs(x="Longitude",y="Latitude")
g

## ---- eval = FALSE------------------------------------------------------------
#  proj_name <- "graphab_example"
#  
#  graphab_project(proj_name = proj_name,
#                  raster = "rast_simul50.tif",
#                  habitat = c(0, 5),
#                  minarea = 200)

## -----------------------------------------------------------------------------
cost <- data.frame(code = 0:5,
                   cost = c(1, 5, 60, 40, 1000, 1))

cost

## ---- eval=FALSE--------------------------------------------------------------
#  
#  graphab_link(proj_name = proj_name,
#               distance = "cost",
#               cost = cost,
#               name = "lkst1",
#               topo = "planar")

## ---- eval=FALSE--------------------------------------------------------------
#  graphab_graph(proj_name = proj_name,
#                linkset = "lkst1",
#                name = "graph")

## ---- eval=FALSE--------------------------------------------------------------
#  # Global metric: PC
#  graphab_metric(proj_name = proj_name,
#                 graph = "graph",
#                 metric = "PC",
#                 dist = 10000,
#                 prob = 0.05,
#                 beta = 1,
#                 cost_conv = TRUE)

## ---- echo = FALSE------------------------------------------------------------
res_g[["PC"]]

## ---- eval=FALSE--------------------------------------------------------------
#  f <- graphab_metric(proj_name = proj_name,
#                 graph = "graph",
#                 metric = "F",
#                 dist = 10000,
#                 prob = 0.05,
#                 beta = 1,
#                 cost_conv = FALSE)

## ---- echo = FALSE------------------------------------------------------------
res_g[["F"]][1]
head(res_g[["F"]][[2]])

## ---- eval=FALSE--------------------------------------------------------------
#  graphab_modul(proj_name = proj_name,
#                graph = "graph",
#                dist = 10000,
#                prob = 0.05,
#                beta = 1)

## -----------------------------------------------------------------------------
# Point data frame
head(pts_pop_simul)

## ---- eval=FALSE--------------------------------------------------------------
#  graphab_pointset(proj_name = proj_name,
#                   linkset = "lkst1",
#                   pointset = pts_pop_simul)

## ---- echo = FALSE------------------------------------------------------------
head(res_g[["PTSG"]])

## ---- eval=FALSE--------------------------------------------------------------
#  get_graphab_linkset(proj_name = proj_name,
#                      linkset = "lkst1")

## ---- echo = FALSE------------------------------------------------------------
head(res_g[["LK"]])

## ---- eval=FALSE--------------------------------------------------------------
#  get_graphab_metric(proj_name = proj_name)

## ---- echo = FALSE------------------------------------------------------------
head(res_g[["MET"]])

## ---- eval=FALSE--------------------------------------------------------------
#  land_graph <- graphab_to_igraph(proj_name = proj_name,
#                                  linkset = "lkst1",
#                                  nodes = "patches",
#                                  weight = "cost",
#                                  fig = TRUE,
#                                  crds = TRUE)
#  
#  crds_patches <- land_graph[[2]]
#  land_graph <- land_graph[[1]]

## ---- echo = FALSE------------------------------------------------------------
crds_patches <- res_g[["CRDS"]]
land_graph <- res_g[["LGRAPH"]]

## -----------------------------------------------------------------------------
plot_graph_lg(land_graph,
              crds = crds_patches,
              mode = "spatial",
              node_size = "Area")

