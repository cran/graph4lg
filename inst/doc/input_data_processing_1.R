## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

library(graph4lg)
library(igraph)


## ---- echo = FALSE, eval = TRUE-----------------------------------------------
# Here, we also rely on a data set created only for the vignettes (`data_tuto`) 
# and containing several objects:

data("data_tuto")

mat_dps <- data_tuto[[1]]
mat_pg <- data_tuto[[2]]
graph_ci <- data_tuto[[3]]
dmc <- data_tuto[[4]]
land_graph <- data_tuto[[5]]
mat_ld <- data_tuto[[6]]

## -----------------------------------------------------------------------------
data_genind <- genepop_to_genind(path = paste0(system.file('extdata', 
                                                           package = 'graph4lg'), 
                                               "/gpop_simul_10_g100_04_20.txt"),
                                 n.loci = 20, pop_names = as.character(1:10))
data_genind

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  genind_to_genepop(x = data_genind,
#                    output = "data_gpop_test.txt")

## -----------------------------------------------------------------------------
loci_names <- paste0("LOCI-", as.character(1:20))

ind_names <- as.character(1:200)
pop_names <- as.character(1:10)
data_paru <- structure_to_genind(path = paste0(system.file('extdata', 
                                                           package = 'graph4lg'), 
                                               "/data_ex_str.txt"),
                                 loci_names = loci_names,
                                 pop_names = pop_names,
                                 ind_names = ind_names)
data_paru

## -----------------------------------------------------------------------------
head(data_ex_gstud)

## -----------------------------------------------------------------------------
gstud_to_genind(x = data_ex_gstud, pop_col = "POP",
                ind_col = "ID")

## -----------------------------------------------------------------------------
gen_div <- pop_gen_index(data_ex_genind)
head(gen_div)

## ----eval=FALSE, echo =TRUE, message = FALSE, warning = FALSE-----------------
#  mat_dps <- mat_gen_dist(x = data_genind, dist = "DPS")

## ---- message = FALSE, warning = FALSE----------------------------------------
mat_dps[1:5, 1:5]

## ---- eval=FALSE, echo =TRUE, message = FALSE, warning = FALSE----------------
#  mat_pg <- mat_gen_dist(x = data_genind, dist = "PG")

## ---- message = FALSE, warning = FALSE----------------------------------------
mat_pg[1:5, 1:5]

## -----------------------------------------------------------------------------
head(pts_pop_simul)

## -----------------------------------------------------------------------------
mat_geo <- mat_geo_dist(data = pts_pop_simul, 
                        ID = "ID", x = "x", y = "y",
                        crds_type = "proj")

## -----------------------------------------------------------------------------
mat_geo[1:5, 1:5]

## -----------------------------------------------------------------------------
city_us <- data.frame(name = c("New York City", "Chicago", 
                               "Los Angeles", "Atlanta"),
                   lat  = c(40.75170,  41.87440,
                            34.05420,  33.75280),
                   lon  = c(-73.99420, -87.63940,
                            -118.24100, -84.39360))

mat_geo_us <- mat_geo_dist(data = city_us,
                           ID = "name", x = "lon", y = "lat",
                           crds_type = "polar")

head(mat_geo_us)


## -----------------------------------------------------------------------------
x <- raster::raster(ncol=10, nrow=10, xmn=0, xmx=100, ymn=0, ymx=100)
raster::values(x) <- sample(c(1,2,3,4), size = 100, replace = TRUE)
pts <- data.frame(ID = 1:4,
                  x = c(10, 90, 10, 90),
                  y = c(90, 10, 90, 10))
cost <- data.frame(code = 1:4,
                   cost = c(1, 10, 100, 1000))
mat_cd <- mat_cost_dist(raster = x,
              pts = pts, cost = cost,
              method = "gdistance")
head(mat_cd)

## ---- eval = FALSE------------------------------------------------------------
#  mat_cost_dist(raster = x,
#                pts = pts, cost = cost,
#                method = "java",
#                parallel.java = 2)
#  

## -----------------------------------------------------------------------------
mat_g1 <- mat_geo_dist(pts_pop_ex, 
                       ID = "ID", x = "x", y = "y", 
                       crds_type = "proj")
head(mat_g1)
row.names(mat_g1)

# Reorder mat_g1
mat_g2 <- reorder_mat(mat_g1, 
                      order = as.character(1:10)[order(as.character(1:10))])
head(mat_g2)
row.names(mat_g2)

## -----------------------------------------------------------------------------
mat_ld <- reorder_mat(mat_ld, order = row.names(mat_geo))

convert_res <- convert_cd(mat_euc = mat_geo, mat_ld = mat_ld, 
                          to_convert = 10000, fig = TRUE, 
                          method = "log-log", pts_col = "grey")
convert_res

## -----------------------------------------------------------------------------
df_dist <- pw_mat_to_df(pw_mat = mat_geo)
head(df_dist)

## -----------------------------------------------------------------------------
mat_dist <- df_to_pw_mat(data = df_dist, 
                         from = "id_1", to = "id_2", value = "value")
mat_dist[1:5, 1:5]

