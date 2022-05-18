# graph4lg 1.6.0

-- Minor changes 
- Update of Graphab software download URL
- Copy of graphab-1.8.jar in the sub-directory plugins within graph4lg_jar for
using the development versions of Graphab functions.
- Add comments regarding development versions of Graphab functions and potential
error when returning metric values already computed previously.


# graph4lg 1.5.0

-- Minor changes 
- Update of pop_rare_gen_index() for allowing the shell/system command to run on both Linux or Windows OS

# graph4lg 1.4.0

-- Major changes: 
- Creation of a new genetic data conversion function (genind_to_structure)
(internal function for the moment)
- Creation of a function for rarefied genetic diversity indices called
pop_rare_gen_index()
(internal function for the moment)
- Creation of a new function for computing corridors in a Graphab project
- Addition of a new contributor to the package
- Modification of graphab_project for allowing users to create multi-habitat
graph projects
- Modification of graphab_metric for allowing users to compute multi-habitat
graph metrics
- Modification of graphab_link for allowing users to include an external
cost surface

-- Minor changes: 
- Update of Graphab software to 2.8
- Creation of a function checking for multihabitat graph projects
- Modification of functions using ggplot2 to avoid aes_string()
- Update of spatstat versions and package names
- Creation of a git repository on Gitlab

-- Bug fixes:

- Modifications of the way the current directory are managed in graphab-based
functions 

# graph4lg 1.2.0

-- Major changes: 
- Creation of 9 new functions linked to Graphab software and internal
- Change of the author email address
- Addition of a CITATION file following the publication of a software paper
 
-- Minor changes: 
- Update of Graphab software to 2.6
- Update of costdist software to 0.4 because of a bug with no data values
- Use of .xml project files for checking arguments
- Update of spatstat versions and package names
 
-- Bug fixes:
- Reset working directory before stopping when it is changed
- Bug due to unique population in file input of genind_to_genepop
- nb argument in graphab_modul was not active


# graph4lg 1.0.1

-- Minor changes: 
- RdMacros packages in Imports instead of Suggests
- RAM allocated to computation when using mat_cost_dist is adjustable. And
 java interface launching is blocked when running this function.
 
-- Bug fixes:
- Order of commands changed in graphab_metric
- Bug due to st_write() fixed

# graph4lg 1.0.0

-- Major changes: 
- Creation of functions to use Graphab java software tool (downloaded by users) 
from R (get_graphab, graphab_project, graphab_link, graphab_graph, graphab_metric, graphab_modul,
graphab_pointset, get_graphab_metric, get_graphab_linkset)
- Creation of mat_cost_dist function integrating costDistance from gdistance
and another java code (downloaded by users)
- Inclusion of geographical distance calculations from polar coordinates in 
mat_geo_dist function
- Compatibility with SNPs markers and 2 allele-digits coded microsatellites
- Creation of compute_node_metric and compute_graph_modul integrating options
previously included in graph_node_compar and graph_modul_compar, and modification
of these two latter functions to integrate the two former
- plot_graph_lg now allows users to change node sizes according to node attribute variables
- Change from proj4string CRS to EPSG when given by users
- Change from rgdal functions to sf functions when the former returned too many
warning messages
- Creation of four vignettes instead of one

- Bug fixes: bug fixes in graph_plan, genepop_to_genind and gstud_to_genind


# graph4lg 0.5.0

- Minor change: change in package dependencies to stick with CRAN requirements

# graph4lg 0.4.0

- Addition of new parameters to gen_graph_topo to allow the construction of k-nearest neighbor graphs
- Removal of a dataset and addition of a new simulated small dataset
- Minor change: modification of function examples to reduce their duration, vignette update, minor changes in functions 'plot_w_hist' and 'scatter_dist_g'
- Bug fixes: modification of pw_mat_to_df function

# graph4lg 0.3.0

- Minor change: weight option in graph_plan
- Bug fixes: vignette modification to solve a Debian error

# graph4lg 0.2.0

- Addition of several functions (pw_mat_to_df, df_to_pw_mat, graph_to_df, graph_plan, set_AB_prob)
- Changes to replace 'class()' instances by 'inherits()' instances in order to make the package compatible with R 4.0.0.
- Minor changes: spelling errors, code lines length

# graph4lg 0.1.1

- Minor changes in the DESCRIPTION file and modifications of the commands writing files in the examples, the vignette and in two functions to conform with CRAN policies.

# graph4lg 0.1.0

- Initial release

