
# Libraries ---------------------------------------------------------------
library(here)
library(quarto)

# Set analyses path -------------------------------------------------------
analyses_folder <- here("analyses")

# Define color palette for 02 and 03
colco <- "darkred"
colli <- "cornflowerblue"

# Define threshold to filter interactions
thr <- 1


# Data transformations ----------------------------------------------------

# To reproduce results from the main text: 
# set log_transform_traits to FALSE and transform_matrix to ""

# To reproduce results from Appendix S1 Section S8, 2 runs are needed:
#   1) log_transform_traits = TRUE and transform_matrix = ""
#   2) log_transform_traits = FALSE and transform_matrix = "log1p"

# log-transform traits?
log_transform_traits <- FALSE

# Transform interaction counts
# Possible values: 
# - "log1p" for log(x + 1),
# - "N2" for N2 preprocessing (see https://doi.org/10.32942/X2TT0B)
# - "" for no transformation
transform_matrix <- ""


# Run analyses ------------------------------------------------------------

# Plot Figure 1 (writes Figure to figures/01_Fig_1_summary/)
source(file.path(analyses_folder, "01_Fig_1_summary.R"))

# Prepare data (writes data to outputs/02_clean_data)
source(file.path(analyses_folder, "02_clean_data.R"))

# Analyze network Peru1
dataset <- "Peru1"
quarto::quarto_render(file.path(analyses_folder, "03_example_network.qmd"),
                      execute_params = list(dataset = dataset,
                                            thr = thr,
                                            colco = colco, 
                                            colli = colli,
                                            log_transform_traits = log_transform_traits))

# Analyze all networks
quarto::quarto_render(file.path(analyses_folder, "04_all_networks.qmd"),
                      execute_params = list(thr = thr, 
                                            colco = colco, 
                                            colli = colli,
                                            log_transform_traits = log_transform_traits,
                                            transform_matrix = transform_matrix))

# Evaluate R2 correction (Appendix S1 Section S6)
# This takes a long time to run (~17 hours) so results are precomputed 
# when rendering the Quarto document
# Change nrep to a smaller number to make tests
nrep <- 300
quarto::quarto_render(file.path(analyses_folder, "05_corrected_R2.qmd"),
                      execute_params = list(nrep = nrep))
