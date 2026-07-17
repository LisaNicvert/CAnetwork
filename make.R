
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

# log-transform traits?
log_transform_traits <- TRUE

# Transform interaction counts using log(x + 1)?
transform_matrix <- FALSE


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
# Change nrep to a smaller number to make tests
nrep <- 300
quarto::quarto_render(file.path(analyses_folder, "05_corrected_R2.qmd"),
                      execute_params = list(nrep = nrep))
