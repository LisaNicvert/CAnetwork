
# Libraries ---------------------------------------------------------------
library(here)
library(quarto)


# Set analyses path -------------------------------------------------------
analyses_folder <- here("analyses")

# Define color palette for 02 and 03
colco <- "darkred"
colli <- "cornflowerblue"

# Run analyses ------------------------------------------------------------

# Prepare data (writes data to outputs/01_clean_data)
source(file.path(analyses_folder, "01_clean_data.R"))

# Analyze network Peru1
dataset <- "Peru1"
quarto::quarto_render(file.path(analyses_folder, "02_example_network.qmd"),
                      execute_params = list(dataset = dataset,
                                            colco = colco, 
                                            colli = colli))

# Analyze all networks
quarto::quarto_render(file.path(analyses_folder, "03_all_networks.qmd"),
                      execute_params = list(colco = colco, 
                                            colli = colli))

# Evaluate R2 correction (Appendix F)
# This takes a long time to run (~10 hours)
# Change nrep to a smaller number to make tests
nrep <- 3
quarto::quarto_render(file.path(analyses_folder, "04_corrected_R2.qmd"),
                      execute_params = list(nrep = nrep))
