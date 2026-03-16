# Libraries ---------------------------------------------------------------
library(here)

devtools::load_all()


# Get data folders ---------------------------------------------------------------
data_folder <- here("data/ANDEAN_frugivory")
out_folder <- here("outputs/01_clean_data")

# ANDEAN cleaning ---------------------------------------------------------
all_files <- list.files(data_folder)

## Read traits tables ----------
plant_traits_file <- file.path(data_folder,
                               grep("Plant_traits", all_files, value = TRUE))
plant_traits_all <- read.csv(plant_traits_file,
                             sep = "\t")

bird_traits_file <- file.path(data_folder,
                                grep("Bird_traits", all_files, value = TRUE))
bird_traits_all <- read.csv(bird_traits_file,
                              sep = "\t")

## Format traits table ----------
plant_traits_all <- plant_traits_all |>
  dplyr::rename("plant_species" = "Species")

bird_traits_all <- bird_traits_all |>
  dplyr::rename("bird_species" = "Species")

## List network files ----------
network_files <- grep("NW_", all_files,
                      value = TRUE)

for (i in 1:length(network_files)) {
  filename_i <- network_files[i]
  print(paste0("Cleaning dataset ", filename_i, " (", 
               i, "/", length(network_files), ") ----------------------"))
  interactions <- read.csv(file.path(data_folder,
                                     filename_i),
                           sep = "\t")

  ## Transform matrix to dataframe ----------
  colnames(interactions)[1] <- "plant_species"
  interactions <- interactions |> 
    tidyr::pivot_longer(cols =2:ncol(interactions),
                        names_to = "bird_species",
                        values_to = "frequency")

  ## Select traits subset ----------
  country <- gsub(filename_i,
                  pattern = "(^NW_)([[:alpha:]]+).*(\\.txt$)",
                  replacement = "\\2")

  bird_traits <- bird_traits_all |>
    filter(Country == country) |>
    filter(bird_species %in% interactions$bird_species)

  plant_traits <- plant_traits_all |>
    filter(Country == country) |>
    filter(plant_species %in% interactions$plant_species)

  ## Create codes ----------
  plant_codes <- create_code(interactions$plant_species,
                             name = "plant_species")
  bird_codes <- create_code(interactions$bird_species,
                              name = "bird_species")

  ## Add codes to data ----------
  ### Add to interactions (remove whole name) -----
  interactions <- interactions |>
    dplyr::left_join(bird_codes, by = "bird_species")

  interactions <- interactions |>
    dplyr::left_join(plant_codes, by = "plant_species")

  interactions <- interactions |>
    dplyr::select(plant_species, plant_species_code,
                  bird_species, bird_species_code,
                  frequency, everything())

  ### Add to traits -----
  bird_traits <- bird_traits |>
    dplyr::left_join(bird_codes, by = "bird_species") |>
    dplyr::relocate(bird_species_code, .before = bird_species) |>
    dplyr::select(-Country)

  plant_traits <- plant_traits |>
    dplyr::left_join(plant_codes, by = "plant_species") |>
    dplyr::relocate(plant_species_code, .before = plant_species) |>
    dplyr::select(-Country)

  subfolder <- gsub(filename_i,
                    pattern = "(^NW_)(.*)(\\.txt$)",
                    replacement = "\\2")
  subfolder <- file.path(out_folder, subfolder)
  
  if (!dir.exists(subfolder)) {
    dir.create(subfolder, recursive = TRUE, showWarnings = FALSE)
  }
  
  write.csv(interactions, file.path(subfolder, "interactions.csv"),
            row.names = FALSE)
  write.csv(plant_traits, file.path(subfolder, "plant_traits.csv"),
            row.names = FALSE)
  write.csv(bird_traits, file.path(subfolder, "bird_traits.csv"),
            row.names = FALSE)
}
