# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2023-12-15
#
# Script Description: code to generate species abundances on an environmental gradient
# (can also be interpreted as consumer-resource species interactions)
# adapted from Stéphane Dray (1/12/2004 Montreal from an idea of J Oksanen)


#' Simulate a community matrix
#'
#' Simulate abundance values based on species traits. The ecological interpretation can be
#' a sites x species matrix or a resource x consumer species interaction matrix.
#'
#' @param nsp Number of species / consumer
#' @param nsite Number of sites / resource species
#' @param le.grad gradient length for traits/environmental variables. The values are
#' comprised between 0 and `le.grad`.
#' @param buffer buffer to allow species / consumer species optima to fall outside the gradient by
#' 0 - `buffer` ; `le.grad` + `buffer`.
#' @param col_abund Optional vector of abundances for species/consumers.
#' If missing, abundances are taken from a uniform law between `min.ab` and `max.ab`.
#' @param row_abund Optional vector of abundances in sites/ for resource species.
#' If missing, abundances are taken from a uniform law between `min.ab` and `max.ab`.
#' @param min.ab minimal abundance value for the uniform law (if `col_abund` and/or `row_abund` are absent)
#' @param max.ab maximal abundance value for the uniform law (if `col_abund` and/or `row_abund` are absent)
#' @param mean.tol Mean value for the species/consumers niche tolerance (standard deviation).
#' The species/consumers niche tolerance are drawn from a normal law of mean `mean.tol`
#' and standard deviation `sd.tol` for each trait axis.
#' @param sd.tol Standard deviation for the species/consumers niche tolerance (standard deviation)
#' The species/consumers niche tolerance are drawn from a normal law of mean `mean.tol`
#' and standard deviation `sd.tol` for each trait axis.
#' @param ratio.grad length of the second gradient as a ratio of the first one.
#' @param col_prefix Prefix for the column names of the matrix. Defaults to "s" for "species".
#' @param row_prefix Prefix for the row names of the matrix. Defaults to "e" for "environment".
#' @param rowvar_prefix Prefix for the column names of the row variable. Defaults to "env" for "environment".
#' @param remove_zeroes If there are unobserved species, should we keep them in the final matrix?
#' @param ninter Total number of observations in the web
#' @param return_intermediate Return intermediate results? If yes, will return the
#' probability of interaction based on matching only, on abundances only and the predicted
#' mix of matching and abundances (before sampling).
#' @param delta Exponent between 0 and 1 to give more (1) or less (0) weight to trait matching relatively to
#' abundance.
#'
#' @return An object of class compas, which is a list with 3 components (4 if `return_intermediate`):
#'
#' + comm: the community matrix of dimension `nsite` x `nsp`. Row names are `row_prefix`
#' followed by 1..`nsite` and column names are `col_prefix` followed by 1..`nsp`.
#' + rowvar: the values of the row variables (resource traits/environment gradient values). It is a
#' nsite x 2 matrix. Row names are `row_prefix` followed by 1..`nsite` and column names are
#' `rowvar_prefix` followed by 1 and 2.
#' + colvar: the values of the column variables (consumer/species traits values). It is a `nsp` x 2 x 2 array.
#' Row names are `row_prefix` followed by 1..`nsp`. Column names in each array dimension are `mean` and `sd` (for
#' the mean and standard deviation of the species niche, respectively). Each array 3rd dimension corresponds
#' to a resource trait or environmental gradient value.
#' + if `return_intermediate`, a fourth component named intermediate is returned. The first element
#' is named `p_matching` and contains the matrix of probabilities of interactions based only on matching.
#' The second `ab_neutral` contains the matrix of count of interactions as predicted by species abundances. The third
#' `p_mix`  contains the matrix of probabilities of interactions taking into account the mix of matching and abundances.
#'
#' @export
"compas2d" <-
function (nsp = 40, le.grad = 100, ratio.grad = 0.8, nsite = 100,
          col_abund = NULL,
          row_abund = NULL,
          col_prefix = "s", row_prefix = "e", rowvar_prefix = "env",
          remove_zeroes = FALSE,
          buffer = 1, min.ab = 1, max.ab = 100,
          ninter = 100,
          delta = 1,
          mean.tol = 2, sd.tol = 10,
          return_intermediate = FALSE)
{

  # Gradient for axis 2 ---
  # The length of the second gradient is a fraction of the length of the
  # first gradient
  gradmin2 <- (1-ratio.grad) / 2 * le.grad
  gradmax2 <- le.grad - (1-ratio.grad) / 2 * le.grad
  # For the columns, we include a buffer that is also on axis 2
  # (if buffer = 0 then gradmin/max2_buffer and gradmin/max2 are equal)
  gradmin2_buffer <- -buffer + gradmin2
  gradmax2_buffer <- buffer + gradmax2

  # Environment/resource species trait gradient ---
  # Initialize empty trait/environment values matrix
  x <- matrix(0, nrow = nsite, ncol = 2)
  # For each site, generate an environmental gradient value at random (uniform) (and sort them)
  x[, 1] <- sort(runif(nsite, min = 0, max = le.grad))
  x[, 2] <- runif(nsite, min = gradmin2, max = gradmax2)

  # Abundances ---
  # Species/consumer species
  if (is.null(col_abund)) {
    # Generate random abundances for each species (uniform law)
    col_abund <- runif(nsp, min = min.ab, max = max.ab)
  }
  # Sites/resource species
  if (is.null(row_abund)) {
    # Generate random abundances for each site (uniform law)
    row_abund <- runif(nsite, min = min.ab, max = max.ab)
  }

  # Species niche/consumer niche ---
  # Initialize an array for species optima and tolerances
  # The array has last dimension 2 (one for each trait)
  p <- array(0, dim = c(nsp, 2, 2))

  # Fill array p

  # -> first dimension
  # Generate environmental optima for each species
  p[, 1, 1] <- runif(nsp, min = 0-buffer, max = le.grad+buffer)
  # Generate random niche widths (sd) for each species
  p[, 2, 1] <- abs(rnorm(nsp, mean = mean.tol, sd = sd.tol))

  # -> second dimension
  # Optima
  p[, 1, 2] <- runif(nsp,
                     min = gradmin2_buffer,
                     max = gradmax2_buffer)
  # Widths (sd)
  p[, 2, 2] <- abs(rnorm(nsp, mean = mean.tol, sd = sd.tol))

  # Probability matrix (only matching) ---
  # Initialize empty community matrix
  p_matching <- matrix(0, nrow = nrow(x), ncol = nsp)

  for (i in 1:nrow(x)) {
    for (j in 1:nsp) {
      # Fill each cell with a "presence probability" or an "interaction probability"
      # from a multivariate normal law depending on:
      # - the matching between resource species trait and consumer species trait
      # - the matching between the site environmental value and the species niche optimum on this gradient
      p_matching[i,j] <- mvtnorm::dmvnorm(x[i,] - p[j, 1, ], sigma = diag(p[j, 2, ]^2))
    }
  }
  # Transform negative values to zero probability
  p_matching <- ifelse(p_matching > 0, p_matching, 0)

  # Make columns a proba distribution
  # - each species (column) is distributed in sites following a probability of choosing this site
  # - each consumer (column) chooses the resource according to a given proba of presence
  p_matching <- sweep(p_matching, 2, STATS = apply(p_matching, 2, sum), FUN = "/")

  # Quick patch (in case there are species with zero obs that became NA at the division step)
  p_matching[is.na(p_matching)] <- 0

  # Theoretical interaction count (based on abundances) ---
  # This code makes sense only for interaction matrices because the abundance of
  # resource species is a limiting factor here.
  # It makes less sense in the context of species x sites association, unless we consider
  # resources in the different sites are limited and a site has a limited number of species.
  prop_row <- row_abund/sum(row_abund) # Get the proportion of each plant (its "availability")
  ab_neutral <- prop_row %*% t(col_abund) # Get the predicted abundance -> the birds
  # choose a plant only based on its availability
  # Each bird abundance is exactly col_abund and the plants abundances are proportional to
  # row_abund

  # Probability matrix (matching and abundance) ---
  # these theoretical abundances of interactions must then be confronted to the probability of interactions
  ab_mix <- ab_neutral*p_matching^delta # We multiply expected interaction number by trait matching
  # Create a vector of probabilities that takes into account the matching
  p_mix <- ab_mix/sum(ab_mix)
  p_mix[is.na(p_mix)] <- 0 # Handle divisions per zero

  # Final community matrix (with counts) ---
  # Then we need to sample observations in this web
  # The method is inspired by Fründ et al 2016
  # We sample ninter interactions from a multinomial distribution, where the probability to
  # draw each interaction depends on the probability taking into account abundance and matching
  ab_obs_vec <- stats::rmultinom(n = 1,
                                 size = ninter,
                                 prob = as.numeric(p_mix))
  ab_obs <- matrix(ab_obs_vec, nrow = nsite, byrow = FALSE) # Reformat to a matrix

  if (remove_zeroes) {
    # Remove species that haven't been observed
    colkeep <- which(colSums(ab_obs) != 0)
    rowkeep <- which(rowSums(ab_obs) != 0)

    ab_obs <- ab_obs[, colkeep]
    ab_obs <- ab_obs[rowkeep, ]

    # Select traits only species that have been observed
    p <- p[colkeep, , ]
    x <- x[rowkeep, ]
  }

  # Format result
  # Rows and columns names
  rnames <- paste0(row_prefix, 1:nrow(ab_obs))
  cnames <- paste0(col_prefix, 1:ncol(ab_obs))

  rownames(ab_obs) <- rnames
  colnames(ab_obs) <- cnames

  rownames(x) <- rnames
  colnames(x) <- paste0(rowvar_prefix, 1:2)

  rownames(p) <- cnames
  colnames(p) <- c("mean", "sd")

  if (return_intermediate) {
    rownames(p_matching) <- rnames
    colnames(p_matching) <- cnames

    rownames(ab_neutral) <- rnames
    colnames(ab_neutral) <- cnames

    rownames(p_mix) <- rnames
    colnames(p_mix) <- cnames

    intermediate <- list(p_matching = p_matching,
                         ab_neutral = ab_neutral,
                         p_mix = p_mix)

    sol <- list(comm = ab_obs,
                rowvar = x,
                colvar = p,
                intermediate = intermediate)
  } else {
    sol <- list(comm = ab_obs,
                rowvar = x,
                colvar = p)
  }


  class(sol) <- "compas"
  sol$call <- match.call()

  return(sol)
}
