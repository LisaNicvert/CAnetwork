# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: function to analyze the outputs of reciprocal scaling (lms)


#' Get mean and variance
#'
#' Return the mean and variances from a reciprocal scaling analysis
#'
#' @param recscal A dataframe (expected to be the output of `reciprocal.coa` function).
#' The first columns must contain the reciprocal scaling scores and the last 3 columns
#' are `Row`, `Col` and `Weight`.
#' @param ax the axes for which to compute mean and variance (the function will use the n-th columns from
#' the `recscal` dataframe).
#'
#' @return A list of 4 matrices, each containing the mean or standard deviation per axes. Each matrix
#' has one column per axis given in `ax` and as many rows as there are grouping levels.
#'
#' + rowsd is the standard deviation per row
#' + rowmean is the mean per row
#' + colsd is the standard deviation per column
#' + colmean is the mean per column
#'
#' Rownames are the row or column groups. Column names are the same as in `recscal`.
#'
#' @export
get_mean_sd <- function(recscal, ax = 1:2) {

  if (length(ax) == 1) {
    warning("farfacwt does not work in one dimension with ade4 version 1.7-22")
  }

  # Get mean and variance for columns
  rowvar <- ade4::varfacwt(recscal[, ax], fac = recscal$Row,
                           wt = recscal$Weight)
  rowmean <- ade4::meanfacwt(recscal[, ax], fac = recscal$Row,
                             wt = recscal$Weight)

  # Get mean and variance for rows
  colvar <- ade4::varfacwt(recscal[, ax], fac = recscal$Col,
                           wt = recscal$Weight)
  colmean <- ade4::meanfacwt(recscal[, ax], fac = recscal$Col,
                             wt = recscal$Weight)


  res <- list(rowsd = sqrt(rowvar),
              rowmean = rowmean,
              colsd = sqrt(colvar),
              colmean = colmean)
  return(res)
}


# Linear models ------------------------------------------------------------


#' Compare 2 linear models
#'
#' This function compares the models `lmsimple` and `lm2` and returns the best
#'
#' @param lmsimple First model (if using method = "LRT", must be nested in lm2)
#' @param lm2 Second model
#' @param method The method to use to compare models: likelihood ratio test
#' (use `method = "LRT"`) (for nested linear models comparison)
#' or AIC (use `method = "AIC"`)
#' @param alpha Threshold to consider for the chi square test
#'
#' @return The best model
#' @export
get_best_model <- function(lmsimple, lm2, method = c("LRT", "AIC"),
                           alpha = 0.05) {

  if (length(method) > 1) {
    method <- method[1]
  }

  if (method == "LRT") {
    # We assume lmsimple is nested in lm2
    lrt <- lmtest::lrtest(lm2, lmsimple)
    pval <- lrt$`Pr(>Chisq)`[2]

    if(pval <= alpha) {
      # Significant difference: the more complex model is better
      return(lm2)
    } else {
      # Non-significant difference: the more complex model is not better
      return(lmsimple)
    }
  } else if (method == "AIC") {
    # Compute models AIC as 2k - 2 log(likelikood)
    AICs <- AIC(lmsimple, lm2)

    AICsimple <- AICs["lmsimple", "AIC"]
    AIC2 <- AICs["lm2", "AIC"]

    # Then we choose the model with the smaller AIC
    if(AICsimple < AIC2) {
      return(lmsimple)
    } else {
      return(lm2)
    }
  }
}

#' Get the predictions of a linear model
#'
#' Returns the model prediction on given data range
#'
#' @param dat_predict The data to predict
#' @param lmpred The model. It must have only one explanatory variable named `mean`.
#' @param by The step to use for the range of `dat_predict` values
#' @param level Confidence level used for the prediction. Defaults to 0.95.
#' @param interval The type of interval to use: the argument is used in `predict.lm` and has the same interpretation.
#' `confidence` gives the confidence interval around the mean of the observations whereas `prediction` gives
#' the confidence interval of predicted values.
#'
#' @return The model's prediction as a dataframe with columns:
#'
#' + `fit`: the predicted value
#' + `lwr`: lower bound of the confidence interval (see arguments `level` and `interval`)
#' + `upr`: upper bound of the confidence interval (see arguments `level` and `interval`)
#' + `x`: the explanatory variable
#'
#' @export
get_pred <- function(dat_predict, lmpred, by = 0.001, level = 0.95,
                     interval = c("confidence", "prediction")) {

  if (length(interval) > 1) {
    int <- interval[1]
  } else {
    int <- interval
  }

  # Get the range of explanatory variable on which to predict values
  newdat <- seq(min(dat_predict), max(dat_predict), by = by)
  newdat <- data.frame(mean = newdat)

  # Predict values over the range of newdat
  pred <- predict(lmpred, interval = int,
                  level = level,
                  newdata = newdat)

  # Return a dataframe
  pred <- as.data.frame(pred)
  pred$x <- newdat$mean

  return(pred)
}

#' Labels of a linear model
#'
#' Get the labels of a linear model to display on a plot. Inspired from https://r-graphics.org/recipe-annotate-facet.
#'
#' @param mod The linear model. It is expected to have 2 or 3 coefficients of the form y = ax + b or y = ax + cx^2 + b.
#' @param a The axis to consider. Used for the subscript of the variables.
#'
#' @return A dataframe with columns `formula` and `r2` Containing respectively the model equation and coefficient of determination.
#' The formulas are written to be parsed later in the plot.
#' The variables' names are s_a^2 for y and m_a for x.
#'
#' @export
lm_labels <- function(mod, a) {

  coef <- stats::coef(mod)
  if (length(coef) == 3) { # quadratic term, slope and intercept
    formula <- sprintf("italic(s[%.0f]) ==  %.2f * italic(m[%.0f]) %+.2f * italic(m[%.0f]^2) %+.2f",
                       a, coef[2], a, coef[3], a, coef[1])
  } else if (length(coef) == 2) { # slope and intercept
    formula <- sprintf("italic(s[%.0f]) ==  %.2f * italic(m[%.0f]) %+.2f",
                       a, coef[2], a, coef[1])
  }

  r2 <- summary(mod)$adj.r.squared
  r2 <- sprintf("italic(R^2) == %.2f", r2)

  res <- data.frame(formula = formula, r2 = r2)
  return(res)
}

# Compare niche measures --------------------------------------------------


#' Get reciprocal scaling niche indices
#'
#' @param recscal reciprocal scaling object
#' @param nfkeep number of axes to keep in the result.
#' This number determines the number of axes coordinates returned as well as
#' the number of axes used in the volumen computation.
#'
#' @return a list with 2 elements: `breadths` and `centroids` (or only `breadths`
#' if `centroids` is `FALSE`).
#' Each element contains a list of 2 dataframes: `rows` with row measures and `cols`
#' with column measures.
#' The elements of `breadths` have `nfkeep` + 2 columns: the niche breadths on each axis
#' (see [get_mean_sd] for details), the area on the first 2 axes and the volume on the first
#' `nfkeep` axes.
#' @export
get_recscal_niche <- function(recscal, nfkeep = 3) {

  niche_mes <- get_mean_sd(recscal, ax = 1:nfkeep)

  # Function to format niche breadths (compute area, volume and select cols)
  fbreadths <- function(df, nfkeep) {
    res <- df

    res$prod <- apply(res, 1, prod)

    res <- res |>
      dplyr::mutate(area = pi*Scor1*Scor2,
                    volume = (pi^(nfkeep/2)/(gamma(nfkeep/2 + 1)))*prod,
                    .before = 1) |>
      dplyr::select(-prod)

    return(res)
  }

  rowsd <- as.data.frame(niche_mes$rowsd)
  rowsd <- fbreadths(rowsd, nfkeep = nfkeep)

  colsd <- as.data.frame(niche_mes$colsd)
  colsd <- fbreadths(colsd, nfkeep = nfkeep)

  rowmean <- as.data.frame(niche_mes$rowmean)
  colmean <- as.data.frame(niche_mes$colmean)

  res <- list(breadths = list(rows = rowsd,
                              cols = colsd),
              centroids = list(rows = rowmean,
                               cols = colmean))
  return(res)
}


#' Get PRN niches
#'
#' @param mat interaction matrix (resources as rows x consumers as columns)
#' @param resource_traits resource traits matrix
#' @param consumer_traits consumer traits matrix
#' @param splits_consumers number of splits for resource niches
#' @param splits_resources number of splits for consumers niches
#' @param centroids return niche centroids
#' @param dimensions the number of dimensions for the convex hull
#' @param ... Additional arguments passed to `fd.niche`
#'
#' @return a list with 2 elements: `breadths` and `centroids` (or only `breadths`
#' if `centroids` is `FALSE`).
#' Each element contains a list of 2 dataframes: `rows` with row measures and `cols`
#' with column measures.
#'
#' @export
get_PRN <- function(mat, resource_traits, consumer_traits,
                    splits_consumers = 100,
                    splits_resources = 100,
                    dimensions,
                    centroids = TRUE,
                    ...) {

  # Consumers
  comm_df <- as.data.frame(matrix_to_df(mat)[, c(2, 1, 3)])
  comm_df <- comm_df |> filter(value != 0)

  PRNj_all <- fd.niche(list(comm_df),
                       list(as.data.frame(resource_traits)),
                       dimensions = dimensions,
                       splits = splits_consumers,
                       ...)
  PRNj <- PRNj_all[[1]][[1]]$contrib.sum

  namej <- names(PRNj)
  newnamej <- stringr::str_replace(pattern = "nw\\d+_", namej, replacement = "")
  names(PRNj) <- newnamej
  PRNj <- PRNj[colnames(mat)]

  # Resources
  comm_df <- as.data.frame(matrix_to_df(mat))
  comm_df <- comm_df |> filter(value != 0)

  PRNi_all <- fd.niche(list(comm_df),
                       list(as.data.frame(consumer_traits)),
                       dimensions = dimensions,
                       splits = splits_resources,
                       ...)

  PRNi <- PRNi_all[[1]][[1]]$contrib.sum

  namei <- names(PRNi)
  newnamei <- stringr::str_replace(pattern = "nw\\d+_", namei, replacement = "")
  names(PRNi) <- newnamei
  PRNi <- PRNi[rownames(mat)]

  if (centroids) {
    PRNcentroidi <- PRNi_all[[1]][[1]]$niche.centroids
    rownames(PRNcentroidi) <- newnamei
    PRNcentroidi <- as.data.frame(PRNcentroidi[rownames(mat),])

    PRNcentroidj <- PRNj_all[[1]][[1]]$niche.centroids
    rownames(PRNcentroidj) <- newnamej
    PRNcentroidj <- as.data.frame(PRNcentroidj[colnames(mat),])

    cols <- data.frame(breadths = PRNj,
                       centroids = PRNcentroidj)
    cols <- data.frame(breadths = PRNj,
                       centroids = PRNcentroidj)

    res <- list(breadths = list(rows = PRNi,
                                cols = PRNj),
                centroids = list(rows = PRNcentroidi,
                                 cols = PRNcentroidj))
  } else {
    res <- list(breadths = list(rows = PRNi,
                                cols = PRNj))
  }

  return(res)
}


#' Get niche breaths measures
#'
#' Returns niche breadths measures computed using different frameworks
#'
#' @param mat interaction matrix
#' @param centroids return niche centroid measures too?
#' @param consumer_traits consumer traits matrix
#' @param resource_traits resource traits matrix
#' @param splits_consumers splits to do for resource traits
#' @param splits_resources splits for consumer traits
#' @param ... additional arguments passed to `fd.niche`
#' @param nf_rec axes to keep for reciprocal scaling
#' @param nf_PRN axes to keep for PRN
#'
#' @return a list with 2 elements: `breadths` and `centroids` (or only `breadths`
#' if `centroids` is `FALSE`).
#' Each element contains a list of 2 dataframes: `rows` with row measures and `cols`
#' with column measures.
#'
#' @export
get_all_niche_measures <- function(mat,
                                   consumer_traits, resource_traits,
                                   nf_rec = 2,
                                   nf_PRN = nf_rec,
                                   splits_consumers = 100,
                                   splits_resources = 100,
                                   centroids = TRUE, ...) {

  mat <- data.frame(mat)

  # Reciprocal scaling ---
  ca <- ade4::dudi.coa(mat, scannf = FALSE, nf = nf_rec)
  rec <- CAnetwork::reciprocal.coa(ca)

  recscal_niche <- get_recscal_niche(rec, nfkeep = nf_rec)

  # d' ---
  di <- bipartite::dfun(mat)
  dj <- bipartite::dfun(t(mat))

  # Partners diversity ---
  # divi <- ade4::divc(data.frame(t(mat)))$diversity
  # divj <- ade4::divc(data.frame(mat))$diversity
  mat01 <- 1*(mat != 0)

  divi <- rowSums(mat01)
  divj <- colSums(mat01)

  # PRN ---
  PRN_niches <- get_PRN(mat = mat,
                        resource_traits = resource_traits,
                        consumer_traits = consumer_traits,
                        splits_consumers = splits_consumers,
                        splits_resources = splits_resources,
                        dimensions = nf_PRN,
                        ...)

  # Traits distance ---
  traitsi <- ade4::divc(data.frame(t(mat)),
                        dis = sqrt(stats::dist(consumer_traits)))$diversity
  traitsj <- ade4::divc(data.frame(mat),
                        dis = sqrt(stats::dist(resource_traits)))$diversity

  if (centroids) {
    # Centroids from traits
    traitsi_centroid <- apply(mat, 1,
                              function(x) meanfacwt(as.matrix(consumer_traits), wt = x))
    traitsi_centroid <- as.data.frame(traitsi_centroid)

    traitsj_centroid <- apply(t(mat), 1,
                              function(x) meanfacwt(as.matrix(resource_traits), wt = x))
    traitsj_centroid <- as.data.frame(traitsj_centroid)
  }

  # Format data ---
  # Niche breadth
  df_row <- data.frame(species = rownames(mat),
                       dprime = di$dprime,
                       d = di$d,
                       PRN = PRN_niches$breadths$rows,
                       spdiv = divi,
                       trdiv = traitsi)
  recscal_nrows <- recscal_niche$breadths$rows |>
    tibble::rownames_to_column("species")
  df_row <- df_row |>
    dplyr::left_join(recscal_nrows, by = "species")

  df_col <- data.frame(species = colnames(mat),
                       dprime = dj$dprime,
                       d = dj$d,
                       PRN = PRN_niches$breadths$cols,
                       spdiv = divj,
                       trdiv = traitsj)
  recscal_ncols <- recscal_niche$breadths$cols |>
    tibble::rownames_to_column("species")
  df_col <- df_col |>
    dplyr::left_join(recscal_ncols, by = "species")

  if (centroids) {
    # Rows ---
    # Format PRN
    centroid_row <- PRN_niches$centroids$rows
    colnames(centroid_row) <- paste0("PRN", 1:ncol(centroid_row))
    centroid_row <- centroid_row |>
      tibble::rownames_to_column("species")

    # Format recscal
    recscal_nrows <- recscal_niche$centroids$rows
    colnames(recscal_nrows) <- paste0("rec", 1:ncol(recscal_nrows))
    recscal_nrows <- recscal_nrows |>
      tibble::rownames_to_column("species")

    # Format traits
    colnames(traitsi_centroid) <- paste0("t", 1:ncol(traitsi_centroid))
    traitsi_centroid <- traitsi_centroid |>
      tibble::rownames_to_column("species")

    # Merge
    centroid_row <- centroid_row |>
      dplyr::left_join(recscal_nrows, by = "species") |>
      dplyr::left_join(traitsi_centroid, by = "species")


    # Cols ---
    # Format PRN
    centroid_col <- PRN_niches$centroids$cols
    colnames(centroid_col) <- paste0("PRN", 1:ncol(centroid_col))
    centroid_col <- centroid_col |>
      tibble::rownames_to_column("species")

    # Format recscal
    recscal_ncols <- recscal_niche$centroids$cols
    colnames(recscal_ncols) <- paste0("rec", 1:ncol(recscal_ncols))
    recscal_ncols <- recscal_ncols |>
      tibble::rownames_to_column("species")

    # Format traits
    colnames(traitsj_centroid) <- paste0("t", 1:ncol(traitsj_centroid))
    traitsj_centroid <- traitsj_centroid |>
      tibble::rownames_to_column("species")

    # Merge
    centroid_row <- centroid_col |>
      dplyr::left_join(recscal_ncols, by = "species") |>
      dplyr::left_join(traitsj_centroid, by = "species")

    # Final data ---
    res <- list(breadths = list(rows = df_row,
                                cols = df_col),
                centroids = list(rows = centroid_row,
                                 cols = centroid_col))
  } else {
    res <- list(breadths = list(rows = df_row,
                                cols = df_col))
  }

  return(res)
}


#' Get correlation table in its long form
#'
#' @param df The dataframe (all columns will be correlated).
#'
#' @return The dataframe with columns `colname`, `rowname` and `value`
#' (the correlation).
#' Some slots in `value` are `NA` and correspond to the repeated values.
#' `colname` and `rowname` are alphabetically ordered factors.
#'
#' @export
get_cor_long <- function(df) {

  cor_df <- cor(df,
                method = "spearman")
  cor_df[lower.tri(cor_df, diag = TRUE)] <- NA
  cor_df <- data.frame(cor_df) |>
    tibble::rownames_to_column()

  cor_long <- data.frame(cor_df) |>
    tidyr::pivot_longer(cols = 2:ncol(cor_df),
                        names_to = "colname")

  # Reorder factors
  cor_long$colname <- factor(cor_long$colname,
                             levels = unique(cor_long$colname))
  cor_long$rowname <- factor(cor_long$rowname,
                             levels = unique(cor_long$rowname))

  return(cor_long)
}
