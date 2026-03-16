# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: reciprocal scaling functions extending reciprocal scaling from Thioulouse & Chessel (1992)

# Different scalings ------------------------------------------------------


#' Reciprocal scaling
#'
#' Perform reciprocal scaling from a count table and optionally
#' 2 other "environment" and "traits" tables.
#'
#' @param Y Count table
#' @param E Environment table (for rows)
#' @param T_ Traits table (for columns)
#'
#' @return A list with R, C and RC results. Each have Scorexx columns, Col, Row and Weight.
#' Each row corresponds to a correspondence in the original table (a nonzero occurrence).
#' Scorexx gives the coordinates of the correspondence in the multivariate space as given
#' by canonical correlation analysis.
#' Row and Col give the row and column this correspondence belongs to.
#' Weight gives its weight (the count of the original cell divided by the sum of the table).
#'
#' @export
recscal <- function(Y, E = NULL, T_ = NULL) {

  # Transform matrix to count table
  P <- Y/sum(Y)
  Pfreq <- as.data.frame(as.table(as.matrix(Y)))
  colnames(Pfreq) <- c("row", "col", "Freq")

  # Remove the cells with no observation
  Pfreq0 <- Pfreq[-which(Pfreq$Freq == 0),]

  # Create indicator tables
  tabR <- ade4::acm.disjonctif(as.data.frame(Pfreq0$row))
  colnames(tabR) <- rownames(Y)
  if (!is.null(E)) {
    tabR <- as.matrix(tabR) %*% as.matrix(E) # duplicate rows of E according to the correspondences of Y
  }

  tabC <- ade4::acm.disjonctif(as.data.frame(Pfreq0$col))
  colnames(tabC) <- colnames(Y)
  if (!is.null(T_)) {
    tabC <- as.matrix(tabC) %*% as.matrix(T_) # duplicate rows of T according to the correspondences of Y
  }

  # Get weights
  wt <- Pfreq0$Freq

  # Center tables
  tabR_scaled <- ade4::scalewt(tabR, wt,
                         scale = FALSE)
  tabC_scaled <- ade4::scalewt(tabC, wt,
                         scale = FALSE)

  res <- stats::cancor(diag(sqrt(wt)) %*% tabR_scaled,
                       diag(sqrt(wt)) %*% tabC_scaled,
                       xcenter = FALSE, ycenter = FALSE)

  # Compute these scores from this coef
  scoreR <- tabR_scaled[, 1:ncol(res$xcoef)]  %*% res$xcoef
  scoreC <- tabC_scaled[, 1:ncol(res$ycoef)]  %*% res$ycoef

  # Get RC score
  mindim <- min(ncol(res$xcoef), ncol(res$ycoef))
  scoreRC <- scoreR[, 1:mindim] + scoreC[, 1:mindim]
  scoreRC_scaled <- ade4::scalewt(scoreRC, wt = wt) # normalisation à 1

  # Summarize the scores ---
  # Columns
  resC <- as.data.frame(scoreC)
  resC <- cbind(resC, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resC) <- c(paste0("Scor", 1:ncol(scoreC)), "Row", "Col", "Weight")
  rownames(resC) <- paste(resC$Row, resC$Col, sep = "_")

  # Rows
  resR <- as.data.frame(scoreR)
  resR <- cbind(resR, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resR) <- c(paste0("Scor", 1:ncol(scoreR)), "Row", "Col", "Weight")
  rownames(resR) <- paste(resR$Row, resR$Col, sep = "_")

  # Rows + Cols
  resRC <- as.data.frame(scoreRC_scaled)
  resRC <- cbind(resRC, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resRC) <- c(paste0("Scor", 1:ncol(scoreRC_scaled)), "Row", "Col", "Weight")
  rownames(resRC) <- paste(resRC$Row, resRC$Col, sep = "_")

  res <- list(R = resR,
              C = resC,
              RC = resRC)
  return(res)
}


#' Reciprocal scaling from multivariate coordinates
#'
#' Perform reciprocal scaling from a multivariate object and
#' the initial table
#'
#' @param Y Count table (converted to relative frequencies inside the function)
#' @param dudi the multivariate object (expects class 'coa', 'pcaiv' or 'dpcaiv'
#' from ade4)
#' @param transpose Whether to invert R and C in the final table (intended for cases
#' where multivariate analysis was performed on a transposed table, like CCA
#' using animal traits for a environment x species table)
#'
#' @return A list with R, C and RC results. Each have Scorexx columns, Col, Row and Weight.
#' Each row corresponds to a correspondence in the original table (a nonzero occurrence).
#' Scorexx gives the coordinates of the correspondence in the multivariate space
#' (it is in essence a duplicated version of the multivariate scores).
#' Row and Col give the row and column this correspondence belongs to.
#' Weight gives its weight (the count of the original cell divided by the sum of the table).
#'
#' @export
recscal_dudi <- function(dudi, Y, transpose = FALSE) {

  # Transform matrix to count table
  P <- Y/sum(Y)
  Pfreq <- as.data.frame(as.table(as.matrix(P)))
  colnames(Pfreq) <- c("row", "col", "Freq")

  # Remove the cells with no observation
  Pfreq0 <- Pfreq[-which(Pfreq$Freq == 0),]

  # Create indicator tables
  if (transpose) {
    mvarR <- dudi$c1
    mvarC <- dudi$l1
  } else {
    mvarR <- dudi$l1
    mvarC <- dudi$c1
  }

  # Get correspondences for rows
  tabR <- ade4::acm.disjonctif(as.data.frame(Pfreq0$row))
  colnames(tabR) <- rownames(Y)
  scoreR <- as.matrix(tabR) %*% as.matrix(mvarR) # duplicate rows of score according to the correspondences of Y

  tabC <- ade4::acm.disjonctif(as.data.frame(Pfreq0$col))
  colnames(tabC) <- colnames(Y)
  scoreC <- as.matrix(tabC) %*% as.matrix(mvarC) # duplicate rows of T according to the correspondences of Y

  # Get weight
  wt <- Pfreq0$Freq

  # Get RC score
  mindim <- min(ncol(scoreC), ncol(scoreR))
  scoreRC <- scoreR[, 1:mindim] + scoreC[, 1:mindim]
  scoreRC_scaled <- ade4::scalewt(scoreRC, wt = wt) # normalisation à 1

  # Summarize the scores ---
  # Columns
  resC <- as.data.frame(scoreC)
  resC <- cbind(resC, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resC) <- c(paste0("Scor", 1:ncol(scoreC)), "Row", "Col", "Weight")
  rownames(resC) <- paste(resC$Row, resC$Col, sep = "_")

  # Rows
  resR <- as.data.frame(scoreR)
  resR <- cbind(resR, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resR) <- c(paste0("Scor", 1:ncol(scoreR)), "Row", "Col", "Weight")
  rownames(resR) <- paste(resR$Row, resR$Col, sep = "_")

  # Rows + Cols
  resRC <- as.data.frame(scoreRC_scaled)
  resRC <- cbind(resRC, Pfreq0$row, Pfreq0$col, Pfreq0$Freq)

  colnames(resRC) <- c(paste0("Scor", 1:ncol(scoreRC_scaled)), "Row", "Col", "Weight")
  rownames(resRC) <- paste(resRC$Row, resRC$Col, sep = "_")

  res <- list(R = resR,
              C = resC,
              RC = resRC)
  return(res)
}


#' Non-reciprocal scaling
#'
#' Preform non-reciprocal scaling and returns the correspondences coordinates in
#' rows and columns multivariate spaces.
#'
#' @param dudi a `dudi` object expected to be of class `coa`, `pcaiv` or `dpcaiv`.
#' @param oritab original data table (matrix or data.frame). It must be a counts matrix
#' (it is converted to relative frequencies inside the function)
#'
#' @return A list with 2 components `H1` and `H2`, containing correspondences respectively
#' in the rows (scaling 1) and the columns space (scaling 2). For detains, see the correspondences
#' private function.
#'
#' @export
nonrecscal <- function(dudi, oritab) {

  # Get correspondences
  if (inherits(oritab, "data.frame")) {
    oritab <- as.matrix(oritab)
  }
  ptab <- oritab/sum(oritab)
  Pfreq <- as.data.frame(as.table(ptab))
  colnames(Pfreq) <- c("row", "col", "Freq")

  nzero <- which(Pfreq$Freq == 0)
  if(length(nzero) != 0) {
    Pfreq0 <- Pfreq[-nzero,]
  } else {
    Pfreq0 <- Pfreq
  }

  # Get H1
  if (inherits(dudi, "coa")) {
    H1_row <- dudi$li
    H1_col <- dudi$c1
  }
  else if (inherits(dudi, "pcaiv")) {
    H1_row <- dudi$ls
    H1_col <- dudi$c1
  } else if (inherits(dudi, "dpcaiv")) {
    H1_row <- dudi$lsR
    H1_col <- dudi$c1
  }

  H1 <- correspondences(Pfreq0 = Pfreq0,
                        Hrow = H1_row, Hcol = H1_col)

  # Get H2
  if (inherits(dudi, "coa")) {
    H2_row <- dudi$l1
    H2_col <- dudi$co
  }
  else if (inherits(dudi, "pcaiv")) {
    H2_row <- dudi$l1
    H2_col <- dudi$co
  } else if (inherits(dudi, "dpcaiv")) {
    H2_row <- dudi$l1
    H2_col <- dudi$lsQ
  }

  H2 <- correspondences(Pfreq0 = Pfreq0,
                        Hrow = H2_row, Hcol = H2_col)

  res <- list("H1" = H1, "H2" = H2)
  return(res)
}


# Correspondences ---------------------------------------------------------


#' Get correspondences
#'
#' Get the correspondences as a mean of the multivariate cordinates, diviede by sqrt(n)
#' .
#'
#' @param Pfreq0 Correspondence table
#' @param Hrow coordinate to use for rows (l x k)
#' @param Hcol coordinate to use for column (c x k)
#'
#' @return A table with k + 3 columns, where the k first columns are the correspondences coordinates
#' and columns `Row`, `Col` and `Weight` are added at the end.
#'
#' @export
correspondences <- function(Pfreq0, Hrow, Hcol) {

  H <- matrix(nrow = nrow(Pfreq0),
              ncol = ncol(Hrow))

  for (k in 1:ncol(Hrow)) { # For each axis
    ind <- 1 # initialize row index
    for (obs in 1:nrow(Pfreq0)) { # For each observation
      i <- Pfreq0$row[obs]
      j <- Pfreq0$col[obs]
      H[ind, k] <- (Hrow[i, k] + Hcol[j, k])/2
      ind <- ind + 1
    }
  }
  H <- cbind(H, Pfreq0[, c("row", "col", "Freq")])
  colnames(H) <- c(paste0("Scor", 1:ncol(Hrow)), "Row", "Col", "Weight")
  rownames(H) <- paste0(H$Row, H$Col)

  return(H)
}


# Get summary statistics --------------------------------------------------


#' Get mean and standard deviation from RC/R/C scores
#'
#' Compute a weighted mean or standard deviation from the R, C or RC scores
#' from reciprocal scaling.
#'
#' @param RCscores The output of reciprocal scaling. A list with components R, C and RC.
#' @param nax Number of axes to take into account
#' @param from Compute the scores from mixed scores (RC) or from the R scores (row scores, used
#' used for the mean/var of columns) and C scores (same in reverse)
#'
#' @return A list of 2 tables: dfli and dfco with the mean and variances per species for each
#' axis up to nax. The column names are "species" (species for which the mean/sd is computed),
#' meanxx (xx is the axis) and sdxx.
#'
#' @export
get_mean_sd_RC <- function(RCscores, nax = 2,
                           from = c("mix", "unique")) {

  from <- from[1]

  if (from == "mix") {
    scoreR <- RCscores$RC
    scoreC <- RCscores$RC
  } else {
    scoreR <- RCscores$C
    scoreC <- RCscores$R
  }
  mrows <- meanfacwt(scoreR[, 1:nax],
                     fac = scoreR$Row,
                     wt = scoreR$Weight)
  mcols <- meanfacwt(scoreC[, 1:nax],
                     fac = scoreC$Col,
                     wt = scoreC$Weight)

  srows <- sqrt(varfacwt(scoreR[, 1:nax],
                         fac = scoreR$Row,
                         wt = scoreR$Weight))
  scols <- sqrt(varfacwt(scoreC[, 1:nax],
                         fac = scoreC$Col,
                         wt = scoreC$Weight))

  dfli <- data.frame(rownames(mrows), mrows, srows)
  colnames(dfli) <- c("species", paste0("mean", 1:nax), paste0("sd", 1:nax))

  dfco <- data.frame(rownames(mcols), mcols, scols)
  colnames(dfco) <- c("species", paste0("mean", 1:nax), paste0("sd", 1:nax))

  res <- list(dfli = dfli,
              dfco = dfco)
  return(res)
}
