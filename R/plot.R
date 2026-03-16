# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2023-12-18
#
# Script Description: plotting functions


# Plot model matrices -----------------------------------------------------


#' Plot matrices for the model
#'
#' Function to plot matrices Y, E (environment/resource traits) and T (species traits/consumer traits)/
#' This function plots the matrices arranged in a square, with the user-given proportions.
#'
#' @param mar horizontal margin between matrices
#' @param marl vertical margin between matrices
#' @param E Plot matrix E (rows variables)?
#' @param T_ Plot matrix T (columns variables)?
#' @param r number of rows of matrix Y
#' @param c number of columns of matrix Y
#' @param l number of columns of matrix E
#' @param k number of columns of matrix T
#' @param Yname name of the matrix Y (center matrix)
#' @param Ename name of the matrix E (row variables matrix)
#' @param Tname name of the matrix T (column variables matrix)
#' @param rname name of the dimension r
#' @param cname name of the dimension c
#' @param horizr display r dimension horizontally? (else, will be vertical)
#' @param plot.margin.x horizontal margin around plot
#' @param plot.margin.y vertical margin around plot
#' @param name_mat display matrices names?
#' @param name_dim display dimension names?
#'
#' @return a ggplot of matrix Y (and optionally E and T),
#' with their dimensions written above and aligned in a square where
#' the fourth-corner is missing.
#' @export
#'
#' @examples plotmat(r = 5, c = 2)
plotmat <- function(mar = 2, marl = 1.2,
                    E = FALSE, T_ = FALSE,
                    r, c, l = NULL, k = NULL,
                    name_mat = TRUE,
                    name_dim = TRUE,
                    Yname = "italic(Y)",
                    Ename = "italic(E)",
                    Tname = "italic(T)",
                    rname = "italic(r)", cname = "italic(c)",
                    horizr = TRUE,
                    plot.margin.x = 1,
                    plot.margin.y = 1) {

  xmin <- 0
  xmax <- c
  ymin <- 0
  ymax <- r

  xmat <- c/2
  ymat <- r/2

  xcol <- c(-marl, c/2)
  ycol <- c(r/2, r+marl)
  vjust <- c(ifelse(horizr, 0.5, 0), 0)
  hjust <- c(ifelse(horizr, 1, 0.5), 0.5)

  labmat <- Yname
  labcol <- c(rname, cname)

  xlim <- c(0, c)
  ylim <- c(0, r)

  if (E) {
    xmin <- c(xmin, c+mar)
    xmax <- c(xmax, c+mar+l)
    ymin <- c(ymin, 0)
    ymax <- c(ymax, r)

    xmat <- c(xmat, mar+c+l/2)
    ymat <- c(ymat, r/2)

    xcol <- c(xcol, c+mar+l/2)
    ycol <- c(ycol, r+marl)
    vjust <- c(vjust, 0)
    hjust <- c(hjust, 0.5)

    labmat <- c(labmat, Ename)
    labcol <- c(labcol, 'italic(l)')

    xlim[2] <- xlim[2] + mar + l
  }
  if (T_) {
    xmin <- c(xmin, 0)
    xmax <- c(xmax, c)
    ymin <- c(ymin, -mar)
    ymax <- c(ymax, -mar-k)

    xmat <- c(xmat, c/2)
    ymat <- c(ymat, -mar-k/2)

    xcol <- c(xcol, -marl)
    ycol <- c(ycol, -mar-k/2)
    vjust <- c(vjust, ifelse(horizr, 0.5, 0))
    hjust <- c(hjust, ifelse(horizr, 1, 0.5))

    labmat <- c(labmat, Tname)
    labcol <- c(labcol, 'italic(k)')

    ylim[1] <- ylim[1] - mar - k
  }

  if (horizr) {
    angle <- 0
  } else {
    angle <- c(90, 0)
    if (E) {
      angle <- c(angle, 0)
    }
    if (T_) {
      angle <- c(angle, 90)
    }
  }

  xlim[1] <- xlim[1] - plot.margin.x*marl
  xlim[2] <- xlim[2] + plot.margin.x*marl

  ylim[1] <- ylim[1] - plot.margin.x*marl
  ylim[2] <- ylim[2] + plot.margin.x*marl

  g <- ggplot() +
    annotate("rect",
             xmin = xmin,
             xmax = xmax,
             ymin = ymin,
             ymax = ymax,
             fill = "white", color = "black") +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    coord_fixed() +
    ggplot2::theme_void()

  if (name_mat) {
    g <- g +
      annotate("text",
               x = xmat,
               y = ymat,
               label = labmat,
               parse = TRUE, size = 12, family = "serif")
  }
  if (name_dim) {
    g <- g +
      annotate("text",
               x = xcol,
               y = ycol,
               vjust = vjust,
               hjust = hjust,
               label = labcol,
               angle = angle,
               parse = TRUE, size = 8, family = "serif")
  }
  g
}

# Multivariate ------------------------------------------------------------


#' Plot eigenvalues
#'
#' Barplot of eigenvalues
#'
#' @param eigenvalues The eigenvalues
#' @param showrank Should the x-axis display the rank of the eigenvalues?
#'
#' @return A ggplot
#' @export
plot_eig <- function(eigenvalues, showrank = FALSE) {

  gg <- ggplot() +
    geom_col(aes(y = eigenvalues, x = factor(1:length(eigenvalues)))) +
    theme_linedraw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    xlab("Axis") +
    ylab("Eigenvalues")

  if (!showrank) {
    gg <- gg +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  gg
}


#' Plot bi, tri or quadriplots
#'
#' Bi, tri or quadriplots from a multivariate analysis
#'
#' @param indiv_row Matrix of individuals coordinates for rows (type `dudi$li`)
#' @param indiv_col Matrix of individuals coordinates for columns (type `dudi$co`)
#' @param indiv_row_lab Labels for row individuals
#' @param indiv_col_lab Labels for columns individuals
#' @param var_row Matrix of variables coordinates for rows (type `dudi$corR`)
#' @param var_col Matrix of variables coordinates for columns (type `dudi$corQ`)
#' @param var_row_lab Labels for row variables
#' @param var_col_lab Labels for columns variables
#' @param row_color Color for the row individuals or variables
#' @param col_color Color for the columns individuals or variables
#' @param eig Eigenvalues vector
#' @param x Axis to use for the x-axis
#' @param y Axis to use for the y-axis
#' @param max.overlaps `max.overlaps` argument for `ggrepel::geom_text_repel`
#' @param mult factor to multiply the vectors (arrows on the plot)
#' @param xlim x-axis limits. Defaults to data range if not specified
#' @param ylim y-axis limits. Defaults to data range if not specified
#' @param grad graduations for the major.grid
#' @param title plot title
#' @param alphapoints Transparency value for points
#' @param labsize text size for point labels
#'
#' @return a ggplot
#'
#' @export
multiplot <- function(indiv_row = NULL, indiv_col = NULL,
                      indiv_row_lab = rownames(indiv_row), indiv_col_lab = rownames(indiv_col),
                      var_row = NULL, var_col = NULL,
                      var_row_lab = rownames(var_row), var_col_lab = rownames(var_col),
                      row_color = "black", col_color = "black",
                      eig = NULL,
                      x = 1, y = 2,
                      xlim = NULL, ylim = NULL,
                      grad = 4,
                      mult = 1,
                      title = NULL,
                      max.overlaps = 20,
                      labsize = 3,
                      alphapoints = 0.8) {

  xlab <- paste0("Axis ", x)
  ylab <- paste0("Axis ", y)

  if (!is.null(eig)) {
    xlab <- paste0(xlab, " (", round(eig[x]/sum(eig)*100, 1), " % variability)")
    ylab <- paste0(ylab, " (", round(eig[y]/sum(eig)*100, 1), " % variability)")
  }


  g <- ggplot() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    theme_linedraw() +
    coord_fixed() +
    xlab(xlab) +
    ylab(ylab)

  # Plot title
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  }

  # Set x axis limits and breaks
  if (!is.null(xlim)) {
    g <- g + scale_x_continuous(breaks = seq(xlim[1] - grad, xlim[2] + grad,
                                             by = grad),
                                limits = xlim)
  } else {
    g <- g + scale_x_continuous(breaks = seq(-50, 50, by = grad))
  }

  # Set y axis limits and breaks
  if (!is.null(ylim)) {
    g <- g + scale_y_continuous(breaks = seq(ylim[1] - grad, ylim[2] + grad,
                                             by = grad),
                                limits = ylim)
  } else {
    g <- g + scale_y_continuous(breaks = seq(-50, 50, by = grad))
  }

  # Plot row individuals
  if (!is.null(indiv_row)) {
    g <- g +
      geom_point(aes(x = indiv_row[, x], y = indiv_row[, y]),
                 col = row_color,
                 alpha = alphapoints)
  }

  # Plot column individuals
  if (!is.null(indiv_col)) {
    g <- g +
      geom_point(aes(x = indiv_col[, x], y = indiv_col[, y]),
                 col = col_color,
                 alpha = alphapoints)

    if (!is.null(indiv_col_lab)) {
      g <- g +
        geom_text_repel(aes(x = indiv_col[, x], y = indiv_col[, y],
                            label = indiv_col_lab),
                        col = col_color,
                        size = labsize,
                        max.overlaps = max.overlaps)
    }

  }

  # Plot row labels over points
  if (!is.null(indiv_row) & !is.null(indiv_row_lab)) {
    g <- g +
      geom_text_repel(aes(x = indiv_row[, x], y = indiv_row[, y],
                          label = indiv_row_lab),
                      col = row_color,
                      size = labsize,
                      max.overlaps = max.overlaps)
  }

  # Plot rows variables
  if (!is.null(var_row)) {
    g <- g +
      geom_segment(aes(x = 0, y = 0,
                       xend = var_row[, x]*mult, yend = var_row[, y]*mult),
                   arrow = arrow(length = grid::unit(0.20, "cm")))
  }

  # Plot columns variables
  if (!is.null(var_col)) {
    g <- g +
      geom_segment(aes(x = 0, y = 0,
                       xend = var_col[, x]*mult, yend = var_col[, y]*mult),
                   arrow = arrow(length = grid::unit(0.20,"cm"))) +
      geom_label(aes(x = var_col[, x]*mult, y = var_col[, y]*mult,
                     label = var_col_lab),
                 col = col_color,
                 vjust = ifelse(var_col[, y] > 0, 0, 1))
  }


  # Plot row variables over arrows
  if (!is.null(var_row)) {
    g <- g +
      geom_label(aes(x = var_row[, x]*mult, y = var_row[, y]*mult,
                     label = var_row_lab),
                 col = row_color,
                 vjust = ifelse(var_row[, y] > 0, 0, 1))
  }

  g
}

#' Plot correlation circle
#'
#' Plot the correlation circle
#'
#' @param cor dataframe containing the coordinates of the variables for each axis in columns
#' @param xax index of the column of the matrix to use for x
#' @param yax index of the column of the matrix to use for y
#' @param xlab Custom x-label
#' @param ylab Custom y-label
#' @param eig Eigenvalues vector
#' @param label label to display on the arrows
#' @param mar margin between arrow tips and label
#' @param col_bg color to use for background elements (circle and y = 0 and x = 0)
#' @param col color for arrows (unique)
#' @param lim limits for x and y: a unique number that gives the distance detween the plot
#' limit and zero (so that the circle is centered)
#' @param clip constrain the labels to be inside the plot? If yes, change to "on".
#' @param cor2 optional second set of variables (for double constrained analyses)
#' @param col2 optional color for the second set of variables
#' @param label2 optional label to display on the arrows for the second set of variables
#' @param lty linetype for the first set of variables
#' @param lty2 linetype for the second set of variables
#' @param labsize Size for arrow labels.
#'
#' @return A ggplot with the correlation circle where variables are represented as arrows
#' inside the circle of radius 1.
#' @export
plot_corcircle <- function(cor,
                           cor2 = NULL,
                           xax = 1, yax = 2,
                           xlab = NULL,
                           ylab = NULL,
                           label = rownames(cor),
                           eig = NULL,
                           mar = 0.01,
                           labsize = 3.88,
                           col_bg = "grey30",
                           col = "black",
                           lty = "solid",
                           col2 = "black",
                           lty2 = "solid",
                           label2 = rownames(cor2),
                           lim = 1,
                           clip = "off") {

  # x and y labels
  if (is.null(xlab)) {
    xlab <- paste0("Axis ", xax)
    if (!is.null(eig)) {
      xlab <- paste0(xlab, " (", round(eig[xax]/sum(eig)*100, 1), " % variability)")
    }
  }

  if (is.null(ylab)) {
    ylab <- paste0("Axis ", yax)
    if (!is.null(eig)) {
      ylab <- paste0(ylab, " (", round(eig[yax]/sum(eig)*100, 1), " % variability)")
    }
  }

  # Create circle
  r <- 1
  angle <- seq(0, 2*pi, by = 0.01)
  circle <- data.frame(x = r*cos(angle),
                       y = r*sin(angle))

  # Extract variables
  varx <- cor[, xax]
  vary <- cor[, yax]

  # Plot
  g <- ggplot() +
    ggplot2::geom_path(data = circle, aes(x = x, y = y),
                       col = col_bg) +
    geom_vline(aes(xintercept = 0), col = col_bg) +
    geom_hline(aes(yintercept = 0), col = col_bg) +
    geom_segment(aes(xend = varx, yend = vary,
                     x = 0, y = 0),
                 linetype = lty,
                 arrow = arrow(length = grid::unit(10, "pt")),
                 col = col) +
    coord_fixed(clip = clip,
                ylim = c(-lim, lim), xlim = c(-lim, lim)) +
    ylab(ylab) +
    xlab(xlab) +
    theme_linedraw()

  if (!is.null(cor2)) { # If there are other arrows to plot
    varx2 <- cor2[, xax]
    vary2 <- cor2[, yax]

    g <- g +
      geom_segment(aes(xend = varx2, yend = vary2,
                       x = 0, y = 0),
                   linetype = lty2,
                   arrow = arrow(length = grid::unit(10, "pt")),
                   col = col2)
  }

  if (!is.null(label)) { # If label not empty, add labels
    g <- g +
      geom_label(aes(x = varx + mar, y = vary, label = label), hjust = 0,
                 col = col, size = labsize)
  }

  if (!is.null(cor2) & !is.null(label2)) { # Same for labels 2
    g <- g +
      geom_label(aes(x = varx2 + mar, y = vary2, label = label2), hjust = 0,
                 col = col2, size = labsize)
  }

  if (clip == "off") {
    # Replace panel.border with panel.background outline, since panel.border
    # is drawn on top of the plot and hence above the labels that are not clipped
    # See https://github.com/tidyverse/ggplot2/issues/2536
    g <- g +
      theme(panel.border = element_blank(),
            panel.background = element_rect(colour = "black"))
  }

  g
}