# Main functions ----------------------------------------------------------


#' Plot traits
#' 
#' Plot trait matrices
#'
#' @param traits traits dataframe (expected to have a `size` column)
#' @param order order to use for the species 
#' @param label species labels
#' @param img species image
#' @param fill fill color (use a unique color to ignore the size)
#' @param type either "resource" (row) or "consumer" (column)
#'
#' @returns a ggplot
#' @export
plot_traits <- function(traits, order, label, img, fill = "grey", type = "consumer") {
  
  # Get number of species
  n <- length(label)
  
  # Plot
  p <- ggplot()
  if (type == "consumer") {
    if (length(fill) == 1) {
      # Constant color -> no size
      p <- p +
        rphylopic::geom_phylopic(img = img,
                      aes(x = rank(order), y = 0), 
                      width = 0.9,
                      fill =  fill)
    } else {
      p <- p +
        rphylopic::geom_phylopic(img = img,
                      aes(x = rank(order), y = 0,
                          width = traits$size / max(traits$size)), 
                      fill =  fill)
    }
    p <- p +
      geom_label(aes(x = rank(order), y = -0.3, label = label), 
                 col = "black", fill = "beige", size = 5, alpha = 0.7, 
                 label.size = NA) +
      scale_width_continuous(range = c(0.4, 1)) +
      labs(tag = expression(bold(v[1]))) +
      xlim(0.5, n + 0.5)  + 
      ylim(-0.5, 0.4)
    plot_tag_position <- c(1.05, 0.5)
  } else if (type == "resource") {
    if (length(fill) == 1) {
      # Constant color -> no size
      p <- p +
        rphylopic::geom_phylopic(img = img,
                      aes(x = 0, y = rank(order)),
                      height = 0.9,
                      fill =  fill)
    } else {
      p <- p +
        rphylopic::geom_phylopic(img = img,
                      aes(x = 0, y = rank(order), 
                          height = traits$size / max(traits$size)), 
                      fill =  fill)
    }
    p <- p +
      geom_label(aes(y = rank(order), x = -0.3, label = label), 
                 col = "black", fill = "beige", size = 5, alpha = 0.7, 
                 label.size = NA) +
      scale_height_continuous(range = c(0.4, 0.9)) +
      labs(tag = expression(bold(u[1]))) +
      xlim(-0.5, 0.4) +
      ylim(0.5, n + 0.5)
    plot_tag_position <- c(0.5, 1.05)
  }
  
  p <- p +
    theme_void() + 
    theme(plot.tag.position = plot_tag_position, 
          plot.tag = element_text(size = 20), 
          legend.position = "none", 
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      linewidth = 3))
  
  p
}


#' Plot traits data
#' 
#' Plot trait data matrices
#'
#' @param traits traits dataframe (expected to have a `size` column)
#' @param label species labels
#' @param img species image
#' @param fill fill color (use a unique color to ignore the size)
#' @param type either "resource" (row) or "consumer" (column)
#' @param reorder should species be reordered by label?
#'
#' @returns a ggplot
#' @export
plot_traits_data <- function(traits, label, img, fill, type = "consumer", reorder = FALSE) {
  # Get number of species
  n <- length(label)
  
  if (reorder) {
    position <- 1:n
  } else {
    position <- label
  }
  # Plot
  p <- ggplot()
  if (type == "consumer") {
    p <- p +
      rphylopic::geom_phylopic(img = img,
                    aes(x = rep(rank(position), 2), 
                        y = rep(c(0, 1), c(n, n)), 
                        width =  c(rep(0.5, n), traits$size / max(traits$size))), 
                    fill =  c(fill, rep("#000000", n))) +
      geom_label(aes(x = rank(position), y = 0.5, label = label), 
                 col = "black", fill = "beige", size = 5, alpha = 0.7, 
                 label.size = NA) +
      scale_width_continuous(range = c(0.4, 1)) +
      labs(tag = expression(bold(Q))) +
      xlim(0.5, n + 0.5)  + 
      ylim(-0.5, 1.5)
    plot_tag_position <- c(1.05, 0.5)
  } else if (type == "resource") {
    p <- p +
      rphylopic::geom_phylopic(img = img,
                    aes(x = rep(c(0, 1), c(n, n)),
                        y = rep(rank(position), 2),
                        height = c(rep(0.5, n), traits$size / max(traits$size))), 
                    fill =  c(fill, rep("#000000", n))) +
      geom_label(aes(y = rank(position), x = 0.5, label = label), 
                 col = "black", fill = "beige", size = 5, alpha = 0.7, 
                 label.size = NA) +
      scale_height_continuous(range = c(0.4, 0.9)) +
      labs(tag = expression(bold(R))) +
      xlim(-0.5, 1.5) +
      ylim(0.5, n + 0.5)
    plot_tag_position <- c(0.5, 1.05)
  }
  
  p <- p +
    theme_void() + 
    theme(plot.tag.position = plot_tag_position, 
          plot.tag = element_text(size = 20), 
          legend.position = "none", 
          panel.border = element_rect(colour = "black", 
                                      fill = NA, 
                                      linewidth = 3))
  
  p
}


#' Plot table
#' 
#' Plot a matrix
#'
#' @param x consumers abundance vector
#' @param y resources abundance vector
#' @param WA display the weighted averaging for one species?
#' @param label_consumer label of consumer species (only useful if WA is TRUE)
#' @param img_consumer image
#'  of consumer species (only useful if WA is TRUE)
#' @param abundance_data abundance data table: expected to have columns
#' 
#' @param abund_table table corresponding to abundance_data (to wide format = interaction matrix)
#' Consumer_sp, Resource_sp, Abundance and coa_c and coa_r of WA is TRUE.
#' @returns a ggplot
#' @export
plot_table <- function(x, y, 
                       abundance_data, abund_table,
                       label_consumer = NULL, 
                       img_consumer = NULL,
                       WA = FALSE){
  # Get resource and consumer species count
  n_consumer <- length(unique(abundance_data$Consumer_sp))
  n_resource <- length(unique(abundance_data$Resource_sp))
  
  p_data <- ggplot(abundance_data, aes(x = x, y = y , size = Abundance))
  
  plot_wa_consumer <- function(p_data, xbase, score, weights, img = img_consumer,
                               label = label_consumer){
    y0 <-  seq( 0.5,  6.5, length = 100)
    x0 <-  xbase - 4 * dnorm(y0, mean = stats::weighted.mean(score, w = weights), 
                             sd = sqrt(varfacwt(score, wt = weights)))
    df <- cbind(yn = c(y0, y0[length(y0)]), xn = c(x0, x0[1])) 
    df_mean <- data.frame(x1 = min(df[,2]),
                          x2 = max(df[,2]), 
                          y1 = stats::weighted.mean(score, w = weights), 
                          y2 = stats::weighted.mean(score, w = weights))
    p_data + geom_polygon(data = df, aes(x = xn, y = yn), 
                          inherit.aes = FALSE, 
                          fill = "burlywood", alpha = 0.4) +
      geom_segment(data = df_mean, 
                   aes(x = x1, xend = x2, y = y1, yend = y2), 
                   col = "burlywood",
                   inherit.aes = FALSE, lwd = 1) +
      geom_point(data = df_mean, 
                 aes(x = x2, y = y2), 
                 col = "burlywood",
                 size = 4) +
      rphylopic::geom_phylopic(img = img, 
                    aes(x = df_mean$x1[1] - 0.3, 
                        y = stats::weighted.mean(score, w = weights)), 
                    width = 0.4, height = NA, fill =  "burlywood") +
      annotate("label", 
               x = df_mean$x1[1]-0.9, 
               y = stats::weighted.mean(score, w = weights), 
               label = label[xbase], 
               col = "black", fill = "beige", size = 5,  alpha = 0.7)
  }
  
  p_data <- p_data + geom_point() + 
    xlim(0.5, n_consumer + 0.5)  +
    ylim(0.5, n_resource + 0.5)  +
    scale_size(range = c(1, 12)) + 
    theme_void() + 
    theme(legend.position = "none", 
          panel.border = element_rect(colour = "black", 
                                      fill= NA, linewidth = 3))
  
  if (WA) {
    highlight_points <- subset(abundance_data, Consumer_sp == "Consumer_H")
    p_data <- p_data + geom_point(data = highlight_points, 
                                  aes(x =  coa_c, y = coa_r, size = Abundance), 
                                  shape = 21, 
                                  color = "burlywood", 
                                  fill = "black",
                                  stroke = 2)
    p_data <- plot_wa_consumer(p_data = p_data, 
                               xbase = 8, 
                               score = 1:6, weights = abund_table[,8])
  }
  
  p_data
  
}


# Helper functions --------------------------------------------------------

#' Patch plot
#' 
#' Arrange plots with patchwork
#'
#' @param p_data data matrix plot (generated with `plot_data()`)
#' @param p_resource resource traits plot (generated with `plot_traits()`)
#' @param p_consumer consumer traits plot (generated with `plot_traits()`)
#' @param n_consumer number of resources for plot height
#' @param n_resource number of consumers for plot width
#' @param traits_dim width/height of the traits matrix (use 1 for ordination, 2 for data)
#'
#' @returns A patchwork ggplot
#' @export
patch_plot <- function(p_data, p_resource, p_consumer, 
                       n_consumer, n_resource, traits_dim = 1) {
  
  patch_1 <-   (p_consumer | plot_spacer() | plot_spacer()) + 
    plot_layout(widths = c(n_consumer, 0.2, traits_dim))
  patch_2 <- (p_data | plot_spacer() | p_resource) + 
    plot_layout(widths = c(n_consumer, 0.2,  traits_dim))
  patch <- patch_1 / plot_spacer() / patch_2
  patch <- patch + plot_layout(heights = c(traits_dim, 0.2, n_resource))  
  patch
}


#' Patch data plot
#' 
#' Arrange data plots with patchwork
#' 
#' @param p_data data matrix plot (generated with `plot_data()`)
#' @param p_resource resource traits plot (generated with `plot_traits_data()`)
#' @param p_consumer consumer traits plot (generated with `plot_traits_data()`)
#' 
#' @returns A patchwork ggplot
#' @export
# patch_plot_data <- function(p_data, p_resource, p_consumer) {
#   patch_1 <-   (p_consumer | plot_spacer() | plot_spacer()) + 
#     plot_layout(widths = c(n_consumer, 0.2, 2))
#   patch_2 <- (p_data | plot_spacer() | p_resource) + 
#     plot_layout(widths = c(n_consumer, 0.2,  2))
#   patch <- patch_1 / plot_spacer() / patch_2
#   patch <- patch  + plot_layout(heights = c(2, 0.2, n_resource))  
#   patch
# }