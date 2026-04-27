# Libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(rphylopic)
library(tidyr)
library(ade4)
library(patchwork)
library(here)

figures_path <- here("figures", "Fig_1_summary")

# Custom functions --------------------------------------------------------
plot_table <- function(x = as.numeric(as.factor(abundance_data$Consumer_sp)), y = as.numeric(as.factor(abundance_data$Resource_sp)), WA = FALSE){
  p_data <- ggplot(abundance_data, aes(x = x, y = y , size = Abundance))
  # xlim(0.5, n_consumer + 0.5) + ylim(0.5, n_resource + 0.5)
  
  plot_wa_consumer <- function(p_data, xbase, score, weights, img = img_consumer){
    y0 <-  seq( 0.5,  6.5, length = 100)
    x0 <-  xbase - 4 * dnorm(y0, mean = weighted.mean(score, w = weights), sd = sqrt(varfacwt(score, wt = weights)))
    df <- cbind(yn = c(y0, y0[length(y0)]), xn = c(x0, x0[1])) 
    df_mean <- data.frame(x1 = min(df[,2]),  x2 = max(df[,2]), y1 = weighted.mean(score, w = weights), y2 = weighted.mean(score, w = weights))
    p_data + geom_polygon(data = df, aes(x = xn, y = yn), inherit.aes = FALSE, fill = "burlywood", alpha = 0.2) +
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), arrow = arrow(type = "closed", length = unit(0.2, "cm")), data = df_mean, inherit.aes = FALSE, lwd = 1, col = "burlywood", ) +
      geom_phylopic(img = img, aes(x = df_mean$x1[1] - 0.3, y = weighted.mean(score, w = weights)), width = 0.4, height = NA, fill =  "burlywood") +
      annotate("label", x = df_mean$x1[1]-0.9, y = weighted.mean(score, w = weights), label = perm_consumer[xbase], col = "black", fill = "beige", size = 5,  alpha = 0.7)
  }
  
  p_data <- p_data + geom_point() + 
    xlim(0.5, n_consumer + 0.5)  +
    ylim(0.5, n_resource + 0.5)  +
    scale_size(range= c(1,12)) + 
    theme_void() + 
    theme(legend.position = "none", 
          panel.border = element_rect(colour = "black", fill= NA, linewidth = 3))# A clean theme
  
  if(WA){
    highlight_points <- subset(abundance_data, Consumer_sp == "Consumer_H")
    p_data <- plot_wa_consumer(p_data, 8, 1:6, abund_table[,8])
    p_data <- p_data + geom_point(data = highlight_points, aes(x =  coa_c, y = coa_r, size = Abundance), 
                                  shape = 21, 
                                  color = "burlywood", 
                                  fill = "black",
                                  stroke = 2)
  }
  
  p_data
  
}

plot_consumer <- function(x = 1:n_consumer, img = img_consumer){
  p_consumer <- ggplot() +
    geom_phylopic(img = img,
                  aes(x = rank(x), y = 0, width =  traits_consumer$size / max(traits_consumer$size)), fill =  col_consumer) + 
    
    geom_label(aes(x = rank(x), y = -0.3, label = perm_consumer), col = "black", fill = "beige", size = 5, alpha = 0.7, 
               label.size = NA) +
    scale_width_continuous(range = c(0.4, 1)) + labs(tag = expression(bold(v[1]))) +
    xlim(0.5, n_consumer + 0.5)  + 
    ylim(-0.5, 0.4)+
    theme_void() + 
    theme(plot.tag.position = c(1.05, 0.5), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3))
  
  p_consumer
}

plot_resource <- function(y = 1:n_resource, img = img_resource){
  p_resource <- ggplot() + 
    geom_phylopic(img = img,
                  aes(x =0, y = rank(y), height = traits_resource$size / max(traits_resource$size)), fill =  col_resource) + 
    geom_label(aes(y = rank(y), x = -0.3, label = perm_resource), col = "black", fill = "beige", size = 5, alpha = 0.7, 
               label.size = NA) +
    scale_height_continuous(range = c(0.4, 0.9)) +
    ylim(0.5, n_resource + 0.5) + labs(tag = expression(bold(u[1]))) +
    xlim(-0.5, 0.4)+
    theme_void() + 
    theme(plot.tag.position = c(0.5, 1.05), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3)) 
  
  p_resource
  
}

plot_consumer_data <- function(x = 1:n_consumer, img = img_consumer){
  p_consumer <- ggplot() +
    geom_phylopic(img = img,
                  aes(x = rep(rank(x), 2), y = rep(c(0,1), c(n_consumer, n_consumer)), width =  c(rep(0.5, n_consumer), traits_consumer$size / max(traits_consumer$size))), fill =  c(col_consumer, rep("#000000", n_consumer))) + 
    
    geom_label(aes(x = rank(x), y = 0.5, label = perm_consumer), col = "black", fill = "beige", size = 5, alpha = 0.7, 
               label.size = NA) +
    scale_width_continuous(range = c(0.4, 1)) +labs(tag = expression(bold(Q))) +
    xlim(0.5, n_consumer + 0.5)  + 
    ylim(-0.5, 1.5)  +
    theme_void() + 
    theme(plot.tag.position = c(1.05, 0.5), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3))
  
  p_consumer
}

plot_resource_data <- function(y = 1:n_resource, img = img_resource){
  p_resource <- ggplot() + 
    geom_phylopic(img = img,
                  aes(x =rep(c(0,1), c(n_resource, n_resource)), y = rep(rank(y), 2), height = c(rep(0.5, n_resource), traits_resource$size / max(traits_resource$size))), fill =  c(col_resource, rep("#000000", n_resource))) + 
    geom_label(aes(y = rank(y), x = 0.5, label = perm_resource), col = "black", size = 5,  alpha = 0.7, 
               label.size = NA, fill = "beige") +
    scale_height_continuous(range = c(0.4, 0.9)) + labs(tag = expression(bold(R))) +
    ylim(0.5, n_resource + 0.5) + 
    xlim(-0.5, 1.5)  +
    theme_void() + 
    theme(plot.tag.position = c(0.5, 1.05), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3)) 
  
  p_resource
  
}


plot_consumer_grey <- function(x = 1:n_consumer, img = img_consumer){
  p_consumer <- ggplot() +
    geom_phylopic(img = img,
                  aes(x = rank(x), y = 0), width = 0.9, fill =  "grey") +
    geom_label(aes(x = rank(x), y = -0.3, label = perm_consumer), col = "black", fill = "beige", size = 5, alpha = 0.7, 
               label.size = NA) +
    
    scale_width_continuous(range = c(0.4, 1)) + labs(tag = expression(bold(v[1]))) +
    xlim(0.5, n_consumer + 0.5)  + 
    ylim(-0.5, 0.4)+
    theme_void() + 
    theme(plot.tag.position = c(1.05, 0.5), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3))
  
  p_consumer
}

plot_resource_grey <- function(y = 1:n_resource, img = img_resource){
  p_resource <- ggplot() + 
    geom_phylopic(img = img,
                  aes(x =0, y = rank(y)), height = 0.9, fill = "grey") +
    geom_label(aes(x = -0.3, y = rank(y), label = perm_resource), col = "black", fill = "beige", size = 5, alpha = 0.7, 
               label.size = NA)+ 
    scale_height_continuous(range = c(0.4, 0.9)) +
    ylim(0.5, n_resource + 0.5) + labs(tag = expression(bold(u[1]))) +
    xlim(-0.5, 0.4)+
    theme_void() + 
    theme(plot.tag.position = c(0.5, 1.05), plot.tag = element_text(size = 20), legend.position = "none", panel.border = element_rect(colour = "black", fill= NA, linewidth = 3)) 
  
  p_resource
  
}

patch_plot <- function(p_data, p_resource, p_consumer){
  patch_1 <-   (p_consumer | plot_spacer() |plot_spacer()) + plot_layout(widths = c(n_consumer, 0.2, 1))
  patch_1 <- patch_1
  patch_2 <- (p_data | plot_spacer() | p_resource) + plot_layout(widths = c(n_consumer, 0.2,  1))
  patch <- patch_1 / plot_spacer() / patch_2
  patch <- patch  + plot_layout(heights = c(1, 0.2, n_resource))  
  patch
  
}


patch_plot_data <- function(p_data, p_resource, p_consumer){
  patch_1 <-   (p_consumer | plot_spacer() |plot_spacer()) + plot_layout(widths = c(n_consumer, 0.2, 2))
  patch_2 <- (p_data | plot_spacer() | p_resource) + plot_layout(widths = c(n_consumer, 0.2,  2))
  patch <- patch_1 / plot_spacer() / patch_2
  patch <- patch  + plot_layout(heights = c(2, 0.2, n_resource))  
  patch
  
}


# Create example data -----------------------------------------------------
n_consumer <- 8
n_resource <- 6

set.seed(69) # for reproducibility

names_consumer <- paste0("Consumer_", LETTERS[1:n_consumer])
names_resource <- paste0("Resource_", LETTERS[1:n_resource])

optim_resource <- seq(12, 38, length = n_resource)
optim_consumer <- seq(15, 35, length.out = n_consumer)
sds_consumer <- runif(length(names_consumer), min = 2, max = 4) # Varying niche breadth

data_list <- list()

# Simulate abundance for each consumer along the gradient of resource
# Distribute optimums across the gradient and vary niche breadth (sd)

for (i in 1:length(names_consumer)) {
  data_list[[i]] <- data.frame(
    Consumer_sp = names_consumer[i],
    Resource_sp = names_resource,
    Abundance = round(dnorm(optim_resource, mean = optim_consumer[i], sd = sds_consumer[i]) * sample(500:1500, 1)),
    Optim_resource = optim_resource,
    Optim_consumer = optim_consumer[i]
    
  )
}

abundance_data_raw <- bind_rows(data_list) # %>% filter(Abundance > 0) # Remove zero abundances for cleaner plot

abund_table <- as.data.frame(pivot_wider(abundance_data_raw[,1:3], 
                                         names_from = Consumer_sp,
                                         values_from = Abundance, id_cols = Resource_sp))
rownames(abund_table) <- abund_table[,1]
abund_table <- abund_table[,-1]

## generate traits
noisy_trait <- function(trait, f = 4) sapply(trait, FUN = function(x) x + runif(1, min = -x / f, max = x / f))
traits_consumer <- data.frame(size = noisy_trait(optim_consumer, f= 0.8), col = noisy_trait(optim_consumer, f = 1))
traits_resource <- data.frame(size = noisy_trait(optim_resource, f = 0.8), col = noisy_trait(optim_resource, f = 0.2))

col_consumer <- colorRampPalette(c("beige", "orange", "maroon"))(n_consumer + 1)[-1]
col_resource <- colorRampPalette(c("pink", "darkmagenta", "brown2"))(n_resource + 1)[-1]

##  add permutations of rows and columns to plot un-ordered table
perm_resource <- sample(n_resource)
perm_consumer <- sample(n_consumer)


# Images for the figures --------------------------------------------------
id_consumer <- "91ffc54e-8a80-498a-8c6a-a93e0cd1339e"
id_resource <- "d59b3e54-faa6-4c5c-9450-0c9d20db7cbd"
img_consumer <- get_phylopic(uuid = id_consumer)
img_resource <- get_phylopic(uuid = id_resource)  
rphylopic::get_attribution(c(id_consumer, id_resource), 
                           text = TRUE)


# Multivariate analysis ---------------------------------------------------
coa_Y <- dudi.coa(abund_table, scannf = FALSE) 
cca_R <- pcaiv(coa_Y, traits_resource, scannf = FALSE) 
cca_Q <- pcaiv(t(coa_Y), traits_consumer, scannf = FALSE) 
dcca_RQ <- dpcaiv(coa_Y, traits_resource, traits_consumer, scannf = FALSE)

res_resource <- data.frame(Resource_sp = names_resource, coa_r = rank(coa_Y$l1[,1]),
                           ccaR_r = rank(cca_R$l1[,1]), ccaQ_r = rank(cca_Q$c1[,1]), 
                           dcca_r = rank(dcca_RQ$l1[,1]), perm_r = perm_resource) 

res_consumer <- data.frame(Consumer_sp = names_consumer, coa_c = rank(coa_Y$c1[,1]),
                           ccaR_c = rank(cca_R$c1[,1]), ccaQ_c = rank(cca_Q$l1[,1]), 
                           dcca_c = rank(dcca_RQ$c1[,1]), perm_c = perm_consumer) 

## build a big table with all results
abundance_data <- abundance_data_raw %>%
  left_join(res_consumer, by = "Consumer_sp") %>% left_join(res_resource, by = "Resource_sp") %>% filter(Abundance > 0)


# Plot --------------------------------------------------------------------

## Plot the original data
g1 <- plot_table() 
g2 <- plot_resource_data()
g3 <- plot_consumer_data()
patch_plot_data(g1, g2, g3)
ggsave(file.path(figures_path, "plot_data0.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot the permuted original data
g1 <- plot_table(x = abundance_data$perm_c, y = abundance_data$perm_r) 
g2 <- plot_resource_data(y = perm_resource)
g3 <- plot_consumer_data(x = perm_consumer)
g_data <- patch_plot_data(g1, g2, g3) + plot_annotation(
    caption = bquote("Original data tables (" * bold(R) * ", " * bold(L) * " and " * bold(Q) * ")"),
    theme = theme(
        plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
    )
)

g_data    
ggsave(file.path(figures_path, "plot_data.pdf"), 
       width = 140, height = 140, unit = "mm")


## Plot COA results
g1 <- plot_table(x = abundance_data$coa_c, y = abundance_data$coa_r, WA = TRUE) 
g2 <- plot_resource_grey(y = unique(abundance_data$coa_r))
g3 <- plot_consumer_grey(x = unique(abundance_data$coa_c))
g_coa <- patch_plot(g1, g2, g3) +  
    plot_annotation(
        caption = bquote("CA of " * bold(L) * " (" * delta[1] * " = " * .(round(coa_Y$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_coa
ggsave(file.path(figures_path, "plot_coa.pdf"), 
       width = 140, height = 140, unit = "mm")


## Plot CCA-R results
g1 <- plot_table(x = abundance_data$ccaR_c, y = abundance_data$ccaR_r) 
g2 <- plot_resource(y = unique(abundance_data$ccaR_r))
g3 <- plot_consumer_grey(x = unique(abundance_data$ccaR_c))
g_ccaR <- patch_plot(g1, g2, g3) + 
    plot_annotation(
        caption = bquote("CCA constrained by " * bold(R) * " (" * delta[1] * " = " * .(round(cca_R$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_ccaR
ggsave(file.path(figures_path, "plot_ccaR.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot CCA-Q results
g1 <- plot_table(x = abundance_data$ccaQ_c, y = abundance_data$ccaQ_r) 
g2 <- plot_resource_grey(y = unique(abundance_data$ccaQ_r))
g3 <- plot_consumer(x = unique(abundance_data$ccaQ_c))
g_ccaQ <- patch_plot(g1, g2, g3) + 
    plot_annotation(
        caption = bquote("CCA constrained by " * bold(Q) * " (" * delta[1] * " = " * .(round(cca_Q$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_ccaQ
ggsave(file.path(figures_path, "plot_ccaQ.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot dc-CA results

g1 <- plot_table(x = abundance_data$dcca_c, y = abundance_data$dcca_r) 
g2 <- plot_resource(y = unique(abundance_data$dcca_r))
g3 <- plot_consumer(x = unique(abundance_data$dcca_c))
g_dcca <- patch_plot(g1, g2, g3) + 
    plot_annotation(
        caption = bquote("dc-CA constrained by " * bold(R) * " and " * bold(Q) * " (" * delta[1] * " = " * .(round(dcca_RQ$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )


g_dcca
ggsave(file.path(figures_path, "plot_dcca.pdf"), 
       width = 140, height = 140, unit = "mm")
