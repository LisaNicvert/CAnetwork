# Libraries ---------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(rphylopic)

library(dplyr)
library(tidyr)

library(ade4)

library(here)

devtools::load_all()

figures_path <- here("figures", "Fig_1_summary")

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

abundance_data_raw <- bind_rows(data_list) # |> filter(Abundance > 0) # Remove zero abundances for cleaner plot

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

# Get the ordering of resource/consumers according to principal coordinates
# (first axis)
res_resource <- data.frame(Resource_sp = names_resource, 
                           coa_r = rank(coa_Y$li[,1]), # WA scores
                           ccaR_r = rank(cca_R$ls[,1]), # LC scores
                           ccaQ_r = rank(cca_Q$co[,1]), # WA scores (permuted table -> resources in co)
                           dcca_r = rank(dcca_RQ$lsR[,1]), # LC scores
                           perm_r = perm_resource)

res_consumer <- data.frame(Consumer_sp = names_consumer, 
                           coa_c = rank(coa_Y$co[,1]), # WA scores
                           ccaR_c = rank(cca_R$co[,1]), # WA scores
                           ccaQ_c = rank(cca_Q$ls[,1]), # LC scores (permuted table -> consumers in ls)
                           dcca_c = rank(dcca_RQ$lsQ[,1]), # LC scores
                           perm_c = perm_consumer) 

## build a big table with all results
abundance_data <- abundance_data_raw |>
  left_join(res_consumer, by = "Consumer_sp") |> 
  left_join(res_resource, by = "Resource_sp") |> 
  filter(Abundance > 0)


# Plot --------------------------------------------------------------------

## Plot the original data ------
g1 <- plot_table(x = as.numeric(as.factor(abundance_data$Consumer_sp)), 
                 y = as.numeric(as.factor(abundance_data$Resource_sp)),
                 abundance_data = abundance_data, abund_table = abund_table)
g2 <- plot_traits_data(traits = traits_resource, label = perm_resource,
                       img = img_resource, fill = col_resource,
                       type = "resource", reorder = TRUE)
g3 <- plot_traits_data(traits = traits_consumer, label = perm_consumer,
                       img = img_consumer, fill = col_consumer,
                       type = "consumer", reorder = TRUE)
patch_plot(g1, g2, g3, 
           n_consumer = n_consumer, n_resource = n_resource,
           traits_dim = 2)
ggsave(file.path(figures_path, "plot_data0.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot the permuted original data ------
g1 <- plot_table(x = abundance_data$perm_c, 
                 y = abundance_data$perm_r,
                 abundance_data = abundance_data, abund_table = abund_table)
g2 <- plot_traits_data(traits = traits_resource, label = perm_resource,
                       img = img_resource, fill = col_resource,
                       type = "resource")
g3 <- plot_traits_data(traits = traits_consumer, label = perm_consumer,
                       img = img_consumer, fill = col_consumer,
                       type = "consumer")

g_data <- patch_plot(g1, g2, g3, 
                     n_consumer = n_consumer, n_resource = n_resource,
                     traits_dim = 2) + 
  plot_annotation(
    caption = bquote("Original data tables (" * bold(R) * ", " * bold(L) * " and " * bold(Q) * ")"),
    theme = theme(
        plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_data    
ggsave(file.path(figures_path, "plot_data.pdf"), 
       width = 140, height = 140, unit = "mm")


## Plot COA results ------
g1 <- plot_table(x = abundance_data$coa_c, 
                 y = abundance_data$coa_r,
                 abundance_data = abundance_data, abund_table = abund_table,
                 label_consumer = perm_consumer,
                 img_consumer = img_consumer,
                 WA = TRUE)
g2 <- plot_traits(traits = traits_resource, order =  unique(abundance_data$coa_r),
                  label = perm_resource, 
                  img = img_resource, type = "resource")
g3 <- plot_traits(traits = traits_consumer, order =  unique(abundance_data$coa_c),
                  label = perm_consumer, 
                  img = img_consumer, type = "consumer")
g_coa <- patch_plot(g1, g2, g3, 
                    n_consumer = n_consumer, n_resource = n_resource) +  
    plot_annotation(
        caption = bquote("CA of " * bold(L) * " (" * delta[1] * " = " * .(round(coa_Y$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_coa
ggsave(file.path(figures_path, "plot_coa.pdf"), 
       width = 140, height = 140, unit = "mm")


## Plot CCA-R results ------
g1 <- plot_table(x = abundance_data$ccaR_c, 
                 y = abundance_data$ccaR_r,
                 abundance_data = abundance_data, abund_table = abund_table) 
g2 <- plot_traits(traits = traits_resource, 
                  order = unique(abundance_data$ccaR_r),
                  label = perm_resource, 
                  img = img_resource, fill = col_resource, type = "resource")
g3 <- plot_traits(traits = traits_consumer, 
                  order = unique(abundance_data$ccaR_c),
                  label = perm_consumer,
                  img = img_consumer, type = "consumer")
g_ccaR <- patch_plot(g1, g2, g3, 
                     n_consumer = n_consumer, n_resource = n_resource) + 
    plot_annotation(
        caption = bquote("CCA constrained by " * bold(R) * " (" * delta[1] * " = " * .(round(cca_R$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_ccaR
ggsave(file.path(figures_path, "plot_ccaR.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot CCA-Q results ------
g1 <- plot_table(x = abundance_data$ccaQ_c, 
                 y = abundance_data$ccaQ_r,
                 abundance_data = abundance_data, abund_table = abund_table) 
g2 <- plot_traits(traits = traits_resource, 
                  order = unique(abundance_data$ccaQ_r),
                  label = perm_resource, 
                  img = img_resource, type = "resource")
g3 <- plot_traits(traits = traits_consumer, 
                  order = unique(abundance_data$ccaQ_c),
                  label = perm_consumer, fill = col_consumer,
                  img = img_consumer, type = "consumer")
g_ccaQ <- patch_plot(g1, g2, g3, 
                     n_consumer = n_consumer, n_resource = n_resource) + 
    plot_annotation(
        caption = bquote("CCA constrained by " * bold(Q) * " (" * delta[1] * " = " * .(round(cca_Q$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )

g_ccaQ
ggsave(file.path(figures_path, "plot_ccaQ.pdf"), 
       width = 140, height = 140, unit = "mm")

## Plot dc-CA results ------
g1 <- plot_table(x = abundance_data$dcca_c, 
                 y = abundance_data$dcca_r,
                 abundance_data = abundance_data, abund_table = abund_table)
g2 <- plot_traits(traits = traits_resource, 
                  order = unique(abundance_data$dcca_r),
                  label = perm_resource, fill = col_resource,
                  img = img_resource, type = "resource")
g3 <- plot_traits(traits = traits_consumer, 
                  order = unique(abundance_data$dcca_c),
                  label = perm_consumer, fill = col_consumer,
                  img = img_consumer, type = "consumer")
g_dcca <- patch_plot(g1, g2, g3, 
                     n_consumer = n_consumer, n_resource = n_resource) + 
    plot_annotation(
        caption = bquote("dc-CA constrained by " * bold(R) * " and " * bold(Q) * " (" * delta[1] * " = " * .(round(dcca_RQ$eig[1], 2)) * ")"),
        theme = theme(
            plot.caption = element_text(hjust = 0.5, vjust = 1, size = 20)
        )
    )


g_dcca
ggsave(file.path(figures_path, "plot_dcca.pdf"), 
       width = 140, height = 140, unit = "mm")
