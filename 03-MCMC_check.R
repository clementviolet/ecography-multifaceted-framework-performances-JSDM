# Library

library(Hmsc)
library(tidyverse)
library(patchwork)


theme_set(theme_minimal())
theme_update(panel.grid.minor = element_blank())
# Load data

model_polychaeta_pa <- readRDS("Data/Model_Output/model_polychaeta_pa.rds")
model_polychaeta_ab <- readRDS("Data/Model_Output/model_polychaeta_ab.rds")
model_polychaeta_phylo_pa <- readRDS("Data/Model_Output/model_polychaeta_phylo_pa.rds")
model_polychaeta_phylo_ab <- readRDS("Data/Model_Output/model_polychaeta_phylo_ab.rds")
model_polychaeta_phylo_traits_pa <- readRDS("Data/Model_Output/model_polychaeta_phylo_traits_pa.rds")
model_polychaeta_phylo_traits_ab <- readRDS("Data/Model_Output/model_polychaeta_phylo_traits_ab.rds")
model_whole_fauna_pa <- readRDS("Data/Model_Output/model_whole_fauna_pa.rds")
model_whole_fauna_ab <- readRDS("Data/Model_Output/model_whole_fauna_ab.rds")

# Function

clean_names <- function(x) {
  x %>%
    str_remove_all("B\\[") %>%
    str_remove_all("G\\[") %>%
    str_remove_all("Omega[:digit:]?\\[") %>%
    str_remove_all("\\]") %>%
    str_remove_all("poly\\(") %>%
    str_replace_all(", degree = 2\\)", "_")
}

models_list <- list(
  "Polychaeta - AB" = model_polychaeta_ab,
  "Polychaeta - PA" = model_polychaeta_pa,
  "Polychaeta Phylogeny - AB" = model_polychaeta_phylo_ab,
  "Polychaeta Phylogeny - PA" = model_polychaeta_phylo_pa,
  "Polychaeta Phylogeny Traits - AB" = model_polychaeta_phylo_traits_ab,
  "Polychaeta Phylogeny Traits - PA" = model_polychaeta_phylo_traits_pa,
  "Whole community - AB" = model_whole_fauna_ab,
  "Whole community - PA" = model_whole_fauna_pa
)

posterior_list <- map(models_list, ~ convertToCodaObject(.x,
  spNamesNumbers = c(TRUE, FALSE),
  covNamesNumbers = c(TRUE, FALSE),
  trNamesNumbers = c(TRUE, FALSE)
))

# Beta Matrix (Environmnetal coefficients)

beta_conv <- map_dfr(posterior_list, function(x) {
  psrf <- gelman.diag(x$Beta, multivariate = FALSE, autoburnin = FALSE)$psrf[, 2]
  ess <- effectiveSize(x$Beta)

  data.frame(
    par = "beta",
    name_par = names(psrf),
    psrf = psrf,
    ess = ess,
    row.names = NULL
  ) %>%
    mutate(name_par_tmp = clean_names(name_par)) %>%
    separate(name_par_tmp, c("Var1", "Var2"), sep = "\\s?,\\s") %>%
    mutate(across(where(is.character), str_trim))
},
.id = "model"
)

beta_conv_plot <- beta_conv %>%
  group_by(model) %>%
  group_map(~., .keep = TRUE) %>%
  map(function(x, ...) {
    data_type <- unique(x$model) %>%
      str_replace("AB", "Abundance") %>%
      str_replace("PA", "Occurence")

    p1 <- ggplot(x, aes(ess)) +
      geom_density() +
      labs(x = "ESS", y = "") +
      scale_x_log10() +
      coord_cartesian(
        expand = FALSE,
        clip = "off"
      ) +
      theme_minimal()

    p2 <- ggplot(x, aes(psrf)) +
      geom_density() +
      geom_segment(aes(
        x = 1.2, xend = 1.2,
        y = 0, yend = Inf
      ),
      linetype = "dashed",
      color = "red"
      ) +
      labs(x = "PSRF", y = "Density") +
      coord_cartesian(
        expand = FALSE,
        clip = "off"
      ) +
      theme_minimal()

    p_final <- p2 + p1 #+ plot_annotation(title = glue::glue("{data_type} {str_to_sentence(unique(x$par))} parameters"))

    return(p_final)
  })

## Rho (Phylogeny coefficient)

rho_conv <- posterior_list[c(
  "Polychaeta Phylogeny - AB", # Keep only models with phylogeny evaluated
  "Polychaeta Phylogeny - PA",
  "Polychaeta Phylogeny Traits - AB",
  "Polychaeta Phylogeny Traits - PA"
)] %>%
  map_dfr(function(x) {
    data.frame(
      par = "rho",
      name_par = "rho",
      Var1 = "rho",
      Var2 = "rho",
      psrf = gelman.diag(x$Rho, autoburnin = FALSE)$psrf[, 2],
      ess = effectiveSize(x$Rho),
      row.names = NULL
    )
  }, .id = "model")

# Gamma (Traits coefficients)

gamma_conv <- map_dfr(posterior_list, function(x) {
  ess.gamma <- effectiveSize(x$Gamma)
  psrf.gamma <- gelman.diag(x$Gamma, multivariate = FALSE, autoburnin = FALSE)$psrf[, 2]

  data.frame(
    par = "gamma",
    name_par = names(psrf.gamma),
    psrf = psrf.gamma,
    ess = ess.gamma,
    row.names = NULL
  ) %>%
    mutate(name_par_tmp = clean_names(name_par)) %>%
    separate(name_par_tmp, c("Var1", "Var2"), sep = "\\s?,\\s") %>%
    mutate(across(where(is.character), str_trim))
}, .id = "model")

gamma_conv_plot <- gamma_conv %>%
  group_by(model) %>%
  group_map(~., .keep = TRUE) %>%
  map(function(x, ...) {
    data_type <- unique(x$model) %>%
      str_replace("AB", "Abundance") %>%
      str_replace("PA", "Occurence")

    p1 <- ggplot(x, aes(ess)) +
      geom_density() +
      labs(x = "ESS", y = "") +
      scale_x_log10() +
      coord_cartesian(
        expand = FALSE,
        clip = "off"
      ) +
      theme_minimal()

    p2 <- ggplot(x, aes(psrf)) +
      geom_density() +
      geom_segment(aes(
        x = 1.2, xend = 1.2,
        y = 0, yend = Inf
      ),
      linetype = "dashed",
      color = "red"
      ) +
      labs(x = "PSRF", y = "Density") +
      coord_cartesian(
        expand = FALSE,
        clip = "off"
      ) +
      theme_minimal()

    p_final <- p2 + p1 #+ plot_annotation(title = glue::glue("{data_type} {str_to_sentence(unique(x$par))} parameters"))

    return(p_final)
  })

# Save object
save(beta_conv, beta_conv_plot, gamma_conv, gamma_conv_plot, rho_conv, file = "convergence_diag.RData")

for (i in 1:length(models_list)) {
  ggsave(
    plot = beta_conv_plot[[i]],
    filename = paste0(
      "fig_supp_conv_beta_",
      str_remove_all(str_remove(names(models_list[i]), " - "), "\\s"),
      ".png"
    ),
    width = 12,
    height = 9,
    device = ragg::agg_png,
    scale = grDevices::dev.size(units = "cm")[1] / 9,
    units = "cm"
  )

  ggsave(
    plot = gamma_conv_plot[[i]],
    filename = paste0(
      "fig_supp_conv_gamma_",
      str_remove_all(str_remove(names(models_list[i]), " - "), "\\s"),
      ".png"
    ),
    width = 12,
    height = 9,
    device = ragg::agg_png,
    scale = grDevices::dev.size(units = "cm")[1] / 9,
    units = "cm"
  )
}

rho_conv_pretty <- rho_conv %>%
  mutate(Model = str_replace(str_replace(model, "AB", "Abundance"), "PA", "Occurence"), .before = par) %>%
  select(-model, -name_par, -Var1, -Var2) %>%
  rename(Parameter = par)

writeLines(
  knitr::kable(rho_conv_pretty),
  "rho_conv_table.txt"
)
