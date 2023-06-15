############## Library

library(tidyverse)
library(magick)
library(patchwork)
library(ggforce)
library(ggdist)
library(Hmsc)
library(FSA)
library(vegan)

library(sf)

theme_set(theme_minimal())
theme_update(panel.grid.minor = element_blank())

############## Data

# Load train models
load("Data/Processed/train_model_polychaeta.RData")
load("Data/Processed/train_model_whole_fauna.RData")

model_polychaeta_ab <- readRDS("Data/Model_Output/model_polychaeta_ab.rds")
model_polychaeta_pa <- readRDS("Data/Model_Output/model_polychaeta_pa.rds")

model_polychaeta_phylo_ab <- readRDS("Data/Model_Output/model_polychaeta_phylo_ab.rds")
model_polychaeta_phylo_pa <- readRDS("Data/Model_Output/model_polychaeta_phylo_pa.rds")

model_polychaeta_phylo_traits_ab <- readRDS("Data/Model_Output/model_polychaeta_phylo_traits_ab.rds")
model_polychaeta_phylo_traits_pa <- readRDS("Data/Model_Output/model_polychaeta_phylo_traits_pa.rds")

model_whole_fauna_ab <- readRDS("Data/Model_Output/model_whole_fauna_ab.rds")
model_whole_fauna_pa <- readRDS("Data/Model_Output/model_whole_fauna_pa.rds")

models <- list(
  "Benchmark - AB" = model_polychaeta_ab,
  "Benchmark - PA" = model_polychaeta_pa,
  "Phylogeny - AB" = model_polychaeta_phylo_ab,
  "Phylogeny - PA" = model_polychaeta_phylo_pa,
  "Phylogeny Traits - AB" = model_polychaeta_phylo_traits_ab,
  "Phylogeny Traits - PA" = model_polychaeta_phylo_traits_pa,
  "Whole community - AB" = model_whole_fauna_ab,
  "Whole community - PA" = model_whole_fauna_pa
)

models_acc <- c(
  "Bench", "Ph",
  "TrPh", "WhC"
)

load("./Data/Raw/trait_poly_estim.Rdata")
trait_list <- read_csv2("./Data/Raw/Trait_list.csv")

# Get names of polychaeta species

polychaeta <- models[[1]]$spNames

# Clean the name of covariables

covNames <- models[[1]]$covNames %>%
  str_replace("\\, degree \\= 2\\)", "_") %>%
  str_remove("poly\\(") %>%
  str_remove_all("\\(|\\)") %>%
  str_remove_all("\\`") %>%
  str_replace("Average", "Fetch") %>%
  str_replace("MO", "Organic Mater") %>%
  str_replace("SAL_sd_w", "Salinity") %>%
  str_replace("TEMP_sd_w", "Temperature") %>%
  str_replace("CURR_mean_w", "Current") %>%
  str_replace("mud", "Mud") %>%
  str_replace("TraskSo", "Trask Index")

# Clean the name of the traits
TrNames <- model_polychaeta_phylo_traits_ab$trNames %>%
  str_remove_all("\\(|\\)")

load("Data/Processed/test_model_polychaeta.RData")
load("Data/Processed/test_model_whole_fauna.RData")

# Load classification data frame
load("Data/Processed/sp_classif.RData")

# Load coordinates of the sampled sites

load("Data/Processed/sites_coordinates.RData")

############## Functions

computeAUC <- function(Y_obs, Y_pred, expected = TRUE) {
      if (!expected) {
        for (i in 1:dim(Y_pred)[3]) {
          Y_pred[, , i] <- vegan::decostand(Y_pred[, , i], method = "pa") # Convert Poisson parameter predictions to P/A
        }
      }

      Y_pred <- apply(Y_pred, MARGIN = c(1, 2), mean, na.rm = TRUE)

      # Vector of output
      AUC <- rep(NA, ncol(Y_obs))

      # Convert obs data to presence/absence
      Y_obs <- vegan::decostand(Y_obs, method = "pa")

      # Compute AUC
      for (i in 1:ncol(Y_obs)) {

        # sel = !is.na(Y[,i])

        if (length(unique(Y_obs[, i])) == 2) { # Make sure that we have at least one presence and one absence for the species in the test set

          AUC[i] <- pROC::auc(Y_obs[, i], Y_pred[, i], levels = c(0, 1), direction = "<")
        }
      }
      return(AUC)
}

computeTjurR2 <- function(Y_obs, Y_pred, expected = TRUE) {
      if (!expected) {
        for (i in 1:dim(Y_pred)[3]) {
          Y_pred[, , i] <- vegan::decostand(Y_pred[, , i], method = "pa") # Convert Poisson parameter predictions to P/A
        }
      }

      Y_pred <- apply(Y_pred, MARGIN = c(1, 2), mean, na.rm = TRUE)

      # Vector of output
      R2 <- rep(NA, ncol(Y_obs))

      # Convert obs data to presence/absence
      Y_obs <- vegan::decostand(Y_obs, method = "pa")

      # Compute TjurR2
      for (i in 1:ncol(Y_obs)) {
        R2[i] <- mean(Y_pred[which(Y_obs[, i] == 1), i]) - mean(Y_pred[which(Y_obs[, i] == 0), i])
      }

      return(R2)
}

computeR2 <- function(Y_obs, Y_pred, expected = TRUE, method = "spearman") {

      # If it is a Normal/Probit model, use the mean of the prediction, if it is a Poisson model, use the median
      if (expected) {
        Y_pred <- apply(Y_pred, MARGIN = c(1, 2), mean, na.rm = TRUE)
      } else {
        Y_pred <- apply(Y_pred, MARGIN = c(1, 2), median, na.rm = TRUE)
      }

      # Vector of output
      R2 <- rep(NA, ncol(Y_obs))

      # Convert obs data to presence/absence
      Y_obs <- vegan::decostand(Y_obs, method = "pa")

      # Compute R2 or Pseudo-R2

      for (i in 1:ncol(Y_obs)) {
        co <- cor(Y_obs[, i], Y_pred[, i], method = method, use = "pairwise")

        R2[i] <- sign(co) * co^2
      }

      return(R2)
}

computeRMSE <- function(Y_obs, Y_pred, expected = TRUE) {

      # If it is a Normal/Probit model, use the mean of the prediction, if it is a Poisson model, use the median
      if (expected) {
        Y_pred <- apply(Y_pred, MARGIN = c(1, 2), mean, na.rm = TRUE)
      } else {
        Y_pred <- apply(Y_pred, MARGIN = c(1, 2), median, na.rm = TRUE)
      }

      # Vector of output
      RMSE <- rep(NA, ncol(Y_obs))

      # Convert obs data to presence/absence
      Y_obs <- vegan::decostand(Y_obs, method = "pa")

      # Compute R2 or Pseudo-R2

      for (i in 1:ncol(Y_obs)) {
        RMSE[i] <- sqrt(mean((Y_obs[, i] - Y_pred[, i])^2, na.rm = TRUE))
      }
      return(RMSE)
}

polynom_classif <- function(x) {
  if (x$Quadratic == 0) {
    if (sign(x$Linear) == 1) {
      out <- c("Constant increase")
    } else if (sign(x$Linear) == -1) {
      out <- c("Constant decrease")
    } else {
      out <- c("Stable")
    }
  } else if (sign(x$deriv_lower) != sign(x$deriv_upper)) {
    if (sign(x$Quadratic) == 1) {
      out <- c("Convex")
    } else {
      out <- c("Concave")
    }
  } else {
    if (sign(x$deriv_upper) == -1 & sign(x$Quadratic * x$curvature) == -1) {
      out <- c("Accelerate decline")
    } else if (sign(x$deriv_upper) == -1 & sign(x$Quadratic * x$curvature) == 1) {
      out <- c("Decelerated decline")
    } else if (sign(x$deriv_upper) == 1 & sign(x$Quadratic * x$curvature) == -1) {
      out <- c("Accelerated increase")
    } else {
      out <- c("Decelerated increase")
    }
  }

  return(out)
}

############## Fig 2

#### Raw predictions

raw_prediction_trainset <- models %>%
  imap(function(x, y) {
    out <- x %>%
      computePredictedValues(expected = FALSE)

    if (str_detect(y, "PA")) {
      out <- out %>%
        apply(MARGIN = c(1, 2), mean, na.rm = TRUE) # If mean pred > 0.5 I consider a presence for Occurence based models
    } else {
      out <- out %>%
        apply(MARGIN = c(1, 2), median, na.rm = TRUE)
    }
    return(out)
  })

raw_prediction_testset <- models %>%
  imap(function(x, y) {
    if (str_detect(y, "PA")) expected <- TRUE else expected <- FALSE

    # Structure of the study design and random factors for predictions

    time_test <- as.matrix(unique(as.numeric(as.character(test_env$annee))))
    rownames(time_test) <- levels(test_env$annee)
    colnames(time_test) <- "annee"

    samplingDesign_test <- data.frame(
      annee = test_env$annee,
      site = as.factor(test_env$site),
      habitat = as.factor(test_env$habitat)
    )

    ## Random levels

    rL1_test <- HmscRandomLevel(sData = time_test)
    rL2_test <- HmscRandomLevel(units = unique(samplingDesign_test$site))
    rL3_test <- HmscRandomLevel(units = unique(samplingDesign_test$habitat))


    if (!str_detect(y, "Whole Community")) {
      test <- test_polychaeta

      out <- Hmsc:::predict.Hmsc(x,
        XData = as.data.frame(test_env[, -(1:3)]),
        studyDesign = samplingDesign_test,
        ranLevels = list(
          annee = rL1_test,
          site = rL2_test,
          habitat = rL3_test
        ),
        expected = expected
      ) %>%
        {
          array(unlist(.),
            dim = c(nrow(.[[1]]), ncol(.[[1]]), length(.)),
            dimnames = list(NULL, colnames(.[[1]]), NULL)
          )
        }
    } else {
      test <- test_whole_fauna

      out <- Hmsc:::predict.Hmsc(x,
        XData = as.data.frame(test_env[, -(1:3)]),
        studyDesign = samplingDesign_test,
        ranLevels = list(
          annee = rL1_test,
          site = rL2_test,
          habitat = rL3_test
        ),
        expected = expected
      ) %>%
        {
          array(unlist(.),
            dim = c(nrow(.[[1]]), ncol(.[[1]]), length(.)),
            dimnames = list(NULL, colnames(.[[1]]), NULL)
          )
        }
    }

    if (expected) {
      out <- out %>%
        apply(MARGIN = c(1, 2), mean, na.rm = TRUE)
      # For probit and normal models, mean of
      # the posterior is use, see
      # ?Hmsc::evaluateModelFit()
    } else {
      out <- out %>%
        apply(MARGIN = c(1, 2), median, na.rm = TRUE)
      # For probit and normal models, median of
      # the posterior is use, see
      # ?Hmsc::evaluateModelFit()
    }

    return(out)
  })

### Computation of the explanatory power

prediction_trainset <- models %>%
  imap_dfr(function(x, y) {
    modelfit <- function(x, expected = TRUE) {
      sp <- x$spNames

      computePredictedValues(x, expected = expected) %>%
        {
          evaluateModelFit(x, .)
        } %>%
        map_dfr(function(x) {
          data.frame(
            species = sp,
            value = x
          )
        }, .id = "measure")
    }

    if (str_detect(y, "PA")) {
      out <- x %>%
        modelfit()
    } else {
      out <- x %>%
        modelfit(expected = FALSE)
    }

    return(out)
  }, .id = "model") %>%
  filter(
    species %in% polychaeta,
    measure %in% c("SR2", "RMSE", "O.AUC", "AUC", "O.TjurR2", "TjurR2")
  )

prediction_testset <- models %>%
  imap_dfr(function(x, y) {

    if (str_detect(y, "PA")) expected <- TRUE else expected <- FALSE # If it's a Poisson model, expected need to be set to FALSE, see ?Hmsc::evaluateModelFit()

    # Structure of the study design and random factors for predictions

    time_test <- as.matrix(unique(as.numeric(as.character(test_env$annee))))
    rownames(time_test) <- levels(test_env$annee)
    colnames(time_test) <- "annee"

    samplingDesign_test <- data.frame(
      annee = test_env$annee, site = as.factor(test_env$site),
      habitat = as.factor(test_env$habitat)
    )

    ## Random levels

    rL1_test <- HmscRandomLevel(sData = time_test)
    rL2_test <- HmscRandomLevel(units = unique(samplingDesign_test$site))
    rL3_test <- HmscRandomLevel(units = unique(samplingDesign_test$habitat))

    if (!str_detect(y, "Whole community")) {
      test <- test_polychaeta

      out <- Hmsc:::predict.Hmsc(x,
        XData = as.data.frame(test_env[, -(1:3)]),
        studyDesign = samplingDesign_test, ranLevels = list(
          annee = rL1_test,
          site = rL2_test, habitat = rL3_test
        ), expected = expected
      ) %>%
        {
          array(unlist(.),
            dim = c(nrow(.[[1]]), ncol(.[[1]]), length(.)),
            dimnames = list(NULL, colnames(.[[1]]), NULL)
          )
        }

      if (expected) {

        # out <- apply(out, MARGIN = c(1, 2), mean, na.rm = TRUE) # For
        # probit and normal models, mean of the posterior is use, see
        # ?Hmsc::evaluateModelFit()

        out <- out %>%
          {
            data.frame(
              species = colnames(.),
              AUC = computeAUC(test, ., expected = expected),
              RMSE = computeRMSE(test, ., expected = expected),
              TjurR2 = computeTjurR2(test, ., expected = expected)
            )
          } %>%
          pivot_longer(-species, names_to = "measure", values_to = "value")
      } else {

        # out <- apply(out, MARGIN = c(1, 2), median, na.rm = TRUE) #
        # For Poisson model, median of the posterior is use, see
        # ?Hmsc::evaluateModelFit()

        out <- out %>%
          {
            data.frame(
              species = colnames(.),
              SR2 = computeR2(test, ., expected = expected),
              RMSE = computeRMSE(test, ., expected = expected),
              O.AUC = computeAUC(test, ., expected = expected),
              O.TjurR2 = computeTjurR2(test, ., expected = expected)
            )
          } %>%
          pivot_longer(-species, names_to = "measure", values_to = "value")
      }
    } else {
      test <- test_whole_fauna

      out <- Hmsc:::predict.Hmsc(x,
        XData = as.data.frame(test_env[, -(1:3)]),
        studyDesign = samplingDesign_test, ranLevels = list(
          annee = rL1_test,
          site = rL2_test, habitat = rL3_test
        ), expected = expected
      ) %>%
        {
          array(unlist(.),
            dim = c(nrow(.[[1]]), ncol(.[[1]]), length(.)),
            dimnames = list(NULL, colnames(.[[1]]), NULL)
          )
        }

      if (expected) {

        # out <- apply(out, MARGIN = c(1, 2), mean, na.rm = TRUE) # For
        # probit and normal models, mean of the posterior is use, see
        # ?Hmsc::evaluateModelFit()

        out <- out %>%
          {
            data.frame(
              species = colnames(.),
              AUC = computeAUC(test, ., expected = expected),
              RMSE = computeRMSE(test, ., expected = expected),
              TjurR2 = computeTjurR2(test, ., expected = expected)
            )
          } %>%
          pivot_longer(-species, names_to = "measure", values_to = "value")
      } else {

        # out <- apply(out, MARGIN = c(1, 2), median, na.rm = TRUE) #
        # For Poisson model, median of the posterior is use, see
        # ?Hmsc::evaluateModelFit()

        out <- out %>%
          {
            data.frame(
              species = colnames(.),
              SR2 = computeR2(test, ., expected = expected),
              RMSE = computeRMSE(test, ., expected = expected),
              O.AUC = computeAUC(test, ., expected = expected),
              O.TjurR2 = computeTjurR2(test, ., expected = expected)
            )
          } %>%
          pivot_longer(-species, names_to = "measure", values_to = "value")
      }
    }

    return(out)
  }, .id = "model") %>%
  filter(species %in% polychaeta)

### Relative change in model performances

# Train

## AUC
fit_train_pa <- prediction_trainset %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter(data_type == "PA", measure == "AUC") %>%
  mutate(
    species = factor(species, levels = polychaeta),
    model = factor(model)
  )

### Extract the baseline model
baseline_train_pa <- fit_train_pa %>%
  filter(model == "Benchmark") %>%
  arrange(species)

### Compute relative change for occurence model
fit_train_pa_relative <- fit_train_pa %>%
  group_by(model) %>%
  group_modify(function(x, y) {
    x %>%
      arrange(species) %>%
      mutate(relative_change = (value - baseline_train_pa$value) / abs(baseline_train_pa$value))
  }) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(!is.na(value)) %>%
  ungroup()

fit_train_pa_relative %>%
    group_by(model) %>%
    summarise(
        mean_AUC = mean(value), 
    )

fit_train_pa_relative %>%
    filter(data_type == "PA", measure == "AUC") %>%
    add_count(model, data_type, measure) %>%
    mutate(
        decrease = if_else(relative_change < 0, 1, 0),
        stable = if_else(relative_change == 0, 1, 0),
        increase = if_else(relative_change > 0, 1, 0)
    ) %>%
    group_by(model) %>%
        summarise(
            n_decrease = sum(decrease),
            n_stable = sum(stable),
            n_increase = sum(increase)
        )

p1 <- ggplot(fit_train_pa_relative, aes(x = model, y = relative_change, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(-.1, .1)) +
  labs(y = "Relative change in AUC - Train dataset") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p2 <- ggplot(
  filter(fit_train_pa_relative, model != "Benchmark"),
  aes(x = relative_change, fill = model)
) +
  geom_boxplot(
    data = filter(fit_train_pa_relative, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_train_pa_relative, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_train_pa_relative, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  coord_flip(xlim = c(-.1, .1)) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

p_train_auc <- p1 + p2 + plot_layout(width = c(2, 0.2), guides = "collect") & theme(
  legend.margin = margin(
    t = 0, r = 0,
    b = 0, l = 0
  ),
  legend.position = c(1.3, 0.5)
)

## RMSE
fit_train_ab <- prediction_trainset %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter(data_type == "AB", measure == "RMSE") %>%
  mutate(
    species = factor(species, levels = polychaeta),
    model = factor(model)
  )

### Extract the baseline model
baseline_train_ab <- fit_train_ab %>%
  filter(model == "Benchmark") %>%
  arrange(species)

fit_train_ab_relative <- fit_train_ab %>%
  group_by(model) %>%
  group_modify(function(x, y) {
    x %>%
      arrange(species) %>%
      mutate(relative_change = (value - baseline_train_ab$value) / abs(baseline_train_ab$value)) %>%
      # Case where we have a RMSE of 0 in the Bench model and in the WhC model
      mutate(
        relative_change = if_else(
          is.nan(relative_change),
          0,
          .$relative_change
        )
      )
  }) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(!is.na(value)) %>%
  ungroup()

fit_train_ab_relative %>%
    filter(data_type == "AB", measure == "RMSE") %>%
    add_count(model, data_type, measure) %>%
    mutate(
        decrease = if_else(relative_change > 0, 1, 0),
        stable = if_else(relative_change == 0, 1, 0),
        increase = if_else(relative_change < 0, 1, 0)
    ) %>%
    group_by(model) %>%
        summarise(
            n_decrease = sum(decrease),
            n_stable = sum(stable),
            n_increase = sum(increase)
        )

p3 <- ggplot(fit_train_ab_relative, aes(x = model, y = relative_change, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_reverse(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(max(fit_train_ab_relative$relative_change), min(fit_train_ab_relative$relative_change))) +
  labs(y = "Relative change in RMSE - Train dataset") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p4 <- ggplot(
  filter(fit_train_ab_relative, model != "Benchmark"),
  aes(x = relative_change, fill = model)
) +
  geom_boxplot(
    data = filter(fit_train_ab_relative, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_train_ab_relative, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_train_ab_relative, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  scale_x_reverse() +
  coord_flip(xlim = c(max(fit_train_ab_relative$relative_change), min(fit_train_ab_relative$relative_change))) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

p_train_rmse <- p3 + p4 + plot_layout(width = c(2, 0.2), guides = "collect") & theme(
  legend.margin = margin(
    t = 0, r = 0,
    b = 0, l = 0
  ),
  legend.position = c(1.3, 0.5)
)

# Test

## AUC

fit_test_pa <- prediction_testset %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter(data_type == "PA", measure == "AUC") %>%
  mutate(
    species = factor(species, levels = polychaeta),
    model = factor(model)
  )

### Extract the baseline model
baseline_test_pa <- fit_test_pa %>%
  filter(model == "Benchmark") %>%
  arrange(species)

### Compute relative change for occurence model

fit_test_pa_relative <- fit_test_pa %>%
  group_by(model) %>%
  group_modify(function(x, y) {
    x %>%
      arrange(species) %>%
      mutate(relative_change = (value - baseline_test_pa$value) / abs(baseline_test_pa$value))
  }) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(!is.na(value)) %>%
  ungroup()

fit_test_pa_relative %>%
    group_by(model) %>%
    summarise(
        mean_AUC = mean(value), 
        mean_relative_change = mean(relative_change)
    )

fit_test_pa_relative %>%
    filter(data_type == "PA", measure == "AUC") %>%
    add_count(model, data_type, measure) %>%
    mutate(
        decrease = if_else(relative_change < 0, 1, 0),
        stable = if_else(relative_change == 0, 1, 0),
        increase = if_else(relative_change > 0, 1, 0)
    ) %>%
    group_by(model) %>%
        summarise(
            n_decrease = sum(decrease),
            n_stable = sum(stable),
            n_increase = sum(increase)
        )

p5 <- ggplot(fit_test_pa_relative, aes(x = model, y = relative_change, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(min(fit_test_pa_relative$relative_change), max(fit_test_pa_relative$relative_change))) +
  labs(y = "Relative change in AUC - Test dataset") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p6 <- ggplot(
  filter(fit_test_pa_relative, model != "Benchmark"),
  aes(x = relative_change, fill = model)
) +
  geom_boxplot(
    data = filter(fit_test_pa_relative, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_test_pa_relative, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_test_pa_relative, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  coord_flip(xlim = c(min(fit_test_pa_relative$relative_change), max(fit_test_pa_relative$relative_change))) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

p_test_auc <- p5 + p6 + plot_layout(width = c(2, 0.2), guides = "collect") & theme(
  legend.margin = margin(
    t = 0, r = 0,
    b = 0, l = 0
  ),
  legend.position = c(1.3, 0.5)
)

## RMSE

fit_test_ab <- prediction_testset %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter(data_type == "AB", measure == "RMSE") %>%
  mutate(
    species = factor(species, levels = polychaeta),
    model = factor(model)
  )

### Extract the baseline model
baseline_test_ab <- fit_test_ab %>%
  filter(model == "Benchmark") %>%
  arrange(species)

### Compute relative change for abundance model
fit_test_ab_relative <- fit_test_ab %>%
  group_by(model) %>%
  group_modify(function(x, y) {
    x %>%
      arrange(species) %>%
      mutate(relative_change = (value - baseline_test_ab$value) / abs(baseline_test_ab$value)) %>%
      # Case where we have a RMSE of 0 in the Bench model and in the WhC model
      mutate(
        relative_change = if_else(
          is.nan(relative_change),
          0,
          .$relative_change
        )
      )
  }) %>%
  ungroup() %>%
  group_by(species) %>%
  filter(!is.na(value)) %>%
  ungroup()


fit_test_ab_relative %>%
    group_by(model) %>%
    summarise(
        mean_rmse = mean(value), 
        mean_relative_change = mean(relative_change)
    )

p7 <- ggplot(fit_test_ab_relative %>% filter(relative_change < Inf), aes(x = model, y = relative_change, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_reverse(
    labels = scales::percent_format(0.4),
    breaks = c(16, 14, 12, 10, 8, 6, 4, 2, 0, -1)
  ) +
  coord_cartesian(ylim = c(17, min(fit_test_ab_relative$relative_change))) +
  labs(y = "Relative change in RMSE - Test dataset") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p8 <- ggplot(
  filter(fit_test_ab_relative, model != "Benchmark", relative_change < Inf),
  aes(x = relative_change, fill = model)
) +
  geom_boxplot(
    data = filter(fit_test_ab_relative, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_test_ab_relative, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(fit_test_ab_relative, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  scale_x_reverse(breaks = c(16, 14, 12, 10, 8, 6, 4, 2, 0, -1)) +
  coord_flip(xlim = c(17, min(fit_test_ab_relative$relative_change, na.rm = TRUE))) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

p_test_rmse <- p7 + p8 + plot_layout(width = c(2, 0.2), guides = "collect") & theme(
    legend.margin = margin(
        t = 0, r = 0,
        b = 0, l = 0
    ),
    legend.position = c(1.3, 0.5)
)

fig2 <- (p1 +
  labs(
    y = "Relavitve change in AUC",
    title = "Train",
    tag = "A"
  ) +
  guides(fill = "none", color = "none") +
  theme(
    plot.title = element_text(hjust = .5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
) +
  p2 +
  (p5 +
    labs(
      y = "",
      title = "Test",
      tag = "B"
    ) +
    guides(fill = "none", color = "none") +
    theme(
      plot.title = element_text(hjust = .5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ) +
  p6 +
  (p3 +
    labs(
      y = "Relavitve change in RMSE",
      tag = "C"
    ) +
    guides(fill = "none", color = "none") +
    scale_x_discrete(labels = models_acc) +
    theme(axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]))
  ) +
  p4 +
  (p7 +
    labs(
      y = "",
      tag = "D"
    ) +
    scale_x_discrete(labels = models_acc) +
    theme(axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]))
  ) +
  p8 +
  plot_layout(ncol = 4, nrow = 2, widths = c(0.5, 0.1, 0.5, 0.1)) &
  theme(
    legend.position = "none",
    legend.margin = margin(
      t = 0, r = 0,
      b = 0, l = 0
    )
  )

############## Fig 3 - Comparison of the model performances to predict the community structure for models fitted with abundance data

# Bray-Curtis/Sorensen dissimilarity between observed community and predict one

diss_train <- raw_prediction_trainset %>%
  imap_dfr(function(x, y){
    if (str_detect(y, "PA")) {
      x <- 1 * (x > 0.5)

      out <- data.frame(
        obs_diss = as.vector(vegdist(decostand(train_polychaeta, method = "pa"))),
          pred_diss = as.vector(vegdist(x, method = "bray")),
          row.names = NULL
        )

    } else {

      out <- data.frame(
        obs_diss = as.vector(vegdist(decostand(train_polychaeta, method = "pa"))),
          pred_diss = as.vector(vegdist(x, method = "bray")),
          row.names = NULL
        )
    }

    return(out)

  }, .id = "model") %>%
  mutate(train_test = "train")

diss_test <- raw_prediction_testset %>%
  imap_dfr(function(x, y){
    if (str_detect(y, "PA")) {
      x <- 1 * (x > 0.5)

      out <- data.frame(
        obs_diss = as.vector(vegdist(decostand(test_polychaeta, method = "pa"))),
          pred_diss = as.vector(vegdist(x, method = "bray")),
          row.names = NULL
        )
    } else {

      out <- data.frame(
        obs_diss = as.vector(vegdist(decostand(test_polychaeta, method = "pa"))),
          pred_diss = as.vector(vegdist(x, method = "bray")),
          row.names = NULL
        )
    }

    return(out)

  }, .id = "model") %>%
  mutate(train_test = "test")

diss <- diss_train %>%
    rbind(diss_test) %>%
    separate(model, into = c('model', 'data_type'), sep = "\\s-\\s")

diss %>%
  group_by(model, data_type, train_test) %>%
  summarise(mean_obs = mean(obs_diss, na.rm =TRUE), mean_pred = mean(pred_diss, na.rm = TRUE))  %>%
  mutate(
    percent = round((mean_pred - mean_obs) / mean_obs * 100, 2),
    diff = mean_pred - mean_obs
  )

# Richness differences between observed community and predict one

richness_train <- 1 * (train_polychaeta > 0) %>%
    rowSums() %>%
    unname()

richness_test <- 1 * (test_polychaeta > 0) %>%
    rowSums() %>%
    unname()

richness_models_train <- raw_prediction_trainset %>%
    map(function(x) {
        rowSums(1 * (x > 0.5)) %>%
            unname()
    }) %>%
    imap_dfr(~ data.frame(
        obs_richness = richness_train,
        pred_richness = .,
        diff_richness = . - richness_train),
        .id = "model") %>%
    mutate(train_test = "train")

richness_models_test <- raw_prediction_testset %>%
    map(function(x) {
        rowSums(1 * (x > 0.5)) %>%
            unname()
    }) %>%
    imap_dfr(~ data.frame(
        obs_richness = richness_test,
        pred_richness = .,
        diff_richness = . - richness_test),
        .id = "model") %>%
    mutate(train_test = "test")

richness_models <- richness_models_train %>%
    rbind(richness_models_test) %>%
    separate(model, into = c('model', 'data_type'), sep = "\\s-\\s")

richness_models %>%
  group_by(model, data_type, train_test) %>%
  summarise(mean_obs = mean(obs_richness), mean_pred = mean(pred_richness))  %>%
  mutate(
    percent = round((mean_pred - mean_obs) / mean_obs * 100, 2),
    diff = mean_pred - mean_obs
  )

# Abundance differences between observed community and predict one

abundance_train <- train_polychaeta %>%
    rowSums() %>%
    unname()

abundance_test <- test_polychaeta %>%
    rowSums() %>%
    unname()

abundance_models_train <- raw_prediction_trainset %>%
    keep(rep(c(TRUE, FALSE), 4)) %>% # Keep only abundance based models
    map(~ unname(rowSums(.x))) %>%
    imap_dfr(~ data.frame(
        obs_abundance = abundance_train,
        pred_abundance = .,
        diff_abundance = . - abundance_train),
        .id = "model") %>%
    mutate(train_test = "train")


abundance_models_test <- raw_prediction_testset %>%
    keep(rep(c(TRUE, FALSE), 4)) %>% # Keep only abundance based models
    map(~unname(rowSums(.x))) %>%
    imap_dfr(~ data.frame(
        obs_abundance = abundance_test,
        pred_abundance = .,
        diff_abundance = . - abundance_test),
        .id = "model") %>%
    mutate(train_test = "test")

abundance_models <- abundance_models_train %>%
  rbind(abundance_models_test) %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s")

abundance_models %>%
  group_by(model, data_type, train_test) %>%
  summarise(mean_obs = mean(obs_abundance), mean_pred = mean(pred_abundance), .groups = "keep") %>%
  mutate(
    percent = round((mean_pred - mean_obs) / mean_obs * 100, 2),
    diff = mean_pred - mean_obs
  )

p1 <- ggplot(
    diss %>% filter(data_type == "AB", !is.na(pred_diss)),
    aes(x = model, y = pred_diss - obs_diss, color = train_test)
) +
stat_pointinterval(aes(group = interaction(model, train_test)),
    position = position_dodge(width = .75),
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Dataset",
    breaks = c("train", "test"),
    labels = c("Train", "Test"),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  labs(
    y = "Difference between predicted and observed Bray-Curtis dissimilarity",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank()
)

p2 <- ggplot(
    richness_models %>% filter(data_type == "AB"),
    aes(x = model, y = diff_richness, color = train_test)
) +
stat_pointinterval(aes(group = interaction(model, train_test)),
    position = position_dodge(width = .75),
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Dataset",
    breaks = c("train", "test"),
    labels = c("Train", "Test"),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  labs(
    y = "Difference between predicted and observed richness",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank()
)

p3 <- ggplot(
    abundance_models,
    aes(x = model, y = diff_abundance, color = train_test)
) +
stat_pointinterval(aes(group = interaction(model, train_test)),
    position = position_dodge(width = .75),
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Dataset",
    breaks = c("train", "test"),
    labels = c("Train", "Test"),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  labs(
    y = "Difference between predicted and observed abundances",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank()
)

fig3 <- (p1 + p2 + p3 & scale_x_discrete(labels = models_acc) &
  theme(
    axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]),
    axis.title.x = element_blank()
  )) +
  plot_layout(guides = "collect")

############## Fig 4 - Change in explained variance related to environmental predictors 

### Computation of the variance partitioning

VP <- models %>%
  map_dfr(function(x) {
    x %>%
      computeVariancePartitioning(
        group = 1,
        groupnames = "Environment"
      ) %>%
      .$vals %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("species") %>%
      filter(species %in% polychaeta) %>%
      pivot_longer(-species, names_to = "variable", values_to = "value")
  }, .id = "model") %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  mutate(species = factor(species, levels = unique(species)))

VP <- prediction_trainset %>%
  separate("model", into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter((data_type == "AB" & measure == "SR2") | (data_type == "PA" & measure == "TjurR2")) %>%
  right_join(VP, by = c("model", "data_type", "species")) %>%
  rename("metric" = "value.x", "explained_var" = "value.y") %>%
  mutate("percent_total_var" = metric * explained_var)

VP %>% filter(variable == "Environment") %>%
  group_by(model, data_type) %>%
  summarise(
    mean_exp_var = mean(explained_var), sd_exp_var = sd(explained_var),
    mean_tot_var = mean(percent_total_var, na.rm = TRUE), sd_tot_var = sd(percent_total_var, na.rm = TRUE)
    ) %>%
  ungroup()

### Relative change in variance partitioning

VP_relatif_explained_var <- VP %>%
  group_by(data_type) %>%
  group_map(function(x, y) {
    benchmark <- x %>%
      mutate(variable = str_remove(variable, ":\\s[:alpha:]*")) %>%
      group_by(model, data_type, species, variable) %>%
      summarise(value = sum(explained_var), .groups = "drop") %>%
      filter(model == "Benchmark") %>%
      pivot_wider(names_from = "variable", values_from = "value")


    out <- x %>%
      mutate(variable = str_remove(variable, ":\\s[:alpha:]*")) %>%
      group_by(model, data_type, species, variable) %>%
      summarise(value = sum(explained_var), .groups = "drop") %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      group_by(model) %>%
      group_map(~ mutate(.x,
        Environment_relat = (Environment - benchmark$Environment) / benchmark$Environment,
        Random_relat = (Random - benchmark$Random) / benchmark$Random
      ), .keep = TRUE)
  }, .keep = TRUE) %>%
  map_dfr(~.x)

VP_relatif_tot_var <- VP %>%
  group_by(data_type) %>%
  group_map(function(x, y) {
    benchmark <- x %>%
      mutate(variable = str_remove(variable, ":\\s[:alpha:]*")) %>%
      group_by(model, data_type, species, variable) %>%
      summarise(value = sum(percent_total_var), .groups = "drop") %>%
      filter(model == "Benchmark") %>%
      pivot_wider(names_from = "variable", values_from = "value")


    out <- x %>%
      mutate(variable = str_remove(variable, ":\\s[:alpha:]*")) %>%
      group_by(model, data_type, species, variable) %>%
      summarise(value = sum(percent_total_var), .groups = "drop") %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      group_by(model) %>%
      group_map(~ mutate(.x,
        Environment_relat = (Environment - benchmark$Environment) / benchmark$Environment,
        Random_relat = (Random - benchmark$Random) / benchmark$Random
      ), .keep = TRUE)
  }, .keep = TRUE) %>%
  map_dfr(~.x)

### Panel of the variance explained by the environment for abundance data

p_vp_ab_env1 <- ggplot(VP_relatif_tot_var %>% filter(data_type == "AB"), aes(x = model, y = Environment_relat, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(
    min(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "AB"], na.rm = TRUE),
    max(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "AB"], na.rm = TRUE)
  )) +
  labs(y = "Relative change in Variance explained by Environment") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p_vp_ab_env2 <- ggplot(
  filter(
    VP_relatif_tot_var, model != "Benchmark",
    data_type == "AB"
  ),
  aes(x = Environment_relat, fill = model)
) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  scale_x_reverse() +
  coord_flip(xlim = c(
    min(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "AB"], na.rm = TRUE),
    max(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "AB"], na.rm = TRUE)
  )) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

### Panel of the variance explained by the random effects for abundance data

p_vp_ab_rnd1 <- ggplot(VP_relatif_tot_var %>% filter(data_type == "AB"), aes(x = model, y = Random_relat, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(
    min(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "AB"], na.rm = TRUE),
    max(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "AB"])
  )) +
  labs(y = "Relative change in Variance explained by random effect") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p_vp_ab_rnd2 <- ggplot(
  filter(
    VP_relatif_tot_var, model != "Benchmark",
    data_type == "AB"
  ),
  aes(x = Random_relat, fill = model)
) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  coord_flip(xlim = c(
    min(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "AB"]),
    max(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "AB"])
  )) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

### Panel of the variance explained by the environment for occurence data

p_vp_pa_env1 <- ggplot(VP_relatif_tot_var %>% filter(data_type == "PA"), aes(x = model, y = Environment_relat, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(
    min(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "PA"], na.rm = TRUE),
    max(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "PA"], na.rm = TRUE)
  )) +
  labs(y = "Relative change in Variance explained by Environment") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p_vp_pa_env2 <- ggplot(
  filter(
    VP_relatif_tot_var, model != "Benchmark",
    data_type == "PA"
  ),
  aes(x = Environment_relat, fill = model)
) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  scale_x_reverse() +
  coord_flip(xlim = c(
    min(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "PA"]),
    max(VP_relatif_tot_var$Environment_relat[VP_relatif_tot_var$data_type == "PA"])
  )) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

### Panel of the variance explained by the random effects for occurence data

p_vp_pa_rnd1 <- ggplot(VP_relatif_tot_var %>% filter(data_type == "PA"), aes(x = model, y = Random_relat, color = model, group = species)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  scale_y_continuous(labels = scales::percent_format(0.4)) +
  coord_cartesian(ylim = c(
    min(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "PA"], na.rm = TRUE),
    max(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "PA"], na.rm = TRUE)
  )) +
  labs(y = "Relative change in Variance explained by random effect") +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "points"), axis.title.x = element_blank())

p_vp_pa_rnd2 <- ggplot(
  filter(
    VP_relatif_tot_var, model != "Benchmark",
    data_type == "PA"
  ),
  aes(x = Random_relat, fill = model)
) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny"),
    aes(y = -2), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Polychaeta Phylogeny Traits"),
    aes(y = -1), outlier.shape = NA, width = 0.5
  ) +
  geom_boxplot(
    data = filter(VP_relatif_tot_var, model == "Whole community"),
    aes(y = 0), outlier.shape = NA, width = 0.5
  ) +
  scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[2:4], guide = "none") +
  coord_flip(xlim = c(
    min(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "PA"]),
    max(VP_relatif_tot_var$Random_relat[VP_relatif_tot_var$data_type == "PA"])
  )) + # Make sure that both plot are on the same y-axis scale
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = unit(c(5.5, 0, 5.5, 0), units = "points"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing = unit(0, units = "points")
  )

### Plot assembly

fig4 <- (p_vp_pa_env1 +
  labs(
    y = "Relative change in variance (%)",
    title = "Environment",
    tag = "A"
  ) +
  guides(fill = "none", color = "none") +
  theme(
    plot.title = element_text(hjust = .5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
) +
  p_vp_pa_env2 +
  (p_vp_pa_rnd1 +
    labs(
      title = "Random effects",
      tag = "B"
    ) +
    guides(fill = "none", color = "none") +
    theme(
      plot.title = element_text(hjust = .5),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ) +
  p_vp_pa_rnd2 +
  (p_vp_ab_env1 +
    labs(
      y = "Relative change in variance (%)",
      tag = "C"
    ) +
    guides(fill = "none", color = "none") +
    scale_x_discrete(labels = models_acc) +
    theme(axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]))
  ) +
  p_vp_ab_env2 +
  (p_vp_ab_rnd1 +
    labs(
      y = NULL,
      tag = "D"
    ) +
    guides(fill = "none", color = "none") +
    scale_x_discrete(labels = models_acc) +
    theme(axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]))
  ) +
  p_vp_ab_rnd2 +
  plot_layout(
    ncol = 4, nrow = 2,
    widths = c(0.5, 0.1, 0.5, 0.1),
    byrow = TRUE
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = .5))
  ) &
  theme(
    legend.position = "none",
    legend.margin = margin(
      t = 0, r = 0,
      b = 0, l = 0
    )
  )

############## Fig 5 - Number & proportion of response curves classify accordingly to Rigal et al. (2020)'s methodology for models fitted with abundance data

### Compute respone curves as Rigal et al. 2020

#### Compute delta and center values as in Rigal et al. 2020
val <- models$`Benchmark - AB`$XData %>%
  summarise(across(everything(), list(min = min, max = max))) %>%
  pivot_longer(everything(), names_to = "Var", values_to = "values") %>%
  mutate(
    Var2 = str_extract(Var, "min|max"),
    Var = str_remove(Var, "_min|_max")
  ) %>%
  pivot_wider(names_from = Var2, values_from = values) %>%
  mutate(
    center = (min + max) / 2,
    delta = (min + center) / 2
  ) %>%
  mutate(Var = unique(str_remove_all(covNames, "_[:digit:]"))[-1])

#### Compute and classify response curves
coeff <- models %>%
  map_dfr(function(x, supportLevel = .95) {
    rowNames <- covNames

    colNames <- x$spNames

    getPostEstimate(x, "Beta") %>%
      {
        .$mean * ((.$support > supportLevel) + (.$support < (1 - supportLevel)) > 0)
      } %>%
      as.data.frame() %>%
      select(any_of(polychaeta)) %>%
      mutate(
        Variable = rowNames,
        .before = `Melinna palmata`
      )
  }, .id = "Model") %>%
  separate(Model, into = c("Model", "Data_type"), sep = "\\s-\\s") %>%
  mutate(
    Model = factor(Model),
    Variable = factor(Variable, levels = covNames)
  ) %>%
  pivot_longer(-c(Model:Variable), names_to = "Species", values_to = "value") %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  pivot_longer(-c(Model:Intercept), names_to = "Variable", values_to = "Values") %>%
  mutate(
    Variable2 = str_remove(Variable, "_[:digit:]"),
    Variable = str_remove(Variable, "\\s"),
    Variable = str_replace_all(Variable, "[:alpha:]*_1", "Linear"),
    Variable = str_replace_all(Variable, "[:alpha:]*_2", "Quadratic")
  ) %>%
  rename(
    Coeff = Variable,
    Variable = Variable2
  ) %>%
  pivot_wider(names_from = "Coeff", values_from = "Values") %>%
  relocate(Variable, .before = Species)

val <- models$`Benchmark - AB`$XData %>%
  summarise(across(everything(), list(min = min, max = max))) %>%
  pivot_longer(everything(), names_to = "Var", values_to = "values") %>%
  mutate(
    Var2 = str_extract(Var, "min|max"),
    Var = str_remove(Var, "_min|_max")
  ) %>%
  pivot_wider(names_from = Var2, values_from = values) %>%
  mutate(
    center = (min + max) / 2,
    delta = (min + center) / 2
  ) %>%
  mutate(Var = unique(str_remove_all(covNames, "_[:digit:]"))[-1])

coeff <- coeff %>%
  left_join(val, by = c("Variable" = "Var")) %>%
  mutate(
    deriv_lower = Linear + 2 * Quadratic * (center - delta),
    deriv_center = Linear + 2 * Quadratic * center,
    deriv_upper = Linear + 2 * Quadratic * (center + delta),
    curvature = ((-12 * Quadratic^2 * (2 * Quadratic * center + Linear)) / (1 + (2 * Quadratic * center + Linear)^2)^(5 / 2))
  ) %>%
  mutate(classif = map_chr(split(., 1:nrow(.)), ~ polynom_classif(.x)))

rigal_class <- coeff %>%
  mutate(classif = factor(classif)) %>%
  group_by(Model, Data_type) %>%
  group_modify(~ count(.x, classif, .drop = FALSE)) %>%
  ungroup() %>%
  group_by(Model, Data_type) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

rigal_class_var <- coeff %>%
  mutate(classif = factor(classif)) %>%
  group_by(Model, Data_type, Variable) %>%
  group_modify(~ count(.x, classif, .drop = FALSE)) %>%
  ungroup() %>%
  group_by(Model, Data_type) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

### Images
accelerate_decline <- image_read("Data/Images/accelerate_decline.png") %>%
    as.raster()

convexe <- image_read("Data/Images/convexe.png") %>%
    as.raster()

accelerate_increase <- image_read("Data/Images/accelerate_increase.png") %>%
    as.raster()

constant_decline <- image_read("Data/Images/constant_decrease.png") %>%
    as.raster()

constant_increase <- image_read("Data/Images/constant_increase.png") %>%
    as.raster()

decelerated_decline <- image_read("Data/Images/decelerate_decline.png") %>%
    as.raster()

concave <- image_read("Data/Images/concave.png") %>%
    as.raster()

decelerated_increase <- image_read("Data/Images/decelerate_increase.png") %>%
    as.raster()


# Abundance

p1 <- ggplot(
  rigal_class %>% filter(
    classif == "Accelerate decline",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(accelerate_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Accelerate decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p2 Convex

p2 <- ggplot(
  rigal_class %>% filter(
    classif == "Convex",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  )

if (nrow(p2$data) > 0) p2 <- p2 + scale_x_discrete(drop = FALSE)
p2 <- p2 +
  annotation_raster(convexe, 0, 1, 300, 500) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Convex") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p3 Accemlerate increase

p3 <- ggplot(
  rigal_class %>% filter(
    classif == "Accelerated increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(accelerate_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Accelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p4 Constant decline

p4 <- ggplot(
  rigal_class %>% filter(
    classif == "Constant decrease",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(constant_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Constant decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p5 Stable

p5 <- ggplot(
  rigal_class %>% filter(
    classif == "Stable",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Stable") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p6 Constant increase

p6 <- ggplot(
  rigal_class %>% filter(
    classif == "Constant increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(constant_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Constant increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p7 Decelerated decline

p7 <- ggplot(
  rigal_class %>% filter(
    classif == "Decelerated decline",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(decelerated_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Decelerated decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p8 Concave

p8 <- ggplot(
  rigal_class %>% filter(
    classif == "Concave",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(concave, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Concave") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p9 Decelerated increase

p9 <- ggplot(
  rigal_class %>% filter(
    classif == "Decelerated increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(decelerated_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Decelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


fig5 <-
  (p1 + p2 + p3 +
    p4 + p5 + p6 +
    p7 + p8 + p9 &
    theme(legend.position = "bottom")) +
  plot_layout(
    ncol = 3, nrow = 3,
    byrow = TRUE,
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = .5))
)

############## Fig 6 - Comparison of residual species-species correlations between WhC model and Bench

residual_corr_raw <- models %>%
  map(function(x) {
    x %>%
      computeAssociations() %>%
      map(function(x, supportLevel = 0) {
        tmp <- x$mean * ((x$support > supportLevel) + (x$support < (1 - supportLevel)) > 0)

        return(tmp)
      }) %>%
      `names<-`(c("year", "site", "habitat")) %>%
      map(function(x) {
        x %>%
          as.data.frame() %>%
          rownames_to_column("species_1") %>%
          filter(species_1 %in% polychaeta) %>%
          column_to_rownames("species_1") %>%
          select(any_of(polychaeta))
      })
  })

# Benchmark model - Abundance

res_corr_bench_ab_df <- residual_corr_raw[[1]] %>%
  map_dfr(function(x) {
    x %>%
      as.data.frame() %>%
      rownames_to_column("species1") %>%
      pivot_longer(-species1, names_to = "species2", values_to = "corr") %>%
      filter(species1 != species2) %>%
      rename(corr_bench = corr)
  }, .id = "random_effect")

# Whole community model - Abundance
res_corr_whole_ab_df <- residual_corr_raw[[7]] %>%
  map_dfr(function(x) {
    x %>%
      as.data.frame() %>%
      rownames_to_column("species1") %>%
      pivot_longer(-species1, names_to = "species2", values_to = "corr") %>%
      filter(species1 != species2) %>%
      rename(corr_whole = corr)
  }, .id = "random_effect")

# Combine both models
res_corr_comp_ab <- res_corr_bench_ab_df %>%
  left_join(res_corr_whole_ab_df, by = c("random_effect", "species1", "species2")) %>%
  drop_na(corr_bench, corr_whole)

# Compute linear R2 between both
r2_df_ab <- res_corr_comp_ab %>%
  group_by(random_effect) %>%
  group_map(function(x, y) {
    cat(str_to_title(unlist(y)), sep = "\n")

    reg <- lm(corr_whole ~ corr_bench, data = x) %>%
      summary() %>%
      .$adj.r.squared

    data.frame(
      random_effect = unlist(y),
      r2 = reg,
      row.names = NULL
    )
  }) %>%
  map_dfr(~.x)

p_res_corr1 <- ggplot(res_corr_comp_ab, aes(x = corr_bench, y = corr_whole)) +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c(name = "Number of neighbors") +
  geom_smooth(
    method = "lm",
    formula = y ~ x
  ) +
  geom_text(data = r2_df_ab, x = -.45, y = .9, aes(label = paste0("R^2:  ", round(r2, 2))), parse = TRUE) +
  facet_row(~random_effect,
    labeller = labeller(random_effect = c(
      "habitat" = "Habitat",
      "site" = "Site",
      "year" = "Year"
    ))
  ) +
  labs(
    x = "Residual Correlations  Bench",
    y = "Residual Correlations WhC"
  ) +
  theme(legend.position = "none")

p_res_corr2 <- ggplot(res_corr_comp_ab, aes(x = abs(corr_whole - corr_bench) * (sign(corr_whole * corr_bench)))) +
  geom_histogram(aes(y = after_stat(density)),
    fill = "transparent",
    col = "black",
    bins = 30
  ) +
  geom_density(col = "red") +
  facet_row(~random_effect,
    labeller = labeller(random_effect = c(
      "habitat" = "Habitat",
      "site" = "Site",
      "year" = "Year"
    ))
  ) +
  labs(
    x = "Index value",
    y = "Density",
  )

fig6 <- p_res_corr1 / p_res_corr2 + plot_annotation(tag_levels = "A")

res_corr_comp_ab %>%
  mutate(index = abs(corr_whole - corr_bench) * (sign(corr_whole * corr_bench))) %>%
  add_count(random_effect) %>%
  filter(abs(index) > 0.25, index < 0) %>%
  add_count(random_effect) %>%
  group_by(random_effect) %>%
  summarise(percent = unique(nn) / unique(n) * 100)


############## Fig S1 - Map of the sampled sites

france <- st_as_sf(maps::map("france", plot = FALSE, fill = TRUE))

world <- map_data("world")

fill_world <- rep("gray90", nrow(world))
fill_world[world$region == "France"] <- "black"

p_site <- ggplot(france) +
  geom_sf() +
  geom_point(
    data = coord_sites,
    aes(
      x = longitude,
      y = latitude,
      color = theme,
      shape = test
    ),
    size = 3,
  ) +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "tl") +
  coord_sf(
    xlim = c(-5.5, -1),
    ylim = c(47.25, 49)
  ) +
  scale_shape_discrete(name = "", labels = c("Train", "Test")) +
  scale_color_manual(
    labels = c("Intertidal seagrass meadows", "Intertidal bare sediments", "Both habitats"),
    values = c("#35B779FF", "#f1a340", "#31688EFF"),
    name = ""
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    # legend.text = element_text(size = 15)
  )

p_france <- ggplot(france) +
  geom_sf(color = NA) +
  geom_sf(
    data = france %>% filter(ID %in% c(
      "Ille-et-Vilaine",
      "Finistere",
      "Cotes-Darmor",
      "Morbihan"
    )),
    color = NA,
    fill = "black",
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(.8)))

p_world <- ggplot(world, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  geom_polygon(fill = fill_world, col = "grey45") +
  coord_map("ortho", orientation = c(41, -5, 0)) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text = element_blank())

layout <- c(
  area(t = 1, l = 1, b = 3, r = 2),
  area(t = 2, l = 3, b = 3, r = 3),
  area(t = 1, l = 3, b = 1, r = 3)
)

fig_supp1 <- (
    p_site +
    p_france +
    p_world
    ) +
  plot_layout(design = layout) &
      theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.box.margin = margin(t = -150)
      )
  
############## Fig S2 - Fuzzy PCA of the species by trait matrix

trait_fuz_estim <- trait_fuz_estim %>%
  tibble::rownames_to_column() %>%
  filter(rowname %in% polychaeta) %>%
  tibble::column_to_rownames()

trait_block <- c(`Maximum size` = 7, `Feeding method` = 8, `Food size` = 2, `Preferred substrate position` = 2, `Living habitat` = 5, `Adult movement capacity (daily)` = 4, Bioturbation = 5, `Sexual differenciation` = 2, `Development mode` = 4, `Reproduction frequency` = 2)

trait_fuz <- ade4::prep.fuzzy.var(trait_fuz_estim, col.blocks = trait_block)

trait_fuz_fpca <- ade4::dudi.fpca(trait_fuz, scannf = FALSE, nf = Inf)

traits_fpca <- rownames_to_column(trait_fuz_fpca$c1[, 1:3], "traits") %>%
  mutate(
    norm_1_2 = sqrt(CS1^2 + CS2^2),
    norm_1_3 = sqrt(CS1^2 + CS3^2),
    traits_1_2 = if_else(norm_1_2 >= 0.4, traits, ""), # Scaling of cleanplot.pca() by P.Legendre
    traits_1_3 = if_else(norm_1_3 >= 0.4, traits, "")
  )

traits_fpca$traits_1_2[traits_fpca$traits_1_2 == ""] <- NA
traits_fpca$traits_1_3[traits_fpca$traits_1_3 == ""] <- NA

p1 <- ggplot(traits_fpca, aes(x = CS1, y = CS2)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(aes(
    x = 0, y = 0,
    xend = CS1, yend = CS2
  ),
  arrow = arrow(length = unit(0.15, units = "cm"), type = "closed")
  ) +
  ggrepel::geom_text_repel(aes(label = traits_1_2),
    segment.size = 0
  ) +
  labs(
    x = paste0(
      "Axis 1 - ",
      round(trait_fuz_fpca$eig / sum(trait_fuz_fpca$eig) * 100, 2)[1], " %"
    ),
    y = paste0(
      "Axis 2 - ",
      round(trait_fuz_fpca$eig / sum(trait_fuz_fpca$eig) * 100, 2)[2], " %"
    )
  ) +
  coord_fixed(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  )

p2 <- ggplot(traits_fpca, aes(x = CS1, y = CS3)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(aes(
    x = 0, y = 0,
    xend = CS1, yend = CS3
  ),
  arrow = arrow(length = unit(0.15, units = "cm"), type = "closed")
  ) +
  ggrepel::geom_text_repel(aes(label = traits_1_3), segment.size = 0) +
  labs(
    x = paste0("Axis 1 - ", round(trait_fuz_fpca$eig / sum(trait_fuz_fpca$eig) * 100, 2)[1], " %"),
    y = paste0("Axis 3 - ", round(trait_fuz_fpca$eig / sum(trait_fuz_fpca$eig) * 100, 2)[3], " %")
  ) +
  coord_fixed(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  )

fig_supp2 <- p1 + p2

############## Fig S3 - Distribution of the richness in the sites of the train/test dataset

p1 <- ggplot(data.frame(richness = richness_train), aes(x = richness)) +
  geom_bar() +
  scale_x_binned(show.limits = TRUE) +
  labs(x = "Richness", y = "Count")

p2 <- ggplot(data.frame(richness = richness_test), aes(x = richness)) +
  geom_bar() +
  scale_x_binned(show.limits = TRUE) +
  labs(x = "Richness", y = "Count")

fig_supp3 <- p1 / p2 + plot_annotation(tag_levels = "A")

############## Fig S4 - Distribution of the total abundance in the sites of the train/test dataset

p1 <- ggplot(data.frame(abundance = abundance_train), aes(x = abundance)) +
  geom_bar() +
  scale_x_binned(n.breaks = 15, show.limits = TRUE) +
  labs(x = "Abundance", y = "Count")

p2 <- ggplot(data.frame(abundance = abundance_test), aes(x = abundance)) +
  geom_bar() +
  scale_x_binned(n.breaks = 15, show.limits = TRUE) +
  labs(x = "Abundance", y = "Count")

fig_supp4 <- p1 / p2 + plot_annotation(tag_levels = "A")

############## Fig S15 - Distribution of PSRF fofr response curve shapes for abundance models 

load("convergence_diag.RData")

convergence_classif_plot <- beta_conv %>%
    filter(Var1 != "(Intercept)") %>%
    separate("model", into = c("Model", "Data_type"), sep = "\\s-\\s") %>%
    mutate(Variable = str_remove(Var1, "_[:digit:]"), Power = str_extract(Var1, "[:digit:]"), Species = Var2) %>%
    mutate(Variable = str_replace_all(Variable, c(
        "Average" = "Fetch",
        "MO" = "Organic Mater",
        "SAL_sd_w" = "Salinity",
        "TEMP_sd_w" = "Temperature",
        "CURR_mean_w" = "Current",
        "mud" = "Mud",
        "\\`Trask\\(So\\)\\`" = "Trask Index"
    ))) %>%
    filter(Variable != "(Intercept)") %>%
        select(Model, Data_type, Variable, Power, Species, psrf, ess) %>%
        left_join(
            select(coeff, Model, Data_type, Variable, Species, classif),
            by = c("Model", "Data_type", "Variable", "Species")
        ) %>%
    filter(!is.na(classif)) %>% # Remove line for non-polychaeta species
    mutate(classif = factor(classif, levels = c(
            "Accelerate decline", "Convex", "Accelerated increase",
            "Constant decrease", "Stable", "Constant increase",
            "Decelerated decline", "Concave", "Decelerated increase"
        )))

fig_supp15 <- ggplot(
    convergence_classif_plot %>% filter(Data_type == "AB"),
    aes(x = classif, y = psrf, fill = Model)
) +
  geom_segment(aes(y = 1.2, yend = 1.2,
                       x = -Inf, xend = Inf),
                   linetype = "dashed",
                   color = "red") +
  geom_boxplot(outlier.shape = NA, color = "black") +
      scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
      ylim(c(1, 2)) +
      labs(x = "", y = "PSRF", title = "Abundance") +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = .5)
      )

############## Fig S16 - Distribution of PSRF for response curve shapes for presence/absence models

fig_supp16 <- ggplot(
    convergence_classif_plot %>% filter(Data_type == "PA"),
    aes(x = classif, y = psrf, fill = Model)
) +
  geom_boxplot(color = "black") +
      scale_fill_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
      ylim(c(1, 1.03)) +
      labs(x = "", y = "PSRF", title = "Presence / Absence") +
      theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = .5)
      )

############## Fig S17 - Comparison of explanatory and predictive power for all models

predictions <- mutate(prediction_trainset, train_test = "train") %>%
  rbind(mutate(prediction_testset, train_test = "test")) %>%
  mutate(measure2 = str_remove(measure, "O.")) %>%
  separate(model, into = c("model", "data_type"), sep = "\\s-\\s") %>%
  filter(
    measure2 %in% c("AUC", "RMSE"),
    !is.na(value)
  ) %>%
  mutate(
    value = if_else(measure2 == "RMSE", log10(value + 1), value),
    train_test = factor(train_test, levels = c("train", "test"))
  )

p1 <- ggplot(
  predictions %>% filter(measure2 == "AUC", train_test == "train"),
  aes(x = model, y = value, fill = model)
) +
  ggdist::stat_slab(
    alpha = .5,
    position = position_nudge(x = .1, y = 0),
    adjust = 1,
    scale = .5
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = .1,
    position = position_dodge(width = .2),
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      x = as.numeric(as.factor(model)) - .15,
      y = value, color = model
    ),
    alpha = .2,
    shape = "|",
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    guide = guide_legend(override.aes = aes(alpha = NA)),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  labs(
    x = "",
    y = "AUC",
    title = "Train",
    tag = "A"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside"
  )

p2 <- ggplot(
  predictions %>% filter(measure2 == "AUC", train_test == "test"),
  aes(x = model, y = value, fill = model)
) +
  ggdist::stat_slab(
    alpha = .5,
    position = position_nudge(x = .1, y = 0),
    adjust = 1,
    scale = .5
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = .1,
    position = position_dodge(width = .2),
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      x = as.numeric(as.factor(model)) - .15,
      y = value, color = model
    ),
    alpha = .2,
    shape = "|",
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    guide = guide_legend(override.aes = aes(alpha = NA)),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  labs(
    x = "",
    y = "",
    title = "Test",
    tag = "B"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.placement = "outside"
  )

p3 <- ggplot(
  predictions %>% filter(measure2 == "RMSE", train_test == "train"),
  aes(x = model, y = value, fill = model)
) +
  ggdist::stat_slab(
    alpha = .5,
    position = position_nudge(x = .1, y = 0),
    adjust = 1,
    scale = .5
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = .1,
    position = position_dodge(width = .2),
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      x = as.numeric(as.factor(model)) - .15,
      y = value, color = model
    ),
    alpha = .2,
    shape = "|",
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    guide = guide_legend(override.aes = aes(alpha = NA)),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  scale_x_discrete(labels = models_acc) +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  ylim(c(0, 4)) +
  labs(
    x = "",
    y = "log10(RMSE)",
    title = "",
    tag = "C"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]),
    axis.ticks.x = element_blank(),
    strip.placement = "outside"
  )

p4 <- ggplot(
  predictions %>% filter(measure2 == "RMSE", train_test == "test"),
  aes(x = model, y = value, fill = model)
) +
  ggdist::stat_slab(
    alpha = .5,
    position = position_nudge(x = .1, y = 0),
    adjust = 1,
    scale = .5
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = .1,
    position = position_dodge(width = .2),
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      x = as.numeric(as.factor(model)) - .15,
      y = value, color = model
    ),
    alpha = .2,
    shape = "|",
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    guide = guide_legend(override.aes = aes(alpha = NA)),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  scale_x_discrete(labels = models_acc) +
  scale_color_manual(values = rcartocolor::carto_pal(12, "Antique")[1:4]) +
  labs(
    x = "",
    y = "",
    title = "",
    tag = "D"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]),
    axis.ticks.x = element_blank(),
    strip.placement = "outside"
  )

fig_supp17 <- (p1 + p2 + p3 + p4 &
  theme(
    plot.title = element_text(hjust = .5)
  )) +
  plot_layout(
    ncol = 2, nrow = 2,
    byrow = TRUE,
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      legend.position = "none"
    )
  )

############## Fig S18 - Mean Predicted abundance vs Mean observed abundance on train dataset 

mean_obs_pred_train <- raw_prediction_trainset[c(1, 3, 5, 7)] %>%
  map_dfr(function(x){

    pred <- data.frame(
      species = colnames(x),
      mean_pred = colMeans(x),
      sd_pred = apply(x, MARGIN = 2, sd),
      row.names = NULL
    ) %>%
    mutate(se_pred = sd_pred / sqrt(180)) # We have 180 obs in the trainning dataset

    obs <- data.frame(
      species = colnames(train_polychaeta),
      mean_obs = colMeans(train_polychaeta),
      sd_obs = apply(train_polychaeta, MARGIN = 2, sd),
      row.names = NULL
      ) %>%
    mutate(se_obs = sd_obs / sqrt(180))

    out <- pred %>%
      left_join(obs, by = "species")

    return(out)
  }, .id = "Model") %>%
  filter(species %in% polychaeta)  %>%
  mutate(Model = str_remove(Model, "\\s-\\sAB"))

mean_obs_pred_test <- raw_prediction_testset[c(1, 3, 5, 7)] %>%
  map_dfr(function(x){

    pred <- data.frame(
      species = colnames(x),
      mean_pred = colMeans(x),
      sd_pred = apply(x, MARGIN = 2, sd),
      row.names = NULL
    ) %>%
    mutate(se_pred = sd_pred / sqrt(180)) # We have 180 obs in the trainning dataset

    obs <- data.frame(
      species = colnames(train_polychaeta),
      mean_obs = colMeans(train_polychaeta),
      sd_obs = apply(train_polychaeta, MARGIN = 2, sd),
      row.names = NULL
      ) %>%
    mutate(se_obs = sd_obs / sqrt(180))

    out <- pred %>%
      left_join(obs, by = "species")

    return(out)
  }, .id = "Model") %>%
  filter(species %in% polychaeta) %>%
  mutate(Model = str_remove(Model, "\\s-\\sAB"))


fig_supp18 <- ggplot(
  mean_obs_pred_train,
  aes(
    x = mean_obs,
    y = mean_pred,
  )) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_errorbarh(aes(
    xmin = mean_obs - se_obs,
    xmax = mean_obs + se_obs
    ), height = .2) +
  geom_errorbar(aes(
      ymin = mean_pred - se_pred,
      ymax = mean_pred + se_pred
    ), width = .2) +
  geom_point(aes(colour = Model)) +
    scale_color_manual(
      name = "Model",
      values = rcartocolor::carto_pal(12, "Antique")[1:4]
    ) +
  labs(
    x = "Mean observed abundance",
    y = "Mean Predicted abundance",
    title = "Train dataset"
    ) +
    coord_fixed(ylim = c(0, 150)) +
    facet_wrap(. ~ Model) +
    theme(plot.title = element_text(hjust = .5))

############## Fig S19 - Mean Predicted abundance vs Mean observed abundance on train dataset 

fig_supp19 <- ggplot(
  mean_obs_pred_test,
  aes(
    x = mean_obs,
    y = mean_pred,
  )) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_errorbarh(aes(
    xmin = mean_obs - se_obs,
    xmax = mean_obs + se_obs
    ), height = .2) +
  geom_errorbar(aes(
      ymin = mean_pred - se_pred,
      ymax = mean_pred + se_pred
    ), width = .2) +
  geom_point(aes(colour = Model)) +
    scale_color_manual(
      name = "Model",
      values = rcartocolor::carto_pal(12, "Antique")[1:4]
    ) +
  labs(
    x = "Mean observed abundance",
    y = "Mean Predicted abundance",
    title = "Test Dataset") +
  coord_fixed(ylim = c(0, 150)) +
    facet_wrap(. ~ Model) +
    theme(plot.title = element_text(hjust = .5))

############## Fig S20 - Relationship between relative change in RMSE for WhC model and mean abundance of species in train dataset

test_abd_occ <- test_polychaeta %>%
  pivot_longer(everything(), names_to = "species", values_to = "abundance") %>%
  group_by(species) %>%
  summarise(mean_abd = mean(abundance), mean_occ = mean(abundance > 0)) %>%
  mutate(species = factor(species, levels = polychaeta))

corr_test_df_whole_fauna <- fit_test_ab_relative %>%
  left_join(test_abd_occ, by = "species") %>%
  filter(model == "Whole community")

cor.test(corr_test_df_whole_fauna$relative_change,
  corr_test_df_whole_fauna$mean_abd,
  method = "kendall"
)

fig_supp20 <- ggplot(corr_test_df_whole_fauna, aes(x = mean_abd, y = relative_change)) +
geom_point() +
geom_smooth(formula = 'y~x', method='loess', se = FALSE) +
labs(x = "Mean abundance", y = "Relative change in RMSE of the WhC model compared to the Benchmark model")

############## Fig S21 - Relationship between relative change in RMSE for WhC model and mean occurence of species in train dataset

cor.test(corr_test_df_whole_fauna$relative_change,
  corr_test_df_whole_fauna$mean_occ,
  method = "kendall"
)

fig_supp21 <- ggplot(corr_test_df_whole_fauna, aes(x = mean_occ, y = relative_change)) +
geom_point() +
geom_smooth(formula = 'y~x', method='loess', se = FALSE) +
labs(x = "Mean occurence", y = "Relative change in RMSE of the WhC model compared to the Benchmark model")


############## Fig S22 - Comparison of the model performances to predict the community structure for models fitted with presence/absence data

p1 <- ggplot(
    diss %>% filter(data_type == "PA", !is.na(pred_diss)),
    aes(x = model, y = pred_diss - obs_diss, color = train_test)
) +
stat_pointinterval(aes(group = interaction(model, train_test)),
    position = position_dodge(width = .75),
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Dataset",
    breaks = c("train", "test"),
    labels = c("Train", "Test"),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  labs(
    y = "Difference between predicted and observed Sørensen dissimilarity",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank()
)

p2 <- ggplot(
    richness_models %>% filter(data_type == "PA"),
    aes(x = model, y = diff_richness, color = train_test)
) +
stat_pointinterval(aes(group = interaction(model, train_test)),
    position = position_dodge(width = .75),
    na.rm = TRUE
  ) +
  scale_color_manual(
    name = "Dataset",
    breaks = c("train", "test"),
    labels = c("Train", "Test"),
    values = rcartocolor::carto_pal(12, "Antique")[1:4]
  ) +
  labs(
    y = "Difference between predicted and observed richness",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank()
)

fig_supp22 <- (p1 + p2 & scale_x_discrete(labels = models_acc) &
  theme(
    axis.text.x = element_text(color = rcartocolor::carto_pal(12, "Antique")[1:4]),
    axis.title.x = element_blank()
)) +
  plot_layout(guides = "collect")

############## Fig S23 - Comparison across alternative models of explained variance

# Arrange dataframe according to benchmark model on abundance data
sp_order_vp_ab <- VP %>%
  filter(
    variable == "Environment",
    data_type == "AB",
    model == "Benchmark"
  ) %>%
  arrange(desc(explained_var)) %>%
  .$species %>%
  as.character()

# Arrange dataframe according to benchmark model on occurence data
sp_order_vp_pa <- VP %>%
  filter(
    variable == "Environment",
    data_type == "PA",
    model == "Benchmark"
  ) %>%
  arrange(desc(explained_var)) %>%
  .$species %>%
  as.character()

VP_ab <- VP %>%
  filter(data_type == "AB") %>%
  mutate(
    variable = str_remove(variable, ":\\s[:alpha:]*"),
    species = factor(species, levels = sp_order_vp_ab)
  ) %>%
  group_by(model, data_type, species, variable) %>%
  summarise(value = sum(explained_var), .groups = "drop") %>%
  {
    ggplot(., aes(x = species, y = value, fill = variable)) +
      geom_col(width = max(diff(.$value)) * 1.05) +
      labs(
        x = "",
        y = "Proportion of explained variance"
      ) +
      scale_fill_manual(
        name = "Variable",
        values = rcartocolor::carto_pal(10, "Prism")[c(5, 10)]
      ) +
      facet_grid(model ~ data_type,
        labeller = labeller(data_type = c("AB" = "Abundance"))
      ) +
      theme(
        axis.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
  }

VP_pa <- VP %>%
  filter(data_type == "PA") %>%
  mutate(
    variable = str_remove(variable, ":\\s[:alpha:]*"),
    species = factor(species, levels = sp_order_vp_pa)
  ) %>%
  group_by(model, data_type, species, variable) %>%
  summarise(value = sum(explained_var), .groups = "drop") %>%
  {
    ggplot(., aes(x = species, y = value, fill = variable)) +
      geom_col(width = max(diff(.$value)) * 1.05) +
      labs(
        x = "",
        y = ""
      ) +
      scale_fill_manual(
        name = "Variable",
        values = rcartocolor::carto_pal(10, "Prism")[c(5, 10)]
      ) +
      facet_grid(model ~ data_type,
        labeller = labeller(data_type = c("PA" = "Occurence"))
      ) +
      theme(
        axis.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
  }


fig_supp23 <- VP_ab + VP_pa + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

############## Fig S24 - Comparison across alternative models of total variance

sp_SR2_na <- VP %>%
  filter(data_type == "AB", is.na(percent_total_var)) %>%
  .$species %>%
  unique()

# Arrange dataframe according to benchmark model on abundance data
sp_order_vp_tot_ab <- VP %>%
  filter(
    !(species %in% sp_SR2_na),
    variable == "Environment",
    data_type == "AB",
    model == "Benchmark"
  ) %>%
  arrange(metric) %>%
  .$species %>%
  as.character()

# Arrange dataframe according to benchmark model on occurence data
sp_order_vp_tot_pa <- VP %>%
  filter(
    !(species %in% sp_SR2_na),
    variable == "Environment",
    data_type == "PA",
    model == "Benchmark"
  ) %>%
  arrange(metric) %>%
  .$species %>%
  as.character()

VP_ab_tot <- VP %>%
  filter(!(species %in% sp_SR2_na), data_type == "AB") %>%
  mutate(
    variable = str_remove(variable, ":\\s[:alpha:]*"),
    species = factor(species, levels = sp_order_vp_tot_ab)
  ) %>%
  group_by(model, data_type, species, variable) %>%
  summarise(value = sum(percent_total_var), .groups = "drop") %>%
  {
    ggplot(., aes(x = species, y = value, fill = variable)) +
      geom_col(width = max(diff(.$value)) * 1.05) +
      labs(
        x = "",
        y = "Proportion of variance"
      ) +
      scale_fill_manual(
        name = "Variable",
        values = rcartocolor::carto_pal(10, "Prism")[c(5, 10)]
      ) +
      facet_grid(model ~ data_type,
        labeller = labeller(data_type = c("AB" = "Abundance"))
      ) +
      theme(
        axis.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
  }

VP_pa_tot <- VP %>%
  filter(!(species %in% sp_SR2_na), data_type == "PA") %>%
  mutate(
    variable = str_remove(variable, ":\\s[:alpha:]*"),
    species = factor(species, levels = sp_order_vp_tot_pa)
  ) %>%
  group_by(model, data_type, species, variable) %>%
  summarise(value = sum(percent_total_var), .groups = "drop") %>%
  {
    ggplot(., aes(x = species, y = value, fill = variable)) +
      geom_col(width = max(diff(.$value)) * 1.6) +
      labs(
        x = "",
        y = ""
      ) +
      scale_fill_manual(
        name = "Variable",
        values = rcartocolor::carto_pal(10, "Prism")[c(5, 10)]
      ) +
      facet_grid(model ~ data_type,
        labeller = labeller(data_type = c("PA" = "Occurence"))
      ) +
      theme(
        axis.text.x = element_blank(),
        strip.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
  }


fig_supp24 <- VP_ab_tot + VP_pa_tot + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

############## Fig S25 - Number & proportion of response curves classify accordingly to Rigal et al. (2020)'s methodology for models fitted with presence/absence data

p1 <- ggplot(
  rigal_class %>% filter(
    classif == "Accelerate decline",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(accelerate_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Accelerate decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p2 Convex

p2 <- ggplot(
  rigal_class %>% filter(
    classif == "Convex",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  )

if (nrow(p2$data) > 0) p2 <- p2 + scale_x_discrete(drop = FALSE)
p2 <- p2 +
  annotation_raster(convexe, 0, 1, 300, 500) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Convex") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p3 Accemlerate increase

p3 <- ggplot(
  rigal_class %>% filter(
    classif == "Accelerated increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(accelerate_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Accelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p4 Constant decline

p4 <- ggplot(
  rigal_class %>% filter(
    classif == "Constant decrease",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(constant_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Constant decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p5 Stable

p5 <- ggplot(
  rigal_class %>% filter(
    classif == "Stable",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Stable") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p6 Constant increase

p6 <- ggplot(
  rigal_class %>% filter(
    classif == "Constant increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(constant_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Constant increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p7 Decelerated decline

p7 <- ggplot(
  rigal_class %>% filter(
    classif == "Decelerated decline",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(decelerated_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Decelerated decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p8 Concave

p8 <- ggplot(
  rigal_class %>% filter(
    classif == "Concave",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(concave, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Concave") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p9 Decelerated increase

p9 <- ggplot(
  rigal_class %>% filter(
    classif == "Decelerated increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Model)
) +
  geom_col() +
  geom_text(aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -.3,
    size = 3
  ) +
  annotation_raster(decelerated_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 600)) +
  scale_fill_manual(
    name = "",
    labels = models_acc,
    values = rcartocolor::carto_pal(12, "Antique")[1:4],
    drop = FALSE
  ) +
  ggtitle("Decelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


fig_supp25 <-
  (p1 + p2 + p3 +
    p4 + p5 + p6 +
    p7 + p8 + p9 &
    theme(legend.position = "bottom")) +
  plot_layout(
    ncol = 3, nrow = 3,
    byrow = TRUE,
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = .5))
  )

############## Fig S26 - Number & proportion of response curves classify accordingly to Rigal et al. (2020)'s methodology for models fitted with abundance data. Each bar is coloured by variable.

p1 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Accelerate decline",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Accelerate decline", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(accelerate_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Accelerate decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p2 Convex

p2 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Convex",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  )

if (nrow(p2$data) > 0) p2 <- p2 + scale_x_discrete(drop = FALSE)
p2 <- p2 +
  annotation_raster(convexe, 0, 1, 300, 500) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Convex") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p3 Accemlerate increase

p3 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Accelerated increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Accelerated increase", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(accelerate_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Accelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p4 Constant decline

p4 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Constant decrease",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Constant decrease", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(constant_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Constant decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p5 Stable

p5 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Stable",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Stable", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Stable") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p6 Constant increase

p6 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Constant increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Constant increase", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(constant_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Constant increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p7 Decelerated decline

p7 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Decelerated decline",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(decelerated_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Decelerated decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p8 Concave

p8 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Concave",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(concave, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Concave") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p9 Decelerated increase

p9 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Decelerated increase",
    Data_type == "AB"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Decelerated increase", Data_type == "AB", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(decelerated_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Decelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

fig_supp26 <- (p1 + p2 + p3 +
  p4 + p5 + p6 +
  p7 + p8 + p9 &
  theme(legend.position = "bottom")) +
  plot_layout(
    ncol = 3, nrow = 3,
    byrow = TRUE,
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = .5))
  )

############## Fig S27 - Number & proportion of response curves classify accordingly to Rigal et al. (2020)'s methodology for models fitted with presence/absence data. Each bar is coloured by variable.

p1 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Accelerate decline",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Accelerate decline", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(accelerate_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Accelerate decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p2 Convex

p2 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Convex",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  )

if (nrow(p2$data) > 0) p2 <- p2 + scale_x_discrete(drop = FALSE)
p2 <- p2 +
  annotation_raster(convexe, 0, 1, 300, 500) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Convex") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p3 Accemlerate increase

p3 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Accelerated increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Accelerated increase", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(accelerate_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Accelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p4 Constant decline

p4 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Constant decrease",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Constant decrease", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(constant_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Constant decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p5 Stable

p5 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Stable",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Stable", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Stable") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p6 Constant increase

p6 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Constant increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Constant increase", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(constant_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Constant increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p7 Decelerated decline

p7 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Decelerated decline",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(decelerated_decline, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Decelerated decline") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

# p8 Concave

p8 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Concave",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(concave, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2)) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Concave") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )


# p9 Decelerated increase

p9 <- ggplot(
  rigal_class_var %>% filter(
    classif == "Decelerated increase",
    Data_type == "PA"
  ),
  aes(x = Model, y = n, fill = Variable)
) +
  geom_col() +
  geom_text(
    data = filter(rigal_class_var, classif == "Decelerated increase", Data_type == "PA", round(percent, 3) > 0),
    aes(label = scales::percent(percent, accuracy = 1)),
    position = position_stack(vjust = .5),
    size = 3
  ) +
  annotation_raster(decelerated_increase, 1, 4, 300, 500) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_discrete(drop = FALSE) +
  ggtitle("Decelerated increase") +
  theme(
    plot.title = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  )

fig_supp27 <- (p1 + p2 + p3 +
  p4 + p5 + p6 +
  p7 + p8 + p9 &
  theme(legend.position = "bottom")) +
  plot_layout(
    ncol = 3, nrow = 3,
    byrow = TRUE,
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = .5))
  )

############## Fig S28 - Relationship between species’ position along the first axis of the fuzzy PCA and the different covariates

trait_axis <- list(AB = model_polychaeta_phylo_traits_ab, PA = model_polychaeta_phylo_traits_pa) %>%
  map_dfr(function(x) {
    getPostEstimate(x, "Gamma") %>%
      .$mean %>%
      `rownames<-`(covNames) %>%
      `colnames<-`(str_remove_all(
        model_polychaeta_phylo_traits_ab$trNames,
        "\\(|\\)"
      )) %>%
      as.data.frame() %>%
      rownames_to_column("Var") %>%
      .[, c(1:3)] %>%
      separate("Var", into = c("Var", "power"), sep = "_") %>%
      replace_na(list(power = "1")) %>%
      pivot_wider(names_from = power, values_from = Axis1) %>%
      rename(coeff_x1 = `1`, coeff_x2 = `2`) %>%
      replace_na(list(coeff_x2 = 0, coeff_x1 = 0)) %>%
      filter(Var != "Intercept") %>%
      mutate(Var = factor(Var, levels = c(
        "Fetch", "Organic Mater", "Salinity",
        "Temperature", "Current", "Mud", "Trask Index"
      ))) %>%
      group_by(Var) %>%
      group_modify(function(x, y) {
        coeff <- colSums(x)

        out <- data.frame(x = seq(range(train_fpca[, 2])[1], range(train_fpca[
          ,
          2
        ])[2], length.out = 1000)) %>%
          mutate(y = coeff[1] + coeff[2] * x + coeff[3] * x^2)

        return(out)
      }) %>%
      ungroup()
  }, .id = "data_type")

fig_supp28 <- ggplot(trait_axis, aes(x = x, y = y, color = Var)) +
  geom_line() +
  scale_color_discrete(name = "Variable") +
  labs(y = "Expected species niche", x = "Axis 1 - Fuzzy-PCA") +
  facet_grid(
    cols = vars(data_type),
    scales = "free_y",
    labeller = labeller(data_type = c(AB = "Abundance", PA = "Presence/Absence"))
  ) +
  coord_flip()

############## Fig S29

# Benchmark model

res_corr_bench_pa_df <- residual_corr_raw[[2]] %>%
  map_dfr(function(x) {
    x %>%
      as.data.frame() %>%
      rownames_to_column("species1") %>%
      pivot_longer(-species1, names_to = "species2", values_to = "corr") %>%
      filter(species1 != species2) %>%
      rename(corr_bench = corr)
  }, .id = "random_effect")

# Whole community model
res_corr_whole_pa_df <- residual_corr_raw[[8]] %>%
  map_dfr(function(x) {
    x %>%
      as.data.frame() %>%
      rownames_to_column("species1") %>%
      pivot_longer(-species1, names_to = "species2", values_to = "corr") %>%
      filter(species1 != species2) %>%
      rename(corr_whole = corr)
  }, .id = "random_effect")

# Combine both models
res_corr_comp_pa <- res_corr_bench_pa_df %>%
  left_join(res_corr_whole_pa_df, by = c("random_effect", "species1", "species2")) %>%
  drop_na(corr_bench, corr_whole)

# Compute linear R2 between both
r2_df_pa <- res_corr_comp_pa %>%
  group_by(random_effect) %>%
  group_map(function(x, y) {
    cat(str_to_title(unlist(y)), sep = "\n")

    reg <- lm(corr_whole ~ corr_bench, data = x) %>%
      summary() %>%
      .$adj.r.squared

    data.frame(
      random_effect = unlist(y),
      r2 = reg,
      row.names = NULL
    )
  }) %>%
  map_dfr(~.x)

p_res_corr1_pa <- ggplot(res_corr_comp_pa, aes(x = corr_bench, y = corr_whole)) +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c(name = "Number of neighbors") +
  geom_smooth(
    method = "lm",
    formula = y ~ x
  ) +
  geom_text(
    data = r2_df_pa, x = -.45, y = .98,
    aes(label = paste0("R^2:  ", round(r2, 2))), parse = TRUE
  ) +
  facet_row(~random_effect,
    labeller = labeller(random_effect = c(
      "habitat" = "Habitat",
      "site" = "Site",
      "year" = "Year"
    ))
  ) +
  labs(
    x = "Residual Correlations Bench",
    y = "Residual Correlations WhC"
  )

p_res_corr2_pa <- ggplot(res_corr_comp_pa, aes(x = abs(corr_whole - corr_bench) * (sign(corr_whole * corr_bench)))) +
  geom_histogram(aes(y = after_stat(density)),
    fill = "transparent",
    col = "black",
    bins = 30
  ) +
  geom_density(col = "red") +
  facet_row(~random_effect,
    labeller = labeller(random_effect = c(
      "habitat" = "Habitat",
      "site" = "Site",
      "year" = "Year"
    ))
  ) +
  labs(
    x = "Index value",
    y = "Density",
  ) +
  theme(legend.position = "none")

fig_supp29 <- p_res_corr1_pa /
  p_res_corr2_pa +
  plot_annotation(tag_levels = "A")

res_corr_comp_pa %>%
      mutate(index = abs(corr_whole - corr_bench) * (sign(corr_whole * corr_bench))) %>%
      add_count(random_effect) %>%
      filter(abs(index) > 0.25, index < 0) %>%
      add_count(random_effect) %>%
      group_by(random_effect) %>%
      summarise(percent = unique(nn) / unique(n) * 100)

############## Saving plots
main_figures <- list(
  fig2,
  fig3,
  fig4,
  fig5,
  fig6
)

supp_figures <- list(
  "fig_supp1" = fig_supp1,
  "fig_supp2" = fig_supp2,
  "fig_supp3" = fig_supp3,
  "fig_supp4" = fig_supp4,
  "fig_supp15" = fig_supp15,
  "fig_supp16" = fig_supp16,
  "fig_supp17" = fig_supp17,
  "fig_supp18" = fig_supp18,
  "fig_supp19" = fig_supp19,
  "fig_supp20" = fig_supp20,
  "fig_supp21" = fig_supp21,
  "fig_supp22" = fig_supp22,
  "fig_supp23" = fig_supp23,
  "fig_supp24" = fig_supp24,
  "fig_supp25" = fig_supp25,
  "fig_supp26" = fig_supp26,
  "fig_supp27" = fig_supp27,
  "fig_supp28" = fig_supp28,
  "fig_supp29" = fig_supp29
)

# Main plot
for (i in 1:length(main_figures)) {
  try({
    ggsave(
      plot = main_figures[[i]],
      filename = paste0("Figs/", "fig", i + 1, ".png"),
      width = 12,
      height = 9,
      device = ragg::agg_png,
      scale = 17.78 / 9, # Same width as on DATARMOR plots
      units = "cm"
    )
  })
}

# Appendix
for (i in 1:length(supp_figures)){
  try({
  ggsave(
    plot = supp_figures[[i]],
    filename = paste0("Figs/", names(supp_figures)[i], ".png"),
    width = 12,
    height = 9,
    device = ragg::agg_png,
    scale = 17.78 / 9,
    units = "cm"
  )
  })
}