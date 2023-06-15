library(tidyverse)
library(Hmsc)

load("./Data/Processed/train_model_polychaeta.RData")

# Create sampling design

time <- as.matrix(unique(as.numeric(as.character(train_env$annee))))
rownames(time) <- levels(train_env$annee)
colnames(time) <- "annee"

samplingDesign <- data.frame(annee = train_env$annee,
                             site = train_env$site,
                             habitat = train_env$habitat)

# Building Hmsc model

## Random levels

rL1 <- HmscRandomLevel(sData = time)
rL2 <- HmscRandomLevel(units = unique(samplingDesign$site))
rL3 <- HmscRandomLevel(units = unique(samplingDesign$habitat))

## Regression formula

XFormula <- ~ poly(Average, degree = 2) + poly(MO, degree = 2) + poly(SAL_sd_w, degree = 2) +
  poly(TEMP_sd_w, degree = 2) + poly(CURR_mean_w, degree = 2) + poly(mud, degree = 2) + poly(`Trask(So)`, degree = 2)

## Fauna

train_polychaeta_pa <- train_polychaeta %>%
  vegan::decostand(method = "pa") %>%
  as.matrix()

# Traits Data

## Species column of traits data frame must match the order of the Y column names

TrData <- train_fpca %>%
  mutate(Species = factor(Species, levels = colnames(train_polychaeta_pa))) %>%
  arrange(Species) %>%
  column_to_rownames("Species")

## Traits regression formula

TrFormula <- ~ Axis1 + Axis2 + Axis3

## Model design

model <- Hmsc(Y = train_polychaeta_pa,
              XData = train_env[, -c(1:3)], # Remove random effect variables
              XFormula = XFormula,
              TrData = TrData,
              TrFormula = TrFormula,
              TrScale = FALSE,
              phyloTree = train_taxo,
              distr = "probit",
              studyDesign = samplingDesign,
              ranLevels = list(annee   = rL1,
                               site    = rL2,
                               habitat = rL3))

# MCMC sampling properties

thin      <- 1000
samples   <- 1000
transient <- round(.5 * thin * samples)
nChains   <- 1

(iter_tot <- thin * samples + transient)

# Running the model

set.seed(42)

model_polychaeta_phylo_traits_pa <- sampleMcmc(model,
    samples = samples,
    thin = thin,
    transient = transient,
    nChains = nChains,
    nParallel = nChains,
    updater = list(GammaEta = FALSE),
    verbose = 1000)

saveRDS(model_polychaeta_phylo_traits_pa, file = paste0("Data/Model_Output/model_polychaeta_phylo_traits_pa_", iter_tot/1E3, "K_", job_index, ".rds"))
