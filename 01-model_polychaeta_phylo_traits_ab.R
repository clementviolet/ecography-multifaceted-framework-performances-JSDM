library(tidyverse)
library(Hmsc)

load("./Data/Processed/train_model_polychaeta.RData")

job_index <- Sys.getenv("PBS_ARRAY_INDEX") %>%
  as.numeric()

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

# Traits Data

## Species column of traits data frame must match the order of the Y column names

TrData <- train_fpca %>%
  mutate(Species = factor(Species, levels = colnames(train_polychaeta))) %>%
  arrange(Species) %>%
  column_to_rownames("Species")

## Traits regression formula

TrFormula <- ~ Axis1 + Axis2 + Axis3

## Model design

model <- Hmsc(Y = as.matrix(train_polychaeta),
              XData = train_env[, -c(1:3)], # Remove random effect variables
              XFormula = XFormula,
              TrData = TrData,
              TrFormula = TrFormula,
              TrScale = FALSE,
              phyloTree = train_taxo,
              distr = "lognormal poisson",
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

set.seed(job_index)

model_polychaeta_phylo_traits_ab <- sampleMcmc(model,
                        samples = samples,
                        thin = thin,
                        transient = transient,
                        nChains = nChains,
                        nParallel = nChains,
                        updater = list(GammaEta = FALSE),
                        verbose = 1000)

saveRDS(model_polychaeta_phylo_traits_ab, file = paste0("Data/Model_Output/model_polychaeta_phylo_traits_ab_", iter_tot/1E3, "K_", job_index, ".rds"))
