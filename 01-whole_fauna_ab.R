library(tidyverse)
library(Hmsc)

load("./Data/Processed/train_model_whole_fauna.RData")

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

## Model design

model <- Hmsc(Y = as.matrix(train_whole_fauna),
              XData = train_env[, -c(1:3)], # Remove random effect variables
              XFormula = XFormula,
              distr = "lognormal poisson",
              studyDesign = samplingDesign,
              ranLevels = list(annee   = rL1,
                               site    = rL2,
                               habitat = rL3))

# MCMC sampling properties

thin <- 500
samples <- 1000
transient <- round(.5 * thin * samples)
nChains   <- 5

iter_tot <- thin * samples + transient

# Running the model

set.seed(42)

model_whole_fauna_ab <- sampleMcmc(model,
    samples = samples,
    thin = thin,
    transient = transient,
    nChains = nChains,
    nParallel = nChains,
    updater = list(GammaEta = FALSE),
    verbose = 1
  )

save(model_whole_fauna_ab, file = "Data/Model_Output/model_whole_fauna_ab.RData")
