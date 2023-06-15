library(tidyverse)
library(Hmsc)

################### Loading the data

models <- c(
    "model_polychaeta_ab",
    "model_polychaeta_pa",
    "model_polychaeta_phylo_ab",
    "model_polychaeta_phylo_pa",
    "model_polychaeta_phylo_traits_ab",
    "model_whole_fauna_ab",
    "model_whole_fauna_pa"
)

for (i in seq_along(models)) {
    for (j in 1:5) {
        assign(paste0(models[i], "_", j), readRDS(paste0("Data/Model_Output/", models[i], "_375K_", j, ".rds")))
    }
}

################### Combine individual chains

model_polychaeta_ab <- c(
    model_polychaeta_ab_1,
    model_polychaeta_ab_2,
    model_polychaeta_ab_3,
    model_polychaeta_ab_4,
    model_polychaeta_ab_5
)

model_polychaeta_phylo_ab <- c(
    model_polychaeta_phylo_ab_1,
    model_polychaeta_phylo_ab_2,
    model_polychaeta_phylo_ab_3,
    model_polychaeta_phylo_ab_4,
    model_polychaeta_phylo_ab_5
)

model_polychaeta_phylo_traits_ab <- c(
    model_polychaeta_phylo_traits_ab_1,
    model_polychaeta_phylo_traits_ab_2,
    model_polychaeta_phylo_traits_ab_3,
    model_polychaeta_phylo_traits_ab_4,
    model_polychaeta_phylo_traits_ab_5
)

model_whole_fauna_ab <- c(
    model_whole_fauna_ab_1,
    model_whole_fauna_ab_2,
    model_whole_fauna_ab_3,
    model_whole_fauna_ab_4,
    model_whole_fauna_ab_5
)

model_polychaeta_pa <- c(
    model_polychaeta_pa_1,
    model_polychaeta_pa_2,
    model_polychaeta_pa_3,
    model_polychaeta_pa_4,
    model_polychaeta_pa_5
)

model_polychaeta_phylo_pa <- c(
    model_polychaeta_phylo_pa_1,
    model_polychaeta_phylo_pa_2,
    model_polychaeta_phylo_pa_3,
    model_polychaeta_phylo_pa_4,
    model_polychaeta_phylo_pa_5
)

model_polychaeta_phylo_traits_pa <- c(
    model_polychaeta_phylo_traits_pa_1,
    model_polychaeta_phylo_traits_pa_2,
    model_polychaeta_phylo_traits_pa_3,
    model_polychaeta_phylo_traits_pa_4,
    model_polychaeta_phylo_traits_pa_5
)

model_whole_fauna_pa <- c(
    model_whole_fauna_pa_1,
    model_whole_fauna_pa_2,
    model_whole_fauna_pa_3,
    model_whole_fauna_pa_4,
    model_whole_fauna_pa_5
)

############################# Save the combined chains

saveRDS(model_polychaeta_ab, "Data/Model_Output/model_polychaeta_ab.rds")
saveRDS(model_polychaeta_pa, "Data/Model_Output/model_polychaeta_pa.rds")

saveRDS(model_polychaeta_phylo_ab, "Data/Model_Output/model_polychaeta_phylo_ab.rds")
saveRDS(model_polychaeta_phylo_pa, "Data/Model_Output/model_polychaeta_phylo_pa.rds")

saveRDS(model_polychaeta_phylo_traits_ab, "Data/Model_Output/model_polychaeta_phylo_traits_ab.rds")
saveRDS(model_polychaeta_phylo_traits_pa, "Data/Model_Output/model_polychaeta_phylo_traits_pa.rds")

saveRDS(model_whole_fauna_ab, "Data/Model_Output/model_whole_fauna_ab.rds")
saveRDS(model_whole_fauna_pa, "Data/Model_Output/model_whole_fauna_pa.rds")