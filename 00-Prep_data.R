library(tidyverse)

library(ade4)
library(adegraphics)

library(ape)

load("./Data/Raw/Faune_cleaned_06-09-18.RData")
load("./Data/Raw/Granulo_06-10-18.Rdata")
load("./Data/Raw/Meteo_01jav_03-09-18.Rdata")
load("./Data/Raw/Previmer_integr_01-01_02-09-18.Rdata")
load("./Data/Raw/fetch_cleaned_10-09-2018.Rdata")
load("./Data/Raw/Distance_cleaned_10-09-18.Rdata")
load("./Data/Raw/trait_poly_estim.Rdata")
trait_list <- read_csv2("./Data/Raw/Trait_list.csv")

# Na imputation for granulometrie

granulo_site <- VIM::kNN(granulo_site, k = 5, numFun = median, useImputedDist = FALSE, imp_var = FALSE)

# Habitat as a factor

granulo_site$habitat <- factor(granulo_site$habitat, ordered = FALSE)
previmer_site$habitat <- factor(previmer_site$habitat, ordered = FALSE)

fetch_site <- fetch_site %>%
  ungroup() %>% # add ungroup() because there is a conflict with old dplyr version
  mutate(habitat = factor(habitat, ordered = FALSE))

# Aggregate fauna data per site

faune_agg <- faune %>%
  group_by_at(vars(-point_name, -replicat, -abondance)) %>%
  summarise_at("abondance", sum) %>%
  ungroup()

all <- suppressWarnings(
  faune_agg %>%
    mutate(habitat = factor(habitat, ordered = FALSE)) %>%
    left_join(meteo_01janv, by = c("habitat", "theme", "site", "method_echant", "annee", "date" = "sampling_date")) %>%
    left_join(granulo_site, by = c("habitat", "theme", "site", "method_echant", "annee", "saison")) %>%
    left_join(previmer_site, by = c("habitat", "theme", "site", "method_echant", "annee", "saison")) %>%
    left_join(fetch_site, by = c("habitat", "theme", "site", "method_echant"))
) %>%
  mutate(annee = as.numeric(annee))

whole_comm <- all %>%
  filter(
    habitat %in% c("Intertidal seagrass beds", "Intertidal bare sediments"),
    abondance > 0
  ) %>%
  group_by(annee, site, habitat, species) %>%
  pivot_wider(
    names_from = species,
    values_from = abondance,
    values_fill = list(abondance = 0)
  ) %>%
  ungroup()

# Get a list of polychaeta species

poly <- sp_classif %>%
  filter(Class == "Polychaeta", !is.na(Species)) %>%
  select(Species) %>%
  unlist(use.names = FALSE) %>%
  as.character()

# Add unknow species but genius known

poly_sp <- sp_classif %>%
  filter(Class == "Polychaeta", is.na(Species)) %>%
  select(species_init) %>%
  unlist(use.names = FALSE) %>%
  as.character()

polychaeta <- c(poly, poly_sp) %>%
  as.factor()

env <- whole_comm %>%
  select(
    site, habitat, annee, Average, MO,
    SAL_sd_w, TEMP_sd_w,
    CURR_mean_w, mud, `Trask(So)`
  ) %>%
  mutate(
    site = as.factor(site),
    annee = as.factor(annee),
    habitat = fct_drop(habitat)
  ) # Remove unseen habitats

# Create training datasets

# Normal datasets

## Fauna

train_whole_fauna <- whole_comm %>%
  filter(!site %in% c("Callot", "Sainte-Marguerite")) %>%
  select(-(habitat:Average)) %>%
  select_if(colSums(.) > 0 & colSums(vegan::decostand(., method = "pa")) >= 4) # Remove species not seen in those sites and seen less than 4

train_polychaeta <- whole_comm %>%
  filter(!site %in% c("Callot", "Sainte-Marguerite")) %>%
  select(any_of(polychaeta)) %>%
  select_if(colSums(.) > 0 & colSums(vegan::decostand(., method = "pa")) >= 4) # Remove species not seen in those sites and seen less than 4

test_whole_fauna <- whole_comm %>%
  filter(site %in% c("Callot", "Sainte-Marguerite")) %>%
  select(any_of(colnames(train_whole_fauna)))

test_polychaeta <- whole_comm %>%
  filter(site %in% c("Callot", "Sainte-Marguerite")) %>%
  select(any_of(colnames(train_polychaeta)))

## Env

train_env <- env %>%
  filter(!site %in% c("Callot", "Sainte-Marguerite")) %>%
  mutate(site = fct_drop(site))

test_env <- env %>%
  filter(site %in% c("Callot", "Sainte-Marguerite"))

## Traits

# Remove species not in training set
trait_fuz_estim <- trait_fuz_estim %>%
  tibble::rownames_to_column() %>%
  filter(rowname %in% names(train_polychaeta)) %>%
  tibble::column_to_rownames()

names(trait_fuz_estim)

print(trait_list)
trait_list %>% count(Trait)

trait_block <- c(`Maximum size` = 7, `Feeding method` = 8, `Food size` = 2, `Preferred substrate position` = 2, `Living habitat` = 5, `Adult movement capacity (daily)` = 4, Bioturbation = 5, `Sexual differenciation` = 2, `Development mode` = 4, `Reproduction frequency` = 2)

trait_fuz <- prep.fuzzy.var(trait_fuz_estim, col.blocks = trait_block)

trait_fuz_fpca <- dudi.fpca(trait_fuz, scannf = FALSE, nf = Inf)

screeplot(trait_fuz_fpca)

summary(trait_fuz_fpca) # If we took the first 3 axes we get 59% of the total variance

train_fpca <- trait_fuz_fpca$li[, 1:3] %>%
  tibble::rownames_to_column("Species")

## Taxonomy

formula_taxo <- ~ Class / Order / Family / Genus / species_init

sp_poly <- sp_classif %>%
  filter(species_init %in% names(train_polychaeta)) %>%
  mutate_all(as.character) %>%
  mutate(
    Genus = if_else(is.na(.$Genus), .$species_init, .$Genus),
    Order = if_else(is.na(.$Order), .$Family, .$Order)
  ) %>%
  mutate_all(as.factor)

train_taxo <- as.phylo(formula_taxo, data = sp_poly, collapse = FALSE)
train_taxo$edge.length <- rep(1, times = length(train_taxo$edge))

plot(train_taxo, cex = 0.3)

## Whole community models, 3 RL, no phylogeny, no trait
train_model_whole_fauna <- c("train_whole_fauna", "train_env")
test_model_whole_fauna <- c("test_whole_fauna", "test_env")

save(list = train_model_whole_fauna, file = "./Data/Processed/train_model_whole_fauna.RData")
save(list = test_model_whole_fauna, file = "./Data/Processed/test_model_whole_fauna.RData")

## Polychaeta models
train_model_polychaeta <- c("train_polychaeta", "train_env", "train_taxo", "train_fpca")
test_model_polychaeta <- c("test_polychaeta", "test_env")

save(list = train_model_polychaeta, file = "./Data/Processed/train_model_polychaeta.RData")
save(list = test_model_polychaeta, file = "./Data/Processed/test_model_polychaeta.RData")

# Map

coord_sites <- read.csv("Data/Raw/Coordinates.csv") %>%
  filter(
    site_name %in% unique(whole_comm$site),
    theme %in% unique(whole_comm$theme)
  ) %>%
  group_by(theme, site_name) %>%
  summarise(
    latitude = mean(latitude),
    longitude = mean(longitude),
    .groups = "drop"
  ) %>%
  group_by(site_name) %>%
  group_map(function(x, y) {
    if (nrow(x) == 2) {
      x <- x %>%
        group_by(site_name) %>%
        summarise(
          latitude = mean(latitude),
          longitude = mean(longitude)
        ) %>%
        mutate(theme = "Both", .before = site_name)
    }

    return(x)
  }, .keep = TRUE) %>%
  map_dfr(~.x) %>%
  mutate(
    theme = str_replace(
      theme, "Herbiers Intertidaux", "Intertidal seagrass meadows"
    )
  ) %>%
  mutate(theme = str_replace(
    theme, "Intertidal Meuble", "Intertidal bare sediments"
  )) %>%
  mutate(
    theme = factor(
      theme,
      levels = c(
        "Intertidal seagrass meadows",
        "Intertidal bare sediments",
        "Both"
      )
    )
  ) %>%
  mutate(
    test = if_else(site_name %in% c("Callot", "Sainte-Marguerite"),
      TRUE,
      FALSE
    ),
    .before = latitude
  )

save(
  coord_sites,
  file = "Data/Processed/sites_coordinates.RData"
)
