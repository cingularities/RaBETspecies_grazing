#Code written by Cynthia L. Norton, University of Arizona, 2021
#code for topoedaphic analysis
#Script written by Cynthia Norton
#Cleaned version: removed unused packages, dead/duplicate code blocks,
#redundant write.csv calls, and commented-out leftovers.

Packages <- c("dplyr", "tidyr", "broom", "tibble", "ggplot2", "tidyverse", "emmeans")
lapply(Packages, library, character.only = TRUE)
setwd("//snre-gaea/projects/RaBET/RaBET_landuse/landuse/")
options(scipen = 100, digits = 4)


###################################################################
### SRER PERCENT COVER
###################################################################
srer_cover <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/srer_landuse_freq_082624.csv') %>%
  dplyr::select(-freq_NaN) %>%
  na.omit()

srer_mesquite   <- srer_cover %>% dplyr::select(-c(cactus, creosote, paloverde, lotebush, bareground, grass)) %>% rename(cover = mesquite)   %>% mutate(species = "mesquite")   %>% na.omit()
srer_cactus     <- srer_cover %>% dplyr::select(-c(mesquite, creosote, paloverde, lotebush, bareground, grass)) %>% rename(cover = cactus)     %>% mutate(species = "cactus")     %>% na.omit()
srer_creosote   <- srer_cover %>% dplyr::select(-c(cactus, mesquite, paloverde, lotebush, bareground, grass))   %>% rename(cover = creosote)   %>% mutate(species = "creosote")   %>% na.omit()
srer_lotebush   <- srer_cover %>% dplyr::select(-c(cactus, creosote, mesquite, paloverde, bareground, grass))   %>% rename(cover = lotebush)   %>% mutate(species = "lotebush")   %>% na.omit()
srer_paloverde  <- srer_cover %>% dplyr::select(-c(cactus, creosote, mesquite, lotebush, bareground, grass))    %>% rename(cover = paloverde)  %>% mutate(species = "paloverde")  %>% na.omit()
srer_bareground <- srer_cover %>% dplyr::select(-c(cactus, creosote, mesquite, paloverde, lotebush, grass))     %>% rename(cover = bareground) %>% mutate(species = "bareground") %>% na.omit()
srer_grass      <- srer_cover %>% dplyr::select(-c(cactus, creosote, mesquite, paloverde, lotebush, bareground)) %>% rename(cover = grass)      %>% mutate(species = "grass")      %>% na.omit()

soil_filter <- function(df) {
  df %>% dplyr::filter(case_when(
    layer == 1 ~ MUSYM == "EbC",
    layer == 2 ~ MUSYM == "An" | MUSYM == "EbC" | MUSYM == "SoB",
    layer == 3 ~ MUSYM == "An" | MUSYM == "SoB",
    layer == 4 ~ MUSYM == "CtB" | MUSYM == "CuC"
  ))
}

srer_cover_rbind <- bind_rows(srer_mesquite, srer_cactus, srer_creosote, srer_lotebush, srer_paloverde) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha >= 300, '1', '0'))) %>%
  soil_filter()

srer_cover_rbind_bg <- bind_rows(srer_bareground, srer_grass) %>%
  na.omit() %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha >= 300, '1', '0'))) %>%
  soil_filter()


###################################################################
### DENSITY
###################################################################
srer_density <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/densities_final_SRER_082624.csv') %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha >= 300, '1', '0'))) %>%
  soil_filter()


###################################################################
### CROWN AREA
###################################################################
srer_mesquite_crownArea  <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/mesquite_crownArea_df_srer_082624.csv')  %>% mutate(species = "mesquite")
srer_cactus_crownArea    <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/cactus_crownArea_df_srer_082624.csv')    %>% mutate(species = "cactus")
srer_creosote_crownArea  <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/creosote_crownArea_df_srer_082624.csv')  %>% mutate(species = "creosote")
srer_lotebush_crownArea  <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/lotebush_crownArea_df_srer_082824.csv')  %>% mutate(species = "lotebush")
srer_paloverde_crownArea <- read.csv('//snre-gaea/projects/RaBET/RaBET_landuse/landuse/paloverde_crownArea_df_srer_082824.csv') %>% mutate(species = "paloverde")

srer_crownarea <- bind_rows(
  srer_mesquite_crownArea, srer_cactus_crownArea, srer_creosote_crownArea,
  srer_lotebush_crownArea, srer_paloverde_crownArea
) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1', '0'))) %>%
  filter(case_when(
    FID_elevat == 0 ~ MUSYM == "EbC",
    FID_elevat == 1 ~ MUSYM == "An" | MUSYM == "EbC" | MUSYM == "SoB",
    FID_elevat == 2 ~ MUSYM == "An" | MUSYM == "SoB",
    FID_elevat == 3 ~ MUSYM == "CtB" | MUSYM == "CuC"
  ))


###################################################################
### SUMMARY TABLE HELPER FUNCTIONS
###################################################################
summarize_cover <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(mean = mean(cover), sd = sd(cover), count = n(),
              se = sd / sqrt(count), upper_limit = mean + se, lower_limit = mean - se,
              .groups = "drop") %>%
    mutate(grouping = paste(group_vars, collapse = "_"))
}

summarize_density <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(mean = mean(density), sd = sd(density), count = n(),
              se = sd / sqrt(count), upper_limit = mean + se, lower_limit = mean - se,
              .groups = "drop") %>%
    mutate(grouping = paste(group_vars, collapse = "_"))
}

summarize_crown <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(mean = mean(crownArea), sd = sd(crownArea), count = n(),
              se = sd / sqrt(count), upper_limit = mean + se, lower_limit = mean - se,
              .groups = "drop") %>%
    mutate(grouping = paste(group_vars, collapse = "_"))
}

cover_summary_all <- bind_rows(
  summarize_cover(srer_cover_rbind, "species"),
  summarize_cover(srer_cover_rbind, c("species", "layer")),
  summarize_cover(srer_cover_rbind, c("MUSYM", "species", "layer"))
)

density_summary_all <- bind_rows(
  summarize_density(srer_density, "species"),
  summarize_density(srer_density, c("species", "layer")),
  summarize_density(srer_density, c("MUSYM", "species", "layer"))
)

crown_summary_all <- bind_rows(
  summarize_crown(srer_crownarea, "species"),
  summarize_crown(srer_crownarea, c("species", "layer")),
  summarize_crown(srer_crownarea, c("MUSYM", "species", "layer"))
)


###################################################################
### PLOTS: SPECIES-LEVEL COVER / DENSITY / CROWN AREA
###################################################################
mean_bareground <- weighted.mean(srer_bareground$cover)
mean_grass      <- weighted.mean(srer_grass$cover)
sd_bareground   <- sd(srer_bareground$cover)
sd_grass        <- sd(srer_grass$cover)

species <- srer_cover_rbind %>%
  group_by(species) %>%
  summarize(mean = mean(cover), sd = sd(cover), count = n(),
            se = sd / sqrt(count), upper_limit = mean + se, lower_limit = mean - se)

mean_wood <- sum(species$mean)

grass <- srer_cover_rbind_bg %>%
  group_by(species, layer, grazing) %>%
  summarize(mean = mean(cover), sd = sd(cover), count = n(),
            se = sd / sqrt(count), upper_limit = mean + sd, lower_limit = mean - sd)

all <- bind_rows(species, grass)

plot_theme <- theme_gray(base_size = 20, base_family = "") +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold", size = 20),
    legend.text  = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text   = element_text(face = "bold")
  )
theme_set(plot_theme)

windows()
ggplot(srer_cover_rbind, aes(x = factor(species), y = cover, fill = species)) +
  scale_fill_manual(values = c("gray", "brown", "lightblue", "lightgreen", "lightpink")) +
  geom_boxplot(alpha = 0.5, size = 1) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, size = 22),
        axis.text.y = element_text(face = "bold", size = 22)) +
  ylab("Mean Cover (%)") + xlab("Species") +
  stat_summary(fun.y = "mean", color = "red", shape = 16, size = 1) +
  geom_line(aes(y = mean_wood, group = 1), size = 1, color = 'red')

# density
species_density <- srer_density %>%
  group_by(species) %>%
  summarize(mean = mean(n), sd = sd(n), count = n(),
            se = sd / sqrt(count), upper_limit = mean + sd, lower_limit = mean - sd)

mean_density <- sum(species_density$mean)

windows()
ggplot(species_density, aes(x = species, y = mean, fill = species)) +
  scale_fill_manual(values = c("gray", "brown", "lightblue", "lightgreen", "lightpink")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, face = "bold", size = 22),
        axis.text.y = element_text(face = "bold", size = 22),
        axis.title  = element_text(face = "bold")) +
  geom_point(shape = 21, size = 6) +
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit),
                width = 0.2, position = position_dodge(0.01), size = 1.2) +
  ylab("Mean Density (stems/ha)")

# crown area
mean_crownarea <- mean(srer_crownarea$crownArea)
sd_crownarea   <- sd(srer_crownarea$crownArea)

species_crownarea <- srer_crownarea %>%
  group_by(species) %>%
  summarize(crownArea_mean = mean(crownArea), sd = sd(crownArea), count = n(),
            se = sd / sqrt(count), upper_limit = crownArea_mean + sd, lower_limit = crownArea_mean - sd) %>%
  mutate(w_mean = crownArea_mean, w_stdev = sd_crownarea, g_mean = mean_grass, g_stdev = sd_grass)

windows()
ggplot(species_crownarea, aes(x = species, y = crownArea_mean, fill = species)) +
  scale_fill_manual(values = c("gray", "brown", "lightblue", "lightgreen", "lightpink")) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit),
                width = 0.2, position = position_dodge(0.01), size = 1.2) +
  geom_line(aes(y = mean_crownarea, group = 1), size = 1, color = 'red') +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, size = 22),
        axis.text.y = element_text(face = "bold", size = 22)) +
  ylab("Mean Crown Area (m²)")


###################################################################
### ANOVA: METRIC ~ ELEVATION, BY SPECIES
###################################################################
anova_cover_species <- srer_cover_rbind %>%
  group_by(species) %>%
  do(tidy(aov(cover ~ factor(layer), data = .))) %>%
  ungroup()

anova_density_species <- srer_density %>%
  group_by(species) %>%
  do(tidy(aov(density ~ factor(layer), data = .))) %>%
  ungroup()

anova_crownarea_species <- srer_crownarea %>%
  group_by(species) %>%
  do(tidy(aov(crownArea ~ factor(layer), data = .))) %>%
  ungroup()

print(anova_cover_species)
print(anova_density_species)
print(anova_crownarea_species)


###################################################################
### GRAZING SUMMARY TABLES (generic helper, reused for every grouping)
### Replaces ~10 near-duplicate copy/pasted blocks in the original script.
###################################################################
summarize_with_grazing <- function(data, value_col, grouping_vars) {
  value <- rlang::sym(value_col)
  
  grazing_test <- data %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(grazing_p = t.test((!!value)[grazing == 0], (!!value)[grazing == 1])$p.value,
              .groups = "drop")
  
  value_summary <- data %>%
    group_by(across(all_of(grouping_vars)), grazing) %>%
    summarise(mean = weighted.mean(!!value, na.rm = TRUE),
              count = n(),
              sd = sd(!!value, na.rm = TRUE),
              .groups = "drop")
  
  value_summary %>% left_join(grazing_test, by = grouping_vars)
}

# ---- Cover, density, crown area summaries at each grouping level ----
cover_grazing_summary_all <- bind_rows(
  summarize_with_grazing(srer_cover_rbind, "cover", "species"),
  summarize_with_grazing(srer_cover_rbind, "cover", c("species", "layer")),
  summarize_with_grazing(srer_cover_rbind, "cover", c("MUSYM", "species", "layer"))
)
write.csv(cover_grazing_summary_all, "cover_grazing_summary_all_110124.csv", row.names = FALSE)

density_grazing_summary_all <- bind_rows(
  summarize_with_grazing(srer_density, "n", "species"),
  summarize_with_grazing(srer_density, "n", c("species", "layer")),
  summarize_with_grazing(srer_density, "n", c("MUSYM", "species", "layer"))
)
write.csv(density_grazing_summary_all, "density_grazing_summary_all_110124.csv", row.names = FALSE)

crownarea_grazing_summary_all <- bind_rows(
  summarize_with_grazing(srer_crownarea, "crownArea", "species"),
  summarize_with_grazing(srer_crownarea, "crownArea", c("species", "layer")),
  summarize_with_grazing(srer_crownarea, "crownArea", c("MUSYM", "species", "layer"))
)
write.csv(crownarea_grazing_summary_all, "crownarea_grazing_summary_all_110124.csv", row.names = FALSE)

# ---- MUSYM x species x layer summaries (used for the bar plots below) ----
species_elev_woody_grazing_soil_cover_summary <-
  summarize_with_grazing(srer_cover_rbind, "cover", c("MUSYM", "species", "layer")) %>%
  rename(cover_mean = mean, cover_count = count, cover_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_soil_cover_summary, "species_elev_woody_grazing_soil_cover_summary_110124.csv", row.names = FALSE)

species_elev_woody_grazing_soil_density_summary <-
  summarize_with_grazing(srer_density, "n", c("MUSYM", "species", "layer")) %>%
  rename(density_mean = mean, density_count = count, density_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_soil_density_summary, "species_elev_woody_grazing_soil_density_summary_110124.csv", row.names = FALSE)

species_elev_woody_grazing_soil_crownarea_summary <-
  summarize_with_grazing(srer_crownarea, "crownArea", c("MUSYM", "species", "layer")) %>%
  rename(crownArea_mean = mean, crownArea_count = count, crownArea_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_soil_crownarea_summary, "species_elev_woody_grazing_soil_crownarea_summary_110124.csv", row.names = FALSE)

# ---- species x layer (no soil) summaries ----
species_elev_woody_grazing_cover_summary     <- summarize_with_grazing(srer_cover_rbind, "cover", c("species", "layer")) %>%
  rename(cover_mean = mean, cover_count = count, cover_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_cover_summary, "species_elev_woody_grazing_cover_summary_110124.csv", row.names = FALSE)

species_elev_woody_grazing_density_summary   <- summarize_with_grazing(srer_density, "n", c("species", "layer")) %>%
  rename(density_mean = mean, density_count = count, density_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_density_summary, "species_elev_woody_grazing_density_summary_110124.csv", row.names = FALSE)

species_elev_woody_grazing_crownarea_summary <- summarize_with_grazing(srer_crownarea, "crownArea", c("species", "layer")) %>%
  rename(crownArea_mean = mean, crownArea_count = count, crownArea_sd = sd, grazing_avg = grazing_p)
write.csv(species_elev_woody_grazing_crownarea_summary, "species_elev_woody_grazing_crownarea_summary_110124.csv", row.names = FALSE)

# ---- layer-only summaries ----
elev_woody_grazing_cover_summary     <- summarize_with_grazing(srer_cover_rbind, "cover", "layer") %>%
  rename(cover_mean = mean, cover_count = count, cover_sd = sd, grazing_avg = grazing_p)
write.csv(elev_woody_grazing_cover_summary, "elev_woody_grazing_cover_summary_110124.csv", row.names = FALSE)

elev_woody_grazing_density_summary   <- summarize_with_grazing(srer_density, "n", "layer") %>%
  rename(density_mean = mean, density_count = count, density_sd = sd, grazing_avg = grazing_p)
write.csv(elev_woody_grazing_density_summary, "elev_woody_grazing_density_summary_110124.csv", row.names = FALSE)

elev_woody_grazing_crownarea_summary <- summarize_with_grazing(srer_crownarea, "crownArea", "layer") %>%
  rename(crownArea_mean = mean, crownArea_count = count, crownArea_sd = sd, grazing_avg = grazing_p)
write.csv(elev_woody_grazing_crownarea_summary, "elev_woody_grazing_crownarea_summary_110124.csv", row.names = FALSE)

# ---- species-only summaries ----
species_woody_grazing_cover_summary     <- summarize_with_grazing(srer_cover_rbind, "cover", "species") %>%
  rename(cover_mean = mean, cover_count = count, cover_sd = sd, grazing_avg = grazing_p)
write.csv(species_woody_grazing_cover_summary, "species_woody_grazing_cover_summary_110124.csv", row.names = FALSE)

species_woody_grazing_density_summary   <- summarize_with_grazing(srer_density, "n", "species") %>%
  rename(density_mean = mean, density_count = count, density_sd = sd, grazing_avg = grazing_p)

species_woody_grazing_crownarea_summary <- summarize_with_grazing(srer_crownarea, "crownArea", "species") %>%
  rename(crownArea_mean = mean, crownArea_count = count, crownArea_sd = sd, grazing_avg = grazing_p)


###################################################################
### SIGNIFICANT GRAZING DIFFERENCES IN CROWN AREA (species x layer x soil)
###################################################################
t_results <- srer_crownarea %>%
  mutate(species = factor(species), grazing = factor(grazing),
         layer = factor(layer), soil = factor(MUSYM)) %>%
  group_by(species, layer, soil) %>%
  summarise(
    p_value = tryCatch(t.test(crownArea ~ grazing)$p.value, error = function(e) NA),
    .groups = "drop"
  )

sig_groups <- t_results %>% filter(!is.na(p_value) & p_value < 0.05)

diff_df <- srer_crownarea %>%
  mutate(species = factor(species), grazing = factor(grazing),
         layer = factor(layer), soil = factor(MUSYM)) %>%
  semi_join(sig_groups, by = c("species", "layer", "soil")) %>%
  group_by(species, layer, soil, grazing) %>%
  summarise(mean_crown = mean(crownArea, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = grazing, values_from = mean_crown, names_prefix = "graz_") %>%
  mutate(diff = graz_1 - graz_0)


###################################################################
### PER-LAYER TABLES FOR BAR PLOTS (cover / density / crown area)
###################################################################
add_plot_fields <- function(df, mean_col, sd_col, p_col, musym_filter) {
  mean_sym <- rlang::sym(mean_col)
  sd_sym   <- rlang::sym(sd_col)
  p_sym    <- rlang::sym(p_col)
  
  df %>%
    mutate(alpha = ifelse(!!p_sym > 0.01, 0.3, 0.8)) %>%
    mutate(across(where(is.numeric), ~ round(., 2))) %>%
    filter(MUSYM %in% musym_filter) %>%
    mutate(landuse = ifelse(grazing < 1, 'protected', 'grazing')) %>%
    mutate(upper_limit = !!mean_sym + !!sd_sym) %>%
    mutate(lower_limit = ifelse(!!mean_sym - !!sd_sym < 0, 0, !!mean_sym - !!sd_sym))
}

musym_by_layer <- list(
  `1` = "EbC",
  `2` = c("An", "EbC", "SoB"),
  `3` = c("An", "SoB"),
  `4` = c("CtB", "CuC")
)

cover_layer1 <- add_plot_fields(filter(species_elev_woody_grazing_soil_cover_summary, layer == 1), "cover_mean", "cover_sd", "grazing_avg", musym_by_layer[["1"]])
cover_layer2 <- add_plot_fields(filter(species_elev_woody_grazing_soil_cover_summary, layer == 2), "cover_mean", "cover_sd", "grazing_avg", musym_by_layer[["2"]])
cover_layer3 <- add_plot_fields(filter(species_elev_woody_grazing_soil_cover_summary, layer == 3), "cover_mean", "cover_sd", "grazing_avg", musym_by_layer[["3"]])
cover_layer4 <- add_plot_fields(filter(species_elev_woody_grazing_soil_cover_summary, layer == 4), "cover_mean", "cover_sd", "grazing_avg", musym_by_layer[["4"]])

density_layer1 <- add_plot_fields(filter(species_elev_woody_grazing_soil_density_summary, layer == 1), "density_mean", "density_sd", "grazing_avg", musym_by_layer[["1"]])
density_layer2 <- add_plot_fields(filter(species_elev_woody_grazing_soil_density_summary, layer == 2), "density_mean", "density_sd", "grazing_avg", musym_by_layer[["2"]])
density_layer3 <- add_plot_fields(filter(species_elev_woody_grazing_soil_density_summary, layer == 3), "density_mean", "density_sd", "grazing_avg", musym_by_layer[["3"]])
density_layer4 <- add_plot_fields(filter(species_elev_woody_grazing_soil_density_summary, layer == 4), "density_mean", "density_sd", "grazing_avg", musym_by_layer[["4"]])

crownarea_layer1 <- add_plot_fields(filter(species_elev_woody_grazing_soil_crownarea_summary, layer == 1), "crownArea_mean", "crownArea_sd", "grazing_avg", musym_by_layer[["1"]])
crownarea_layer2 <- add_plot_fields(filter(species_elev_woody_grazing_soil_crownarea_summary, layer == 2), "crownArea_mean", "crownArea_sd", "grazing_avg", musym_by_layer[["2"]])
crownarea_layer3 <- add_plot_fields(filter(species_elev_woody_grazing_soil_crownarea_summary, layer == 3), "crownArea_mean", "crownArea_sd", "grazing_avg", musym_by_layer[["3"]])
crownarea_layer4 <- add_plot_fields(filter(species_elev_woody_grazing_soil_crownarea_summary, layer == 4), "crownArea_mean", "crownArea_sd", "grazing_avg", musym_by_layer[["4"]])


###################################################################
### BAR PLOTS: COVER / DENSITY / CROWN AREA BY ELEVATION & SOIL
### Same bars/colors/alpha/diff-labels as before, but reorganized so
### each elevation is its own row (matching the difference-plot layout).
###################################################################
library(patchwork)

faceted_theme <- theme_bw(base_size = 16) +
  theme(axis.text.x  = element_text(angle = 90, face = "bold", size = 16),
        axis.text.y  = element_text(face = "bold", size = 16),
        axis.title   = element_text(face = "bold", size = 16),
        strip.text   = element_text(face = "bold", size = 16),
        legend.text  = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16),
        plot.title   = element_text(face = "bold", size = 16))

grazing_labels <- c("0" = "Protected", "1" = "Grazed")

# Helper: for each species x MUSYM x elevation cell, compute Grazed - Protected
# and a y-position (just above the taller bar/error bar) to place the label.
compute_diff_labels <- function(df, value_col) {
  value_sym <- rlang::sym(value_col)
  
  protected <- df %>%
    filter(grazing == 0) %>%
    select(species, MUSYM, layer_label, val0 = !!value_sym, up0 = upper_limit)
  
  grazed <- df %>%
    filter(grazing == 1) %>%
    select(species, MUSYM, layer_label, val1 = !!value_sym, up1 = upper_limit)
  
  inner_join(protected, grazed, by = c("species", "MUSYM", "layer_label")) %>%
    mutate(
      diff  = val1 - val0,
      label = paste0(ifelse(diff >= 0, "+", ""), round(diff, 2)),
      y_pos = pmax(up0, up1, na.rm = TRUE)
    )
}

# Shared plotting helper: one row per elevation, columns = soil types present
# in that elevation (same idea as plot_grazing_soil_diff, applied to bars).
plot_bars_by_elevation <- function(df, diff_labels, value_col, y_label, plot_title) {
  value_sym <- rlang::sym(value_col)
  
  # Pull the "(unit)" portion out of y_label (e.g. "Mean Cover (%)" -> "%")
  # to append to the plot title instead of repeating it on every row.
  unit_match  <- regmatches(y_label, regexpr("\\(([^)]+)\\)", y_label))
  unit_label  <- if (length(unit_match) > 0) gsub("[()]", "", unit_match) else NA
  title_final <- if (!is.na(unit_label)) paste0(plot_title, " (", unit_label, ")") else plot_title
  
  elevations   <- levels(droplevels(df$layer_label))
  y_max        <- max(df$upper_limit, na.rm = TRUE) * 1.15
  alpha_limits <- range(df$alpha, na.rm = TRUE)
  
  row_plots <- lapply(seq_along(elevations), function(i) {
    elev       <- elevations[i]
    is_last    <- i == length(elevations)
    df_sub     <- df %>% filter(layer_label == elev)
    labels_sub <- diff_labels %>% filter(layer_label == elev)
    
    ggplot(df_sub, aes(x = species, y = !!value_sym, fill = factor(grazing))) +
      geom_col(aes(alpha = alpha), position = "identity") +
      geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit, linetype = as.factor(grazing)),
                    width = 0.4, position = position_dodge(width = 1), size = 0.5) +
      geom_text(data = labels_sub, aes(x = species, y = y_pos, label = label),
                inherit.aes = FALSE, vjust = -0.5, fontface = "bold", size = 3) +
      facet_wrap(~ MUSYM, nrow = 1) +
      coord_cartesian(ylim = c(0, y_max)) +
      scale_alpha_continuous(range = c(0.2, 0.8), limits = alpha_limits) +
      scale_fill_discrete(name = "Grazing", labels = grazing_labels) +
      scale_linetype_discrete(name = "Grazing", labels = grazing_labels) +
      labs(x = NULL, y = as.character(elev)) +
      faceted_theme +
      theme(
        plot.title    = element_blank(),
        axis.title.y  = element_text(face = "bold", size = 13, angle = 90, vjust = 0.5),
        axis.text.x   = if (is_last) element_text(angle = 90, face = "bold", size = 14) else element_blank(),
        axis.ticks.x  = if (is_last) element_line() else element_blank(),
        panel.spacing = unit(3, "pt"),
        plot.margin   = margin(t = 2, r = 4, b = 2, l = 4)
      )
  })
  
  wrap_plots(row_plots, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title_final,
      subtitle = "Labels show Grazed − Protected difference",
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.subtitle = element_text(size = 13))
    ) &
    theme(plot.margin = margin(t = 1, r = 4, b = 1, l = 4))
}

# ---- Cover ----
cover_all_layers <- bind_rows(cover_layer1, cover_layer2, cover_layer3, cover_layer4) %>%
  mutate(layer_label = factor(paste0("Elevation ", layer),
                              levels = c("Elevation 1", "Elevation 2", "Elevation 3", "Elevation 4")))
cover_diff_labels <- compute_diff_labels(cover_all_layers, "cover_mean")

windows()
print(plot_bars_by_elevation(cover_all_layers, cover_diff_labels, "cover_mean",
                             "Mean Cover (%)", "Percent Cover by Elevation & Soil Type"))

# ---- Density ----
density_all_layers <- bind_rows(density_layer1, density_layer2, density_layer3, density_layer4) %>%
  mutate(layer_label = factor(paste0("Elevation ", layer),
                              levels = c("Elevation 1", "Elevation 2", "Elevation 3", "Elevation 4")))
density_diff_labels <- compute_diff_labels(density_all_layers, "density_mean")

windows()
print(plot_bars_by_elevation(density_all_layers, density_diff_labels, "density_mean",
                             "Mean Density (stems/ha)", "Density by Elevation & Soil Type"))

# ---- Crown Area ----
crown_all_layers <- bind_rows(crownarea_layer1, crownarea_layer2, crownarea_layer3, crownarea_layer4) %>%
  mutate(layer_label = factor(paste0("Elevation ", layer),
                              levels = c("Elevation 1", "Elevation 2", "Elevation 3", "Elevation 4")))
crown_diff_labels <- compute_diff_labels(crown_all_layers, "crownArea_mean")

windows()
print(plot_bars_by_elevation(crown_all_layers, crown_diff_labels, "crownArea_mean",
                             "Mean Crown Area (m²)", "Crown Area by Elevation & Soil Type"))



###################################################################
### DIFFERENCE PLOT: GRAZING IMPACT ON CROWN AREA, WITH SOIL BREAKDOWN
### High Grazing minus Low/No Grazing, within elevation x soil (MUSYM) x species.
### Fit separately per elevation layer since each layer only contains a
### specific subset of soil types (see soil_filter()).
###################################################################
run_emmeans_grazing_soil <- function(data, response_col) {
  purrr::map_dfr(sort(unique(data$layer)), function(lyr) {
    layer_data <- data %>%
      filter(layer == lyr) %>%
      mutate(
        grazing = factor(grazing, levels = c(0, 1), labels = c("Low/No Grazing", "High Grazing")),
        species = factor(species),
        MUSYM   = factor(MUSYM)
      )
    
    n_musym <- nlevels(layer_data$MUSYM)
    
    tryCatch({
      if (n_musym > 1) {
        # Normal case: multiple soil types at this elevation, MUSYM is a real model term
        fit <- lm(reformulate("grazing * species * MUSYM", response = response_col), data = layer_data)
        emm <- emmeans(fit, ~ grazing | species + MUSYM)
        
        contrast(emm, method = "revpairwise", by = c("species", "MUSYM")) %>%
          as.data.frame() %>%
          mutate(layer = lyr)
      } else {
        # Only ONE soil type at this elevation (e.g. Elevation 1 = EbC only).
        # A single-level factor in the model formula throws
        # "contrasts can be applied only to factors with 2 or more levels",
        # so drop MUSYM from the model and re-attach its (constant) value after.
        fit <- lm(reformulate("grazing * species", response = response_col), data = layer_data)
        emm <- emmeans(fit, ~ grazing | species)
        
        contrast(emm, method = "revpairwise", by = "species") %>%
          as.data.frame() %>%
          mutate(layer = lyr, MUSYM = as.character(levels(layer_data$MUSYM)))
      }
    }, error = function(e) {
      message("run_emmeans_grazing_soil: skipped layer ", lyr, " for '", response_col,
              "' due to error: ", conditionMessage(e))
      NULL
    })
  }) %>%
    filter(!is.na(estimate)) %>%
    mutate(
      layer_label = factor(paste0("Elevation ", layer),
                           levels = c("Elevation 1", "Elevation 2", "Elevation 3", "Elevation 4")),
      sig = case_when(
        p.value < 0.001 ~ "p < 0.001",
        p.value < 0.01  ~ "p < 0.01",
        p.value < 0.05  ~ "p < 0.05",
        TRUE            ~ "n.s."
      )
    )
}

crownarea_soil_contrast_df <- run_emmeans_grazing_soil(srer_crownarea, "crownArea")
write.csv(crownarea_soil_contrast_df, "crownarea_emmeans_grazing_contrast_by_soil_110124.csv", row.names = FALSE)

# Shared plotting helper: one row per elevation, columns = soil types present
# in that elevation. facet_wrap() alone can't guarantee this (it packs panels
# left-to-right regardless of elevation, so elevations with few soil types end
# up sharing a row with the next elevation) -- so each elevation is built as
# its own mini-plot faceted only by soil type, then stacked with patchwork.
library(patchwork)

plot_grazing_soil_diff <- function(contrast_df, y_label, plot_title) {
  # Color = direction of the effect, but ONLY for significant points;
  # non-significant points are grey regardless of direction.
  plot_colors <- c("Higher under Grazing" = "#D55E00", "Higher under Protected" = "#0072B2",
                   "n.s." = "gray70")
  # Shape = whether the difference is statistically significant (filled vs. open point)
  sig_shapes <- c("Significant (p < 0.05)" = 16, "Not significant" = 1)
  
  # Pull the "(unit)" portion out of y_label (e.g. "Difference in Cover (%), High − Low" -> "%")
  # so it can be appended to the elevation label on every row's y-axis.
  unit_match <- regmatches(y_label, regexpr("\\(([^)]+)\\)", y_label))
  unit_label <- if (length(unit_match) > 0) gsub("[()]", "", unit_match) else NA
  
  contrast_df <- contrast_df %>%
    mutate(
      direction   = ifelse(estimate >= 0, "Higher under Grazing", "Higher under Protected"),
      significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not significant"),
      plot_color  = ifelse(significant == "Significant (p < 0.05)", direction, "n.s.")
    )
  
  elevations <- levels(droplevels(contrast_df$layer_label))
  
  # Common y-range across all rows so bars are visually comparable
  y_range <- range(c(contrast_df$estimate - contrast_df$SE * 1.96,
                     contrast_df$estimate + contrast_df$SE * 1.96), na.rm = TRUE)
  
  row_plots <- lapply(seq_along(elevations), function(i) {
    elev    <- elevations[i]
    is_last <- i == length(elevations)
    df_sub  <- contrast_df %>% filter(layer_label == elev)
    
    row_y_label <- if (!is.na(unit_label)) paste0(elev, " (", unit_label, ")") else as.character(elev)
    
    ggplot(df_sub, aes(x = species, y = estimate, color = plot_color, shape = significant)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_errorbar(aes(ymin = estimate - SE * 1.96, ymax = estimate + SE * 1.96), width = 0.3, size = 1) +
      geom_point(size = 3.5, stroke = 1.3) +
      geom_text(aes(label = round(estimate, 2)), vjust = -0.9, size = 3.3, fontface = "bold", show.legend = FALSE) +
      facet_wrap(~ MUSYM, nrow = 1) +
      coord_cartesian(ylim = y_range) +
      scale_color_manual(values = plot_colors, drop = FALSE,
                         breaks = c("Higher under Grazing", "Higher under Protected", "n.s."),
                         labels = c("Higher under Grazing", "Higher under Protected", "Not significant")) +
      scale_shape_manual(values = sig_shapes, drop = FALSE) +
      labs(x = NULL, y = row_y_label,
           color = "Effect", shape = "Significance") +
      faceted_theme +
      theme(
        plot.title    = element_blank(),
        axis.title.y  = element_text(face = "bold", size = 13, angle = 90, vjust = 0.5),
        axis.text.x   = if (is_last) element_text(angle = 90, face = "bold", size = 14) else element_blank(),
        axis.ticks.x  = if (is_last) element_line() else element_blank(),
        panel.spacing = unit(3, "pt"),
        plot.margin   = margin(t = 2, r = 4, b = 2, l = 4)
      )
  })
  
  wrap_plots(row_plots, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = plot_title,
      subtitle = paste0(y_label, " — color = which side is higher, shape = significance"),
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.subtitle = element_text(size = 13))
    ) &
    theme(plot.margin = margin(t = 1, r = 4, b = 1, l = 4))
}


windows()
print(plot_grazing_soil_diff(crownarea_soil_contrast_df, "Difference in Crown Area (m²), High − Low",
                             "Grazing Impact on Crown Area (High Grazing vs. Protected)"))

# ---- Density: same soil-broken-down difference analysis ----
density_soil_contrast_df <- run_emmeans_grazing_soil(srer_density, "n")
write.csv(density_soil_contrast_df, "density_emmeans_grazing_contrast_by_soil_110124.csv", row.names = FALSE)

windows()
print(plot_grazing_soil_diff(density_soil_contrast_df, "Difference in Density (stems/ha), High − Low",
                             "Grazing Impact on Density (High Grazing vs. Protected)"))

# ---- Percent cover: same soil-broken-down difference analysis ----
cover_soil_contrast_df <- run_emmeans_grazing_soil(srer_cover_rbind, "cover")
write.csv(cover_soil_contrast_df, "cover_emmeans_grazing_contrast_by_soil_110124.csv", row.names = FALSE)

windows()
print(plot_grazing_soil_diff(cover_soil_contrast_df, "Difference in Cover (%), High − Low",
                             "Grazing Impact on Percent Cover (High Grazing vs. Protected)"))



###################################################################
### COMBINED FOREST PLOT: GRAZING EFFECT ACROSS ALL THREE METRICS
### One row-group per elevation, one column per metric (crown area /
### density / cover), species on the y-axis -- mirrors a standard
### multi-panel Delta-metric coefficient plot.
###################################################################
metric_levels <- c("Δ Crown Area (m²)", "Δ Density (stems/ha)", "Δ Cover (%)")

prep_forest_df <- function(df, metric_name) {
  df %>%
    mutate(
      metric      = metric_name,
      layer_label = factor(paste0("Elevation ", layer),
                           levels = c("Elevation 1", "Elevation 2", "Elevation 3", "Elevation 4")),
      direction   = ifelse(estimate >= 0, "Higher under Grazing", "Higher under Protected"),
      significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not significant"),
      plot_color  = ifelse(significant == "Significant (p < 0.05)", direction, "n.s.")
    )
}

forest_df <- bind_rows(
  prep_forest_df(crownarea_contrast_df, "Δ Crown Area (m²)"),
  prep_forest_df(density_contrast_df,   "Δ Density (stems/ha)"),
  prep_forest_df(cover_contrast_df,     "Δ Cover (%)")
) %>%
  mutate(metric = factor(metric, levels = metric_levels))


write.csv(forest_df, "grazing_effect_combined_forest_data_110124.csv", row.names = FALSE)

forest_colors <- c("Higher under Grazing" = "#D55E00", "Higher under Protected" = "#0072B2", "n.s." = "gray70")
forest_shapes <- c("Significant (p < 0.05)" = 16, "Not significant" = 1)

windows()
print(
  ggplot(forest_df, aes(x = estimate, y = species, color = plot_color, shape = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_errorbar(aes(xmin = estimate - SE * 1.96, xmax = estimate + SE * 1.96),
                  width = 0.3, orientation = "y") +
    geom_point(size = 3, stroke = 1.2) +
    facet_grid(layer_label ~ metric, scales = "free_x") +
    scale_color_manual(values = forest_colors, drop = FALSE,
                       breaks = c("Higher under Grazing", "Higher under Protected", "n.s."),
                       labels = c("Higher under Grazing", "Higher under Protected", "Not significant")) +
    scale_shape_manual(values = forest_shapes, drop = FALSE) +
    labs(title = "Grazing Effect Across Metrics, by Elevation & Species",
         subtitle = "High Grazing − Protected estimate (95% CI)",
         x = "Estimate", y = NULL, color = "Effect", shape = "Significance") +
    faceted_theme +
    theme(strip.text.y = element_text(angle = 0),
          panel.spacing = unit(4, "pt"))
)




###################################################################
### EXPORT: ALL RESULT TABLES TO A SINGLE EXCEL WORKBOOK
### (means, High-Low differences, p-values -- one tab per table)
###################################################################
library(openxlsx)

# Round all numeric columns for cleaner display in Excel
round_numeric <- function(df, digits = 3) {
  df %>% mutate(across(where(is.numeric), ~ round(., digits)))
}

# Write one data frame to a formatted worksheet: bold header, frozen top row,
# auto column widths.
add_result_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeDataTable(
    wb, sheet_name, round_numeric(df),
    tableStyle = "TableStyleLight9",
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
  )
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
}

wb <- createWorkbook()

# --- Crown area ---
add_result_sheet(wb, "CrownArea_EMMeans",          crownarea_emm_df)
add_result_sheet(wb, "CrownArea_Diff_Species",      crownarea_contrast_df)
add_result_sheet(wb, "CrownArea_Diff_Species_Soil", crownarea_soil_contrast_df)

# --- Density ---
add_result_sheet(wb, "Density_EMMeans",          density_emm_df)
add_result_sheet(wb, "Density_Diff_Species",      density_contrast_df)
add_result_sheet(wb, "Density_Diff_Species_Soil", density_soil_contrast_df)

# --- Percent cover ---
add_result_sheet(wb, "Cover_EMMeans",          cover_emm_df)
add_result_sheet(wb, "Cover_Diff_Species",      cover_contrast_df)
add_result_sheet(wb, "Cover_Diff_Species_Soil", cover_soil_contrast_df)

# --- Combined (all 3 metrics, one table) ---
add_result_sheet(wb, "Combined_Forest_Data", forest_df)

saveWorkbook(wb, "grazing_impact_analysis_summary_110124.xlsx", overwrite = TRUE)
