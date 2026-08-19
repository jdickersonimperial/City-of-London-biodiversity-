#!/usr/bin/env Rscript
#
# gamm_seasonal_analysis_v3.R
# City of London Acoustic Biodiversity Study
#
# Different filename from v2, and writes to a NEW output folder
# (GAMM_Outputs_v3) so nothing in GAMM_Outputs_v2 gets overwritten.
#
# I have not been able to run this against your real data (no access
# to the HPC filesystem), so treat it as a careful draft rather than
# a guaranteed pass. The one thing I can promise: a single failing
# model will no longer take the rest of the script down with it -
# every model fit is wrapped so it logs and moves on instead of
# halting (see safe_fit() in Section 6).
#
# Changes from v2, in the order they're discussed in chat:
#  - New output_dir (v3) and an Errors_Log_v3.txt that records any
#    model that fails, instead of Rscript halting at the first error
#  - Section 3B: a permanent site-coverage check (this is what caught
#    Cleary having a results folder but nothing surviving the
#    confidence threshold, rather than never being processed at all)
#  - Section 7: Activity is refit as a bam() negative-binomial model
#    with an AR1 term, after confirming severe overdispersion
#    (ratio ~649) in the original Poisson fit. The original Poisson
#    fit is kept for comparison, clearly labelled. An outlier check
#    on total_activity is logged before any of this, since a ratio
#    that large is often a handful of extreme rows rather than a
#    smooth pattern - worth checking before trusting either fit.
#  - Section 7B: NDVI decomposed into between-site and within-site
#    components (the pooled coefficient looked like it was averaging
#    over opposite-signed within-site relationships)
#  - Section 7C: Height as a covariate (loaded in v2 but never used),
#    plus a smooth NDVI x Veg_Density interaction for Richness
#  - Section 9B: multiple-comparisons correction across the species
#    NDVI tests, pulled live from the fitted models rather than typed
#    in by hand
#  - Section 11: roof-type vs the site-level residuals. The ANOVA +
#    boxplot version is the one to trust. The GAMM version is kept
#    but simplified to additive-only (no interaction) - the
#    interaction version is what crashed last time (site-constant
#    predictor x site-constant predictor, on top of a per-site random
#    intercept, with only ~13 effective sites - too little
#    information to estimate that many parameters). Even the additive
#    version is exploratory; don't be surprised if it's unstable.
#  - Section 12: the seasonal timeseries plot no longer lets loess
#    interpolate across large data gaps (this is what produced the
#    impossible ~29-species peak in the Charterhouse panel), and adds
#    a proper quantitative test of "do high-vegetation sites show
#    bigger seasonal swings" instead of an eyeball comparison across
#    facets with free_y-scales
#  - Section 13: adds the fitted-effect plot for the activity
#    interaction, which didn't exist in v2 for any interaction model
#
# Still outstanding, unchanged from v2:
#  - 120_fenchurch stays excluded (no Area figure to compute
#    Veg_Density from)
#  - Cleary's confidence-threshold issue is a data/validation problem
#    upstream of this script - see the diagnostic snippet from
#    earlier in this chat. Nothing in this script needs to change
#    once that's sorted; the site-coverage check in 3B will just stop
#    flagging it and it'll flow through the existing merges
#    automatically once it has data - Veg_Volume/Area for Cleary are
#    already in veg_data below.

library(data.table)
# library(tidyverse) replaced with its individual pieces below - the
# tidyverse meta-package's startup banner (via crayon/cli) is what's
# crashing R on this node with an illegal-operation fault, before any
# of this script's own code runs. This avoids loading that code path
# at all, without changing what functions are available.
library(dplyr)
library(stringr)
library(lubridate)
library(readxl)
library(ggplot2)
library(ggrepel)
library(vegan)
library(mgcv)
library(nlme)
library(patchwork)

# ---------------------------------------------------------
# 1. FILE PATHS & SETUP
# ---------------------------------------------------------
results_dir     <- "/rds/general/user/jd1322/home/BirdNET_Results_v3"
thresholds_file <- "/rds/general/user/jd1322/home/BirdNET_Validation_Stratified/threshold_results.csv"
ecology_file    <- "/rds/general/user/jd1322/home/All_London_Sites_Ecology_Data_2024_2026.xlsx"
heights_file    <- "/rds/general/user/jd1322/home/site_heights.csv"
output_dir      <- "/rds/general/user/jd1322/home/GAMM_Outputs_v5"   # NEW - v4 outputs untouched
dir.create(output_dir, showWarnings = FALSE)

error_log <- file.path(output_dir, "Errors_Log_v3.txt")
if (file.exists(error_log)) file.remove(error_log)

# ---------------------------------------------------------
# 2. LOAD & FILTER BIRD DATA
# ---------------------------------------------------------
cat("Loading and filtering data based on validated thresholds...\n")
thresholds <- data.table::fread(thresholds_file)[threshold_for_0.9 != "not achieved"]
thresholds[, threshold := as.numeric(threshold_for_0.9)]

csv_files <- list.files(path = results_dir, pattern = "_detections_.*\\.csv$", recursive = TRUE, full.names = TRUE)
raw_data <- rbindlist(lapply(csv_files, data.table::fread), fill = TRUE)
raw_data <- raw_data[, .(site, date = ymd(date), species = common_name, confidence)]
filtered_birds <- raw_data[thresholds, on = .(species), nomatch = NULL][confidence >= threshold]

# ---------------------------------------------------------
# 3. RICHNESS, TOTAL ACTIVITY & SITE NAME STANDARDISATION
# ---------------------------------------------------------
filtered_birds[, slot_5d := floor_date(date, "5 days")]

name_map <- c(
  "120fenchurchst" = "120_fenchurch", "79charterhouse" = "charterhouse",
  "cannon_bridge"  = "cannon", "christchurch"   = "chrstchurch",
  "inner_temple"   = "temple_shape_test", "stdunstan"      = "st_dunstans",
  "walbrook_wharf" = "walbrook"
)
filtered_birds[, site_match := tolower(site)]
filtered_birds[, site_match := ifelse(site_match %in% names(name_map), name_map[site_match], site_match)]

slot_richness <- filtered_birds[, .(richness = uniqueN(species)), by = .(site_match, slot_5d)]
slot_activity <- filtered_birds[, .(total_activity = .N), by = .(site_match, slot_5d)]

cat("\nSpecies name check - confirm these match what you expect:\n")
print(unique(filtered_birds$species)[grepl("swift|peregrine|redstart|sparrow",
                                            unique(filtered_birds$species), ignore.case = TRUE)])

species_of_interest <- c("Peregrine Falcon", "Black Redstart", "House Sparrow", "Common Swift")

# ---------------------------------------------------------
# 3B. SITE COVERAGE CHECK
# Would have caught Cleary automatically instead of by chance -
# distinguishes "no raw detections at all" from "has detections but
# none survive the confidence threshold" from "excluded later for a
# different reason (e.g. missing vegetation Area)".
# ---------------------------------------------------------
expected_sites <- c("chrstchurch", "st_dunstans", "120_fenchurch", "temple_shape_test",
                     "cleary", "barber_surgeons", "cannon", "aldgate_school",
                     "nomura", "wood_street", "guildhall", "weil",
                     "charterhouse", "mansion_house", "walbrook")

sites_with_raw_folders <- tolower(basename(list.dirs(results_dir, recursive = FALSE)))
sites_after_filtering  <- unique(filtered_birds$site_match)

coverage_report <- capture.output({
  cat("--- Site coverage check ---\n")
  cat("Expected sites with no BirdNET_Results_v3 folder at all:\n")
  print(setdiff(expected_sites, sites_with_raw_folders))
  cat("Expected sites with a folder but nothing surviving the confidence threshold:\n")
  print(setdiff(expected_sites, sites_after_filtering))
})
cat(paste(coverage_report, collapse = "\n"), "\n")
writeLines(coverage_report, file.path(output_dir, "Site_Coverage_Check.txt"))

# ---------------------------------------------------------
# 4. ECOLOGY, HEIGHT & VEGETATION DENSITY
# ---------------------------------------------------------
cat("\nMerging ecology and vegetation data...\n")
ecology_clean <- as.data.table(openxlsx::read.xlsx(ecology_file, detectDates = TRUE))
# readxl::read_excel() crashed here (illegal-operation fault inside
# vctrs's column-name repair, called via tibble construction) -
# openxlsx::read.xlsx() reads straight into a data.frame with no
# tibble/vctrs involved anywhere, avoiding that code path entirely.
# detectDates=TRUE gets Date back out as an actual Date rather than
# an Excel serial number - if the Date column comes back as numeric
# instead, use as.Date(Date, origin = "1899-12-30") to convert it.
ecology_clean <- ecology_clean[, .(site = layer, date = ymd(Date), NDVI, Cloud_Percent)][Cloud_Percent <= 15]
ecology_clean[, slot_5d := floor_date(date, "5 days")]
ecology_clean[, site_match := tolower(site)]
ecology_clean <- ecology_clean[, .(NDVI = mean(NDVI, na.rm = TRUE)), by = .(site_match, slot_5d)]

heights_df <- data.table::fread(heights_file)
heights_df[, site_match := tolower(Site)]

veg_data <- data.table(
  site_match = c("chrstchurch", "st_dunstans", "120_fenchurch", "temple_shape_test",
                 "cleary", "barber_surgeons", "cannon", "aldgate_school",
                 "nomura", "wood_street", "guildhall", "weil",
                 "charterhouse", "mansion_house", "walbrook"),
  Veg_Volume = c(9553.20, 5107.77, NA, 70866.42, 2852.32, 7625.14, 970.94,
                 68.89, 214.86, 113.72, 0.71, 0.57, 0, 0, 0),
  Area = c(2052.92, 1211.44, NA, 18281.81, 1137.27, 3454.42, 2961.31,
           226.58, 1965.57, 1456.34, 102.74, 264.79, 94.88, 96.83, 417.41)
)

# Roof-type classification, for the Section 11 comparison against the
# site-level residuals. Reconstructed from your GIS classification -
# double check this matches. Cleary is included so it's ready to use
# once its detection data comes through.
roof_type_lookup <- data.table(
  site_match = c("nomura", "wood_street", "guildhall",
                 "120_fenchurch", "cannon", "aldgate_school", "weil",
                 "chrstchurch", "st_dunstans", "temple_shape_test", "cleary", "barber_surgeons",
                 "charterhouse", "mansion_house", "walbrook"),
  roof_type = c("extensive", "extensive", "extensive",
                "intensive", "intensive", "intensive", "intensive",
                "ground_garden", "ground_garden", "ground_garden", "ground_garden", "ground_garden",
                "none", "none", "none")
)

# ---------------------------------------------------------
# 5. BUILD MODEL_DATA
# ---------------------------------------------------------
model_data <- merge(ecology_clean, slot_richness, by = c("site_match", "slot_5d"), all.x = FALSE)
model_data <- merge(model_data, slot_activity, by = c("site_match", "slot_5d"), all.x = TRUE)
model_data[is.na(total_activity), total_activity := 0]
model_data <- merge(model_data, heights_df[, .(site_match, Height)], by = "site_match", all.x = FALSE)
model_data <- merge(model_data, veg_data, by = "site_match", all.x = FALSE)

model_data <- model_data[!is.na(Veg_Volume)]
model_data[, DOY := yday(slot_5d)]
model_data[, Veg_Density := Veg_Volume / Area]
model_data <- model_data[complete.cases(model_data[, .(NDVI, Veg_Density, DOY, site_match, slot_5d)])]
model_data[, site_match := factor(site_match)]

dup_check <- model_data[, .N, by = .(site_match, slot_5d)][N > 1]
if (nrow(dup_check) > 0) {
  cat("\nWARNING: duplicate site/slot rows still present -\n")
  print(dup_check)
}

ndvi_veg_cor <- cor(model_data$NDVI, model_data$Veg_Density)
cat("\nNDVI vs Veg_Density correlation:", round(ndvi_veg_cor, 3), "\n")
writeLines(paste("NDVI vs Veg_Density correlation:", round(ndvi_veg_cor, 3)),
           file.path(output_dir, "NDVI_VegDensity_Collinearity.txt"))

# ---------------------------------------------------------
# 5B. TEMPORAL ORDERING
# Needed for the AR1 term in the new negative-binomial activity
# models (Section 7) and for the gap-aware timeseries plot
# (Section 12). gap_days/segment flag large breaks in a site's
# recording history so a smooth doesn't get drawn across them.
# ---------------------------------------------------------
model_data <- model_data[order(site_match, slot_5d)]
model_data[, gap_days := c(0, diff(as.numeric(slot_5d))), by = site_match]
model_data[, ar_start := c(TRUE, rep(FALSE, .N - 1)), by = site_match]
model_data[, segment := cumsum(gap_days > 60), by = site_match]   # 60 days is a guess - tune to your real deployment gaps if you know them

# ---------------------------------------------------------
# 6. REUSABLE GAMM TEMPLATES + DEFENSIVE FITTING HELPER
# safe_fit() means one failing model logs an error and returns NULL
# instead of halting the whole script - every downstream use checks
# for NULL before trying to use the result.
# ---------------------------------------------------------
fit_gamm <- function(response, data, family = gaussian(), k = 8, use_correlation = TRUE) {
  form <- as.formula(paste0(response, " ~ NDVI + Veg_Density + s(DOY, bs = 'cc', k = ", k, ")"))
  if (use_correlation) {
    gamm(form,
         random = list(site_match = ~1),
         correlation = corCAR1(form = ~ as.numeric(slot_5d) | site_match),
         knots = list(DOY = c(0, 366)),
         family = family,
         data = data)
  } else {
    gamm(form,
         random = list(site_match = ~1),
         knots = list(DOY = c(0, 366)),
         family = family,
         data = data)
  }
}

fit_gamm_custom <- function(response, predictors, data, family = gaussian(), k = 8, use_correlation = TRUE) {
  form <- as.formula(paste0(response, " ~ ", predictors, " + s(DOY, bs = 'cc', k = ", k, ")"))
  if (use_correlation) {
    gamm(form,
         random = list(site_match = ~1),
         correlation = corCAR1(form = ~ as.numeric(slot_5d) | site_match),
         knots = list(DOY = c(0, 366)),
         family = family,
         data = data)
  } else {
    gamm(form,
         random = list(site_match = ~1),
         knots = list(DOY = c(0, 366)),
         family = family,
         data = data)
  }
}

safe_fit <- function(label, fit_expr) {
  tryCatch({
    result <- fit_expr
    cat("OK -", label, "\n")
    result
  }, error = function(e) {
    msg <- paste0("[", format(Sys.time()), "] FAILED - ", label, ": ", conditionMessage(e))
    cat(msg, "\n")
    write(msg, file = error_log, append = TRUE)
    NULL
  })
}

# ---------------------------------------------------------
# 7. RICHNESS & TOTAL ACTIVITY MODELS
# ---------------------------------------------------------
cat("\nRunning Richness GAMM (k=10, WITHOUT corCAR1 - tested and confirmed unnecessary:",
    "estimated phi = 5.3e-08, LRT p = 0.9996 - see test_all_autocorrelation.R output)...\n")
m_richness_gam <- safe_fit("Richness GAMM", fit_gamm("richness", model_data, k = 10, use_correlation = FALSE))
if (!is.null(m_richness_gam)) {
  capture.output(summary(m_richness_gam$gam), file = file.path(output_dir, "GAMM_Richness.txt"))
  capture.output(gam.check(m_richness_gam$gam), file = file.path(output_dir, "GAM_Check_Richness.txt"))
}

cat("Running Total Activity GAMM (original Poisson, kept for comparison)...\n")
m_activity_gam <- safe_fit("Activity GAMM (Poisson)", fit_gamm("total_activity", model_data, family = poisson(), k = 10))
if (!is.null(m_activity_gam)) {
  capture.output(summary(m_activity_gam$gam), file = file.path(output_dir, "GAMM_Activity_Poisson.txt"))
  overdisp <- sum(residuals(m_activity_gam$gam, type = "pearson")^2) / m_activity_gam$gam$df.residual
  cat("Poisson overdispersion ratio:", round(overdisp, 2), "\n")
  writeLines(paste("Poisson overdispersion ratio:", round(overdisp, 2),
                    "\n(this is why a negative-binomial refit is used below instead of trusting this one)"),
             file.path(output_dir, "Activity_Overdispersion_Ratio.txt"))
}

# Outlier check, logged before refitting anything - a ratio in the
# hundreds is often a handful of extreme rows rather than a smooth
# pattern, worth knowing either way.
outlier_report <- capture.output({
  cat("Summary of total_activity:\n"); print(summary(model_data$total_activity))
  cat("\nTop 15 highest activity site/slot rows:\n")
  print(model_data[order(-total_activity)][1:15, .(site_match, slot_5d, total_activity)])
})
writeLines(outlier_report, file.path(output_dir, "Activity_Outlier_Check.txt"))
cat("\nSaved Activity_Outlier_Check.txt - worth a look before trusting any Activity model.\n")

# Negative-binomial refit via bam(), with its own AR1 term, in place
# of the corCAR1/gamm() approach used elsewhere. A small rho grid
# search picks the autocorrelation parameter by AIC rather than
# guessing a fixed value.
cat("\nRefitting Activity as negative-binomial (bam), searching for rho...\n")
rho_grid <- seq(0, 0.9, by = 0.1)
rho_search <- sapply(rho_grid, function(r) {
  m <- safe_fit(paste("Activity NB rho =", r),
                bam(total_activity ~ NDVI + Veg_Density + s(DOY, bs = "cc", k = 10) + s(site_match, bs = "re"),
                    knots = list(DOY = c(0, 366)), family = nb(), rho = r,
                    AR.start = model_data$ar_start, data = model_data))
  if (is.null(m)) NA_real_ else AIC(m)
})
best_rho <- if (all(is.na(rho_search))) {
  cat("All rho values failed to fit - falling back to rho = 0\n"); 0
} else {
  rho_grid[which.min(rho_search)]
}
cat("Best rho by AIC:", best_rho, "\n")

m_activity_nb <- safe_fit("Activity NB final fit",
  bam(total_activity ~ NDVI + Veg_Density + s(DOY, bs = "cc", k = 10) + s(site_match, bs = "re"),
      knots = list(DOY = c(0, 366)), family = nb(), rho = best_rho,
      AR.start = model_data$ar_start, data = model_data))
if (!is.null(m_activity_nb)) {
  capture.output(summary(m_activity_nb), file = file.path(output_dir, "GAMM_Activity_NB.txt"))
}

# Interaction models - Poisson kept for comparison, NB is the one to trust
model_data[, NDVI_c := scale(NDVI, scale = FALSE)]
model_data[, Veg_Density_c := scale(Veg_Density, scale = FALSE)]

m_act_gam_int <- safe_fit("Activity Poisson interaction",
  fit_gamm_custom("total_activity", "NDVI_c * Veg_Density_c", model_data, family = poisson(), k = 10))
if (!is.null(m_act_gam_int)) {
  capture.output(summary(m_act_gam_int$gam), file = file.path(output_dir, "GAMM_Activity_Poisson_Interaction.txt"))
}

m_act_nb_int <- safe_fit("Activity NB interaction",
  bam(total_activity ~ NDVI_c * Veg_Density_c + s(DOY, bs = "cc", k = 10) + s(site_match, bs = "re"),
      knots = list(DOY = c(0, 366)), family = nb(), rho = best_rho,
      AR.start = model_data$ar_start, data = model_data))
if (!is.null(m_act_nb_int)) {
  capture.output(summary(m_act_nb_int), file = file.path(output_dir, "GAMM_Activity_NB_Interaction.txt"))
}

# ---------------------------------------------------------
# 7B. NDVI DECOMPOSED INTO BETWEEN-SITE AND WITHIN-SITE COMPONENTS
# Site_NDVI_Trend.csv showed within-site NDVI-richness correlations
# from -0.50 to +0.52 across sites - a single pooled NDVI coefficient
# is averaging over opposite-signed relationships.
# ---------------------------------------------------------
model_data[, NDVI_site_mean := mean(NDVI), by = site_match]
model_data[, NDVI_within := NDVI - NDVI_site_mean]

m_richness_decomp <- safe_fit("Richness NDVI decomposed",
  gamm(richness ~ NDVI_site_mean + NDVI_within + Veg_Density + s(DOY, bs = "cc", k = 10),
       random = list(site_match = ~1),
       knots = list(DOY = c(0, 366)), data = model_data))
if (!is.null(m_richness_decomp)) {
  capture.output(summary(m_richness_decomp$gam), file = file.path(output_dir, "GAMM_Richness_NDVI_Decomposed.txt"))
}

m_richness_randslope <- safe_fit("Richness NDVI random slope (heavier model - may not converge)",
  gamm(richness ~ NDVI_site_mean + NDVI_within + Veg_Density + s(DOY, bs = "cc", k = 10),
       random = list(site_match = ~1 + NDVI_within),
       knots = list(DOY = c(0, 366)), data = model_data))
if (!is.null(m_richness_randslope)) {
  capture.output(summary(m_richness_randslope$gam), file = file.path(output_dir, "GAMM_Richness_NDVI_RandomSlope.txt"))
}

# ---------------------------------------------------------
# 7C. HEIGHT COVARIATE + SMOOTH NDVI x Veg_Density INTERACTION
# Height was loaded and merged in v2 but never actually used.
# ---------------------------------------------------------
m_richness_height <- safe_fit("Richness with Height",
  fit_gamm_custom("richness", "NDVI + Veg_Density + Height", model_data, k = 10, use_correlation = FALSE))
if (!is.null(m_richness_height)) {
  capture.output(summary(m_richness_height$gam), file = file.path(output_dir, "GAMM_Richness_with_Height.txt"))
}

m_richness_te <- safe_fit("Richness tensor interaction",
  gamm(richness ~ te(NDVI, Veg_Density, k = c(5, 5)) + s(DOY, bs = "cc", k = 10),
       random = list(site_match = ~1),
       knots = list(DOY = c(0, 366)), data = model_data))
# Note: the p=0.007 result already in the write-up was fit WITH
# corCAR1 in place - re-check GAMM_Richness_TensorInteraction.txt
# after this change and update the write-up if the coefficient/p-value
# shifted meaningfully (expect it to be similar, given corCAR1 turned
# out to contribute almost nothing to the base Richness model, but
# this hasn't been independently confirmed for this specific variant).
if (!is.null(m_richness_te)) {
  capture.output(summary(m_richness_te$gam), file = file.path(output_dir, "GAMM_Richness_TensorInteraction.txt"))
}

# ---------------------------------------------------------
# 8. COMMUNITY COMPOSITION (PCoA)
# ---------------------------------------------------------
cat("\n--- COMMUNITY COMPOSITION (PCoA) ---\n")

sp_counts <- filtered_birds[, .(abundance = .N), by = .(site_match, slot_5d, species)]
sp_wide   <- dcast(sp_counts, site_match + slot_5d ~ species, value.var = "abundance", fill = 0)

aligned_data <- merge(model_data, sp_wide, by = c("site_match", "slot_5d"), all.x = FALSE)

metadata_cols <- names(model_data)
species_cols  <- setdiff(names(aligned_data), metadata_cols)

species_matrix_full <- as.matrix(aligned_data[, ..species_cols])
species_occurrence  <- colSums(species_matrix_full > 0)
keep_species <- names(species_occurrence[species_occurrence >= 3])
cat(length(species_cols) - length(keep_species), "rare species (<3 occurrences) excluded from ordination\n")
species_matrix <- species_matrix_full[, keep_species]

cat("Calculating Bray-Curtis distance and running PCoA...\n")
bray_dist <- vegdist(species_matrix, method = "bray")
pcoa_res  <- cmdscale(bray_dist, k = 2, eig = TRUE, add = TRUE)

explained_var <- round(pcoa_res$eig[1:2] / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
aligned_data[, PCoA1 := pcoa_res$points[, 1]]
aligned_data[, PCoA2 := pcoa_res$points[, 2]]

cat("Running GAMMs on PCoA Axes (PCoA1 at k=10 - v2 flagged it as too constrained)...\n")
m_pcoa1_gam <- safe_fit("PCoA1 GAMM", fit_gamm("PCoA1", aligned_data, k = 10))
m_pcoa2_gam <- safe_fit("PCoA2 GAMM", fit_gamm("PCoA2", aligned_data, k = 10))

if (!is.null(m_pcoa1_gam) && !is.null(m_pcoa2_gam)) {
  sink(file.path(output_dir, "GAMM_PCoA_Results.txt"))
  cat("PCoA1 (", explained_var[1], "% variance)\n"); print(summary(m_pcoa1_gam$gam))
  cat("\nPCoA2 (", explained_var[2], "% variance)\n"); print(summary(m_pcoa2_gam$gam))
  sink()
}
if (!is.null(m_pcoa1_gam)) capture.output(gam.check(m_pcoa1_gam$gam), file = file.path(output_dir, "GAM_Check_PCoA1.txt"))
if (!is.null(m_pcoa2_gam)) capture.output(gam.check(m_pcoa2_gam$gam), file = file.path(output_dir, "GAM_Check_PCoA2.txt"))

# ---------------------------------------------------------
# 9. SPECIES-SPECIFIC MODELS (5-day, auto-falls-back to monthly)
# Same logic as v2, but each fit is defensive now and models are kept
# in a list so Section 9B can pull p-values back out programmatically.
# ---------------------------------------------------------
filtered_birds[, month_yr := floor_date(date, "month")]
model_data[, month_yr := floor_date(slot_5d, "month")]
monthly_eco <- model_data[, .(NDVI = mean(NDVI, na.rm = TRUE),
                               Veg_Density = mean(Veg_Density, na.rm = TRUE),
                               DOY = mean(DOY, na.rm = TRUE)),
                           by = .(site_match, month_yr)]

species_models <- list()

for (sp in species_of_interest) {
  cat("\n==============================\nAnalysing:", sp, "\n==============================\n")

  sp_slot <- filtered_birds[species == sp, .(n_detections = .N), by = .(site_match, slot_5d)]
  sp_data <- merge(model_data[, .(site_match, slot_5d, NDVI, Veg_Density, DOY)],
                    sp_slot, by = c("site_match", "slot_5d"), all.x = TRUE)
  sp_data[is.na(n_detections), n_detections := 0]
  sp_data[, detected := as.integer(n_detections > 0)]

  cat("5-day slots - detections:", sum(sp_data$n_detections),
      "| slots with detection:", sum(sp_data$detected), "/", nrow(sp_data), "\n")

  if (sum(sp_data$detected) >= 10) {
    m_sp <- safe_fit(paste(sp, "5-day GAMM"), fit_gamm("detected", sp_data, family = binomial()))
    if (!is.null(m_sp)) {
      capture.output(summary(m_sp$gam), file = file.path(output_dir, paste0("GAMM_", gsub(" ", "_", sp), "_5day.txt")))
      species_models[[sp]] <- m_sp
      cat("5-day model saved.\n")
    }
    next
  }

  cat("Too few 5-day detections - trying monthly aggregation instead...\n")
  sp_month <- filtered_birds[species == sp, .(n_detections = .N), by = .(site_match, month_yr)]
  sp_month_data <- merge(monthly_eco, sp_month, by = c("site_match", "month_yr"), all.x = TRUE)
  sp_month_data[is.na(n_detections), n_detections := 0]
  sp_month_data[, detected := as.integer(n_detections > 0)]

  cat("Monthly - months with detection:", sum(sp_month_data$detected), "/", nrow(sp_month_data), "\n")

  if (sum(sp_month_data$detected) >= 10) {
    m_sp_month <- safe_fit(paste(sp, "monthly GAMM"),
      gamm(detected ~ NDVI + Veg_Density + s(DOY, bs = "cc", k = 5),
           random = list(site_match = ~1),
           knots = list(DOY = c(0, 366)),
           family = binomial(),
           data = sp_month_data))
    if (!is.null(m_sp_month)) {
      capture.output(summary(m_sp_month$gam), file = file.path(output_dir, paste0("GAMM_", gsub(" ", "_", sp), "_monthly.txt")))
      species_models[[sp]] <- m_sp_month
      cat("Monthly model saved.\n")
    }
  } else {
    cat("Still too few detections even monthly - reporting raw counts only.\n")
    print(sp_data[, .(total_detections = sum(n_detections), slots_detected = sum(detected)), by = site_match])
    fwrite(sp_data[, .(total_detections = sum(n_detections), slots_detected = sum(detected)), by = site_match],
           file.path(output_dir, paste0("RawCounts_", gsub(" ", "_", sp), ".csv")))
  }
}

# ---------------------------------------------------------
# 9B. MULTIPLE-COMPARISONS CORRECTION ACROSS SPECIES NDVI TESTS
# Pulled live from whichever species models actually fit, rather than
# typed in by hand.
# ---------------------------------------------------------
species_ndvi_p <- sapply(species_models, function(m) {
  pt <- summary(m$gam)$p.table
  if ("NDVI" %in% rownames(pt)) pt["NDVI", "Pr(>|t|)"] else NA_real_
})
species_ndvi_p <- species_ndvi_p[!is.na(species_ndvi_p)]
if (length(species_ndvi_p) > 0) {
  bh_adjusted <- round(p.adjust(species_ndvi_p, method = "BH"), 3)
  cat("\nBH-adjusted species NDVI p-values:\n"); print(bh_adjusted)
  capture.output(print(bh_adjusted), file = file.path(output_dir, "Species_NDVI_BH_Adjusted.txt"))
}

# ---------------------------------------------------------
# 10. PANEL REGRESSION PLOTS
# ---------------------------------------------------------
cat("\nGenerating panel regression plots...\n")

make_panel <- function(data, yvar, xvar, cyclic = FALSE) {
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.35, color = "grey40") +
    theme_classic() +
    labs(x = xvar, y = yvar)
  if (cyclic) {
    p <- p + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc"),
                          method.args = list(knots = list(x = c(0, 366))),
                          color = "darkred")
  } else {
    p <- p + geom_smooth(method = "lm", color = "darkred")
  }
  p
}

panel_plot <- (make_panel(aligned_data, "PCoA1", "NDVI") |
               make_panel(aligned_data, "PCoA1", "Veg_Density") |
               make_panel(aligned_data, "PCoA1", "DOY", cyclic = TRUE)) /
              (make_panel(aligned_data, "PCoA2", "NDVI") |
               make_panel(aligned_data, "PCoA2", "Veg_Density") |
               make_panel(aligned_data, "PCoA2", "DOY", cyclic = TRUE))

ggsave(file.path(output_dir, "PCoA_Panel_Regressions.png"), panel_plot, width = 13, height = 8, dpi = 300)

# ---------------------------------------------------------
# 11. SITE-LEVEL COMPARISON + ROOF TYPE
# ---------------------------------------------------------
cat("\nGenerating site-level trend comparison...\n")

site_means <- model_data[, .(
  mean_richness = mean(richness, na.rm = TRUE),
  Veg_Density = mean(Veg_Density),
  n_obs = .N
), by = site_match]

fit_between <- lm(mean_richness ~ Veg_Density, data = site_means)
site_means[, predicted := predict(fit_between)]
site_means[, residual := mean_richness - predicted]
fwrite(site_means[order(-residual)], file.path(output_dir, "Site_VegDensity_Residuals.csv"))

veg_residual_plot <- ggplot(site_means, aes(x = Veg_Density, y = mean_richness, label = site_match)) +
  geom_point(size = 3, color = "#2C5F2D") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  geom_text_repel(size = 3) +
  theme_classic() +
  labs(title = "Which sites over/under-perform relative to vegetation density?",
       subtitle = "Above the line: richer than expected. Below: poorer than expected.")
ggsave(file.path(output_dir, "Site_VegDensity_Residuals.png"), veg_residual_plot, width = 9, height = 7, dpi = 300)

site_ndvi_trend <- model_data[, .(
  cor_ndvi = cor(NDVI, richness, use = "complete.obs"), n_obs = .N
), by = site_match][order(-cor_ndvi)]
fwrite(site_ndvi_trend, file.path(output_dir, "Site_NDVI_Trend.csv"))

# Roof type vs the residuals - the ANOVA + boxplot is the version to
# trust. The additive GAMM is exploratory only.
site_means <- merge(site_means, roof_type_lookup, by = "site_match", all.x = TRUE)
site_means[, roof_type := factor(roof_type)]

roof_anova <- safe_fit("Roof type ANOVA on richness residuals", aov(residual ~ roof_type, data = site_means))
if (!is.null(roof_anova)) {
  capture.output(summary(roof_anova), file = file.path(output_dir, "RoofType_Residual_ANOVA.txt"))
}

roof_boxplot <- ggplot(site_means, aes(x = roof_type, y = residual, color = roof_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2.5) +
  geom_text_repel(aes(label = site_match), size = 3, show.legend = FALSE) +
  theme_classic() +
  labs(y = "Richness residual (observed - expected from Veg_Density)",
       title = "Does over/underperformance track roof type?")
ggsave(file.path(output_dir, "RoofType_Residual_Boxplot.png"), roof_boxplot, width = 8, height = 6, dpi = 300)

# Exploratory only, additive not interaction (the interaction version
# is what crashed - a site-constant predictor crossed with another
# site-constant predictor, on top of a per-site random intercept,
# with only ~13 effective sites, is too little information for that
# many parameters). Even this simplified version may not converge
# cleanly - that's expected, not a sign something else is broken.
model_data <- merge(model_data, roof_type_lookup, by = "site_match", all.x = TRUE)
model_data[, roof_type := factor(roof_type)]
m_richness_roof <- safe_fit("Richness + roof_type (additive, exploratory)",
  fit_gamm_custom("richness", "NDVI + Veg_Density + roof_type", model_data, k = 10, use_correlation = FALSE))
# IMPORTANT: this is the model behind the headline roof-type result in
# the write-up (roof_type coefficients, R-sq = 0.635). Re-check
# GAMM_Richness_RoofType.txt against those numbers after this change
# and update the write-up if they moved - this hasn't been
# independently re-verified for this specific model variant, only for
# the simpler base Richness model.
if (!is.null(m_richness_roof)) {
  capture.output(summary(m_richness_roof$gam), file = file.path(output_dir, "GAMM_Richness_RoofType.txt"))
}

# ---------------------------------------------------------
# 12. SEASONAL TIME SERIES BY SITE (gap-aware) + AMPLITUDE TEST
# ---------------------------------------------------------
cat("Generating seasonal amplitude plot, ordered by vegetation density...\n")

site_veg_order <- unique(model_data[, .(site_match, Veg_Density)])[order(-Veg_Density)]
model_data[, site_match := factor(site_match, levels = as.character(site_veg_order$site_match))]

# Fixed version: split into segments wherever there's a >60-day gap
# so loess doesn't interpolate across it (this is what produced the
# impossible ~29-species peak for Charterhouse), and free both axes
# so each panel shows its own real date range.
seasonal_ts <- ggplot(model_data, aes(x = slot_5d, y = richness, group = interaction(site_match, segment))) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", span = 0.3, color = "darkgreen", se = FALSE) +
  facet_wrap(~ site_match, ncol = 3, scales = "free") +
  theme_bw() +
  labs(x = "Date", y = "Species Richness",
       title = "Smoothed richness over time, ordered high \u2192 low vegetation density")
ggsave(file.path(output_dir, "Smoothed_Timeseries_OrderedByVeg.png"), seasonal_ts, width = 14, height = 10, dpi = 300)

# Formal test of "do high-vegetation sites show bigger seasonal
# swings", replacing the free_y-scale eyeball comparison.
site_amplitude <- model_data[, {
  amp_result <- tryCatch({
    lo <- loess(richness ~ as.numeric(slot_5d), span = 0.3)
    fv <- predict(lo)
    list(amplitude = max(fv) - min(fv), sd_fitted = sd(fv))
  }, error = function(e) list(amplitude = NA_real_, sd_fitted = NA_real_))
  .(amplitude = amp_result$amplitude, sd_fitted = amp_result$sd_fitted, Veg_Density = unique(Veg_Density))
}, by = site_match]
site_amplitude <- site_amplitude[!is.na(amplitude)]

fit_amplitude <- safe_fit("Seasonal amplitude vs Veg_Density", lm(amplitude ~ Veg_Density, data = site_amplitude))
if (!is.null(fit_amplitude)) {
  capture.output(summary(fit_amplitude), file = file.path(output_dir, "Seasonal_Amplitude_vs_VegDensity.txt"))
}

amplitude_plot <- ggplot(site_amplitude, aes(x = Veg_Density, y = amplitude, label = site_match)) +
  geom_point(size = 3, color = "#2C5F2D") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  geom_text_repel(size = 3) +
  theme_classic() +
  labs(title = "Does seasonal richness amplitude scale with vegetation density?",
       y = "Seasonal amplitude (max - min fitted richness)")
ggsave(file.path(output_dir, "Seasonal_Amplitude_vs_VegDensity.png"), amplitude_plot, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------
# 13. GAM FITTED-EFFECT PLOTS
# ---------------------------------------------------------
cat("\nGenerating GAM fitted-effect plots...\n")

# For gamm()-based models (Richness, PCoA) - unchanged from v2
plot_gam_effect <- function(model, data, xvar, response, x_label, family = gaussian()) {
  grid <- data.frame(
    NDVI = mean(data$NDVI, na.rm = TRUE),
    Veg_Density = mean(data$Veg_Density, na.rm = TRUE),
    DOY = mean(data$DOY, na.rm = TRUE)
  )
  x_seq <- if (xvar == "DOY") seq(1, 366, length.out = 100) else
    seq(min(data[[xvar]], na.rm = TRUE), max(data[[xvar]], na.rm = TRUE), length.out = 100)
  grid <- grid[rep(1, length(x_seq)), ]
  grid[[xvar]] <- x_seq

  pred <- predict(model$gam, newdata = grid, se.fit = TRUE)
  linkinv <- family$linkinv
  grid$fit   <- linkinv(pred$fit)
  grid$lower <- linkinv(pred$fit - 1.96 * pred$se.fit)
  grid$upper <- linkinv(pred$fit + 1.96 * pred$se.fit)

  ggplot() +
    geom_point(data = data, aes(x = .data[[xvar]], y = .data[[response]]), alpha = 0.15, color = "grey40") +
    geom_ribbon(data = grid, aes(x = .data[[xvar]], ymin = lower, ymax = upper), fill = "#97BC62", alpha = 0.35) +
    geom_line(data = grid, aes(x = .data[[xvar]], y = fit), color = "#2C5F2D", linewidth = 1.1) +
    theme_classic() +
    labs(x = x_label, y = response)
}

# For bam()-based models (Activity NB) - separate function because
# bam objects aren't a $gam-wrapped list, and the random intercept
# term needs excluding to get a population-level curve
plot_gam_effect_bam <- function(model, data, xvar, response, x_label) {
  grid <- data.frame(
    NDVI = mean(data$NDVI, na.rm = TRUE),
    Veg_Density = mean(data$Veg_Density, na.rm = TRUE),
    DOY = mean(data$DOY, na.rm = TRUE),
    site_match = data$site_match[1]
  )
  x_seq <- if (xvar == "DOY") seq(1, 366, length.out = 100) else
    seq(min(data[[xvar]], na.rm = TRUE), max(data[[xvar]], na.rm = TRUE), length.out = 100)
  grid <- grid[rep(1, length(x_seq)), ]
  grid[[xvar]] <- x_seq

  pred <- predict(model, newdata = grid, se.fit = TRUE, exclude = "s(site_match)")
  grid$fit   <- exp(pred$fit)
  grid$lower <- exp(pred$fit - 1.96 * pred$se.fit)
  grid$upper <- exp(pred$fit + 1.96 * pred$se.fit)

  ggplot() +
    geom_point(data = data, aes(x = .data[[xvar]], y = .data[[response]]), alpha = 0.15, color = "grey40") +
    geom_ribbon(data = grid, aes(x = .data[[xvar]], ymin = lower, ymax = upper), fill = "#97BC62", alpha = 0.35) +
    geom_line(data = grid, aes(x = .data[[xvar]], y = fit), color = "#2C5F2D", linewidth = 1.1) +
    theme_classic() +
    labs(x = x_label, y = response)
}

if (!is.null(m_richness_gam)) {
  richness_effects <- (plot_gam_effect(m_richness_gam, model_data, "NDVI", "richness", "NDVI") |
                        plot_gam_effect(m_richness_gam, model_data, "Veg_Density", "richness", "Vegetation Density (m3/m2)") |
                        plot_gam_effect(m_richness_gam, model_data, "DOY", "richness", "Day of Year")) +
    plot_annotation(title = "GAM Fitted Effects \u2014 Species Richness")
  ggsave(file.path(output_dir, "GAM_Effects_Richness.png"), richness_effects, width = 13, height = 4.5, dpi = 300)
}

if (!is.null(m_activity_nb)) {
  activity_effects_nb <- (plot_gam_effect_bam(m_activity_nb, model_data, "NDVI", "total_activity", "NDVI") |
                           plot_gam_effect_bam(m_activity_nb, model_data, "Veg_Density", "total_activity", "Vegetation Density (m3/m2)") |
                           plot_gam_effect_bam(m_activity_nb, model_data, "DOY", "total_activity", "Day of Year")) +
    plot_annotation(title = "GAM Fitted Effects (Negative Binomial) \u2014 Total Acoustic Activity")
  ggsave(file.path(output_dir, "GAM_Effects_Activity_NB.png"), activity_effects_nb, width = 13, height = 4.5, dpi = 300)
}
if (!is.null(m_activity_gam)) {
  activity_effects_poisson <- (plot_gam_effect(m_activity_gam, model_data, "NDVI", "total_activity", "NDVI", poisson()) |
                                plot_gam_effect(m_activity_gam, model_data, "Veg_Density", "total_activity", "Vegetation Density (m3/m2)", poisson()) |
                                plot_gam_effect(m_activity_gam, model_data, "DOY", "total_activity", "Day of Year", poisson())) +
    plot_annotation(title = "GAM Fitted Effects (Poisson - overdispersed, for comparison only) \u2014 Total Acoustic Activity")
  ggsave(file.path(output_dir, "GAM_Effects_Activity_Poisson_ForComparison.png"), activity_effects_poisson, width = 13, height = 4.5, dpi = 300)
}

if (!is.null(m_pcoa1_gam)) {
  pcoa1_effects <- (plot_gam_effect(m_pcoa1_gam, aligned_data, "NDVI", "PCoA1", "NDVI") |
                     plot_gam_effect(m_pcoa1_gam, aligned_data, "Veg_Density", "PCoA1", "Vegetation Density (m3/m2)") |
                     plot_gam_effect(m_pcoa1_gam, aligned_data, "DOY", "PCoA1", "Day of Year")) +
    plot_annotation(title = "GAM Fitted Effects \u2014 PCoA Axis 1")
  ggsave(file.path(output_dir, "GAM_Effects_PCoA1.png"), pcoa1_effects, width = 13, height = 4.5, dpi = 300)
}

if (!is.null(m_pcoa2_gam)) {
  pcoa2_effects <- (plot_gam_effect(m_pcoa2_gam, aligned_data, "NDVI", "PCoA2", "NDVI") |
                     plot_gam_effect(m_pcoa2_gam, aligned_data, "Veg_Density", "PCoA2", "Vegetation Density (m3/m2)") |
                     plot_gam_effect(m_pcoa2_gam, aligned_data, "DOY", "PCoA2", "Day of Year")) +
    plot_annotation(title = "GAM Fitted Effects \u2014 PCoA Axis 2")
  ggsave(file.path(output_dir, "GAM_Effects_PCoA2.png"), pcoa2_effects, width = 13, height = 4.5, dpi = 300)
}

# The interaction figure that didn't exist in v2 for any interaction
# model - shows the crossover directly instead of just a coefficient.
#
# Uses mean Veg_Density per roof-type group as the three representative
# levels, rather than 10th/50th/90th percentiles - Veg_Density is
# heavily right-skewed, so percentiles put "Low" and "Median" almost
# on top of each other (visible in the original plot), while roof type
# is the grouping that actually structures this variable elsewhere in
# the analysis and gives genuinely distinct, meaningfully-labelled
# levels instead.
#
# Also adds exclude.too.far() masking (the same tool used for the
# tensor interaction plot, applied here to three lines instead of a
# full grid) - each line only draws over the NDVI range actually
# observed near real data at that Veg_Density level, rather than
# extrapolating across the full NDVI range regardless of whether sites
# at that vegetation density ever reached it.
if (!is.null(m_act_nb_int)) {
  veg_by_roof <- model_data[, .(mean_veg = mean(Veg_Density)), by = roof_type][order(mean_veg)]
  cat("Mean Veg_Density by roof type, used as the three interaction plot levels:\n")
  print(veg_by_roof)

  ndvi_seq <- seq(min(model_data$NDVI), max(model_data$NDVI), length.out = 100)

  pred_grid <- do.call(rbind, lapply(seq_len(nrow(veg_by_roof)), function(i) {
    data.frame(
      NDVI_c = ndvi_seq - mean(model_data$NDVI),
      Veg_Density_c = veg_by_roof$mean_veg[i] - mean(model_data$Veg_Density),
      DOY = mean(model_data$DOY),
      site_match = model_data$site_match[1],
      NDVI = ndvi_seq,
      Veg_Density = veg_by_roof$mean_veg[i],
      roof_type_label = as.character(veg_by_roof$roof_type[i])
    )
  }))
  pred_grid$veg_level <- factor(pred_grid$roof_type_label, levels = veg_by_roof$roof_type)

  pred_grid$too_far <- exclude.too.far(pred_grid$NDVI, pred_grid$Veg_Density,
                                        model_data$NDVI, model_data$Veg_Density, dist = 0.1)
  pred_grid <- pred_grid[!pred_grid$too_far, ]

  pred <- predict(m_act_nb_int, newdata = pred_grid, se.fit = TRUE, exclude = "s(site_match)")
  pred_grid$fit   <- exp(pred$fit)
  pred_grid$lower <- exp(pred$fit - 1.96 * pred$se.fit)
  pred_grid$upper <- exp(pred$fit + 1.96 * pred$se.fit)

  interaction_plot <- ggplot(pred_grid, aes(x = NDVI, y = fit, color = veg_level, fill = veg_level)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    theme_classic() +
    labs(y = "total_activity (fitted)", color = "Roof type", fill = "Roof type",
         title = "NDVI's effect on activity depends on vegetation density (negative-binomial model)",
         subtitle = "Lines only shown where NDVI x Veg_Density combinations were actually observed for that roof type")
  ggsave(file.path(output_dir, "GAM_Effects_Activity_NB_Interaction.png"), interaction_plot, width = 8, height = 5.5, dpi = 300)
}

# ---------------------------------------------------------
# 14. QUICK EXTRA: CANNON'S EARLY CRASH-AND-REBOUND
# Not a model fit, just pulling the raw rows so you can eyeball
# whether it's a recorder fault, a short/duplicated batch of files,
# or a genuine disturbance event.
# ---------------------------------------------------------
cat("\nCannon - earliest raw detections:\n")
print(filtered_birds[site_match == "cannon"][order(date)][1:15])
cat("\nCannon - earliest model_data rows:\n")
print(model_data[site_match == "cannon"][order(slot_5d)][1:10])

cat("\nAll done. Outputs saved to:", output_dir, "\n")
cat("Check", error_log, "for anything that failed along the way.\n")
