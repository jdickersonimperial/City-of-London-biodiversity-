# =========================================================================
# extract_statistics.R
# Produces everything the Results text and figure captions currently lack:
#   1. Summary numbers for the Results opening paragraph and captions
#   2. Confidence intervals for every parametric term
#   3. Back-transformed (multiplicative) effects for log-link models
#   4. Percentile-based effect sizes (10th -> 90th percentile)
#   5. Concurvity - how far NDVI, vegetation density and s(DOY) overlap
#   6. k.check - are the smooth basis dimensions adequate?
#   7. Sensitivity refits: with a recording-effort offset, and with height
#   8. Spatial autocorrelation across sites (Moran's I, permutation test)
#
# Everything is printed to the console AND written to
#   tables_final/statistics_digest.txt
#
#   source("figure_setup.R")
#   source("extract_statistics.R")
# =========================================================================

if (!exists("lab")) stop("Run source('figure_setup.R') first.")
stopifnot(exists("model_data"))

D  <- add_roof_label(as.data.table(model_data))
SC <- col_of(D, "site"); RC <- col_of(D, "roof")
NC <- col_of(D, "ndvi"); VC <- col_of(D, "veg"); YC <- col_of(D, "doy")
RICH <- col_of(D, "richness"); ACT <- col_of(D, "activity")

digest_file <- file.path(tab_dir, "statistics_digest.txt")
con <- file(digest_file, open = "wt")
sink(con, split = TRUE)   # print to console and file simultaneously

rule <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p)))

# =========================================================================
rule("1. SUMMARY NUMBERS (for the Results opening paragraph and captions)")
# =========================================================================
cat("Sites analysed:              ", uniqueN(D[[SC]]), "\n")
cat("5-day windows:               ", nrow(D), "\n")
cat("Windows per site (min/med/max):",
    paste(D[, .N, by = c(SC)][, c(min(N), as.integer(median(N)), max(N))], collapse = " / "), "\n")
cat("Total detections:            ", format(sum(D[[ACT]], na.rm = TRUE), big.mark = ","), "\n")
cat("Sites per roof type:\n"); print(D[, .(sites = uniqueN(get(SC))), by = roof_label])
cat("\nNDVI:            ", sprintf("%.3f to %.3f (mean %.3f)",
    min(D[[NC]], na.rm = TRUE), max(D[[NC]], na.rm = TRUE), mean(D[[NC]], na.rm = TRUE)), "\n")
cat("Vegetation density:", sprintf("%.2f to %.2f (mean %.2f); %d unique values",
    min(D[[VC]], na.rm = TRUE), max(D[[VC]], na.rm = TRUE), mean(D[[VC]], na.rm = TRUE),
    uniqueN(round(D[[VC]], 4))), "\n")
cat("Richness:        ", sprintf("%d to %d (mean %.1f)",
    min(D[[RICH]], na.rm = TRUE), max(D[[RICH]], na.rm = TRUE), mean(D[[RICH]], na.rm = TRUE)), "\n")
cat("Activity:        ", sprintf("%d to %d (median %.0f, mean %.1f)",
    min(D[[ACT]], na.rm = TRUE), max(D[[ACT]], na.rm = TRUE),
    median(D[[ACT]], na.rm = TRUE), mean(D[[ACT]], na.rm = TRUE)), "\n")
WC <- col_of(D, "window")
if (!is.na(WC)) cat("Distinct NDVI/window dates:  ", uniqueN(D[[WC]]), "\n")
cat("\nNOTE: vegetation density has only one value per site, so the tensor\n")
cat("smooth te(NDVI, Veg_Density) is effectively evaluated at",
    uniqueN(round(D[[VC]], 4)), "discrete density\n")
cat("levels. Say this in the Fig 6 caption - it explains the banding.\n")

# =========================================================================
rule("2-3. PARAMETRIC TERMS WITH 95% CI AND BACK-TRANSFORMED EFFECTS")
# =========================================================================
report_terms <- function(m, label) {
  g <- as_gam(m); if (is.null(g)) return(invisible(NULL))
  s <- summary(g); lnk <- g$family$link; z <- qnorm(0.975)
  cat("\n---", label, sprintf("[%s, %s link]", g$family$family, lnk), "---\n")
  if (!is.null(s$p.table) && nrow(s$p.table)) {
    pt <- as.data.frame(s$p.table)
    for (i in seq_len(nrow(pt))) {
      e <- pt[i, 1]; se <- pt[i, 2]; lo <- e - z * se; hi <- e + z * se
      cat(sprintf("  %-34s b = %8.4f  SE = %6.4f  95%% CI [%8.4f, %8.4f]  %s = %6.3f  P = %s\n",
                  rownames(pt)[i], e, se, lo, hi, colnames(pt)[3], pt[i, 3], fmt_p(pt[i, 4])))
      if (lnk %in% c("log", "logit"))
        cat(sprintf("  %-34s %s = %7.4f  95%% CI [%7.4f, %7.4f]\n", "",
                    if (lnk == "log") "exp(b)" else "odds ratio",
                    exp(e), exp(lo), exp(hi)))
    }
  }
  if (!is.null(s$s.table) && nrow(s$s.table)) {
    stb <- as.data.frame(s$s.table)
    for (i in seq_len(nrow(stb)))
      cat(sprintf("  %-34s edf = %5.2f  Ref.df = %5.2f  %s = %7.2f  P = %s\n",
                  rownames(stb)[i], stb[i, 1], stb[i, 2], colnames(stb)[3],
                  stb[i, 3], fmt_p(stb[i, 4])))
  }
  cat(sprintf("  Model: n = %s | adj. R2 = %s | deviance explained = %s\n",
              if (!is.null(s$n)) s$n else nrow(g$model),
              if (!is.null(s$r.sq)) sprintf("%.3f", s$r.sq) else "NA",
              if (!is.null(s$dev.expl)) sprintf("%.1f%%", 100 * s$dev.expl) else "NA"))
}
all_keys <- c("richness_main", "activity_main", "richness_roof", "activity_roof",
              "richness_tensor", "activity_interact", "richness_within",
              "activity_within", "pcoa1", "pcoa2")
for (k in all_keys) {
  m <- get_model(k, quiet = TRUE)
  if (!is.null(m)) report_terms(m, k)
}
if (exists("species_models"))
  for (sp in names(species_models)) report_terms(species_models[[sp]], paste("focal:", sp))

# =========================================================================
rule("4. PERCENTILE-BASED EFFECT SIZES (10th -> 90th percentile)")
# =========================================================================
cat("Use these in the Results instead of 'per unit increase'. A one-unit\n")
cat("change in vegetation density spans most of the observed range, so a\n")
cat("per-unit statement is not interpretable.\n")

pct_effect <- function(m, data, xvar, label) {
  g <- as_gam(m); if (is.null(g)) return(invisible(NULL))
  q <- quantile(data[[xvar]], c(0.10, 0.90), na.rm = TRUE)
  nd <- newdata_from(data, vary = setNames(list(as.numeric(q)), xvar))
  pr <- predict_ci(m, nd)
  lnk <- g$family$link
  cat(sprintf("\n  %s | %s: %.3f -> %.3f\n", label, xvar, q[1], q[2]))
  cat(sprintf("    fitted: %.2f [%.2f, %.2f]  ->  %.2f [%.2f, %.2f]\n",
              pr$fit[1], pr$lo[1], pr$hi[1], pr$fit[2], pr$lo[2], pr$hi[2]))
  if (lnk == "log")
    cat(sprintf("    ratio: x%.2f (a %+.0f%% change in %s)\n",
                pr$fit[2] / pr$fit[1], 100 * (pr$fit[2] / pr$fit[1] - 1), label))
  else
    cat(sprintf("    difference: %+.2f\n", pr$fit[2] - pr$fit[1]))
}
for (k in c("richness_main", "activity_main")) {
  m <- get_model(k, quiet = TRUE); if (is.null(m)) next
  for (xv in c(VC, NC)) pct_effect(m, D, xv, k)
}

# =========================================================================
rule("5. CONCURVITY - do NDVI, vegetation density and s(DOY) overlap?")
# =========================================================================
cat("Interpretation: 'worst' > 0.8 means a term is largely reproducible from\n")
cat("the others, so its individual effect is poorly identified. If NDVI shows\n")
cat("high concurvity with s(DOY), that is a substantive part of the\n")
cat("explanation for the null NDVI main effect - report it and discuss it.\n")
for (k in all_keys) {
  m <- get_model(k, quiet = TRUE); if (is.null(m)) next
  cc <- tryCatch(concurvity(as_gam(m), full = TRUE), error = function(e) NULL)
  if (is.null(cc)) next
  cat("\n---", k, "(full) ---\n"); print(round(cc, 3))
  cw <- tryCatch(concurvity(as_gam(m), full = FALSE), error = function(e) NULL)
  if (!is.null(cw)) { cat("  pairwise 'worst':\n"); print(round(cw$worst, 3)) }
}

# =========================================================================
rule("6. k.check - ARE THE SMOOTH BASIS DIMENSIONS ADEQUATE?")
# =========================================================================
cat("Rule of thumb: k-index well below 1 with a low p-value means k is too\n")
cat("small and the smooth is oversmoothed. Report this in the Methods.\n")
for (k in all_keys) {
  m <- get_model(k, quiet = TRUE); if (is.null(m)) next
  kc <- tryCatch(k.check(as_gam(m)), error = function(e) NULL)
  if (is.null(kc)) next
  cat("\n---", k, "---\n"); print(round(kc, 4))
}

# =========================================================================
rule("7. SENSITIVITY REFITS: RECORDING EFFORT, AND BUILDING HEIGHT")
# =========================================================================
# These rebuild the base models with extra terms. Check that the formulas
# below match your v5 specification before quoting the output.
RE_NAME <- SC
cyc_knots <- list(c(0.5, 366.5)); names(cyc_knots) <- YC
base_rhs  <- sprintf("%s + %s + s(%s, bs = 'cc', k = 10)", NC, VC, YC)

# --- 7a. recording effort ------------------------------------------------
EFC <- col_of(D, "effort")
if (is.na(EFC) && exists("recording_effort")) {
  eff <- as.data.table(recording_effort)
  ec  <- col_of(eff, "effort"); ew <- col_of(eff, "window"); es <- col_of(eff, "site")
  if (!is.na(ec) && !is.na(ew) && !is.na(es) && !is.na(WC)) {
    D <- merge(D, eff[, c(es, ew, ec), with = FALSE],
               by.x = c(SC, WC), by.y = c(es, ew), all.x = TRUE)
    EFC <- ec
    cat("Merged recording effort from recording_effort into the modelling table.\n")
  }
}
if (!is.na(EFC)) {
  cat("\nEffort column:", EFC, "\n")
  cat("Effort by roof type (is it balanced?):\n")
  print(D[, .(mean_effort = mean(get(EFC), na.rm = TRUE),
              sd_effort   = sd(get(EFC), na.rm = TRUE)), by = roof_label])
  ke <- tryCatch(kruskal.test(D[[EFC]] ~ D$roof_label), error = function(e) NULL)
  if (!is.null(ke)) cat(sprintf("Kruskal-Wallis on effort ~ roof type: chi2 = %.2f, df = %d, P = %s\n",
                                ke$statistic, ke$parameter, fmt_p(ke$p.value)))

  D2 <- D[!is.na(get(EFC)) & get(EFC) > 0]
  # Activity: log-effort as an offset turns a count into a rate
  f_act <- as.formula(sprintf("%s ~ %s + s(%s, bs = 're') + offset(log(%s))",
                              ACT, base_rhs, RE_NAME, EFC))
  m_act_off <- tryCatch(bam(f_act, family = nb(), data = D2, knots = cyc_knots,
                            discrete = TRUE), error = function(e) {cat("  bam failed:", conditionMessage(e), "\n"); NULL})
  if (!is.null(m_act_off)) report_terms(m_act_off, "ACTIVITY with log(effort) offset")

  # Richness: effort as a covariate (an offset is not appropriate for a
  # Gaussian richness model). If the vegetation-density coefficient is
  # unchanged, effort is not driving your result - which is what you want.
  f_rich <- as.formula(sprintf("%s ~ %s + log(%s)", RICH, base_rhs, EFC))
  m_rich_eff <- tryCatch(gamm(f_rich, random = setNames(list(~1), RE_NAME),
                              data = D2, knots = cyc_knots),
                         error = function(e) {cat("  gamm failed:", conditionMessage(e), "\n"); NULL})
  if (!is.null(m_rich_eff)) report_terms(m_rich_eff, "RICHNESS with log(effort) covariate")
} else {
  cat("\nNo recording-effort column found. This is the single most important\n")
  cat("sensitivity check you are missing. Build it from the SongMeter summary\n")
  cat("logs as effort per site per 5-day window, add it to model_data, then\n")
  cat("re-run this script.\n")
}

# --- 7b. building height ------------------------------------------------
HC <- col_of(D, "height")
if (!is.na(HC)) {
  cat("\nHeight column:", HC, "\n")
  cat("Height by roof type (this is the confound):\n")
  print(D[, .(mean_height = mean(get(HC), na.rm = TRUE),
              min = min(get(HC), na.rm = TRUE), max = max(get(HC), na.rm = TRUE)),
           by = roof_label])
  cat(sprintf("Correlation, height vs vegetation density: r = %.3f\n",
              cor(D[[HC]], D[[VC]], use = "complete.obs")))

  # Base model (no height) - the comparison partner for both AIC and LRT
  m_base <- get_model("richness_main", quiet = TRUE)

  f_rh <- as.formula(sprintf("%s ~ %s + %s", RICH, base_rhs, HC))
  m_rh <- tryCatch(gamm(f_rh, random = setNames(list(~1), RE_NAME),
                        data = D, knots = cyc_knots), error = function(e) NULL)
  if (!is.null(m_rh)) report_terms(m_rh, "RICHNESS + height")

  # --- AIC and LRT: does height improve the model at all? -----------------
  # For gamm() objects, AIC/anova must compare the $lme component, not the
  # object itself - comparing $gam components or the whole gamm list gives
  # meaningless or silently wrong numbers. Both models must be fitted by ML
  # (not REML) for the comparison to be valid, since REML likelihoods are
  # not comparable across models with different fixed-effect structures.
  # gamm() defaults to REML, so refit both with method = "ML" for this
  # specific comparison only - do not use these ML refits for reporting
  # coefficients elsewhere, use the REML fits for that.
  if (!is.null(m_base) && !is.null(m_rh)) {
    cat("\n--- Does building height improve the richness model? ---\n")
    m_base_ml <- tryCatch(gamm(formula(as_gam(m_base)),
                               random = setNames(list(~1), RE_NAME),
                               data = D, knots = cyc_knots, method = "ML"),
                          error = function(e) {cat("  base ML refit failed:", conditionMessage(e), "\n"); NULL})
    m_rh_ml   <- tryCatch(gamm(f_rh, random = setNames(list(~1), RE_NAME),
                               data = D, knots = cyc_knots, method = "ML"),
                          error = function(e) {cat("  +height ML refit failed:", conditionMessage(e), "\n"); NULL})
    if (!is.null(m_base_ml) && !is.null(m_rh_ml)) {
      a0 <- AIC(as_lme(m_base_ml)); a1 <- AIC(as_lme(m_rh_ml))
      lrt <- tryCatch(anova(as_lme(m_base_ml), as_lme(m_rh_ml)), error = function(e) NULL)
      cat(sprintf("  AIC without height: %.2f\n", a0))
      cat(sprintf("  AIC with height:    %.2f  (delta AIC = %+.2f)\n", a1, a1 - a0))
      cat("  Rule of thumb: delta AIC < 2 means the extra term is not earning\n")
      cat("  its complexity; delta AIC > 2 favours including it; > 10 is strong.\n")
      if (!is.null(lrt) && nrow(lrt) == 2) {
        cat(sprintf("  Likelihood ratio test: L.Ratio = %.3f, df = %d, P = %s\n",
                    lrt$L.Ratio[2], lrt$df[2] - lrt$df[1], fmt_p(lrt$`p-value`[2])))
      }
      cat("\n  Report this as: 'Adding building height did not improve model fit\n")
      cat("  (delta AIC = ", sprintf("%+.2f", a1 - a0), ", LRT P = ",
          if (!is.null(lrt)) fmt_p(lrt$`p-value`[2]) else "NA",
          "), and its coefficient was non-significant (see above) - so the\n", sep = "")
      cat("  roof-type effect does not appear to be a height artefact.'\n")
      cat("  Adjust the wording if delta AIC or the LRT instead favour height -\n")
      cat("  in that case say so plainly rather than keeping this sentence.\n")
    }
  } else {
    cat("\n  Cannot compute delta AIC - richness_main or the +height model\n")
    cat("  did not fit. Check MODELS$richness_main in figure_setup.R.\n")
  }

  if (!is.na(RC)) {
    f_rhr <- as.formula(sprintf("%s ~ %s + %s + %s", RICH, base_rhs, RC, HC))
    m_rhr <- tryCatch(gamm(f_rhr, random = setNames(list(~1), RE_NAME),
                           data = D, knots = cyc_knots), error = function(e) NULL)
    if (!is.null(m_rhr)) {
      report_terms(m_rhr, "RICHNESS + roof type + height")
      cat("\n>> KEY QUESTION: do the roof-type contrasts survive with height in\n")
      cat("   the model? If they do, your claim is much stronger than it is now.\n")
      cat("   If they shrink, say so - that is a genuinely interesting result.\n")
    }

    # Same AIC/LRT logic for roof type + height vs roof type alone
    m_roof <- get_model("richness_roof", quiet = TRUE)
    if (!is.null(m_roof) && !is.null(m_rhr)) {
      cat("\n--- Does height improve the model further, beyond roof type? ---\n")
      m_roof_ml <- tryCatch(gamm(formula(as_gam(m_roof)),
                                 random = setNames(list(~1), RE_NAME),
                                 data = D, knots = cyc_knots, method = "ML"),
                            error = function(e) NULL)
      m_rhr_ml  <- tryCatch(gamm(f_rhr, random = setNames(list(~1), RE_NAME),
                                 data = D, knots = cyc_knots, method = "ML"),
                            error = function(e) NULL)
      if (!is.null(m_roof_ml) && !is.null(m_rhr_ml)) {
        a0 <- AIC(as_lme(m_roof_ml)); a1 <- AIC(as_lme(m_rhr_ml))
        lrt2 <- tryCatch(anova(as_lme(m_roof_ml), as_lme(m_rhr_ml)), error = function(e) NULL)
        cat(sprintf("  AIC, roof type alone:        %.2f\n", a0))
        cat(sprintf("  AIC, roof type + height:      %.2f  (delta AIC = %+.2f)\n", a1, a1 - a0))
        if (!is.null(lrt2) && nrow(lrt2) == 2)
          cat(sprintf("  Likelihood ratio test: L.Ratio = %.3f, df = %d, P = %s\n",
                      lrt2$L.Ratio[2], lrt2$df[2] - lrt2$df[1], fmt_p(lrt2$`p-value`[2])))
      }
    }
  }
} else {
  cat("\nNo height column in model_data. You have this from the GIS work -\n")
  cat("merge it in, then re-run. It is the largest scientific gap in the draft.\n")
}

# =========================================================================
rule("8. SPATIAL AUTOCORRELATION ACROSS SITES (Moran's I)")
# =========================================================================
morans_I <- function(v, lat, lon, nperm = 4999, seed = 1) {
  ok <- is.finite(v) & is.finite(lat) & is.finite(lon)
  v <- v[ok]; lat <- lat[ok]; lon <- lon[ok]; n <- length(v)
  if (n < 5) { cat("  too few sites for a meaningful test.\n"); return(invisible(NULL)) }
  dm <- as.matrix(dist(cbind(lon, lat)))
  W <- 1 / dm; diag(W) <- 0; W[!is.finite(W)] <- 0; W <- W / sum(W)
  I_of <- function(x) { z <- x - mean(x); n * sum(W * outer(z, z)) / sum(z^2) }
  obs <- I_of(v)
  set.seed(seed)
  perm <- replicate(nperm, I_of(sample(v)))
  p <- (1 + sum(abs(perm) >= abs(obs))) / (nperm + 1)
  cat(sprintf("  Moran's I = %+.3f (expected %.3f under no structure), P = %.4f (%d permutations)\n",
              obs, -1 / (n - 1), p, nperm))
  invisible(list(I = obs, p = p))
}
lat_c <- col_of(D, "lat"); lon_c <- col_of(D, "lon")
site_xy <- NULL
if (!is.na(lat_c) && !is.na(lon_c)) {
  site_xy <- D[, .(lat = first(get(lat_c)), lon = first(get(lon_c))), by = c(SC)]
} else {
  si <- find_obj(c("site_info", "sites", "site_table", "site_coords"), quiet = TRUE)
  if (!is.null(si)) {
    si <- as.data.table(si)
    la <- col_of(si, "lat"); lo <- col_of(si, "lon"); ss <- col_of(si, "site")
    if (!is.na(la) && !is.na(lo) && !is.na(ss))
      site_xy <- si[, .(lat = get(la), lon = get(lo)), by = c(ss)]
  }
}
if (!is.null(site_xy)) {
  m <- get_model("richness_main", quiet = TRUE)
  if (!is.null(m)) {
    g <- as_gam(m)
    r <- data.table(site = as.character(g$model[[RE_NAME]]), res = residuals(g, type = "response"))
    if (all(is.na(r$site))) r <- data.table(site = as.character(D[[SC]]), res = residuals(g, type = "response"))
    rs <- r[, .(res = mean(res, na.rm = TRUE)), by = site]
    setnames(site_xy, 1, "site"); site_xy[, site := as.character(site)]
    mm <- merge(rs, site_xy, by = "site")
    cat("Site-level mean residuals from the richness model, n =", nrow(mm), "sites:\n")
    morans_I(mm$res, mm$lat, mm$lon)
    cat("\nAll 13 sites lie within roughly 3 km2, so if this is significant you\n")
    cat("must say so; if it is not, that is a sentence worth having in the\n")
    cat("Discussion because a reader will wonder.\n")
  }
} else {
  cat("No site coordinates found. Add Latitude/Longitude columns (you have them\n")
  cat("in the QGIS site layer) and re-run - this is a cheap, high-value check.\n")
}

sink(); close(con)
cat("\nWritten to:", normalizePath(digest_file), "\n")
