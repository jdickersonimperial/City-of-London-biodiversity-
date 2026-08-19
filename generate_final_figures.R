# =========================================================================
# generate_final_figures.R   (v3 - Fig 7 now two-panel, richness + activity)
# Run AFTER your v5 analysis script, in the same session:
#
#   source("figure_setup.R")
#   source("generate_final_figures.R")
#
# CHANGES FROM v2 (this version):
#   - Fig 7 is now two panels: (a) species richness, (b) total activity,
#     matching Fig 8's layout and encoding. Panel (a) uses the same model
#     and logic as before (m_richness_roof); panel (b) is new and uses
#     a_roof_full (registered below as a candidate for MODELS$activity_roof
#     if not already present in figure_setup.R).
#   - REMINDER, not fixed here: lab("richness_pred") in figure_setup.R still
#     reads "Predicted species richness". Change it to "Species richness"
#     in figure_setup.R for consistency with Fig 8 - a one-line edit there,
#     not duplicated here, so the label stays centralised as intended.
#   - Both Fig 7 panels evaluate NDVI and vegetation density at the DATASET
#     mean (via newdata_from's default behaviour), not at each roof type's
#     own mean. Since ground garden sites sit above the dataset mean on
#     both covariates, their model-based point can therefore lie above the
#     raw mean of the plotted site points - this is an "adjusted for other
#     covariates" estimate, not a bug, but say so in the caption. An
#     alternative (each type evaluated at its own mean NDVI/vegetation) is
#     included below as a commented-out block if you and your supervisor
#     prefer that reading instead - it is a modelling choice, not a display
#     one, so it's left off by default rather than switched silently.
#
# FIGURE NUMBERING - differs from your current draft, see the renumbering
# map in figure_legends_and_placement.md.
#   Fig 1  Map of study sites (QGIS - not produced here)
#   Fig 2  Sampling effort
#   Fig 3  Seasonal NDVI trajectories by site
#   Fig 4  Richness: main effects (a-c)
#   Fig 5  Activity: main effects (a-c)
#   Fig 6  NDVI x vegetation density interaction (a tensor, b activity)
#   Fig 7  Predicted richness and activity by roof type (a-b)
#   Fig 8  Seasonal curves by roof type (a richness, b activity)
#   Fig 9  Peak timing posteriors (a richness, b activity)
#   Fig 10 Within-site NDVI x seasonal amplitude
#   Fig 11 PCoA partial effects (a-f)
#   Fig 12 Focal species
#   Fig S1 Predictor collinearity
#   Fig S2 NDVI trajectories, faceted
#   Fig S3+ gam.check diagnostics
# =========================================================================

if (!exists("lab")) stop("Run source('figure_setup.R') first.")

# Guard against a common mistake: sourcing this script after theme_fig was
# set in an earlier session, or after some other code called theme_set()
# in between. Re-assert it here so every figure in this run uses the
# intended theme regardless of what happened earlier in the session.
theme_set(theme_fig)

# Register the real activity+roof-type model object as a candidate, in case
# figure_setup.R's MODELS$activity_roof list doesn't already include it.
# Harmless to prepend even if it's already there or doesn't exist yet -
# get_model() just checks each candidate in order and uses the first that
# resolves to a real object.
MODELS$activity_roof <- unique(c("a_roof_full", MODELS$activity_roof))

for (nm in c("model_data", "aligned_data", "site_means", "recording_effort")) {
  if (exists(nm)) assign(nm, add_roof_label(get(nm)))
}
stopifnot(exists("model_data"))
D  <- as.data.table(model_data)
SC <- col_of(D, "site"); RC <- col_of(D, "roof")
NC <- col_of(D, "ndvi"); VC <- col_of(D, "veg"); YC <- col_of(D, "doy")
RICH <- col_of(D, "richness"); ACT <- col_of(D, "activity")
WC <- col_of(D, "window")

cat("\nColumns detected: site =", SC, "| roof =", RC, "| NDVI =", NC,
    "| density =", VC, "| DOY =", YC, "\n")

# =========================================================================
# FIG 2 - SAMPLING EFFORT
# =========================================================================
cat("\nFig 2: sampling effort\n")
if (exists("recording_effort")) {
  eff <- as.data.table(recording_effort)
  EC  <- col_of(eff, "effort"); EW <- col_of(eff, "window")
  eff <- merge(eff, unique(D[, c(SC, "roof_label"), with = FALSE]), by = SC, all.x = TRUE)
  eff <- eff[!is.na(roof_label)]
  p <- ggplot(eff, aes(.data[[EW]], reorder(.data[[SC]], as.numeric(roof_label)),
                       fill = .data[[EC]])) +
    geom_tile() +
    scale_fill_viridis_c(name = "Recording\neffort", direction = -1) +
    labs(x = lab("date"), y = NULL)
} else {
  cat("  recording_effort not found - showing windows with data instead.\n")
  cat("  The caption MUST say this is coverage, not measured effort.\n")
  p <- ggplot(D, aes(.data[[WC]], reorder(.data[[SC]], as.numeric(roof_label)),
                     fill = roof_label)) +
    geom_tile() + scale_roof(line = FALSE) + labs(x = lab("date"), y = NULL)
}
save_fig(p, "Fig02_Sampling_effort", 10, 6)

# =========================================================================
# FIG 3 - SEASONAL NDVI TRAJECTORIES
# =========================================================================
cat("\nFig 3: NDVI trajectories\n")
p3 <- ggplot(D, aes(.data[[YC]], .data[[NC]], group = .data[[SC]],
                    colour = roof_label, linetype = roof_label)) +
  geom_line_w(alpha = 0.9, lw = 0.7) +
  geom_point(aes(shape = roof_label), size = 1.2, alpha = 0.6) +
  scale_roof() +
  scale_shape_manual(values = roof_shapes, name = lab("roof"), drop = FALSE) +
  labs(x = lab("DOY"), y = lab("NDVI"))
save_fig(p3, "Fig03_NDVI_trajectories", 9, 6)
save_fig(p3 + facet_wrap(as.formula(paste("~", SC)), ncol = 4) +
           theme(legend.position = "bottom"),
         "FigS02_NDVI_trajectories_faceted", 12, 8)

# =========================================================================
# FIGS 4, 5 - MAIN EFFECT PANELS
# =========================================================================
plot_effect <- function(m, data, xvar, yvar, log_y = FALSE, ylab_key = yvar) {
  xs <- if (xvar == YC) seq(1, 366, length.out = 300) else
    seq(min(data[[xvar]], na.rm = TRUE), max(data[[xvar]], na.rm = TRUE), length.out = 300)
  nd <- newdata_from(data, vary = setNames(list(xs), xvar))
  pr <- predict_ci(m, nd)
  p <- ggplot() +
    geom_point(data = data, aes(.data[[xvar]], .data[[yvar]]),
               alpha = 0.15, colour = "black", size = 1) +
    geom_ribbon(data = pr, aes(.data[[xvar]], ymin = lo, ymax = hi),
                fill = RIBBON_COL, alpha = 0.35) +
    geom_line_w(data = pr, mapping = aes(.data[[xvar]], fit),
                colour = FIT_COL, lw = 1.1) +
    labs(x = lab(xvar), y = lab(ylab_key))
  if (log_y) p <- p + scale_y_log10()
  p
}

m_rich <- get_model("richness_main")
if (!is.null(m_rich)) {
  cat("\nFig 4: richness main effects\n")
  save_fig(tag_panels(
    plot_effect(m_rich, D, NC, RICH) | plot_effect(m_rich, D, VC, RICH) |
      plot_effect(m_rich, D, YC, RICH)),
    "Fig04_Richness_effects", 13, 4.5)
}

m_act <- get_model("activity_main")
if (!is.null(m_act)) {
  cat("\nFig 5: activity main effects\n")
  save_fig(tag_panels(
    plot_effect(m_act, D, NC, ACT, ylab_key = "total_activity_short") |
      plot_effect(m_act, D, VC, ACT, ylab_key = "total_activity_short") |
      plot_effect(m_act, D, YC, ACT, ylab_key = "total_activity_short")),
    "Fig05_Activity_effects", 13, 4.5)
}

# =========================================================================
# FIG 6 - NDVI x VEGETATION DENSITY INTERACTION
# (a) tensor smooth on richness, extrapolation masked
# (b) NDVI effect on activity at low / median / high vegetation density
# =========================================================================
cat("\nFig 6: NDVI x vegetation density interaction\n")
p6a <- NULL; p6b <- NULL

m_te <- get_model("richness_tensor")
if (!is.null(m_te)) {
  gr <- expand.grid(
    a = seq(min(D[[NC]], na.rm = TRUE), max(D[[NC]], na.rm = TRUE), length.out = 120),
    b = seq(min(D[[VC]], na.rm = TRUE), max(D[[VC]], na.rm = TRUE), length.out = 120))
  nd <- newdata_from(D, vary = setNames(list(gr$a, gr$b), c(NC, VC)))
  pr <- predict_ci(m_te, nd)
  tf <- mask_too_far(pr, D, NC, VC, dist = 0.08)
  cat("  masked", round(100 * mean(tf)), "% of the grid as extrapolation",
      "-> quote this in the caption\n")
  cat("  vegetation density has", uniqueN(round(D[[VC]], 4)),
      "unique values (one per site) -> explains the banding\n")
  pr <- pr[!tf, ]
  p6a <- ggplot(pr, aes(.data[[NC]], .data[[VC]], fill = fit)) +
    geom_raster(interpolate = TRUE) +
    geom_contour_w(mapping = aes(z = fit), colour = "white", lw = 0.3, alpha = 0.7) +
    geom_point(data = D, mapping = aes(.data[[NC]], .data[[VC]]), inherit.aes = FALSE,
               colour = "black", size = 0.5, alpha = 0.5) +
    scale_fill_viridis_c(name = "Fitted\nspecies\nrichness") +
    labs(x = lab("NDVI"), y = lab("Veg_Density")) +
    coord_cartesian(expand = FALSE)
}

m_ai <- get_model("activity_interact")
if (!is.null(m_ai)) {
  qs   <- quantile(D[[VC]], c(0.10, 0.50, 0.90), na.rm = TRUE)
  qlab <- sprintf(c("Low (10th pct, %.2f)", "Median (%.2f)", "High (90th pct, %.2f)"), qs)
  xs   <- seq(min(D[[NC]], na.rm = TRUE), max(D[[NC]], na.rm = TRUE), length.out = 300)
  pr <- rbindlist(lapply(seq_along(qs), function(i) {
    nd <- newdata_from(D, vary = setNames(list(xs, rep(qs[i], length(xs))), c(NC, VC)))
    out <- as.data.table(predict_ci(m_ai, nd)); out[, veg_level := qlab[i]]; out
  }))
  pr[, veg_level := factor(veg_level, levels = qlab)]
  p6b <- ggplot(pr, aes(.data[[NC]], fit, colour = veg_level,
                        fill = veg_level, linetype = veg_level)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line_w(lw = 1.1) +
    scale_colour_manual(values = c("#D55E00", "#009E73", "#0072B2"),
                        name = "Vegetation density") +
    scale_fill_manual(values = c("#D55E00", "#009E73", "#0072B2"),
                      name = "Vegetation density") +
    scale_linetype_manual(values = c("dashed", "solid", "longdash"),
                          name = "Vegetation density") +
    labs(x = lab("NDVI"), y = lab("total_activity_short"))
}

if (!is.null(p6a) && !is.null(p6b)) {
  save_fig(tag_panels(p6a | p6b), "Fig06_NDVI_VegDensity_interaction", 13, 5.5)
} else if (!is.null(p6a)) {
  save_fig(p6a, "Fig06a_Tensor_richness", 7, 5.5)
} else if (!is.null(p6b)) {
  save_fig(p6b, "Fig06b_Activity_interaction", 7.5, 5.5)
} else {
  cat("  Neither interaction model found. Run ls(pattern = 'tensor|interact')\n")
  cat("  and add the real object names to MODELS in figure_setup.R.\n")
}

# =========================================================================
# FIG 7 - PREDICTED RICHNESS AND ACTIVITY BY ROOF TYPE (two panels)
# Site means behind the predictions so n = 3-5 per group is visible.
#
# =========================================================================
build_roof_panel <- function(model_key, response_col, ylab_key, logy = FALSE) {
  m <- get_model(model_key)
  if (is.null(m)) return(NULL)
  lv <- model_levels(m, "roof")
  if (is.null(lv)) lv <- unique(as.character(D[[RC]]))

  # newdata_from(..., roof = lv[1]) sets roof_type (and roof_type_ordered,
  # if the model has one) for a single row; replicate to one row per level,
  # then overwrite both columns so they actually vary across rows rather
  # than all repeating lv[1]. NDVI, vegetation density and DOY are left at
  # newdata_from's default - the dataset mean - for every row; see the
  # note at the top of this file on what that does to the ground-garden
  # estimate specifically.
  nd <- newdata_from(D, vary = list(), m = m, roof = lv[1])
  nd <- nd[rep(1, length(lv)), , drop = FALSE]
  nd[[RC]] <- factor(lv, levels = lv)
  ro <- col_of(D, "roof_ord")
  if (!is.na(ro)) {
    lvo <- model_levels(m, "roof_ord"); if (is.null(lvo)) lvo <- lv
    nd[[ro]] <- factor(lv, levels = lvo, ordered = TRUE)
  }

  pr <- as.data.table(predict_ci(m, nd))
  pr[, roof_label := factor(roof_labels[as.character(get(RC))], levels = roof_display_order)]

  site_pts <- D[, .(mean_val = mean(get(response_col), na.rm = TRUE)),
               by = c(SC, "roof_label")]

  p <- ggplot() +
    geom_point(data = site_pts, aes(roof_label, mean_val, shape = roof_label),
               colour = "black", size = 2.2, alpha = 0.7,
               position = position_jitter(width = 0.12, height = 0)) +
    geom_pointrange_w(data = pr,
                      mapping = aes(roof_label, fit, ymin = lo, ymax = hi,
                                    colour = roof_label),
                      size = 0.9, lw = 1.1) +
    scale_roof(fill = FALSE, line = FALSE, legend = FALSE) +
    scale_shape_manual(values = roof_shapes, guide = "none") +
    labs(x = lab("roof"), y = lab(ylab_key))
  if (logy) p <- p + scale_y_log10()
  p
}

cat("\nFig 7: roof type (richness + activity)\n")
p7a <- build_roof_panel("richness_roof", RICH, "richness_pred", logy = FALSE)
p7b <- build_roof_panel("activity_roof", ACT,  "total_activity_short", logy = TRUE)

if (!is.null(p7a) && !is.null(p7b)) {
  save_fig(tag_panels(p7a | p7b), "Fig07_RoofType_levels", 13, 5.5)
} else if (!is.null(p7a)) {
  cat("  activity_roof model not found - saving richness panel only.\n")
  cat("  Fit a_roof_full (or your equivalent activity+roof model) and re-run.\n")
  save_fig(p7a, "Fig07a_RoofType_richness", 7.5, 5.5)
} else {
  cat("  richness_roof model not found - skipping Fig 7 entirely.\n")
}

# --- Alternative reading (commented out): evaluate each roof type at its
# OWN mean NDVI and vegetation density, rather than the dataset mean. This
# changes what the model-based point represents - "a typical site of this
# type" rather than "this type, if it had average-for-the-whole-dataset
# vegetation" - and is a modelling choice to agree with your supervisor
# before switching to it, not a formatting change. Left inactive by
# default; the two build_roof_panel() calls above are what Fig07 actually
# uses unless you replace them with this version.
#
# build_roof_panel_owncov <- function(model_key, response_col, ylab_key, logy = FALSE) {
#   m <- get_model(model_key); if (is.null(m)) return(NULL)
#   lv <- model_levels(m, "roof"); if (is.null(lv)) lv <- unique(as.character(D[[RC]]))
#   means_by_roof <- D[, .(ndvi_m = mean(get(NC), na.rm = TRUE),
#                          veg_m  = mean(get(VC), na.rm = TRUE)), by = c(RC)]
#   nd <- newdata_from(D, vary = list(), m = m, roof = lv[1])
#   nd <- nd[rep(1, length(lv)), , drop = FALSE]
#   nd[[RC]] <- factor(lv, levels = lv)
#   nd[[NC]] <- means_by_roof[[which(names(means_by_roof) == "ndvi_m")]][match(lv, means_by_roof[[RC]])]
#   nd[[VC]] <- means_by_roof[[which(names(means_by_roof) == "veg_m")]][match(lv, means_by_roof[[RC]])]
#   ro <- col_of(D, "roof_ord")
#   if (!is.na(ro)) {
#     lvo <- model_levels(m, "roof_ord"); if (is.null(lvo)) lvo <- lv
#     nd[[ro]] <- factor(lv, levels = lvo, ordered = TRUE)
#   }
#   pr <- as.data.table(predict_ci(m, nd))
#   pr[, roof_label := factor(roof_labels[as.character(get(RC))], levels = roof_display_order)]
#   site_pts <- D[, .(mean_val = mean(get(response_col), na.rm = TRUE)), by = c(SC, "roof_label")]
#   p <- ggplot() +
#     geom_point(data = site_pts, aes(roof_label, mean_val, shape = roof_label),
#                colour = "black", size = 2.2, alpha = 0.7,
#                position = position_jitter(width = 0.12, height = 0)) +
#     geom_pointrange_w(data = pr, aes(roof_label, fit, ymin = lo, ymax = hi,
#                                      colour = roof_label), size = 0.9, lw = 1.1) +
#     scale_roof(fill = FALSE, line = FALSE, legend = FALSE) +
#     scale_shape_manual(values = roof_shapes, guide = "none") +
#     labs(x = lab("roof"), y = lab(ylab_key))
#   if (logy) p <- p + scale_y_log10()
#   p
# }

# =========================================================================
# FIG 8 - SEASONAL CURVES BY ROOF TYPE
# =========================================================================
plot_seasonal <- function(m, data, ylab_key, log_y = FALSE) {
  lv <- model_levels(m, "roof"); if (is.null(lv)) lv <- unique(as.character(data[[RC]]))
  doy <- 1:366
  pr <- rbindlist(lapply(lv, function(rl) {
    nd <- newdata_from(data, vary = setNames(list(doy), YC), m = m, roof = rl)
    as.data.table(predict_ci(m, nd))
  }))
  pr[, roof_label := factor(roof_labels[as.character(get(RC))], levels = roof_display_order)]
  p <- ggplot(pr, aes(.data[[YC]], fit, colour = roof_label,
                      fill = roof_label, linetype = roof_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line_w(lw = 1.1) + scale_roof() +
    labs(x = lab("DOY"), y = lab(ylab_key))
  if (log_y) p <- p + scale_y_log10(breaks = c(10, 100, 1000, 10000),
                                    labels = c("10", "100", "1,000", "10,000"))
  p
}
m_rs <- get_model("richness_shape"); m_as <- get_model("activity_shape")
if (!is.null(m_rs) || !is.null(m_as)) {
  cat("\nFig 8: seasonal curves by roof type\n")
  ps <- list()
  if (!is.null(m_rs)) ps$a <- plot_seasonal(m_rs, D, "richness") +
    theme(legend.position = "none")
  if (!is.null(m_as)) ps$b <- plot_seasonal(m_as, D, "total_activity_log", log_y = TRUE) +
    theme(legend.position = "none")
  if (HAVE_COWPLOT) {
    # Build the legend from a throwaway plot carrying the full roof_label
    # factor levels, so it always shows all four categories even if one
    # panel happens not to display every level in its plotted range.
    leg_src <- ggplot(data.frame(x = 1, roof_label = factor(roof_display_order,
                                                             levels = roof_display_order)),
                      aes(x, x, colour = roof_label, linetype = roof_label)) +
      geom_line_w(lw = 1.1) + scale_roof() +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal")
    shared_legend <- cowplot::get_legend(leg_src)
    fig8 <- tag_panels(wrap_plots(ps, ncol = 2))
    fig8 <- cowplot::plot_grid(fig8, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
    save_fig(fig8, "Fig08_Seasonal_byRoof", 13, 6)
  } else {
    # Fallback: single collected legend on the right via patchwork, at least
    # functional even without cowplot
    ps$a <- ps$a + theme(legend.position = "right")
    save_fig(tag_panels(wrap_plots(ps, ncol = 2, guides = "collect")) &
               theme(legend.position = "right"),
             "Fig08_Seasonal_byRoof", 14, 5.5)
  }
}


qbin <- function(x) cut(x, quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                        include.lowest = TRUE, labels = c("low", "mid", "high"))
rng_lab <- function(x, b) {
  r <- tapply(x, b, range, na.rm = TRUE)
  sprintf("%.2f-%.2f", vapply(r, `[`, numeric(1), 1), vapply(r, `[`, numeric(1), 2))
}
if (!"veg_bin"  %in% names(D)) D[, veg_bin  := qbin(get(VC))]
if (!"ndvi_bin" %in% names(D)) D[, ndvi_bin := qbin(get(NC))]
if (!"veg_ord"  %in% names(D)) D[, veg_ord  := as.ordered(veg_bin)]
if (!"ndvi_ord" %in% names(D)) D[, ndvi_ord := as.ordered(ndvi_bin)]
VLAB <- rng_lab(D[[VC]], D$veg_bin)
NLAB <- rng_lab(D[[NC]], D$ndvi_bin)
cat("Vegetation density tertiles (m3/m2):", paste(VLAB, collapse = " | "), "\n")
cat("NDVI tertiles:                      ", paste(NLAB, collapse = " | "), "\n")

# --- reuse existing model objects; fit only if genuinely missing ---------
need_fit <- function(nm) !exists(nm, envir = .GlobalEnv)
if (need_fit("m_veg") || need_fit("m_ndvi") || need_fit("a_veg") || need_fit("a_ndvi")) {
  cat("One or more of m_veg / m_ndvi / a_veg / a_ndvi not found - fitting\n")
  cat("them now with the same specification used previously (gam() with\n")
  cat("s(site, bs = 're'), not gamm()/bam()) so results stay comparable.\n")
  D_df <- as.data.frame(D)
  fitv <- function(resp, byvar, fam) gam(as.formula(sprintf(
    "%s ~ %s + s(%s, bs='cc', k=10) + s(%s, by=%s, bs='cc', k=10) + s(%s, bs='re')",
    resp, byvar, YC, YC, byvar, SC)), data = D_df, knots = cyc_knots,
    method = "REML", family = fam)
  if (need_fit("m_veg"))  m_veg  <- fitv(RICH, "veg_ord",  gaussian())
  if (need_fit("m_ndvi")) m_ndvi <- fitv(RICH, "ndvi_ord", gaussian())
  if (need_fit("a_veg"))  a_veg  <- fitv(ACT,  "veg_ord",  nb())
  if (need_fit("a_ndvi")) a_ndvi <- fitv(ACT,  "ndvi_ord", nb())
}

# --- panel builder, reusing predict_ci()/lab()/theme_fig from figure_setup.R
BIN_COLS <- if (identical(PALETTE, "safe")) c("#E69F00", "#56B4E9", "#009E73") else
  c("#D95F02", "#7570B3", "#1B9E77")
BIN_LTYS <- c("solid", "dashed", "longdash")

bin_panel <- function(m, byvar, level_labels, ylab_key, logy, title) {
  fac_levels <- levels(D[[byvar]])
  grid <- expand.grid(doy_seq = 1:366, lvl = fac_levels, stringsAsFactors = FALSE)
  names(grid)[1] <- YC
  grid[[byvar]] <- factor(grid$lvl, levels = fac_levels, ordered = TRUE)
  grid[[SC]] <- D[[SC]][1]
  pr <- as.data.table(predict_ci(m, grid))
  pr[, grp := factor(level_labels[match(get(byvar), fac_levels)], levels = level_labels)]
  p <- ggplot(pr, aes(.data[[YC]], fit, colour = grp, fill = grp, linetype = grp)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line_w(lw = 1.1) +
    scale_colour_manual(values = BIN_COLS, name = NULL) +
    scale_fill_manual(values = BIN_COLS, name = NULL) +
    scale_linetype_manual(values = BIN_LTYS, name = NULL) +
    labs(x = lab("DOY"), y = lab(ylab_key), title = title) +
    theme(legend.position = "bottom")
  if (logy) p <- p + scale_y_log10()
  p
}

cat("\nFig 8B: seasonal effects of vegetation density and NDVI (tertile bins)\n")
p1 <- bin_panel(m_veg,  "veg_ord",  VLAB, "richness",             FALSE, "Vegetation density")
p2 <- bin_panel(m_ndvi, "ndvi_ord", NLAB, "richness",             FALSE, "NDVI")
p3 <- bin_panel(a_veg,  "veg_ord",  VLAB, "total_activity_short", TRUE,  "Vegetation density")
p4 <- bin_panel(a_ndvi, "ndvi_ord", NLAB, "total_activity_short", TRUE,  "NDVI")

fig8b <- tag_panels(wrap_plots(list(p1, p2, p3, p4), ncol = 2))
save_fig(fig8b, "Fig08B_Seasonal_VegNDVI_bins", 13, 10)

# =========================================================================
# FIG 9 - PEAK TIMING POSTERIORS (keeps all 1000 draws)
# =========================================================================
simulate_peak_draws <- function(m, data, n_sim = 1000, seed = 42) {
  g <- as_gam(m); set.seed(seed); doy <- 1:366
  bs <- tryCatch(rmvn(n_sim, coef(g), vcov(g, unconditional = TRUE)),
                 error = function(e) tryCatch(rmvn(n_sim, coef(g), vcov(g)),
                                              error = function(e2) NULL))
  if (is.null(bs)) { cat("  rmvn failed - skipping\n"); return(NULL) }
  lv <- model_levels(m, "roof"); if (is.null(lv)) lv <- unique(as.character(data[[RC]]))
  ex <- site_smooths(g)
  rbindlist(lapply(lv, function(rl) {
    nd <- newdata_from(data, vary = setNames(list(doy), YC), m = m, roof = rl)
    Xp <- if (length(ex)) predict(g, newdata = nd, type = "lpmatrix", exclude = ex)
          else            predict(g, newdata = nd, type = "lpmatrix")
    data.table(roof = rl, peak_day = doy[apply(Xp %*% t(bs), 2, which.max)])
  }))
}
peak_plot <- function(draws, ttl) {
  draws[, roof_label := factor(roof_labels[roof], levels = roof_display_order)]
  s <- draws[, .(med = median(peak_day),
                 lo = quantile(peak_day, 0.025),
                 hi = quantile(peak_day, 0.975)), by = roof_label]
  cat("  ", ttl, "peak days (median [95% CI]):\n")
  for (i in seq_len(nrow(s)))
    cat(sprintf("     %-14s %3.0f [%3.0f, %3.0f]\n", as.character(s$roof_label[i]),
                s$med[i], s$lo[i], s$hi[i]))
  ggplot(draws, aes(peak_day, roof_label, fill = roof_label)) +
    geom_violin(scale = "width", alpha = 0.45, colour = NA) +
    geom_linerange_w(data = s, mapping = aes(xmin = lo, xmax = hi, y = roof_label),
                     inherit.aes = FALSE, lw = 0.7) +
    geom_point(data = s, mapping = aes(med, roof_label, shape = roof_label),
               inherit.aes = FALSE, size = 2.6) +
    scale_fill_manual(values = roof_cols, guide = "none") +
    scale_shape_manual(values = roof_shapes, guide = "none") +
    scale_x_continuous(limits = c(1, 366), breaks = seq(0, 360, 60)) +
    labs(x = lab("peak_day"), y = NULL, title = ttl)
}
pk <- list()
if (!is.null(m_rs)) {
  d <- simulate_peak_draws(m_rs, D)
  if (!is.null(d)) { pk$a <- peak_plot(d, "Species richness")
                     fwrite(d, file.path(tab_dir, "peak_draws_richness.csv")) }
}
if (!is.null(m_as)) {
  d <- simulate_peak_draws(m_as, D)
  if (!is.null(d)) { pk$b <- peak_plot(d, "Total activity")
                     fwrite(d, file.path(tab_dir, "peak_draws_activity.csv")) }
}
if (length(pk)) {
  cat("\nFig 9: peak timing posteriors\n")
  save_fig(tag_panels(wrap_plots(pk, ncol = 1)), "Fig09_Peak_timing", 8, 7)
}


# =========================================================================
# FIG 11 - PCoA PARTIAL EFFECTS
# =========================================================================
m_p1 <- get_model("pcoa1"); m_p2 <- get_model("pcoa2")
AD <- if (exists("aligned_data")) as.data.table(aligned_data) else D
if (!is.null(m_p1)) {
  cat("\nFig 11: PCoA partial effects\n")
  panels <- list()
  for (xv in c(NC, VC, YC)) panels[[length(panels) + 1]] <- plot_effect(m_p1, AD, xv, "PCoA1")
  if (!is.null(m_p2))
    for (xv in c(NC, VC, YC)) panels[[length(panels) + 1]] <- plot_effect(m_p2, AD, xv, "PCoA2")
  save_fig(tag_panels(wrap_plots(panels, ncol = 3)),
           "Fig11_PCoA_effects", 13, if (!is.null(m_p2)) 8.5 else 4.5)
}

# =========================================================================

# =========================================================================
# FIG S1 - PREDICTOR COLLINEARITY
# =========================================================================
cat("\nFig S1: predictor collinearity\n")
pv <- c(NC, VC, YC)
for (k in c("height", "area")) { cc <- col_of(D, k); if (!is.na(cc)) pv <- c(pv, cc) }
pv <- unique(pv[!is.na(pv)])
cm  <- cor(as.data.frame(D)[, pv, drop = FALSE], use = "pairwise.complete.obs")
cmd <- as.data.table(as.table(cm)); setnames(cmd, c("V1", "V2", "r"))
cmd[, `:=`(V1 = factor(V1, levels = pv), V2 = factor(V2, levels = rev(pv)))]
pS1 <- ggplot(cmd, aes(V1, V2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(pS1, "FigS01_Predictor_collinearity", 7, 6)

# =========================================================================
# FIG S3+ - gam.check DIAGNOSTICS
# =========================================================================
cat("\nFig S3+: model diagnostics\n")
diag_keys <- c(Richness = "richness_main", Activity = "activity_main",
               RichnessRoof = "richness_roof", Tensor = "richness_tensor",
               PCoA1 = "pcoa1", PCoA2 = "pcoa2")
i <- 3
for (nm in names(diag_keys)) {
  m <- get_model(diag_keys[[nm]], quiet = TRUE); if (is.null(m)) next
  png(file.path(fig_dir, sprintf("FigS%02d_Diagnostics_%s.png", i, nm)),
      width = 1800, height = 1600, res = 200)
  par(mfrow = c(2, 2)); gam.check(as_gam(m)); dev.off()
  cat("  saved: diagnostics", nm, "\n"); i <- i + 1
}

cat("\nDone. Figures in:", normalizePath(fig_dir), "\n")
cat("Modelling factors were never altered - only display labels were added,\n")
cat("so every coefficient reported in the text remains valid.\n")
