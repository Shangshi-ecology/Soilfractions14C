# =============================================================================
# Figures S2 & S3 — Stoichiometric and stable-isotope depth profiles of paired
#                   POM and MAOM density fractions.
#
#   Supplementary Figure S2: free POM vs MAOM
#   Supplementary Figure S3: occluded POM vs M
AOM
#
## Authors:   Shangshi Liu, Jonathan Sanderman, Mark A. Bradford
# Citation:  [add manuscript citation when accepted]
#
# These scripts reproduce Supplementary Figures S2 and S3 of the manuscript.
# They apply the pairing logic and statistical model of Figure1_DepthProfile.R
# to three non-radiocarbon properties — δ¹³C, δ¹⁵N and the C:N ratio — to test
# whether the measured fractions carry the compositional signatures expected
# if POM and MAOM form through distinct processes.
#
# Each figure is a 3 x 2 composite:
#   Row 1 (a, b):  δ¹³C  depth profile  +  paired POM − MAOM difference
#   Row 2 (c, d):  δ¹⁵N  depth profile  +  paired difference
#   Row 3 (e, f):  C:N   depth profile  +  paired difference
#
# Sample-inclusion criteria are set in the `Inclusion criteria` block below
# (section 0).  They follow Figure1_DepthProfile.R with one deliberate
# exception, documented there: no sampling-year restriction is applied.
#
# Output files:
#   FigureS2_Stoichiometry_FreePOM.pdf / .png       6-panel composite
#   FigureS2_Stoichiometry_FreePOM_table.csv        numerical summary
#   FigureS3_Stoichiometry_OccludedPOM.pdf / .png   6-panel composite
#   FigureS3_Stoichiometry_OccludedPOM_table.csv    numerical summary
#
# Dependencies:
#   R >= 4.3
#   tidyverse, lme4, lmerTest, patchwork
#
# Input data:
#   ISRaD v2.9.9 flat-fraction product, downloaded from
#   https://www.soilradiocarbon.org/database-1
#   and placed in `./ISRaD_database_files/` with the filename
#   `ISRaD_extra_flat_fraction_v 2.9.9.2025-08-14.csv`.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(lme4)
  library(lmerTest)
})

# -----------------------------------------------------------------------------
# 0. Inclusion criteria
# -----------------------------------------------------------------------------
# Layers thicker than MAX_THICK_CM are excluded, as in Figure1_DepthProfile.R,
# so that no single value averages across several of the 10-cm depth
# increments used for binning.
#
# No sampling-year restriction is applied here.  In the radiocarbon analysis,
# samples collected before 2000 are excluded because the non-monotonic
# atmospheric bomb-¹⁴C curve makes Δ¹⁴C depend on collection date.  δ¹³C, δ¹⁵N
# and C:N carry no such time-dependent atmospheric signal, so restricting the
# sampling window would discard informative observations without removing any
# comparable artefact.  The full record is therefore retained for these
# properties.
MAX_THICK_CM <- 30      # maximum layer thickness retained (cm)

# Plausibility limits on the reported C:N ratio, applied per pool.  Values
# above these thresholds are treated as reporting or unit errors rather than
# real measurements and set to NA.
CN_MAX_MAOM <- 100
CN_MAX_POM  <- 200

# -----------------------------------------------------------------------------
# 1. Load the ISRaD flat-fraction product
# -----------------------------------------------------------------------------
ISRAD_CSV <- "ISRaD_database_files/ISRaD_extra_flat_fraction_v 2.9.9.2025-08-14.csv"
stopifnot(file.exists(ISRAD_CSV))
raw <- read_csv(ISRAD_CSV, show_col_types = FALSE, guess_max = 50000)

# -----------------------------------------------------------------------------
# 2. Map ISRaD density-fraction labels onto the three operational pools
# -----------------------------------------------------------------------------
# ISRaD v2.9.9 records density-fractionation results with a controlled
# vocabulary in `frc_property`: "free light", "occluded light", "heavy",
# and a small number of "sand" (coarse-fraction) rows.  We map the first
# three onto the operational pools used throughout the manuscript; the
# "sand" rows are returned as NA so they are excluded from pairing.
classify_pool <- function(prop) {
  case_when(
    prop == "free light"     ~ "fPOM",
    prop == "occluded light" ~ "oPOM",
    prop == "heavy"          ~ "MAOM",
    TRUE                     ~ NA_character_
  )
}

# -----------------------------------------------------------------------------
# 3. Filter and prepare the analysis frame
# -----------------------------------------------------------------------------
# We retain density-fractionation measurements with known depth bounds, and
# exclude organic horizons, above-ground litter, and layers thicker than
# MAX_THICK_CM.
#
# `lyr_top >= 0` excludes above-ground litter: a layer straddling the
# organic-mineral interface has lyr_top < 0 but can still have lyr_mid >= 0,
# so the test is applied to the upper bound rather than to the midpoint.
# `lyr_all_org_neg` flags layers recorded as entirely organic; it is stored
# as the character value "yes" or as NA.
d <- raw |>
  filter(frc_scheme == "density",
         is.na(lyr_all_org_neg) | lyr_all_org_neg != "yes",
         !is.na(lyr_top), !is.na(lyr_bot)) |>
  mutate(lyr_mid    = (lyr_top + lyr_bot) / 2,
         thickness  = lyr_bot - lyr_top,
         pool       = classify_pool(frc_property),
         entry_name = as.character(entry_name),
         site_id    = paste0(entry_name, "::", site_name),
         profile_id = paste0(site_id,   "::", pro_name)) |>
  filter(lyr_top   >= 0,
         thickness <= MAX_THICK_CM,
         !is.na(pool))

# The C:N ratio is taken from the curated `frc_c_to_n` column rather than
# computed from `frc_c_perc` and `frc_n_tot`.  Nitrogen concentrations in
# ISRaD are not unit-standardised across contributing studies, so ratios
# derived from them are unreliable for a large share of heavy-fraction rows;
# `frc_c_to_n` is the contributor-reported, internally consistent ratio.
# Non-positive and non-finite values, and values beyond the per-pool
# plausibility limits set in section 0, are treated as missing.
d <- d |>
  mutate(cn = as.numeric(frc_c_to_n),
         cn = case_when(
           !is.finite(cn) | cn <= 0                        ~ NA_real_,
           pool == "MAOM"              & cn > CN_MAX_MAOM  ~ NA_real_,
           pool %in% c("fPOM", "oPOM") & cn > CN_MAX_POM   ~ NA_real_,
           TRUE                                             ~ cn))

# -----------------------------------------------------------------------------
# 4. Depth bins and the "All paired" overall row
# -----------------------------------------------------------------------------
# Layer midpoints are assigned to 10-cm bins from 0 to 100 cm; observations
# below 100 cm are pooled into a single ">100 cm" bin because of sparse data
# at depth.  The `OVERALL` row pools every observation without depth
# stratification and is drawn at the bottom of each panel.
bin_edges <- c(seq(0, 100, by = 10), Inf)
bin_core  <- c(paste0(seq(0, 90, 10), "-", seq(10, 100, 10)), ">100")
OVERALL   <- "All"
bin_labs  <- c(bin_core, OVERALL)

d <- d |>
  mutate(bin = cut(lyr_mid, breaks = bin_edges, labels = bin_core,
                    right = FALSE, include.lowest = TRUE)) |>
  filter(!is.na(bin)) |>
  mutate(bin = factor(bin, levels = bin_labs))

# -----------------------------------------------------------------------------
# 5. Maximal three-level LMM (singular fits accepted, hard-error backstop)
# -----------------------------------------------------------------------------
# Identical statistical machinery to Figure1_DepthProfile.R.  For every
# stratum we fit
#
#     y ~ 1 + (1 | entry / site / profile)
#
# Random effects nested as profiles within sites within studies are the
# maximal structure justified by the design.  Following Barr et al. (2013)
# and Bates et al. (2015), this structure is retained even when lme4 reports
# a singular fit.  A sequence of simpler structures (entry/site → entry →
# OLS) is entered only as a hard-error backstop, when lmer() cannot return a
# model at all.  Confidence intervals are Wald 95 % (mean ± 1.96 × SE) with
# the matching Wald-z p-value, so the CI bound and the significance label are
# algebraically consistent.  Strata with one contributing study (k = 1) carry
# no between-study variance information; the sample mean is reported with no
# CI and no p-value.
fit_nested <- function(y, entry, site, profile) {
  ok <- !is.na(y);  y <- y[ok]; entry <- entry[ok]
  site <- site[ok]; profile <- profile[ok]
  n <- length(y);   k <- length(unique(entry))

  if (n < 2)
    return(list(mean = if (n == 1) y else NA, lo = NA, hi = NA, p = NA,
                 model = "n<2 (none)", n = n, k = k))
  if (k < 2)
    return(list(mean = mean(y), lo = NA, hi = NA, p = NA,
                 model = "k=1 (no inference)", n = n, k = k))

  dat <- tibble(y = y, entry = factor(entry), site = factor(site),
                profile = factor(profile))
  fit_lmer <- function(formula_str) {
    tryCatch(suppressMessages(suppressWarnings(
      lmer(as.formula(formula_str), data = dat,
           control = lmerControl(
             check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
             optimizer = "bobyqa",
             optCtrl   = list(maxfun = 1e5))))),
      error = function(e) NULL)
  }
  m <- fit_lmer("y ~ 1 + (1 | entry/site/profile)"); used <- "lmm3"
  if (is.null(m)) { m <- fit_lmer("y ~ 1 + (1 | entry/site)"); used <- "lmm2" }
  if (is.null(m)) { m <- fit_lmer("y ~ 1 + (1 | entry)");      used <- "lmm1" }
  if (is.null(m)) {
    m  <- lm(y ~ 1, data = dat)
    co <- summary(m)$coefficients
    mu <- co[1, 1]; se <- co[1, 2]
    return(list(mean = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se,
                 p = 2 * (1 - pnorm(abs(mu / se))),
                 model = "ols", n = n, k = k))
  }
  co <- summary(m)$coefficients
  mu <- co[1, 1]; se <- co[1, 2]
  list(mean = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se,
       p = 2 * (1 - pnorm(abs(mu / se))), model = used, n = n, k = k)
}

# -----------------------------------------------------------------------------
# 6. Within-aliquot pairing
# -----------------------------------------------------------------------------
# Fractions are matched on the hierarchical key
#     entry / site / profile / layer / frc_input,
# where `frc_input` identifies the physical soil aliquot from which each
# density fraction was separated.  Where a study reports a pool as multiple
# density sub-cuts or as analytical replicates, the rows for that pool at that
# aliquot are averaged before pivoting, so the wide table holds exactly one
# value per (aliquot, pool).  Differences are then computed within each
# aliquot, separately for the two POM types.
pair_key <- c("entry_name", "site_name", "pro_name", "lyr_name", "frc_input")

build_pair_table <- function(var_col) {
  agg <- d |>
    group_by(across(all_of(c(pair_key, "pool")))) |>
    summarise(value   = mean(.data[[var_col]], na.rm = TRUE),
              lyr_mid = mean(lyr_mid, na.rm = TRUE),
              .groups = "drop")
  agg |>
    pivot_wider(id_cols = all_of(pair_key),
                names_from = pool, values_from = value) |>
    left_join(agg |> group_by(across(all_of(pair_key))) |>
                summarise(lyr_mid = mean(lyr_mid, na.rm = TRUE),
                          .groups = "drop"),
              by = pair_key) |>
    mutate(bin = cut(lyr_mid, breaks = bin_edges, labels = bin_core,
                      right = FALSE, include.lowest = TRUE)) |>
    filter(!is.na(bin)) |>
    mutate(bin     = factor(bin, levels = bin_labs),
           entry   = entry_name,
           site    = paste0(entry_name, "::", site_name),
           profile = paste0(site,       "::", pro_name))
}

# Paired long-format table — used both for the LMM per-bin profile means and
# for the raw-point overlay in the profile panels.  Restricting the profile
# panels to paired observations keeps their point density consistent with the
# paired-difference panels beside them.
paired_long <- function(wide, pom_col) {
  sub <- wide |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM))
  bind_rows(
    sub |> transmute(value = .data[[pom_col]], bin,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = pom_col),
    sub |> transmute(value = MAOM, bin,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = "MAOM")
  )
}

# -----------------------------------------------------------------------------
# 7. Per-stratum summary functions
# -----------------------------------------------------------------------------
# `pool_summary` returns the LMM-estimated mean of each pool per depth bin,
# plus an "All" row pooled across depth (left-hand profile panels).
# `diff_summary` returns the LMM-estimated within-aliquot difference per depth
# bin and an "All" row (right-hand difference panels).
pool_summary <- function(long_df, pool_filter) {
  sub <- long_df |> filter(pool == pool_filter)
  per_bin <- sub |>
    group_by(bin) |>
    summarise(s = list(fit_nested(value, entry_name, site_id, profile_id)),
              .groups = "drop")
  all_row <- tibble(bin = factor(OVERALL, levels = bin_labs),
                    s   = list(fit_nested(sub$value, sub$entry_name,
                                           sub$site_id, sub$profile_id)))
  bind_rows(per_bin, all_row) |>
    mutate(pool = pool_filter,
           mean = map_dbl(s, "mean"), lo = map_dbl(s, "lo"),
           hi   = map_dbl(s, "hi"),   p  = map_dbl(s, "p"),
           model = map_chr(s, "model"),
           k = map_int(s, "k"), n = map_int(s, "n")) |>
    select(-s)
}

diff_summary <- function(wide, pom_col) {
  wd <- wide |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM)) |>
    mutate(diff = .data[[pom_col]] - MAOM)
  per_bin <- wd |>
    group_by(bin) |>
    summarise(s = list(fit_nested(diff, entry, site, profile)),
              .groups = "drop")
  all_row <- tibble(bin = factor(OVERALL, levels = bin_labs),
                    s   = list(fit_nested(wd$diff, wd$entry,
                                           wd$site, wd$profile)))
  bind_rows(per_bin, all_row) |>
    mutate(mean = map_dbl(s, "mean"), lo = map_dbl(s, "lo"),
           hi   = map_dbl(s, "hi"),   p  = map_dbl(s, "p"),
           model = map_chr(s, "model"),
           k = map_int(s, "k"), n = map_int(s, "n"),
           sig = case_when(is.na(p) ~ "",
                            p < 0.001 ~ "***", p < 0.01 ~ "**",
                            p < 0.05 ~ "*",   TRUE ~ "ns")) |>
    select(-s)
}

# -----------------------------------------------------------------------------
# 8. Plot styling
# -----------------------------------------------------------------------------
# Colour palette follows Figures 1 and 2: red = free POM, orange = occluded
# POM, blue = MAOM, dark grey = within-aliquot paired difference.  Marker fill
# on the difference panels encodes LMM significance (filled when p < 0.05,
# hollow otherwise).
C_FPOM <- "#D62728"
C_OPOM <- "#D55E00"
C_MAOM <- "#0072B2"
C_DIFF <- "#444444"
C_ZERO <- "#8B7355"

# One entry per figure row: the column to plot, its axis labels and the panel
# title.
#
# Axis ranges are NOT fixed here.  Both the profile and the difference panels
# take their x-range from the data actually plotted — every raw observation
# and every confidence bound — so that no point and no interval is ever drawn
# outside the panel or hidden by a scale limit.  A fixed range would silently
# omit extreme values while still letting them influence the plotted means.
# `xrange` may be set to a numeric pair to override the automatic range for a
# row; leave it NULL for the data-derived range.
VAR_SPEC <- list(
  list(var = "frc_13c", title = "δ¹³C",
       lab_prof = "δ¹³C (‰)",     lab_diff = "POM − MAOM δ¹³C (‰)",
       xrange = NULL),
  list(var = "frc_15n", title = "δ¹⁵N",
       lab_prof = "δ¹⁵N (‰)",     lab_diff = "POM − MAOM δ¹⁵N (‰)",
       xrange = NULL),
  list(var = "cn",      title = "C : N",
       lab_prof = "C : N (mass)", lab_diff = "POM − MAOM  C:N",
       xrange = NULL)
)

# Helper: the range spanned by every value that will be drawn, with a small
# margin so markers are not clipped by the panel edge.
data_range <- function(..., pad = 0.04) {
  v <- c(...)
  v <- v[is.finite(v)]
  r <- range(v)
  r + c(-1, 1) * max(pad * diff(r), 1e-8)
}

theme_nature <- function() {
  theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(axis.line   = element_line(linewidth = 0.4),
          axis.ticks  = element_line(linewidth = 0.35),
          axis.title  = element_text(size = 7),
          axis.text   = element_text(size = 6.5),
          plot.title  = element_text(size = 7.5, color = "grey25",
                                       face = "plain"),
          plot.tag    = element_text(face = "bold", size = 9),
          panel.grid.major.x = element_line(color = "grey94",
                                              linewidth = 0.25),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 6.5),
          legend.key.size = unit(0.3, "cm"))
}

# Alternating row-stripe background.  The y-axis is reversed so position 1
# (= "All") sits at the bottom; every other row is shaded to help the eye
# track across panels.
stripe_df <- tibble(pos = seq_along(bin_labs)) |> filter(pos %% 2 == 1)

# Thin horizontal separator drawn between ">100" (y = 2) and the pooled
# "All" row (y = 1).
SEP_Y <- 1.5

# -----------------------------------------------------------------------------
# 9. Profile panels (a, c, e) — pool means per depth bin
# -----------------------------------------------------------------------------
plot_profile <- function(pom_pool, pom_color, raw_long, summary_df,
                          xrange, panel_tag, var_title, x_lab) {
  pal     <- setNames(c(pom_color, C_MAOM), c(pom_pool, "MAOM"))
  raw_sub <- raw_long   |> filter(pool %in% c(pom_pool, "MAOM"))
  ps_sub  <- summary_df |> filter(pool %in% c(pom_pool, "MAOM"))

  # A zero reference line is drawn only where zero falls inside the panel.
  # It is meaningful for δ¹⁵N and C:N but not for δ¹³C, whose values are all
  # strongly negative.
  zero_layer <- if (0 >= xrange[1] && 0 <= xrange[2]) {
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.85)
  } else NULL

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE, fill = "#F4F4F2", alpha = 0.5) +
    zero_layer +
    geom_hline(yintercept = SEP_Y, color = "#888", linewidth = 0.3) +
    geom_jitter(data = raw_sub |> filter(bin != OVERALL),
                aes(x = value, y = bin, color = pool),
                width = 0, height = 0.18, size = 0.35, alpha = 0.30) +
    geom_linerange(data = ps_sub,
                    aes(xmin = lo, xmax = hi, y = bin, color = pool),
                    linewidth = 0.6, na.rm = TRUE) +
    geom_point(data = ps_sub,
               aes(x = mean, y = bin, fill = pool),
               shape = 21, color = "black", size = 1.6, stroke = 0.25,
               na.rm = TRUE) +
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual( values = pal, guide = "none") +
    scale_y_discrete(limits = rev(bin_labs)) +
    coord_cartesian(xlim = xrange, clip = "on") +
    labs(x = x_lab, y = "Depth (cm)",
         tag = panel_tag, title = var_title) +
    theme_nature()
}

# -----------------------------------------------------------------------------
# 10. Difference panels (b, d, f) — POM − MAOM per depth bin
# -----------------------------------------------------------------------------
# The x-axis range is set from the union of CI bounds and raw points, with
# space reserved on the right for the count annotation.  The annotation lies
# outside the data range, so the limit is applied through coord_cartesian()
# rather than through a scale limit, which would drop the text.
plot_diff <- function(diff_df, wide_pair, pom_col, panel_tag, x_lab) {
  raw <- wide_pair |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM)) |>
    mutate(value = .data[[pom_col]] - MAOM)

  ds <- diff_df |>
    mutate(sig_zero = !is.na(p) & p < 0.05,
           fill_col = ifelse(sig_zero, C_DIFF, "white"),
           label    = sprintf("k=%d,n=%d %s", k, n, sig))

  # The range spans every raw paired difference and every confidence bound,
  # so nothing plotted falls outside the panel.
  finite  <- ds |> filter(!is.na(hi), !is.na(lo))
  rng     <- data_range(finite$hi, finite$lo, raw$value)
  xmin    <- rng[1]; xmax <- rng[2]
  xpad    <- max((xmax - xmin) * 0.07, 0.5)
  x_annot <- xmax + xpad

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE, fill = "#F4F4F2", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.85) +
    geom_hline(yintercept = SEP_Y, color = "#888", linewidth = 0.3) +
    geom_jitter(data = raw, aes(x = value, y = bin),
                width = 0, height = 0.18,
                color = C_DIFF, size = 0.35, alpha = 0.32) +
    geom_linerange(data = ds,
                    aes(xmin = lo, xmax = hi, y = bin),
                    color = C_DIFF, linewidth = 0.6, na.rm = TRUE) +
    geom_point(data = ds,
               aes(x = mean, y = bin, fill = I(fill_col)),
               shape = 21, color = "black", size = 1.6, stroke = 0.35,
               na.rm = TRUE) +
    geom_text(data = ds, aes(x = x_annot, y = bin, label = label),
              hjust = 0, size = 1.8, color = "grey25", na.rm = TRUE) +
    scale_y_discrete(limits = rev(bin_labs)) +
    coord_cartesian(xlim = c(xmin - xpad, x_annot + (xmax - xmin) * 0.55),
                     clip = "off") +
    labs(x = x_lab, y = NULL, tag = panel_tag) +
    theme_nature() +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

# -----------------------------------------------------------------------------
# 11. Legend strip
# -----------------------------------------------------------------------------
# A small helper that renders a stand-alone legend strip; the underlying data
# points are placed off-canvas via NA coordinates so only the legend appears.
make_legend <- function(pom_label, pom_color) {
  items  <- c(pom_label, "MAOM", "POM − MAOM (paired)", "ns (CI crosses 0)")
  colors <- setNames(c(pom_color, C_MAOM, C_DIFF, "white"), items)
  df     <- tibble(what = factor(items, levels = items),
                   x = NA_real_, y = NA_real_)
  ggplot(df, aes(x, y, fill = what)) +
    geom_point(shape = 21, color = "black", size = 2.4, stroke = 0.3,
                na.rm = TRUE) +
    scale_fill_manual(values = colors, name = NULL) +
    guides(fill = guide_legend(nrow = 1,
                                override.aes = list(size = 2.4))) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void(base_size = 7, base_family = "Helvetica") +
    theme(legend.position    = "bottom",
          legend.text        = element_text(size = 6.5),
          legend.key.size    = unit(0.3, "cm"),
          legend.box.spacing = unit(0, "pt"),
          plot.margin        = margin(0, 0, 0, 0))
}

# -----------------------------------------------------------------------------
# 12. Assemble the 6-panel composite
# -----------------------------------------------------------------------------
# Rows: δ¹³C, δ¹⁵N, C:N.  Columns: left = depth profile of both pools;
# right = paired within-aliquot difference.  The numerical summary of every
# cell is written alongside the figure so each drawn element is traceable to
# a CSV row.
build_six_panel <- function(pom_pool, pom_color, pom_label, out_prefix) {
  panel_letters <- list(c("a", "b"), c("c", "d"), c("e", "f"))
  rows  <- list()
  table <- list()

  for (i in seq_along(VAR_SPEC)) {
    spec   <- VAR_SPEC[[i]]
    wide   <- build_pair_table(spec$var)
    paired <- paired_long(wide, pom_pool)
    ps     <- bind_rows(pool_summary(paired, pom_pool),
                        pool_summary(paired, "MAOM"))
    ds     <- diff_summary(wide, pom_pool)

    # Profile-panel range covers every raw observation and every confidence
    # bound drawn in that row, unless VAR_SPEC overrides it.
    xr <- if (is.null(spec$xrange)) {
      data_range(paired$value, ps$lo, ps$hi, ps$mean)
    } else spec$xrange

    rows[[i]] <-
      plot_profile(pom_pool, pom_color, paired, ps, xr,
                   panel_letters[[i]][1], spec$title, spec$lab_prof) |
      plot_diff(ds, wide, pom_pool,
                panel_letters[[i]][2], spec$lab_diff)

    # Reproducibility check: nothing drawn should fall outside its panel.
    n_out <- sum(paired$value < xr[1] | paired$value > xr[2], na.rm = TRUE)
    if (n_out > 0)
      warning(sprintf("%s / %s: %d observations outside the panel range",
                      out_prefix, spec$title, n_out), call. = FALSE)

    table[[i]] <- bind_rows(
      ps |> mutate(variable = spec$title, panel = "profile"),
      ds |> mutate(variable = spec$title, panel = "difference",
                   pool = paste0(pom_pool, "-MAOM")))
  }

  bind_rows(table) |> write_csv(paste0(out_prefix, "_table.csv"))

  fig <- rows[[1]] / rows[[2]] / rows[[3]] +
    plot_layout(heights = c(1, 1, 1)) &
    theme(plot.margin = margin(2, 4, 2, 2))

  final <- wrap_elements(full = fig) /
           wrap_elements(full = make_legend(pom_label, pom_color)) +
           plot_layout(heights = c(30, 1))

  ggsave(paste0(out_prefix, ".pdf"), final,
         width = 6, height = 7.5, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), final,
         width = 6, height = 7.5, dpi = 600, bg = "white")
  message("Wrote ", out_prefix, ".pdf, .png and _table.csv")
}

build_six_panel("fPOM", C_FPOM, "Free POM",
                "FigureS2_Stoichiometry_FreePOM")

build_six_panel("oPOM", C_OPOM, "Occluded POM",
                "FigureS3_Stoichiometry_OccludedPOM")
