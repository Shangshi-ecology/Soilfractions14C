# =============================================================================
# Figure 1 — Δ¹⁴C of paired POM and MAOM density fractions, by depth.
#
# This script reproduces Figure 1 of the manuscript. It loads the ISRaD
# flat-fraction product, pairs free-light and occluded-light POM with the
# corresponding heavy (MAOM) fraction separated from the same physical
# soil aliquot, and fits a three-level linear mixed-effects model to the
# within-aliquot paired difference (POM − MAOM) at each depth bin.
#
# Output files:
#   Figure1_DepthProfile.pdf / .png        4-panel composite figure
#   Figure1_DepthProfile_table.csv         numerical summary of every cell
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
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 1. Load the ISRaD flat-fraction product
# -----------------------------------------------------------------------------
# The flat-fraction table contains one row per measured density fraction with
# the hierarchical identifiers that allow fractions separated from the same
# physical aliquot to be linked.  Required columns include the study/site/
# profile/layer identifiers, the layer depth bounds, the fraction property
# label, and Δ¹⁴C.
ISRAD_CSV <- "ISRaD_database_files/ISRaD_extra_flat_fraction_v 2.9.9.2025-08-14.csv"
stopifnot(file.exists(ISRAD_CSV))
raw <- read_csv(ISRAD_CSV, show_col_types = FALSE, guess_max = 50000)

# -----------------------------------------------------------------------------
# 2. Map ISRaD density-fraction labels onto the three operational pools
# -----------------------------------------------------------------------------
# ISRaD v2.9.9 records density-fractionation results with a controlled
# vocabulary in `frc_property`.  Density-fractionation rows take one of
# four values: "free light", "occluded light", "heavy", and "sand".
# We use the first three as our operational POM/MAOM pools; the very small
# number of "sand" rows are coarse-fraction products that do not map onto
# any of these pools and are returned as NA so they are excluded from
# pairing.  Throughout the manuscript, "fPOM", "oPOM" and "MAOM" refer to
# the free light, occluded light, and heavy density fractions exclusively.
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
# We retain density-fractionation measurements with reported Δ¹⁴C and known
# depth bounds, and exclude above-ground litter (layer midpoint < 0 cm).
d <- raw |>
  filter(frc_scheme == "density",
         !is.na(frc_14c),
         !is.na(lyr_top), !is.na(lyr_bot)) |>
  mutate(lyr_mid    = (lyr_top + lyr_bot) / 2,
         pool       = classify_pool(frc_property),
         entry_name = as.character(entry_name),
         site_id    = paste0(entry_name, "::", site_name),
         profile_id = paste0(site_id,   "::", pro_name)) |>
  filter(lyr_mid >= 0, !is.na(pool))

# -----------------------------------------------------------------------------
# 4. Depth bins and the "All paired" overall row
# -----------------------------------------------------------------------------
# Layer midpoints are assigned to 10-cm bins from 0 to 100 cm; observations
# below 100 cm are pooled into a single ">100 cm" bin because of sparse data
# at depth.  The `OVERALL` row pools every observation without depth
# stratification.  
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
# For every stratum (depth bin, or the pooled "All" row), we fit:
#
#     y ~ 1 + (1 | entry / site / profile)
#
# Random effects nested as profiles within sites within studies are the
# maximal structure justified by the design.  Following Barr et al. (2013)
# and Bates et al. (2015), we keep this structure even when lme4 reports a
# singular fit; the `check.conv.singular = ignore` setting implements this.
# A sequence of simpler RE structures (entry/site → entry → OLS) is only
# used as a hard-error backstop, entered when lmer() cannot return a model
# at all (PIRLS failure or n too small for the requested structure).  The
# `model` column of the output CSV records which level was actually used
# at each cell; in practice almost every cell is fit with the full 3-level
# LMM.
#
# Confidence intervals are Wald 95 % (mean ± 1.96 × SE) and the matching
# Wald-z p-value is reported, so the CI bound and the significance label
# are algebraically consistent (CI excludes 0 ⇔ p < 0.05).  Strata with
# only one contributing study (k = 1) carry no between-study variance
# information; we report the sample mean only with no CI and no p.
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
# density fraction was separated.  This is the standard within-sample test
# for paired-fraction comparisons.  Where a study reports a pool as
# multiple density sub-cuts (e.g. heavy at > 1.6, > 1.8, > 2.0 g cm-3) or
# as analytical replicates, the rows for that pool at that aliquot are
# averaged before pivoting so the wide table has exactly one value per
# (aliquot, pool).  Differences are then computed within each aliquot as
# Δ¹⁴C(POM) − Δ¹⁴C(MAOM), separately for the two POM types.
pair_key <- c("entry_name", "site_name", "pro_name", "lyr_name", "frc_input")

agg <- d |>
  group_by(across(all_of(c(pair_key, "pool")))) |>
  summarise(frc_14c = mean(frc_14c, na.rm = TRUE),
            lyr_mid = mean(lyr_mid, na.rm = TRUE),
            .groups = "drop")

wide <- agg |>
  pivot_wider(id_cols = all_of(pair_key),
              names_from = pool, values_from = frc_14c) |>
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

# Paired long-format tables — used both for the LMM-based per-bin profile
# means and for the raw-point overlay in the profile panels.  Restricting
# the profile-panel input to paired observations keeps the dot density on
# (a, c) consistent with the paired-difference panels on (b, d).
paired_long <- function(wd, pom_col) {
  sub <- wd |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM))
  bind_rows(
    sub |> transmute(frc_14c    = .data[[pom_col]], bin,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = pom_col),
    sub |> transmute(frc_14c    = MAOM, bin,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = "MAOM")
  )
}

# -----------------------------------------------------------------------------
# 7. Per-stratum summary functions
# -----------------------------------------------------------------------------
# `pool_summary` returns the LMM-estimated mean Δ¹⁴C of each pool per depth
# bin and an "All" row pooled across depth (used for the left-hand profile
# panels).  `diff_summary` returns the LMM-estimated within-aliquot
# difference per depth bin and an "All" row (used for the right-hand diff
# panels).  Both keep "All" as a factor level appended to the end of
# `bin_labs` so the discrete y-axis can place it at the bottom of the
# figure via `limits = rev(bin_labs)`.
pool_summary <- function(long_df, pool_filter) {
  sub <- long_df |> filter(pool == pool_filter)
  per_bin <- sub |>
    group_by(bin) |>
    summarise(s = list(fit_nested(frc_14c, entry_name, site_id, profile_id)),
              .groups = "drop")
  all_row <- tibble(bin = factor(OVERALL, levels = bin_labs),
                    s   = list(fit_nested(sub$frc_14c, sub$entry_name,
                                           sub$site_id, sub$profile_id)))
  bind_rows(per_bin, all_row) |>
    mutate(pool = pool_filter,
           mean = map_dbl(s, "mean"), lo = map_dbl(s, "lo"),
           hi   = map_dbl(s, "hi"),   p  = map_dbl(s, "p"),
           model = map_chr(s, "model"),
           k = map_int(s, "k"), n = map_int(s, "n")) |>
    select(-s)
}

diff_summary <- function(wd, pom_col) {
  wd2 <- wd |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM)) |>
    mutate(diff = .data[[pom_col]] - MAOM)
  per_bin <- wd2 |>
    group_by(bin) |>
    summarise(s = list(fit_nested(diff, entry, site, profile)),
              .groups = "drop")
  all_row <- tibble(bin = factor(OVERALL, levels = bin_labs),
                    s   = list(fit_nested(wd2$diff, wd2$entry,
                                           wd2$site, wd2$profile)))
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

paired_fp       <- paired_long(wide, "fPOM")
paired_op       <- paired_long(wide, "oPOM")
pool_summary_fp <- bind_rows(pool_summary(paired_fp, "fPOM"),
                              pool_summary(paired_fp, "MAOM"))
pool_summary_op <- bind_rows(pool_summary(paired_op, "oPOM"),
                              pool_summary(paired_op, "MAOM"))
free_d <- diff_summary(wide, "fPOM")
occ_d  <- diff_summary(wide, "oPOM")

# Long-format paired differences — used as the raw-point overlay on the
# diff panels (b, d).  Carries the depth-bin factor so it lines up with
# the LMM summaries.
wide_free <- wide |> filter(!is.na(fPOM), !is.na(MAOM)) |>
  mutate(value = fPOM - MAOM)
wide_occ  <- wide |> filter(!is.na(oPOM), !is.na(MAOM)) |>
  mutate(value = oPOM - MAOM)

# Numerical summary written to disk so every line in the figure is
# traceable to a CSV row.  The `model` column records which LMM level
# (or OLS) was actually used at each cell.
bind_rows(
  pool_summary_fp |> mutate(panel = "fPOM/MAOM profile"),
  pool_summary_op |> mutate(panel = "oPOM/MAOM profile"),
  free_d  |> mutate(pool = "fPOM-MAOM", panel = "free diff"),
  occ_d   |> mutate(pool = "oPOM-MAOM", panel = "occluded diff")
) |> write_csv("Figure1_DepthProfile_table.csv")

# -----------------------------------------------------------------------------
# 8. Plot styling
# -----------------------------------------------------------------------------
# Colour palette follows a colour-blind-safe scheme: red = free POM,
# orange = occluded POM, blue = MAOM, dark grey = within-aliquot paired
# difference.  Marker fill on the diff panels encodes LMM significance
# (filled when p < 0.05, hollow otherwise).
C_FPOM <- "#D62728"
C_OPOM <- "#D55E00"
C_MAOM <- "#0072B2"
C_DIFF <- "#444444"
C_ZERO <- "#8B7355"
XRANGE_PROFILE <- c(-1050, 250)

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
# (= "All") sits at the bottom; we shade every other row to help the eye
# track across panels.
stripe_df <- tibble(pos  = seq_along(bin_labs)) |>
  filter(pos %% 2 == 1)

# Thin horizontal separator drawn between ">100" (y = 2) and the pooled
# "All" row (y = 1).
SEP_Y <- 1.5

# -----------------------------------------------------------------------------
# 9. Profile panels (a, c) — Δ¹⁴C means of POM and MAOM per depth bin
# -----------------------------------------------------------------------------
# Raw paired observations are jittered vertically within their depth bin
# (excluded for the "All" row, which is summary-only).  LMM summary markers
# and CI bars are drawn on top in the pool colour.
plot_profile <- function(pools, raw_long, summary_df,
                          panel_tag, panel_title) {
  raw_sub <- raw_long   |> filter(pool %in% pools)
  ps_sub  <- summary_df |> filter(pool %in% pools)
  pal     <- c(fPOM = C_FPOM, oPOM = C_OPOM, MAOM = C_MAOM)

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE,
              fill = "#F4F4F2", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.85) +
    geom_hline(yintercept = SEP_Y, color = "#888", linewidth = 0.3) +
    geom_jitter(data = raw_sub |> filter(bin != OVERALL),
                aes(x = frc_14c, y = bin, color = pool),
                width = 0, height = 0.18,
                size = 0.35, alpha = 0.30) +
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
    coord_cartesian(xlim = XRANGE_PROFILE, clip = "off") +
    labs(x = expression(Delta^14*"C (‰)"), y = "Depth (cm)",
         tag = panel_tag, title = panel_title) +
    theme_nature()
}

# -----------------------------------------------------------------------------
# 10. Difference panels (b, d) — POM − MAOM per depth bin
# -----------------------------------------------------------------------------
# Raw paired-difference observations are jittered within their depth bin
# (the "All" row is summary-only).  Marker fill on the LMM summary point
# encodes significance: filled when p < 0.05, hollow otherwise.  Counts
# (k = number of studies, n = number of paired observations) and the
# significance star are annotated to the right of each CI bar.
plot_diff <- function(diff_summary, wide_pair_df, panel_tag, panel_title) {
  ds <- diff_summary |>
    mutate(sig_zero = !is.na(p) & p < 0.05,
           fill_col = ifelse(sig_zero, C_DIFF, "white"),
           label    = sprintf("k=%d, n=%d %s", k, n, sig))

  finite  <- ds |> filter(!is.na(hi), !is.na(lo))
  xmax    <- max(c(finite$hi, wide_pair_df$value, 200), na.rm = TRUE)
  xmin    <- min(c(finite$lo, wide_pair_df$value, -200), na.rm = TRUE)
  xpad    <- max(60, 0.05 * (xmax - xmin))
  x_annot <- xmax + xpad

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE,
              fill = "#F4F4F2", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.85) +
    geom_hline(yintercept = SEP_Y, color = "#888", linewidth = 0.3) +
    geom_jitter(data = wide_pair_df |> filter(bin != OVERALL),
                aes(x = value, y = bin),
                width = 0, height = 0.18,
                color = C_DIFF, size = 0.35, alpha = 0.32) +
    geom_linerange(data = ds,
                    aes(xmin = lo, xmax = hi, y = bin),
                    color = C_DIFF, linewidth = 0.6, na.rm = TRUE) +
    geom_point(data = ds,
               aes(x = mean, y = bin, fill = I(fill_col)),
               shape = 21, color = "black",
               size = 1.6, stroke = 0.35, na.rm = TRUE) +
    geom_text(data = ds, aes(x = x_annot, y = bin, label = label),
              hjust = 0, size = 1.95, color = "grey25", na.rm = TRUE) +
    scale_y_discrete(limits = rev(bin_labs)) +
    coord_cartesian(xlim = c(xmin - xpad, x_annot + 240), clip = "off") +
    labs(x = expression("POM − MAOM  "*Delta^14*"C (‰)"),
         y = "Depth (cm)",
         tag = panel_tag, title = panel_title) +
    theme_nature()
}

# -----------------------------------------------------------------------------
# 11. Legend strip
# -----------------------------------------------------------------------------
# A small helper that renders a stand-alone legend strip; the underlying
# data points are placed off-canvas via NA coordinates so only the legend
# appears.
make_legend_panel <- function(items, colors) {
  df <- tibble(what = factor(items, levels = items),
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
# 12. Assemble the 4-panel composite
# -----------------------------------------------------------------------------
# Rows: (a, b) = free POM vs MAOM; (c, d) = occluded POM vs MAOM.
# Columns: left = depth profile of both pools; right = paired difference.
pa <- plot_profile(c("fPOM", "MAOM"), paired_fp, pool_summary_fp,
                    "a", "Free POM vs. MAOM (paired)")
pb <- plot_diff(   free_d, wide_free, "b", "Free POM − MAOM")
pc <- plot_profile(c("oPOM", "MAOM"), paired_op, pool_summary_op,
                    "c", "Occluded POM vs. MAOM (paired)")
pd <- plot_diff(   occ_d,  wide_occ,  "d", "Occluded POM − MAOM")

leg <- make_legend_panel(
  items  = c("Free POM", "Occluded POM", "MAOM",
              "POM − MAOM (paired)", "ns (CI crosses 0)"),
  colors = c("Free POM" = C_FPOM, "Occluded POM" = C_OPOM,
              "MAOM" = C_MAOM, "POM − MAOM (paired)" = C_DIFF,
              "ns (CI crosses 0)" = "white"))

fig <- (pa | pb) / (pc | pd) +
  plot_layout(widths = c(1, 1.05)) &
  theme(plot.margin = margin(2, 4, 2, 2))

final <- wrap_elements(full = fig) /
         wrap_elements(full = leg) +
         plot_layout(heights = c(20, 1))

ggsave("Figure1_DepthProfile.pdf", final,
       width = 6, height = 6, device = cairo_pdf)
ggsave("Figure1_DepthProfile.png", final,
       width = 6, height = 6, dpi = 600, bg = "white")
message("Wrote Figure1_DepthProfile.pdf, .png and _table.csv")
