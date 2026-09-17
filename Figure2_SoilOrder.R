# =============================================================================
# Figure 2 — Δ¹⁴C of paired POM and MAOM density fractions, by USDA soil order.
#
# Authors:   Shangshi Liu, Jonathan Sanderman, Mark A. Bradford
# Citation:  [add manuscript citation when accepted]
# Licence:   MIT (code) / CC-BY 4.0 (figure outputs)
#
# This script reproduces Figure 2 of the manuscript. It uses the same data
# subset, pairing logic, and statistical model as Figure 1, but aggregates
# the within-aliquot paired difference (POM − MAOM) by USDA soil order
# (with depths pooled).  Rows are ordered top-to-bottom by descending
# free-POM − MAOM gap.
#
# Sample-inclusion criteria are set in the `Inclusion criteria` block below
# (section 0), matching Figure1_DepthProfile.R.
#
# Output files:
#   Figure2_SoilOrder.pdf / .png           4-panel composite figure
#   Figure2_SoilOrder_table.csv            numerical summary of every cell
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
# 0. Inclusion criteria  (identical to Figure1_DepthProfile.R)
# -----------------------------------------------------------------------------
# Samples are restricted to those collected in or after YEAR_MIN. Layers thicker than
# MAX_THICK_CM are excluded so that no single Δ¹⁴C value averages over an
# excessive depth range.
YEAR_MIN     <- 2000    # earliest sampling year retained
MAX_THICK_CM <- 30      # maximum layer thickness retained (cm)

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
# 3. Filter, prepare hierarchical identifiers, and attach soil order
# -----------------------------------------------------------------------------
# We retain density-fractionation measurements with reported Δ¹⁴C and known
# depth bounds, and exclude organic horizons, above-ground litter, samples
# collected before YEAR_MIN, and layers thicker than MAX_THICK_CM.
#

obs_year_col <- if ("lyr_obs_date_y.x" %in% names(raw)) "lyr_obs_date_y.x" else
                                                        "lyr_obs_date_y"
raw <- raw |> mutate(obs_year = .data[[obs_year_col]])

d <- raw |>
  filter(frc_scheme == "density",
         is.na(lyr_all_org_neg) | lyr_all_org_neg != "yes",
         !is.na(frc_14c),
         !is.na(lyr_top), !is.na(lyr_bot)) |>
  mutate(lyr_mid    = (lyr_top + lyr_bot) / 2,
         thickness  = lyr_bot - lyr_top,
         pool       = classify_pool(frc_property),
         entry_name = as.character(entry_name),
         site_id    = paste0(entry_name, "::", site_name),
         profile_id = paste0(site_id,   "::", pro_name)) |>
  filter(lyr_top   >= 0,
         obs_year  >= YEAR_MIN,
         thickness <= MAX_THICK_CM,
         !is.na(pool))

# USDA soil order is read from `pro_usda_soil_order`, with the legacy
# `pro_soilOrder_USDA` column used as a fallback when the primary column is
# missing.  Source values are trimmed and title-cased; the historical
# misspelling "Spodisols" is re-mapped to the canonical "Spodosols".
soil_col <- ifelse("pro_usda_soil_order" %in% names(d),
                    "pro_usda_soil_order", "pro_soilOrder_USDA")
d <- d |>
  mutate(order = coalesce(.data[[soil_col]], .data[["pro_soilOrder_USDA"]])) |>
  filter(!is.na(order)) |>
  mutate(order = str_to_title(str_trim(order)),
         order = recode(order, Spodisols = "Spodosols"))

# -----------------------------------------------------------------------------
# 4. Maximal three-level LMM (singular fits accepted, hard-error backstop)
# -----------------------------------------------------------------------------
# Identical statistical machinery to Figure 1 — see the comment block in
# Figure1_DepthProfile.R for full motivation.  Briefly:
#   diff ~ 1 + (1 | entry / site / profile)   maximal RE structure
# Following Barr et al. (2013) and Bates et al. (2015), singular fits are
# accepted as-is (`check.conv.singular = ignore`).  Simpler RE structures
# are only entered as a hard-error backstop when lmer() cannot return a
# model at all (PIRLS failure, n too small for the requested structure).
# Wald 95 % CI and matching Wald-z p-value; k = 1 strata report mean only.
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
# 5. Within-aliquot pairing
# -----------------------------------------------------------------------------
# Pair on entry / site / profile / layer / frc_input so that only fractions
# derived from the same physical aliquot are compared.  Soil order is
# carried through alongside the pair so it can be used as the grouping
# variable in the per-order summaries.
pair_key <- c("entry_name", "site_name", "pro_name", "lyr_name", "frc_input")

agg <- d |>
  group_by(across(all_of(c(pair_key, "pool")))) |>
  summarise(frc_14c = mean(frc_14c, na.rm = TRUE),
            lyr_mid = mean(lyr_mid, na.rm = TRUE),
            order   = first(order),
            .groups = "drop")

wide <- agg |>
  pivot_wider(id_cols = all_of(pair_key),
              names_from = pool, values_from = frc_14c) |>
  left_join(agg |> group_by(across(all_of(pair_key))) |>
              summarise(lyr_mid = mean(lyr_mid, na.rm = TRUE),
                        order   = first(order),
                        .groups = "drop"),
            by = pair_key) |>
  filter(!is.na(order)) |>
  mutate(entry   = entry_name,
         site    = paste0(entry_name, "::", site_name),
         profile = paste0(site,       "::", pro_name))

# -----------------------------------------------------------------------------
# 6. Per-soil-order summary functions
# -----------------------------------------------------------------------------
# `pool_summary` returns the LMM-estimated mean Δ¹⁴C of each fraction per
# soil order (used for the left-hand profile panels).
# `diff_summary` returns the LMM-estimated within-aliquot difference per
# soil order (used for the right-hand diff panels).
paired_long <- function(wd, pom_col) {
  sub <- wd |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM))
  bind_rows(
    sub |> transmute(frc_14c    = .data[[pom_col]],
                      order, lyr_mid,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = pom_col),
    sub |> transmute(frc_14c    = MAOM,
                      order, lyr_mid,
                      entry_name = entry, site_id = site,
                      profile_id = profile, pool = "MAOM")
  )
}

pool_summary <- function(long_df) {
  long_df |>
    filter(!is.na(order)) |>
    group_by(pool, order) |>
    summarise(s = list(fit_nested(frc_14c, entry_name, site_id, profile_id)),
              .groups = "drop") |>
    mutate(mean = map_dbl(s, "mean"), lo = map_dbl(s, "lo"),
           hi   = map_dbl(s, "hi"),   p  = map_dbl(s, "p"),
           model = map_chr(s, "model"),
           k = map_int(s, "k"), n = map_int(s, "n")) |>
    select(-s)
}

diff_summary <- function(wd, pom_col) {
  wd |> filter(!is.na(.data[[pom_col]]), !is.na(MAOM)) |>
    mutate(diff = .data[[pom_col]] - MAOM) |>
    filter(!is.na(order)) |>
    group_by(order) |>
    summarise(s = list(fit_nested(diff, entry, site, profile)),
              .groups = "drop") |>
    mutate(mean = map_dbl(s, "mean"), lo = map_dbl(s, "lo"),
           hi   = map_dbl(s, "hi"),   p  = map_dbl(s, "p"),
           model = map_chr(s, "model"),
           k = map_int(s, "k"), n = map_int(s, "n"),
           sig = case_when(is.na(p) ~ "",
                            p < 0.001 ~ "***", p < 0.01 ~ "**",
                            p < 0.05 ~ "*",   TRUE ~ "ns")) |>
    select(-s)
}

paired_fp <- paired_long(wide, "fPOM")
paired_op <- paired_long(wide, "oPOM")
pool_fp   <- pool_summary(paired_fp)
pool_op   <- pool_summary(paired_op)
free_d    <- diff_summary(wide, "fPOM")
occ_d     <- diff_summary(wide, "oPOM")

# Ordering rule: every soil order with at least one paired observation in
# either panel is retained.  Rows are sorted top-to-bottom by descending
# mean free-POM − MAOM gap; soil orders that appear only in the occluded-
# POM panel are placed at the bottom of the figure.
all_orders <- union(free_d$order, occ_d$order)
all_orders <- all_orders[!is.na(all_orders)]
ord_rank   <- free_d |> filter(order %in% all_orders) |>
  arrange(desc(mean)) |> pull(order)
ORDERED_ORDERS <- c(ord_rank, setdiff(all_orders, ord_rank))

# Numerical summary written to disk so every line in the figure is
# traceable to a CSV row.
bind_rows(
  pool_fp |> mutate(panel = "fPOM/MAOM profile"),
  pool_op |> mutate(panel = "oPOM/MAOM profile"),
  free_d  |> mutate(pool = "fPOM-MAOM", panel = "free diff"),
  occ_d   |> mutate(pool = "oPOM-MAOM", panel = "occluded diff")
) |> write_csv("Figure2_SoilOrder_table.csv")

# -----------------------------------------------------------------------------
# 7. Plot styling
# -----------------------------------------------------------------------------
C_FPOM <- "#D62728"
C_OPOM <- "#D55E00"
C_MAOM <- "#0072B2"
C_DIFF <- "#444444"
C_ZERO <- "#8B7355"

# Single goldenrod for all raw-observation points (no depth coding).
RAW_POINT_COL  <- "#F2C14E"

# Axis ranges are derived from the data rather than fixed, so that every raw
# observation and every confidence bound is drawn inside its panel.  A fixed
# or capped range would omit extreme values from the figure while still
# letting them influence the plotted means, which cannot be reconciled by a
# reader inspecting the points.
#
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
          legend.title    = element_text(size = 6.5),
          legend.text     = element_text(size = 6.5),
          legend.key.size = unit(0.3, "cm"))
}

# Alternating row-stripe background so the eye can follow rows across the
# four panels in the composite.
stripe_df <- tibble(pos = seq_len(length(ORDERED_ORDERS)) - 1) |>
  filter(pos %% 2 == 1)

# -----------------------------------------------------------------------------
# 8. Profile panels (a, c) — Δ¹⁴C of each fraction per soil order
# -----------------------------------------------------------------------------
plot_profile <- function(pools, raw_long, summary_df,
                          panel_tag, panel_title) {
  raw_sub <- raw_long |> filter(pool %in% pools, order %in% ORDERED_ORDERS)
  ps_sub  <- summary_df |> filter(pool %in% pools, order %in% ORDERED_ORDERS)
  xrange  <- data_range(raw_sub$frc_14c, ps_sub$lo, ps_sub$hi, ps_sub$mean)
  pal     <- c(fPOM = C_FPOM, oPOM = C_OPOM, MAOM = C_MAOM)
  dy_map  <- setNames(c(-0.15, 0.15), pools)

  # Vertical jitter offset per pool so the two pools are visually separated
  # within each soil-order row.  An additional small jitter is added to the
  # raw points so they do not stack on top of each other.
  raw_sub <- raw_sub |>
    mutate(dy = dy_map[pool],
           yi = match(order, ORDERED_ORDERS) - 1 + dy +
                runif(n(), -0.06, 0.06))
  ps_sub  <- ps_sub  |>
    mutate(dy = dy_map[pool],
           yi = match(order, ORDERED_ORDERS) - 1 + dy)

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE, fill = "#F4F4F2", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.95) +
    geom_point(data = raw_sub, aes(x = frc_14c, y = yi),
                color = RAW_POINT_COL, size = 0.5, alpha = 0.4) +
    geom_linerange(data = ps_sub,
                    aes(xmin = lo, xmax = hi, y = yi),
                    color = "black", linewidth = 0.5, na.rm = TRUE) +
    geom_point(data = ps_sub,
               aes(x = mean, y = yi, fill = pool),
               shape = 21, color = "black", size = 1.8, stroke = 0.3,
               na.rm = TRUE) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(breaks = seq_along(ORDERED_ORDERS) - 1,
                        labels = ORDERED_ORDERS,
                        trans  = "reverse",
                        expand = expansion(add = 0.5)) +
    coord_cartesian(xlim = xrange, clip = "on") +
    labs(x = expression(Delta^14*"C (‰)"), y = NULL,
         tag = panel_tag, title = panel_title) +
    theme_nature() + theme(legend.position = "none")
}

# -----------------------------------------------------------------------------
# 9. Difference panels (b, d) — POM − MAOM per soil order
# -----------------------------------------------------------------------------
# `label_nudge` is an optional named vector of per-soil-order vertical offsets
# for the "k=, n=" annotation, in y data units.  The y axis is reversed, so a
# negative value raises that row's label.  It is used where a wide confidence
# interval would otherwise run underneath the text.  Soil orders not named in
# the vector are left in place.
plot_diff <- function(diff_summary, wide_pair_df, pom_col,
                       panel_tag, panel_title,
                       label_nudge = NULL) {
  raw <- wide_pair_df |>
    filter(!is.na(.data[[pom_col]]), !is.na(MAOM),
           order %in% ORDERED_ORDERS) |>
    mutate(value = .data[[pom_col]] - MAOM,
           yi    = match(order, ORDERED_ORDERS) - 1 + runif(n(), -0.18, 0.18))

  ds <- diff_summary |> filter(order %in% ORDERED_ORDERS) |>
    mutate(yi = match(order, ORDERED_ORDERS) - 1,
           sig_zero = !is.na(p) & p < 0.05,
           fill_col = ifelse(sig_zero, C_DIFF, "white"),
           label    = sprintf("k=%d, n=%d %s", k, n, sig),
           y_lab    = yi + if (is.null(label_nudge)) 0 else
                      coalesce(unname(label_nudge[as.character(order)]), 0))

  # X-axis range spans every CI bound and every raw paired difference, so
  # nothing plotted falls outside the panel.
  finite <- ds |> filter(!is.na(lo), !is.na(hi))
  rng  <- data_range(finite$hi, finite$lo, raw$value)
  xmin <- rng[1]; xmax <- rng[2]
  xpad <- max(60, 0.05 * (xmax - xmin))
  x_annot <- xmax + xpad

  ggplot() +
    geom_rect(data = stripe_df,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = pos - 0.5, ymax = pos + 0.5),
              inherit.aes = FALSE, fill = "#F4F4F2", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = C_ZERO,
               linewidth = 0.4, alpha = 0.95) +
    geom_point(data = raw, aes(x = value, y = yi),
                color = RAW_POINT_COL, size = 0.5, alpha = 0.4) +
    geom_linerange(data = ds,
                    aes(xmin = lo, xmax = hi, y = yi),
                    color = C_DIFF, linewidth = 0.5, na.rm = TRUE) +
    geom_point(data = ds,
               aes(x = mean, y = yi, fill = I(fill_col)),
               shape = 21, color = "black", size = 1.8, stroke = 0.35,
               na.rm = TRUE) +
    geom_text(data = ds, aes(x = x_annot, y = y_lab, label = label),
              hjust = 0, size = 1.95, color = "grey25", na.rm = TRUE) +
    scale_y_continuous(breaks = seq_along(ORDERED_ORDERS) - 1,
                        labels = ORDERED_ORDERS,
                        trans  = "reverse",
                        expand = expansion(add = 0.5)) +
    coord_cartesian(xlim = c(xmin - xpad, x_annot + 230), clip = "off") +
    labs(x = expression("POM − MAOM  "*Delta^14*"C (‰)"), y = NULL,
         tag = panel_tag, title = panel_title) +
    theme_nature() + theme(legend.position = "none")
}

# -----------------------------------------------------------------------------
# 10. Legend strip — pool/marker fill convention shared across panels
# -----------------------------------------------------------------------------
# A small helper that renders a stand-alone legend strip; the underlying
# data points are placed off-canvas so only the legend appears.
make_legend_panel <- function(items, colors, key_name, fill = TRUE) {
  df <- tibble(what = factor(items, levels = items),
               x = NA_real_, y = NA_real_)
  p <- if (fill) {
    ggplot(df, aes(x, y, fill = what)) +
      geom_point(shape = 21, color = "black", size = 2.4, stroke = 0.3,
                  na.rm = TRUE) +
      scale_fill_manual(values = colors, name = key_name)
  } else {
    ggplot(df, aes(x, y, color = what)) +
      geom_point(size = 2.4, na.rm = TRUE) +
      scale_color_manual(values = colors, name = key_name)
  }
  p +
    guides(fill  = if (fill) guide_legend(nrow = 1,
              override.aes = list(size = 2.4)) else "none",
            color = if (!fill) guide_legend(nrow = 1,
              override.aes = list(size = 2.4)) else "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void(base_size = 7, base_family = "Helvetica") +
    theme(legend.position    = "bottom",
          legend.title       = element_text(size = 6.5),
          legend.text        = element_text(size = 6.5),
          legend.key.size    = unit(0.3, "cm"),
          legend.box.spacing = unit(0, "pt"),
          plot.margin        = margin(0, 0, 0, 0))
}

# -----------------------------------------------------------------------------
# 11. Assemble the 4-panel composite
# -----------------------------------------------------------------------------
# Rows: (a, b) = free POM vs MAOM; (c, d) = occluded POM vs MAOM.
# Columns: left = soil-order profile of both fractions; right = paired
# difference per soil order.
pa <- plot_profile(c("fPOM", "MAOM"), paired_fp, pool_fp,
                    "a", "Free POM vs. MAOM")
pb <- plot_diff(   free_d, wide, "fPOM",
                    "b", "Free POM − MAOM",
                    label_nudge = c(Gelisols = -0.30))
pc <- plot_profile(c("oPOM", "MAOM"), paired_op, pool_op,
                    "c", "Occluded POM vs. MAOM")
pd <- plot_diff(   occ_d,  wide, "oPOM",
                    "d", "Occluded POM − MAOM")

leg_pool <- make_legend_panel(
  items  = c("Free POM", "Occluded POM", "MAOM",
              "POM − MAOM (paired)", "ns (CI crosses 0)"),
  colors = c("Free POM" = C_FPOM, "Occluded POM" = C_OPOM,
              "MAOM" = C_MAOM, "POM − MAOM (paired)" = C_DIFF,
              "ns (CI crosses 0)" = "white"),
  key_name = NULL, fill = TRUE)

fig <- (pa | pb) / (pc | pd) +
  plot_layout(widths = c(1, 1.05)) &
  theme(plot.margin = margin(2, 4, 2, 2))

final <- wrap_elements(full = fig) /
         wrap_elements(full = leg_pool) +
         plot_layout(heights = c(28, 1))

ggsave("Figure2_SoilOrder.pdf", final,
       width = 6.5, height = 6, device = cairo_pdf)
ggsave("Figure2_SoilOrder.png", final,
       width = 6.5, height = 6, dpi = 600, bg = "white")
message("Wrote Figure2_SoilOrder.pdf, .png and _table.csv")
