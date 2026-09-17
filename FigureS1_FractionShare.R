# =============================================================================
# Figure S1 — Share of soil organic carbon held by each density fraction,
#             top 30 cm.
#
## Authors:   Shangshi Liu, Jonathan Sanderman, Mark A. Bradford
# Citation:  [add manuscript citation when accepted]

#
# This script reproduces Supplementary Figure 1 of the manuscript. It loads
# the ISRaD flat-fraction product, retains soil samples for which all three
# density fractions were measured on the same physical aliquot, and expresses
# the carbon content of each fraction as its share of total layer organic
# carbon within three 10-cm depth bins.  Unlike Figures 1 and 2, the summaries
# here are descriptive (median, interquartile range and arithmetic mean); no
# mixed-effects model is fitted, because the quantity of interest is the
# composition of the measured samples rather than a global mean.
#
# Sample-inclusion criteria are set in the `Inclusion criteria` block below
# (section 0), matching Figure1_DepthProfile.R and Figure2_SoilOrder.R.
#
# Output files:
#   FigureS1_FractionShare.pdf / .png      single-panel figure
#   FigureS1_FractionShare_table.csv       numerical summary of every cell
#
# Dependencies:
#   R >= 4.3
#   tidyverse
#
# Input data:
#   ISRaD v2.9.9 flat-fraction product, downloaded from
#   https://www.soilradiocarbon.org/database-1
#   and placed in `./ISRaD_database_files/` with the filename
#   `ISRaD_extra_flat_fraction_v 2.9.9.2025-08-14.csv`.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# 0. Inclusion criteria  (identical to Figure1_DepthProfile.R)
# -----------------------------------------------------------------------------
# Samples are restricted to those collected in or after YEAR_MIN and to layers
# no thicker than MAX_THICK_CM, so that the composition reported here describes
# the same set of samples as the radiocarbon analyses in Figures 1 and 2.
# MAX_DEPTH_CM sets the depth range shown; carbon-share data become sparse
# below 30 cm.
YEAR_MIN     <- 2000    # earliest sampling year retained
MAX_THICK_CM <- 30      # maximum layer thickness retained (cm)
MAX_DEPTH_CM <- 30      # deepest layer midpoint shown (cm)

# Mass-balance window.  The three shares are reported independently by the
# original studies and need not sum to exactly 100 %; samples whose three
# fractions sum outside this window are treated as incomplete recoveries and
# excluded, and the remainder are renormalised to 100 %.
MASS_BAL <- c(90, 110)

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
# "sand" rows are returned as NA so they are excluded.
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
# `frc_c_perc` is the fraction's proportion of total layer carbon, so the
# three pools of a complete sample sum to 100 %.  We retain density-
# fractionation measurements with a reported share and known depth bounds, and
# exclude organic horizons, above-ground litter, samples collected before
# YEAR_MIN, and layers thicker than MAX_THICK_CM.
#

obs_year_col <- if ("lyr_obs_date_y.x" %in% names(raw)) "lyr_obs_date_y.x" else
                                                        "lyr_obs_date_y"
raw <- raw |> mutate(obs_year = .data[[obs_year_col]])

d <- raw |>
  filter(frc_scheme == "density",
         is.na(lyr_all_org_neg) | lyr_all_org_neg != "yes",
         !is.na(frc_c_perc),
         !is.na(lyr_top), !is.na(lyr_bot)) |>
  mutate(lyr_mid    = (lyr_top + lyr_bot) / 2,
         thickness  = lyr_bot - lyr_top,
         pool       = classify_pool(frc_property),
         entry_name = as.character(entry_name)) |>
  filter(lyr_top   >= 0,
         obs_year  >= YEAR_MIN,
         thickness <= MAX_THICK_CM,
         lyr_mid   <  MAX_DEPTH_CM,
         !is.na(pool))

# -----------------------------------------------------------------------------
# 4. Depth bins
# -----------------------------------------------------------------------------
# Three 10-cm bins spanning the top 30 cm.  There is no pooled row: the
# composition of the profile changes systematically with depth, so a single
# depth-averaged figure would obscure the pattern the panel is meant to show.
bin_edges <- seq(0, MAX_DEPTH_CM, by = 10)
bin_labs  <- c("0-10", "10-20", "20-30")   # top -> bottom

# -----------------------------------------------------------------------------
# 5. Complete-sample carbon shares
# -----------------------------------------------------------------------------
# Fractions are matched on the hierarchical key
#     entry / site / profile / layer / frc_input,
# where `frc_input` identifies the physical soil aliquot from which each
# density fraction was separated.  Where a study reports a pool as several
# density sub-cuts (e.g. heavy at > 1.6 and > 2.0 g cm-3), the shares for that
# pool at that aliquot are SUMMED, so that multi-cut schemes are recovered as
# a single value per pool.  Only aliquots with all three pools present are
# retained, the mass balance is checked against MASS_BAL, and the three shares
# are renormalised to 100 %.
pair_key <- c("entry_name", "site_name", "pro_name", "lyr_name", "frc_input")

pool_sum <- d |>
  group_by(across(all_of(c(pair_key, "pool")))) |>
  summarise(share   = sum(frc_c_perc, na.rm = TRUE),
            lyr_mid = mean(lyr_mid,   na.rm = TRUE),
            .groups = "drop")

wide <- pool_sum |>
  pivot_wider(id_cols = all_of(pair_key),
              names_from = pool, values_from = share) |>
  filter(!is.na(fPOM), !is.na(oPOM), !is.na(MAOM)) |>
  mutate(total = fPOM + oPOM + MAOM) |>
  filter(total >= MASS_BAL[1], total <= MASS_BAL[2]) |>
  mutate(across(c(fPOM, oPOM, MAOM), ~ 100 * .x / total)) |>
  left_join(pool_sum |> group_by(across(all_of(pair_key))) |>
              summarise(lyr_mid = mean(lyr_mid, na.rm = TRUE), .groups = "drop"),
            by = pair_key)

long <- wide |>
  pivot_longer(c(fPOM, oPOM, MAOM), names_to = "pool", values_to = "share") |>
  mutate(bin = cut(lyr_mid, breaks = bin_edges, labels = bin_labs,
                    right = FALSE, include.lowest = TRUE)) |>
  filter(!is.na(bin)) |>
  mutate(bin  = factor(bin,  levels = bin_labs),
         pool = factor(pool, levels = c("fPOM", "oPOM", "MAOM")))

# -----------------------------------------------------------------------------
# 6. Per-bin summary
# -----------------------------------------------------------------------------
# Descriptive statistics only: the arithmetic mean is overplotted on the box,
# and k (number of contributing studies) and n (number of complete samples)
# are annotated to the right of each depth row.
summ <- long |>
  group_by(bin, pool) |>
  summarise(mean   = mean(share),
            median = median(share),
            q25    = quantile(share, 0.25),
            q75    = quantile(share, 0.75),
            n      = dplyr::n(),
            k      = dplyr::n_distinct(entry_name),
            .groups = "drop")

# Numerical summary written to disk so every element of the figure is
# traceable to a CSV row.
write_csv(summ, "FigureS1_FractionShare_table.csv")

# -----------------------------------------------------------------------------
# 7. Plot styling
# -----------------------------------------------------------------------------
# Colour palette matches Figures 1 and 2: red = free POM, orange = occluded
# POM, blue = MAOM.  Occluded POM additionally carries a triangular marker so
# the three fractions remain distinguishable in greyscale.
C_FPOM <- "#D62728"
C_OPOM <- "#D55E00"
C_MAOM <- "#0072B2"

pal <- c(fPOM = C_FPOM, oPOM = C_OPOM, MAOM = C_MAOM)
shp <- c(fPOM = 21,     oPOM = 24,     MAOM = 21)

theme_nature <- function() {
  theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(axis.line   = element_line(linewidth = 0.4),
          axis.ticks  = element_line(linewidth = 0.35),
          axis.title  = element_text(size = 7),
          axis.text   = element_text(size = 6.5),
          plot.title  = element_text(size = 7.5, color = "grey25",
                                       face = "plain"),
          panel.grid.major.x = element_line(color = "grey94",
                                              linewidth = 0.25),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 6.5),
          legend.key.size = unit(0.3, "cm"))
}

# -----------------------------------------------------------------------------
# 8. Vertical positions
# -----------------------------------------------------------------------------
# Depth rows are drawn on a numeric y axis so the three fractions can be
# offset within each row.  Position 0 is the deepest bin, so that depth
# increases downwards once the labels are reversed.  DODGE separates the
# fractions vertically; without it the three boxes would overlap.
ypos_of <- function(b) length(bin_labs) - match(as.character(b), bin_labs)
DODGE   <- c(fPOM = +0.26, oPOM = 0.0, MAOM = -0.26)

long <- long |> mutate(y0 = ypos_of(bin), y = y0 + DODGE[as.character(pool)])
summ <- summ |> mutate(y0 = ypos_of(bin), y = y0 + DODGE[as.character(pool)])

# Alternating row-stripe background so the eye can follow each depth row
# across the full width of the panel.
stripe_df <- tibble(y0 = seq_along(bin_labs) - 1) |> filter(y0 %% 2 == 1)

# Counts are annotated just beyond the right-hand end of the 0-100 % axis.
# The x limit is therefore set by coord_cartesian() rather than by
# scale_x_continuous(), because a scale limit would drop the annotation.
n_lab   <- summ |> filter(pool == "MAOM") |>
  transmute(y0, label = sprintf("k=%d, n=%d", k, n))
x_annot <- 101

# -----------------------------------------------------------------------------
# 9. Assemble the figure
# -----------------------------------------------------------------------------
# Horizontal boxplots (orientation = "y") show the median, interquartile range
# and 1.5 x IQR whiskers of each fraction within each depth bin.  Outlier
# points are suppressed because every observation is already drawn as a
# jittered point over the box.  The filled marker is the arithmetic mean.
p <- ggplot() +
  geom_rect(data = stripe_df,
            aes(xmin = -Inf, xmax = Inf,
                ymin = y0 - 0.5, ymax = y0 + 0.5),
            inherit.aes = FALSE, fill = "#F4F4F2", alpha = 0.5) +
  geom_boxplot(data = long,
               aes(x = share, y = y, group = interaction(bin, pool),
                   color = pool, fill = pool),
               width = 0.22, linewidth = 0.4, outlier.shape = NA,
               alpha = 0.22, orientation = "y") +
  geom_jitter(data = long,
              aes(x = share, y = y, color = pool),
              shape = 16, width = 0, height = 0.09,
              size = 0.5, alpha = 0.40) +
  geom_point(data = summ,
             aes(x = mean, y = y, fill = pool, shape = pool),
             color = "black", size = 1.9, stroke = 0.3) +
  geom_text(data = n_lab,
            aes(x = x_annot, y = y0, label = label),
            hjust = 0, size = 1.8, color = "grey45") +
  scale_color_manual(values = pal, guide = "none") +
  scale_shape_manual(values = shp,
                      breaks = c("fPOM", "oPOM", "MAOM"),
                      labels = c("Free POM", "Occluded POM", "MAOM")) +
  scale_fill_manual(values = pal,
                     breaks = c("fPOM", "oPOM", "MAOM"),
                     labels = c("Free POM", "Occluded POM", "MAOM")) +
  scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(0)) +
  scale_y_continuous(breaks = seq_along(bin_labs) - 1,
                      labels = rev(bin_labs),
                      expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  labs(x = "Share of total soil organic carbon (%)", y = "Depth (cm)") +
  theme_nature() +
  theme(legend.position = "bottom",
        plot.margin     = margin(t = 4, r = 48, b = 4, l = 4)) +
  guides(fill  = "none",
         shape = guide_legend(nrow = 1,
                   override.aes = list(size = 2.2,
                                        fill  = c(C_FPOM, C_OPOM, C_MAOM),
                                        color = "black")))

ggsave("FigureS1_FractionShare.pdf", p,
       width = 5.8, height = 3.4, device = cairo_pdf)
ggsave("FigureS1_FractionShare.png", p,
       width = 5.8, height = 3.4, dpi = 600, bg = "white")
message("Wrote FigureS1_FractionShare.pdf, .png and _table.csv")

cat("\nMean share of total soil organic carbon (%), top 30 cm:\n")
print(summ |>
        mutate(across(c(mean, median), ~ round(.x, 1))) |>
        select(bin, pool, mean, median, k, n) |>
        arrange(bin, pool))
