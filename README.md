# Soilfractions14C

R code and Source Data for:

> Liu, S., Sanderman, J. & Bradford, M. A. Comparable soil carbon persistence through physical and mineral protection. 


GitHub: https://github.com/Shangshi-ecology/Soilfractions14C

Contact: Shangshi Liu (liushangshi@gmail.com)

## Overview

The scripts compare the radiocarbon (Δ¹⁴C) content of density-separated soil organic matter fractions, using pairs of fractions separated from the same soil aliquot in the International Soil Radiocarbon Database (ISRaD). The three fractions are free particulate organic matter (free POM, ISRaD "free light"), occluded POM ("occluded light") and mineral-associated organic matter (MAOM, "heavy"). Every figure in the paper, and every number shown on it, is reproduced from the public ISRaD file.

## Input data

ISRaD flat-fraction product, version 2.9.9 (released 2025-08-14), available from https://www.soilradiocarbon.org/database-1.

Save the file with its original name in a folder called `ISRaD_database_files/`, placed next to the code folder:

```
project_folder/
├── ISRaD_database_files/
│   └── ISRaD_extra_flat_fraction_v 2.9.9.2025-08-14.csv
└── Soilfractions14C/            (this folder)
```

## Scripts

| Script | Reproduces | Output files |
|---|---|---|
| `Figure1_DepthProfile.R` | Fig. 1: Δ¹⁴C of paired fractions by depth | `Figure1_DepthProfile.pdf/.png`, `_table.csv` |
| `Figure2_SoilOrder.R` | Fig. 2: Δ¹⁴C of paired fractions by USDA soil order | `Figure2_SoilOrder.pdf/.png`, `_table.csv` |
| `FigureS1_FractionShare.R` | Extended Data Fig. 1: share of soil organic carbon in each fraction, top 30 cm | `FigureS1_FractionShare.pdf/.png`, `_table.csv` |
| `FigureS2S3_Stoichiometry.R` | Extended Data Figs. 2 and 3: δ¹³C, δ¹⁵N and C:N of free and occluded POM vs MAOM | `FigureS2_Stoichiometry_FreePOM.*`, `FigureS3_Stoichiometry_OccludedPOM.*` |
| `FigureS4_SoilOrder_Top30.R` | Extended Data Fig. 4: as Fig. 2, top 30 cm only | `FigureS4_SoilOrder_Top30.pdf/.png`, `_table.csv` |

The Extended Data scripts keep the file names `FigureS1`–`FigureS4` from an earlier version of the manuscript, in which these figures were Supplementary Figures.

Each `_table.csv` file lists every plotted summary value (mean, 95% confidence interval, P value, number of studies k, number of paired observations n, and the model used), so each element of a figure can be traced to a row.

## How to run

Requirements: R ≥ 4.3 with the packages `tidyverse`, `lme4`, `lmerTest` and `patchwork`.

Run each script with the project folder (the folder that contains `ISRaD_database_files/`) as the working directory, for example:

```
cd project_folder
Rscript Soilfractions14C/Figure1_DepthProfile.R
```

In RStudio, set the working directory to the project folder (Session > Set Working Directory) and click **Source**. Outputs are written to the project folder. Each script runs independently and takes less than 20 seconds.

## Analysis

**Inclusion criteria** are set at the top of each script.
- All analyses:
  - density fractionation only;
  - layers with known depth and no thicker than 30 cm;
  - above-ground litter layers (upper depth < 0 cm) and profiles whose mineral interface is unknown are excluded.
- Figs. 1 and 2 and Extended Data Figs. 1 and 4: samples collected in or after 2000, with a reported value for the analysed variable.
- δ¹³C, δ¹⁵N and C:N (Extended Data Figs. 2 and 3): no sampling-year restriction, because these variables carry no bomb-¹⁴C signal. C:N values above 100 (MAOM) or 200 (POM) are treated as reporting errors.
- Extended Data Fig. 1:
  - layer midpoints shallower than 30 cm;
  - only samples in which all three fractions were measured on the same aliquot, and whose fraction carbon shares sum to 90–110% (shares are then rescaled to 100%).

**Pairing.** Fractions are matched on study / site / profile / layer / aliquot (`entry_name / site_name / pro_name / lyr_name / frc_input`), so that only fractions separated from the same physical aliquot are compared. Replicates or density sub-cuts of the same fraction within an aliquot are averaged before pairing.

**Statistics.** For each depth bin (10 cm; layers with midpoints ≥ 100 cm are pooled into one bin), for all depths combined, or soil order, the paired difference (POM − MAOM) is modelled as `y ~ 1 + (1 | study/site/profile)` with `lme4`. Singular fits are accepted. Simpler random-effect structures, and finally ordinary least squares, are used only when the full model cannot be fitted.

The model used for each stratum is recorded in the `model` column of the tables and the Source Data:

| Code | Model |
|---|---|
| `lmm3` | study/site/profile random intercepts |
| `lmm2` | study/site random intercepts |
| `lmm1` | study random intercepts |
| `ols` | ordinary least squares |
| `k=1 (no inference)`, `n<2 (none)` | Mean only, no test |

Means are shown with 95% confidence intervals (mean ± 1.96 SE). P values come from two-sided Wald z-tests of the intercept against zero, without adjustment for multiple comparisons. Strata with a single study (k = 1) are shown without a confidence interval or P value. Extended Data Fig. 1 is descriptive only.

## Source Data

The Source Data workbooks, one per figure, are provided with the paper and in the `Source_Data/` folder of this archive:

| Workbook | Figure |
|---|---|
| `Source_Data_Fig1.xlsx` | Fig. 1 |
| `Source_Data_Fig2.xlsx` | Fig. 2 |
| `Source_Data_Extended_Data_Fig1.xlsx` | Extended Data Fig. 1 |
| `Source_Data_Extended_Data_Fig2.xlsx` | Extended Data Fig. 2 |
| `Source_Data_Extended_Data_Fig3.xlsx` | Extended Data Fig. 3 |
| `Source_Data_Extended_Data_Fig4.xlsx` | Extended Data Fig. 4 |

Each workbook has a `README` sheet that defines its columns, plus two kinds of data sheet:

- **`obs` sheets:** one row per paired sample (the points in the figure). Each row gives the identifiers (study, site, profile, layer, aliquot), the layer midpoint, the depth bin or soil order, the POM and MAOM values, and their difference.
- **`summary` sheets:** one row per depth bin or soil order (the circles and error bars). Each row gives the model mean, 95% confidence interval, k, n and model. For the difference panels it also gives the standard error, the z statistic and the P value.

The summary values are identical to those in the `_table.csv` file written by the corresponding script. SE and z are recovered from the confidence interval (CI = mean ± 1.96 SE).

## Licence

ISRaD data are subject to the ISRaD terms of use; please cite ISRaD (Lawrence et al. 2020, *Earth Syst. Sci. Data* 12, 61–76; von Fromm et al. 2026, *Earth Syst. Sci. Data* 18, 5793–5815) when using them.
