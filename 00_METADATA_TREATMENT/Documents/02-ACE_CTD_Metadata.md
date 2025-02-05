02-ACE CTD Metadata
================
Lisa-Marie Delpech
1/26/2022

- [Load packages and data](#load-packages-and-data)
- [Format data](#format-data)
- [Grouping variables](#grouping-variables)
  - [Qualitative depth](#qualitative-depth)
  - [Format water mass](#format-water-mass)
- [Missing values](#missing-values)
  - [Explore missing values](#explore-missing-values)
  - [Imputation](#imputation)
  - [Inspect imputed data sets](#inspect-imputed-data-sets)
  - [Inspect effect of imputation in the TS
    space](#inspect-effect-of-imputation-in-the-ts-space)
- [Subset to metagenome samples](#subset-to-metagenome-samples)
- [PCA](#pca)
  - [All variables and all samples](#all-variables-and-all-samples)
  - [ANTA and APLR](#anta-and-aplr)
  - [Surface samples](#surface-samples)
  - [With physical variables](#with-physical-variables)
- [Match Latent Variables from Landwehr *et al.* (2021) with CTD
  samples](#match-latent-variables-from-landwehr-et-al-2021-with-ctd-samples)
- [Hierarchical clustering](#hierarchical-clustering)
- [Merged CTD and UDW metadata files for metagenomic
  samples](#merged-ctd-and-udw-metadata-files-for-metagenomic-samples)
- [CTD and UDW metadata table for metagenomic
  samples](#ctd-and-udw-metadata-table-for-metagenomic-samples)

# Load packages and data

``` r
set.seed(3)
```

``` r
require(tidyverse)
require(vegan)
require(viridis)
require(ggsci)
require(ggrepel)
require(rcompanion)
require(caret)
require(VIM)
require(naniar)
require(mice)
require(missForest)
require(FactoMineR)
require(Amelia)
require(RANN)
require(cowplot)
require(lubridate)
require(ggdendro)
require(factoextra)
require(ggOceanMaps)
require(ggOceanMapsData)
```

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'ggOceanMapsData'

``` r
source("../HOME_FUNCTIONS/ggPCA.R")
source("../HOME_FUNCTIONS/ggRDA.R")
```

``` r
meta <- read.table(file = "./Metadata/meta_CTD_ess_V4_4.csv", header=TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE)
str(meta)
```

    ## 'data.frame':    345 obs. of  93 variables:
    ##  $ Event_number                         : int  75 75 75 75 75 75 123 123 123 123 ...
    ##  $ CTD_cast_number                      : int  1 1 1 1 1 1 3 3 3 3 ...
    ##  $ ACE_station_number                   : int  2 2 2 2 2 2 8 8 8 8 ...
    ##  $ TM_station_number                    : int  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ ACE_seq_name                         : chr  "ACE1_75_11A_15m" "ACE1_75_11B_15m" "ACE1_75_13A_150m" "ACE1_75_13B_150m" ...
    ##  $ Size_fraction                        : chr  ">3 µm" "0.2-3 µm" ">3 µm" "0.2-3 µm" ...
    ##  $ Depth_m                              : int  15 15 150 150 60 60 15 15 15 150 ...
    ##  $ Depth_q                              : chr  "ZZZ" "ZZZ" "ZZZ" "ZZZ" ...
    ##  $ DNA_vol_µl                           : int  98 98 98 98 98 98 99 NA 98 NA ...
    ##  $ DNA_Conc_ng.µl                       : num  8.07 11.7 2.17 3.57 12.25 ...
    ##  $ DNA_qty_ng                           : num  791 1147 213 350 1200 ...
    ##  $ RNA_vol_µl                           : num  58 57 58 56 58 55 58 NA 55 NA ...
    ##  $ RNA_Conc_ng.µl                       : num  31.5 44.9 8 5.26 201 ...
    ##  $ RNA_qty_ng                           : num  1827 2559 464 295 11658 ...
    ##  $ Event_start_date                     : chr  "2016-12-28T06:55:00Z+00" "2016-12-28T06:55:00Z+00" "2016-12-28T06:55:00Z+00" "2016-12-28T06:55:00Z+00" ...
    ##  $ Event_end_date                       : chr  "2016-12-28T08:00:00Z+00" "2016-12-28T08:00:00Z+00" "2016-12-28T08:00:00Z+00" "2016-12-28T08:00:00Z+00" ...
    ##  $ Latitude                             : num  -46.9 -46.9 -46.9 -46.9 -46.9 ...
    ##  $ Longitude                            : num  38.1 38.1 38.1 38.1 38.1 ...
    ##  $ Temperature_ITS90                    : num  NA NA NA NA NA ...
    ##  $ Salinity_PSS78                       : num  NA NA NA NA NA ...
    ##  $ Density_kg.m2                        : num  NA NA NA NA NA ...
    ##  $ dO2_µmol.kg                          : num  NA NA NA NA NA ...
    ##  $ Fluorescence_mg.m3                   : num  NA NA NA NA NA NA NA NA NA 0.043 ...
    ##  $ PAR_µmol.m2.s                        : num  NA NA NA NA NA NA NA NA NA 0.0195 ...
    ##  $ NOx_µmol.L                           : num  19.7 19.7 24 24 20.4 ...
    ##  $ Silicic_acid_µmol.L                  : num  2.95 2.95 12.64 12.64 3.92 ...
    ##  $ Nitrite_µmol.L                       : num  0.19 0.19 0.17 0.17 0.18 0.18 0.24 0.24 0.24 0.19 ...
    ##  $ Phosphate_µmol.L                     : num  1.34 1.34 1.45 1.45 1.33 1.33 1.42 1.42 1.42 1.46 ...
    ##  $ Ammonium_µmol.L                      : num  0.28 0.28 0.87 0.87 0.44 0.44 0.61 0.61 0.61 0.33 ...
    ##  $ Nitrate_µmol.L                       : num  19.5 19.5 23.8 23.8 20.3 ...
    ##  $ Bsi_µM                               : num  1.09 1.09 1.17 1.17 3.19 ...
    ##  $ dCd_nmol.kg                          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ dFe_nmol.kg                          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ dZn_nmol.kg                          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ dCu_nmol.kg                          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ dNi_nmol.kg                          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Fe_Ligands.nM                        : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ log_Kcond_FeL                        : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ TPZT_µM.C                            : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Humics_µg.SRFAeq.L                   : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ POC_µM                               : num  NA NA 1.46 1.46 7.7 7.7 3.67 3.67 3.67 1.2 ...
    ##  $ d13C                                 : num  NA NA -25 -25 -24.4 ...
    ##  $ PON.µmol                             : num  NA NA 0.32 0.32 1.37 1.37 0.65 0.65 0.65 0.28 ...
    ##  $ d15N                                 : num  NA NA 3.37 3.37 2.39 2.39 0.95 0.95 0.95 6.19 ...
    ##  $ Bacteria_HDNA.cell.mL                : int  300000 300000 255000 255000 404000 404000 118000 118000 118000 535000 ...
    ##  $ Bacteria_LDNA.cell.mL                : int  249000 249000 161000 161000 268000 268000 111000 111000 111000 433000 ...
    ##  $ Tot_bacteria_cell.mL                 : int  549000 549000 416000 416000 671000 671000 229000 229000 229000 968000 ...
    ##  $ Synechos_cell.mL                     : int  222 222 NA NA 54 54 NA NA NA NA ...
    ##  $ Picoeuk1_cell.mL                     : int  46 46 NA NA 409 409 NA NA NA NA ...
    ##  $ Picoeuk2_cell.mL                     : int  73 73 NA NA 744 744 NA NA NA NA ...
    ##  $ Nanoeuk_cell.mL                      : int  138 138 NA NA 546 546 NA NA NA NA ...
    ##  $ Cryptomonas_cell.mL                  : int  0 0 NA NA 0 0 0 0 0 NA ...
    ##  $ TEP_µg.XG.eq.L                       : num  10.8 10.8 15.6 15.6 28.5 ...
    ##  $ CSP_µg.BSA.eq.L                      : num  NA NA 7.07 7.07 7.68 7.68 6.42 6.42 6.42 NA ...
    ##  $ Chlorophyll_chemtax                  : num  0.0116 0.0116 0.0063 0.0063 0.0136 0.0136 0.056 0.056 0.056 0 ...
    ##  $ Crypto1_chemtax                      : num  0.0006 0.0006 0.0006 0.0006 0.0037 0.0037 0.0256 0.0256 0.0256 0.0039 ...
    ##  $ Cyano2_chemtax                       : num  0.0041 0.0041 0 0 0.0001 0.0001 0.0183 0.0183 0.0183 0.0061 ...
    ##  $ DiatA_chemtax                        : num  0.0756 0.0756 0.0506 0.0506 0.3182 ...
    ##  $ DiatB_chemtax                        : num  0.0178 0.0178 0.0342 0.0342 0.12 0.12 0.0296 0.0296 0.0296 0.0048 ...
    ##  $ DinoA_chemtax                        : num  0.0075 0.0075 0 0 0.0158 0.0158 0.0001 0.0001 0.0001 0 ...
    ##  $ Hapto8_chemtax                       : num  0.0565 0.0565 0.0176 0.0176 0.1054 ...
    ##  $ Hapto67_chemtax                      : num  0.1001 0.1001 0.0174 0.0174 0.1184 ...
    ##  $ Pras3_chemtax                        : num  0.0082 0.0082 0.0006 0.0006 0.0104 ...
    ##  $ Pelago_chemtax                       : num  0.011 0.011 0.001 0.001 0.0204 0.0204 0.0023 0.0023 0.0023 0.0026 ...
    ##  $ Fv_Fm                                : num  0.158 0.158 0.342 0.342 0.355 0.355 0.433 0.433 0.433 0.283 ...
    ##  $ Uabs_Glorys_mean                     : num  0.281 0.281 0.281 0.281 0.281 ...
    ##  $ Ekin_Glorys_mean                     : num  0.0398 0.0398 0.0398 0.0398 0.0398 ...
    ##  $ Vorticity_Glorys_mean                : num  -0.205 -0.205 -0.205 -0.205 -0.205 ...
    ##  $ OW_Glorys_mean                       : num  0.921 0.921 0.921 0.921 0.921 ...
    ##  $ EulerDiverg_Glorys_mean              : num  0.0922 0.0922 0.0922 0.0922 0.0922 ...
    ##  $ LagrDiverg_Glorys_05daysbackward_mean: num  0.0422 0.0422 0.0422 0.0422 0.0422 ...
    ##  $ LagrDiverg_Glorys_10daysbackward_mean: num  -0.00717 -0.00717 -0.00717 -0.00717 -0.00717 ...
    ##  $ LagrDiverg_Glorys_15daysbackward_mean: num  -0.00851 -0.00851 -0.00851 -0.00851 -0.00851 ...
    ##  $ LagrDiverg_Glorys_20daysbackward_mean: num  -0.00971 -0.00971 -0.00971 -0.00971 -0.00971 ...
    ##  $ LagrDiverg_Glorys_30daysbackward_mean: num  -0.00823 -0.00823 -0.00823 -0.00823 -0.00823 ...
    ##  $ LagrDiverg_Glorys_60daysbackward_mean: num  -0.00543 -0.00543 -0.00543 -0.00543 -0.00543 ...
    ##  $ RetentionTime_Glorys_mean            : num  1.09 1.09 1.09 1.09 1.09 ...
    ##  $ Ftle_Glorys_05daysBackward_mean      : num  0.0858 0.0858 0.0858 0.0858 0.0858 ...
    ##  $ Ftle_Glorys_10daysBackward_mean      : num  0.168 0.168 0.168 0.168 0.168 ...
    ##  $ Ftle_Glorys_15daysBackward_mean      : num  0.148 0.148 0.148 0.148 0.148 ...
    ##  $ Ftle_Glorys_20daysBackward_mean      : num  0.138 0.138 0.138 0.138 0.138 ...
    ##  $ Ftle_Glorys_30daysBackward_mean      : num  0.0821 0.0821 0.0821 0.0821 0.0821 ...
    ##  $ Ftle_Glorys_60daysBackward_mean      : num  0.0574 0.0574 0.0574 0.0574 0.0574 ...
    ##  $ Betw_Glorys_05daysCentered_mean      : num  4.86 4.86 4.86 4.86 4.86 ...
    ##  $ Betw_Glorys_10daysCentered_mean      : num  6.03 6.03 6.03 6.03 6.03 ...
    ##  $ Betw_Glorys_15daysCentered_mean      : num  12.2 12.2 12.2 12.2 12.2 ...
    ##  $ Betw_Glorys_20daysCentered_mean      : num  19 19 19 19 19 ...
    ##  $ Betw_Glorys_30daysCentered_mean      : num  41.1 41.1 41.1 41.1 41.1 ...
    ##  $ Betw_Glorys_60daysCentered_mean      : num  230 230 230 230 230 ...
    ##  $ Longhurst_Prov                       : chr  "SANT" "SANT" "SANT" "SANT" ...
    ##  $ MertzGlacier                         : chr  "FAUX" "FAUX" "FAUX" "FAUX" ...
    ##  $ Water_mass_detailed                  : chr  "SASW" "SASW" "SASW" "SASW" ...
    ##  $ Water_mass                           : chr  "SASW" "SASW" "SASW" "SASW" ...

``` r
meta.udw <- read.table(file = "./Metadata/meta_UDW_ess_V4_4.csv", header=TRUE, sep = ";", dec = ".", stringsAsFactors = TRUE)
str(meta.udw) 
```

    ## 'data.frame':    303 obs. of  87 variables:
    ##  $ ACE_seq_name                         : Factor w/ 301 levels "ACE1_11_3A_4m",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Event_number                         : int  11 11 115 115 13 13 150 150 154 154 ...
    ##  $ Size_fraction                        : Factor w/ 5 levels "","<0.2 µm",">3 µm",..: 3 4 3 4 3 4 3 4 3 4 ...
    ##  $ Depth_m                              : int  4 4 4 4 4 4 4 4 4 4 ...
    ##  $ Depth_q                              : Factor w/ 1 level "SRF": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Water_mass                           : Factor w/ 3 levels "","AASW","SASW": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ RNA_vol_µl                           : num  60 NA 60 NA 60 NA 62 NA 62 NA ...
    ##  $ RNA_Conc_ng.l                        : num  0 NA 0 NA 0 NA 0 NA 0 NA ...
    ##  $ RNA_qty_ng                           : num  0 NA 0 NA 0 NA 0 NA 0 NA ...
    ##  $ Event_start_date                     : Factor w/ 146 levels "2016-12-23T15:00:00Z+00",..: 3 3 12 12 4 4 13 13 14 14 ...
    ##  $ Event_end_date                       : Factor w/ 146 levels "2016-12-23T15:00:00Z+00",..: 3 3 12 12 4 4 13 13 14 14 ...
    ##  $ Latitude                             : num  -43.2 -43.2 -46.2 -46.2 -44.2 ...
    ##  $ Longitude                            : num  32.4 32.4 46.9 46.9 33.9 ...
    ##  $ SOG                                  : num  7.31 7.31 7.48 7.48 7.06 ...
    ##  $ Distance                             : num  703607 703607 530419 530419 510244 ...
    ##  $ Temperature_degC                     : num  11.18 11.18 6.02 6.02 9.54 ...
    ##  $ Salinity                             : num  NA NA 33.8 33.8 NA ...
    ##  $ sigma0                               : num  NA NA 26.6 26.6 NA ...
    ##  $ MLD_m                                : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Sea_ice                              : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ RR_200m                              : num  0 0 0.123 0.123 0.465 ...
    ##  $ Snowfall_mm.h                        : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Total_precipitation_mm.h             : num  0 0 0.153 0.153 0.867 ...
    ##  $ Nitrate_uM                           : num  6.46 6.46 NA NA 13.99 ...
    ##  $ Nitrite_uM                           : num  0.102 0.102 NA NA 0.225 ...
    ##  $ Phosphate_µM                         : num  0.783 0.783 NA NA 1.269 ...
    ##  $ Silicate_µM                          : num  0.934 0.934 NA NA 3.15 ...
    ##  $ Ammonium_µM                          : num  NA NA NA NA NA NA NA NA 1.57 1.57 ...
    ##  $ POC_µM                               : num  13.1 13.1 NA NA 15.3 ...
    ##  $ PON_µM                               : num  1.73 1.73 NA NA 1.42 ...
    ##  $ NCP                                  : num  NA NA 12.4 12.4 NA ...
    ##  $ Chla_cor_µg.l                        : num  1.427 1.427 NA NA 0.649 ...
    ##  $ Total_Chla_Merge_mg.m3               : num  1.06 1.06 NA NA 0.662 ...
    ##  $ Fv_Fm_night                          : num  NA NA 0.248 0.248 NA ...
    ##  $ Fv_Fm_day                            : num  NA NA NA NA NA ...
    ##  $ Chloro                               : num  0.1181 0.1181 NA NA 0.0959 ...
    ##  $ Crypto1                              : num  0.0175 0.0175 NA NA 0.0195 0.0195 0.0008 0.0008 0.0117 0.0117 ...
    ##  $ Cyano2                               : num  0.0216 0.0216 NA NA 0.0266 0.0266 0.0021 0.0021 0.0014 0.0014 ...
    ##  $ DiatA                                : num  0.0201 0.0201 NA NA 0.0421 0.0421 0.043 0.043 0.0177 0.0177 ...
    ##  $ DiatB                                : num  0.0864 0.0864 NA NA 0.0114 0.0114 0 0 0.0138 0.0138 ...
    ##  $ DinoA                                : num  0.0907 0.0907 NA NA 0.0649 0.0649 0.0019 0.0019 0.0046 0.0046 ...
    ##  $ Hapto8                               : num  0.2454 0.2454 NA NA 0.0946 ...
    ##  $ Hapto67                              : num  0.369 0.369 NA NA 0.212 ...
    ##  $ Pras3                                : num  0.0236 0.0236 NA NA 0.0393 0.0393 0.0166 0.0166 0.0397 0.0397 ...
    ##  $ Pelago                               : num  0.0419 0.0419 NA NA 0.0483 0.0483 0.0204 0.0204 0.0336 0.0336 ...
    ##  $ HDNA_bacteria_cells.ml               : int  379636 379636 NA NA 730569 730569 191936 191936 444875 444875 ...
    ##  $ LDNA_bacteria_cells.ml               : int  340911 340911 NA NA 664419 664419 159909 159909 411845 411845 ...
    ##  $ Tot_bacteria_cells.ml                : int  720547 720547 NA NA 1394989 1394989 351845 351845 856720 856720 ...
    ##  $ Synechococcus_cells.ml               : int  85805 85805 NA NA 14422 14422 1689 1689 4461 4461 ...
    ##  $ Picoeukaryotes1_cells.ml             : int  8389 8389 NA NA 7053 7053 2030 2030 6879 6879 ...
    ##  $ Picoeukaryotes2_cells.ml             : int  8027 8027 NA NA 2328 2328 1731 1731 4188 4188 ...
    ##  $ Nanoeukaryotes_cells.ml              : int  1965 1965 NA NA 491 491 156 156 445 445 ...
    ##  $ Cryptomonas_cells.ml                 : int  204 204 NA NA 91 91 6 6 62 62 ...
    ##  $ Picoeukaryotes_cells.ml              : int  16416 16416 NA NA 9381 9381 3761 3761 11067 11067 ...
    ##  $ TEP                                  : num  51.2 51.2 NA NA 25.9 ...
    ##  $ CSP                                  : num  24.2 24.2 NA NA 5.83 5.83 4.59 4.59 4.29 4.29 ...
    ##  $ Total_DMSP_nM                        : num  68.3 68.3 NA NA 42 42 15.8 15.8 21.2 21.2 ...
    ##  $ DMS_nM                               : num  1.74 1.74 NA NA 0.9 0.9 0.28 0.28 0.47 0.47 ...
    ##  $ CDOM_ppb                             : num  -0.26 -0.26 NA NA -0.363 ...
    ##  $ Uabs_Glorys_mean                     : num  0.0785 0.0785 0.0866 0.0866 0.1478 ...
    ##  $ Ekin_Glorys_mean                     : num  0.00386 0.00386 0.00379 0.00379 0.01196 ...
    ##  $ Vorticity_Glorys_mean                : num  -0.705 -0.705 -0.174 -0.174 -0.907 ...
    ##  $ OW_Glorys_mean                       : num  0.0279 0.0279 0.0602 0.0602 -0.0833 ...
    ##  $ EulerDiverg_Glorys_mean              : num  0.0643 0.0643 -0.023 -0.023 0.0789 ...
    ##  $ LagrDiverg_Glorys_05daysbackward_mean: num  0.00665 0.00665 -0.02524 -0.02524 0.10595 ...
    ##  $ LagrDiverg_Glorys_10daysbackward_mean: num  -0.0178 -0.0178 -0.0311 -0.0311 0.0727 ...
    ##  $ LagrDiverg_Glorys_15daysbackward_mean: num  -0.03 -0.03 -0.0351 -0.0351 0.0563 ...
    ##  $ LagrDiverg_Glorys_20daysbackward_mean: num  -0.0177 -0.0177 -0.031 -0.031 0.0475 ...
    ##  $ LagrDiverg_Glorys_30daysbackward_mean: num  -0.0139 -0.0139 -0.0108 -0.0108 0.0295 ...
    ##  $ LagrDiverg_Glorys_60daysbackward_mean: num  -0.01215 -0.01215 -0.00586 -0.00586 0.01537 ...
    ##  $ RetentionTime_Glorys_mean            : num  1.02 1.02 1 1 1.35 ...
    ##  $ Ftle_Glorys_05daysBackward_mean      : num  0.294 0.294 0.153 0.153 0.239 ...
    ##  $ Ftle_Glorys_10daysBackward_mean      : num  0.2063 0.2063 0.0802 0.0802 0.1602 ...
    ##  $ Ftle_Glorys_15daysBackward_mean      : num  0.17 0.17 0.12 0.12 0.133 ...
    ##  $ Ftle_Glorys_20daysBackward_mean      : num  0.135 0.135 0.106 0.106 0.11 ...
    ##  $ Ftle_Glorys_30daysBackward_mean      : num  0.0969 0.0969 0.0717 0.0717 0.0789 ...
    ##  $ Ftle_Glorys_60daysBackward_mean      : num  0.0525 0.0525 0.0495 0.0495 0.0467 ...
    ##  $ Betw_Glorys_05daysCentered_mean      : num  4.57 4.57 2.42 2.42 5.48 ...
    ##  $ Betw_Glorys_10daysCentered_mean      : num  10.3 10.3 3.7 3.7 14.2 ...
    ##  $ Betw_Glorys_15daysCentered_mean      : num  20.38 20.38 6.44 6.44 23.06 ...
    ##  $ Betw_Glorys_20daysCentered_mean      : num  32.62 32.62 9.38 9.38 31.4 ...
    ##  $ Betw_Glorys_30daysCentered_mean      : num  58.1 58.1 26.3 26.3 48.4 ...
    ##  $ Betw_Glorys_60daysCentered_mean      : num  167.5 167.5 54.4 54.4 72.5 ...
    ##  $ DNA_Conc_ng_l_new                    : num  41.7 42.9 6.97 12.95 30.2 ...
    ##  $ DNA_vol_l_new                        : int  96 101 101 100 104 103 100 98 98 98 ...
    ##  $ DNA_qty_ng_new                       : num  4003 4333 704 1295 3141 ...
    ##  $ Longhurst_Prov                       : Factor w/ 6 levels "ANTA","APLR",..: 6 6 5 5 5 5 5 5 5 5 ...

``` r
metagenomes <- read.table("./Metadata/ACEsamples_With_Ace_seq_name.tsv", header = TRUE) # File listing samples with metagenomes
```

# Format data

The aim here is to have a matrix that can be roughly investigated
through a PCA so as to have a better overview of the station
characteristics and how they cluster. This requires getting rid of some
columns (especially the sequencing variables, that vary according to
filter), and to pool samples by event number and depth, regardless of
the filter.

``` r
meta <- meta %>% 
  mutate(MertzGlacier.new = case_when(MertzGlacier == "FAUX" ~ "FALSE",
                                      TRUE ~ "TRUE")) %>% 
  mutate(MertzGlacier = as.factor(MertzGlacier.new)) %>% 
  select(-MertzGlacier.new)

meta <- meta %>% 
  mutate_if(.predicate = is.character, .funs = as.factor)
```

``` r
# samples_id <- meta %>%
#   select(c("Event_number","CTD_cast_number","ACE_station_number","TM_station_number","ACE_seq_name","Latitude","Longitude"))

meta.filtered <- meta %>% 
  select(-c("CTD_cast_number","TM_station_number","DNA_vol_µl","DNA_Conc_ng.µl","DNA_qty_ng","RNA_vol_µl","RNA_Conc_ng.µl","RNA_qty_ng","Event_start_date","Event_end_date"))

meta.filtered.phys <- meta.filtered

meta.filtered <- meta.filtered %>% select(-c(contains("Glorys"))) # Leave out physical and dynamical properties
glimpse(meta.filtered)
```

    ## Rows: 345
    ## Columns: 59
    ## $ Event_number          <int> 75, 75, 75, 75, 75, 75, 123, 123, 123, 123, 123,…
    ## $ ACE_station_number    <int> 2, 2, 2, 2, 2, 2, 8, 8, 8, 8, 8, 8, 8, 16, 16, 1…
    ## $ ACE_seq_name          <fct> ACE1_75_11A_15m, ACE1_75_11B_15m, ACE1_75_13A_15…
    ## $ Size_fraction         <fct> >3 µm, 0.2-3 µm, >3 µm, 0.2-3 µm, >3 µm, 0.2-3 µ…
    ## $ Depth_m               <int> 15, 15, 150, 150, 60, 60, 15, 15, 15, 150, 150, …
    ## $ Depth_q               <fct> ZZZ, ZZZ, ZZZ, ZZZ, DCM, DCM, ZZZ, ZZZ, ZZZ, ZZZ…
    ## $ Latitude              <dbl> -46.91799, -46.91799, -46.91799, -46.91799, -46.…
    ## $ Longitude             <dbl> 38.11201, 38.11201, 38.11201, 38.11201, 38.11201…
    ## $ Temperature_ITS90     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 4.7296, 4.72…
    ## $ Salinity_PSS78        <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 33.8687, 33.…
    ## $ Density_kg.m2         <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 1027.501, 10…
    ## $ dO2_µmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 284.191, 284…
    ## $ Fluorescence_mg.m3    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.0430, 0.04…
    ## $ PAR_µmol.m2.s         <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.0195, 0.01…
    ## $ NOx_µmol.L            <dbl> 19.67, 19.67, 23.96, 23.96, 20.44, 20.44, 21.61,…
    ## $ Silicic_acid_µmol.L   <dbl> 2.95, 2.95, 12.64, 12.64, 3.92, 3.92, 7.04, 7.04…
    ## $ Nitrite_µmol.L        <dbl> 0.19, 0.19, 0.17, 0.17, 0.18, 0.18, 0.24, 0.24, …
    ## $ Phosphate_µmol.L      <dbl> 1.34, 1.34, 1.45, 1.45, 1.33, 1.33, 1.42, 1.42, …
    ## $ Ammonium_µmol.L       <dbl> 0.28, 0.28, 0.87, 0.87, 0.44, 0.44, 0.61, 0.61, …
    ## $ Nitrate_µmol.L        <dbl> 19.51, 19.51, 23.82, 23.82, 20.29, 20.29, 21.40,…
    ## $ Bsi_µM                <dbl> 1.086, 1.086, 1.166, 1.166, 3.186, 3.186, 0.191,…
    ## $ dCd_nmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ dFe_nmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ dZn_nmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ dCu_nmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ dNi_nmol.kg           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ Fe_Ligands.nM         <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ log_Kcond_FeL         <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ TPZT_µM.C             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ Humics_µg.SRFAeq.L    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ POC_µM                <dbl> NA, NA, 1.46, 1.46, 7.70, 7.70, 3.67, 3.67, 3.67…
    ## $ d13C                  <dbl> NA, NA, -25.00, -25.00, -24.36, -24.36, -26.81, …
    ## $ PON.µmol              <dbl> NA, NA, 0.32, 0.32, 1.37, 1.37, 0.65, 0.65, 0.65…
    ## $ d15N                  <dbl> NA, NA, 3.37, 3.37, 2.39, 2.39, 0.95, 0.95, 0.95…
    ## $ Bacteria_HDNA.cell.mL <int> 300000, 300000, 255000, 255000, 404000, 404000, …
    ## $ Bacteria_LDNA.cell.mL <int> 249000, 249000, 161000, 161000, 268000, 268000, …
    ## $ Tot_bacteria_cell.mL  <int> 549000, 549000, 416000, 416000, 671000, 671000, …
    ## $ Synechos_cell.mL      <int> 222, 222, NA, NA, 54, 54, NA, NA, NA, NA, NA, 20…
    ## $ Picoeuk1_cell.mL      <int> 46, 46, NA, NA, 409, 409, NA, NA, NA, NA, NA, 70…
    ## $ Picoeuk2_cell.mL      <int> 73, 73, NA, NA, 744, 744, NA, NA, NA, NA, NA, 21…
    ## $ Nanoeuk_cell.mL       <int> 138, 138, NA, NA, 546, 546, NA, NA, NA, NA, NA, …
    ## $ Cryptomonas_cell.mL   <int> 0, 0, NA, NA, 0, 0, 0, 0, 0, NA, NA, 2, 2, NA, N…
    ## $ TEP_µg.XG.eq.L        <dbl> 10.84, 10.84, 15.58, 15.58, 28.50, 28.50, 13.33,…
    ## $ CSP_µg.BSA.eq.L       <dbl> NA, NA, 7.07, 7.07, 7.68, 7.68, 6.42, 6.42, 6.42…
    ## $ Chlorophyll_chemtax   <dbl> 0.0116, 0.0116, 0.0063, 0.0063, 0.0136, 0.0136, …
    ## $ Crypto1_chemtax       <dbl> 0.0006, 0.0006, 0.0006, 0.0006, 0.0037, 0.0037, …
    ## $ Cyano2_chemtax        <dbl> 0.0041, 0.0041, 0.0000, 0.0000, 0.0001, 0.0001, …
    ## $ DiatA_chemtax         <dbl> 0.0756, 0.0756, 0.0506, 0.0506, 0.3182, 0.3182, …
    ## $ DiatB_chemtax         <dbl> 0.0178, 0.0178, 0.0342, 0.0342, 0.1200, 0.1200, …
    ## $ DinoA_chemtax         <dbl> 0.0075, 0.0075, 0.0000, 0.0000, 0.0158, 0.0158, …
    ## $ Hapto8_chemtax        <dbl> 0.0565, 0.0565, 0.0176, 0.0176, 0.1054, 0.1054, …
    ## $ Hapto67_chemtax       <dbl> 0.1001, 0.1001, 0.0174, 0.0174, 0.1184, 0.1184, …
    ## $ Pras3_chemtax         <dbl> 0.0082, 0.0082, 0.0006, 0.0006, 0.0104, 0.0104, …
    ## $ Pelago_chemtax        <dbl> 0.0110, 0.0110, 0.0010, 0.0010, 0.0204, 0.0204, …
    ## $ Fv_Fm                 <dbl> 0.158, 0.158, 0.342, 0.342, 0.355, 0.355, 0.433,…
    ## $ Longhurst_Prov        <fct> SANT, SANT, SANT, SANT, SANT, SANT, SANT, SANT, …
    ## $ MertzGlacier          <fct> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…
    ## $ Water_mass_detailed   <fct> SASW, SASW, SASW, SASW, SASW, SASW, SASW, SASW, …
    ## $ Water_mass            <fct> SASW, SASW, SASW, SASW, SASW, SASW, SASW, SASW, …

Here we group by event number and depth, hence the size fractions (which
affect the biology) will be pooled, as the physics and chemistry are
repeated values for the different filters (we could have taken values
from samples of SF \<0.2 µm as well). To create a metadata file that
matches metagenome samples, we will need to include the size fraction
column.

``` r
meta.sum <- meta.filtered %>%
  group_by(Event_number,Depth_m,Depth_q,Longhurst_Prov,MertzGlacier,Water_mass) %>% 
  summarise_if(.predicate = is.numeric, .funs = mean)

meta.sum <- meta.sum %>% 
  mutate(Sample_name=paste(Event_number,Depth_m, sep = "_")) %>% 
  column_to_rownames(var = "Sample_name")

head(meta.sum)
```

    ##         Event_number Depth_m Depth_q Longhurst_Prov MertzGlacier Water_mass
    ## 75_15             75      15     ZZZ           SANT        FALSE       SASW
    ## 75_60             75      60     DCM           SANT        FALSE       SASW
    ## 75_150            75     150     ZZZ           SANT        FALSE       SASW
    ## 123_15           123      15     ZZZ           SANT        FALSE       SASW
    ## 123_60           123      60     DCM           SANT        FALSE       SASW
    ## 123_150          123     150     ZZZ           SANT        FALSE       SASW
    ##         ACE_station_number  Latitude Longitude Temperature_ITS90 Salinity_PSS78
    ## 75_15                    2 -46.91799  38.11201                NA             NA
    ## 75_60                    2 -46.91799  38.11201                NA             NA
    ## 75_150                   2 -46.91799  38.11201                NA             NA
    ## 123_15                   8 -46.17711  51.00348                NA             NA
    ## 123_60                   8 -46.17711  51.00348            5.4717        33.8171
    ## 123_150                  8 -46.17711  51.00348            4.7296        33.8687
    ##         Density_kg.m2 dO2_µmol.kg Fluorescence_mg.m3 PAR_µmol.m2.s NOx_µmol.L
    ## 75_15              NA          NA                 NA            NA      19.67
    ## 75_60              NA          NA                 NA            NA      20.44
    ## 75_150             NA          NA                 NA            NA      23.96
    ## 123_15             NA          NA                 NA            NA      21.61
    ## 123_60       1026.960     289.703             0.1478        2.6775      21.85
    ## 123_150      1027.501     284.191             0.0430        0.0195      23.62
    ##         Silicic_acid_µmol.L Nitrite_µmol.L Phosphate_µmol.L Ammonium_µmol.L
    ## 75_15                  2.95           0.19             1.34            0.28
    ## 75_60                  3.92           0.18             1.33            0.44
    ## 75_150                12.64           0.17             1.45            0.87
    ## 123_15                 7.04           0.24             1.42            0.61
    ## 123_60                 7.02           0.23             1.50            0.44
    ## 123_150               10.55           0.19             1.46            0.33
    ##         Nitrate_µmol.L Bsi_µM dCd_nmol.kg dFe_nmol.kg dZn_nmol.kg dCu_nmol.kg
    ## 75_15            19.51  1.086          NA          NA          NA          NA
    ## 75_60            20.29  3.186          NA          NA          NA          NA
    ## 75_150           23.82  1.166          NA          NA          NA          NA
    ## 123_15           21.40  0.191          NA          NA          NA          NA
    ## 123_60           21.65  0.364          NA          NA          NA          NA
    ## 123_150          23.45  0.517          NA          NA          NA          NA
    ##         dNi_nmol.kg Fe_Ligands.nM log_Kcond_FeL TPZT_µM.C Humics_µg.SRFAeq.L
    ## 75_15            NA            NA            NA        NA                 NA
    ## 75_60            NA            NA            NA        NA                 NA
    ## 75_150           NA            NA            NA        NA                 NA
    ## 123_15           NA            NA            NA        NA                 NA
    ## 123_60           NA            NA            NA        NA                 NA
    ## 123_150          NA            NA            NA        NA                 NA
    ##         POC_µM   d13C PON.µmol d15N Bacteria_HDNA.cell.mL Bacteria_LDNA.cell.mL
    ## 75_15       NA     NA       NA   NA                300000                249000
    ## 75_60     7.70 -24.36     1.37 2.39                404000                268000
    ## 75_150    1.46 -25.00     0.32 3.37                255000                161000
    ## 123_15    3.67 -26.81     0.65 0.95                118000                111000
    ## 123_60    3.33 -26.66     0.62 6.45                432000                381000
    ## 123_150   1.20 -25.40     0.28 6.19                535000                433000
    ##         Tot_bacteria_cell.mL Synechos_cell.mL Picoeuk1_cell.mL Picoeuk2_cell.mL
    ## 75_15                 549000              222               46               73
    ## 75_60                 671000               54              409              744
    ## 75_150                416000               NA               NA               NA
    ## 123_15                229000               NA               NA               NA
    ## 123_60                814000             2013             7081             2144
    ## 123_150               968000               NA               NA               NA
    ##         Nanoeuk_cell.mL Cryptomonas_cell.mL TEP_µg.XG.eq.L CSP_µg.BSA.eq.L
    ## 75_15               138                   0          10.84              NA
    ## 75_60               546                   0          28.50            7.68
    ## 75_150               NA                  NA          15.58            7.07
    ## 123_15               NA                   0          13.33            6.42
    ## 123_60              996                   2          18.10            9.47
    ## 123_150              NA                  NA             NA              NA
    ##         Chlorophyll_chemtax Crypto1_chemtax Cyano2_chemtax DiatA_chemtax
    ## 75_15                0.0116          0.0006         0.0041        0.0756
    ## 75_60                0.0136          0.0037         0.0001        0.3182
    ## 75_150               0.0063          0.0006         0.0000        0.0506
    ## 123_15               0.0560          0.0256         0.0183        0.0081
    ## 123_60               0.0487          0.0300         0.0184        0.0499
    ## 123_150              0.0000          0.0039         0.0061        0.0127
    ##         DiatB_chemtax DinoA_chemtax Hapto8_chemtax Hapto67_chemtax
    ## 75_15          0.0178        0.0075         0.0565          0.1001
    ## 75_60          0.1200        0.0158         0.1054          0.1184
    ## 75_150         0.0342        0.0000         0.0176          0.0174
    ## 123_15         0.0296        0.0001         0.0472          0.0478
    ## 123_60         0.0195        0.0002         0.0705          0.1052
    ## 123_150        0.0048        0.0000         0.0197          0.0222
    ##         Pras3_chemtax Pelago_chemtax Fv_Fm
    ## 75_15          0.0082         0.0110 0.158
    ## 75_60          0.0104         0.0204 0.355
    ## 75_150         0.0006         0.0010 0.342
    ## 123_15         0.1544         0.0023 0.433
    ## 123_60         0.1230         0.0279 0.425
    ## 123_150        0.0268         0.0026 0.283

``` r
dim(meta.sum)
```

    ## [1] 150  56

``` r
meta.group <- meta.sum %>% 
  select(c("Event_number","Depth_m","Depth_q","Longhurst_Prov","MertzGlacier","Water_mass")) # Store grouping variables

meta.sum <- meta.sum %>%
  select(-c("Event_number","ACE_station_number","Depth_q","Longhurst_Prov","MertzGlacier","Water_mass")) # Remove grouping variables
```

``` r
meta.phys <- meta.filtered.phys %>%
  group_by(Event_number,Depth_m,Depth_q,Longhurst_Prov,MertzGlacier,Water_mass) %>% 
  summarise_if(.predicate = is.numeric, .funs = mean, na.rm = TRUE) %>% 
  mutate(Sample_name=paste(Event_number,Depth_m, sep = "_")) %>% 
  column_to_rownames(var = "Sample_name")

dim(meta.phys)
```

    ## [1] 150  80

``` r
meta.phys <- meta.phys %>% 
  select(-c("Event_number","ACE_station_number","Depth_q","Longhurst_Prov","MertzGlacier","Water_mass"))
```

# Grouping variables

## Qualitative depth

Modify the depth qualitative variables to another variable that
discriminates the `ZZZ` depth better.

``` r
meta.group <- meta.group %>% 
  mutate(Depth_q_new = as.factor(case_when(Depth_m <= 15 ~ "SRF",
                                           Depth_m == 150 ~ "150m",
                                           Depth_m == 1000 ~ "1000m",
                                           Depth_m >= 3460 ~ "Deep",
                                           Depth_q == "DCM" ~ "DCM",
                                           TRUE & Depth_m < 150 ~ "ZZZ < 150",
                                           TRUE & Depth_m > 150 ~ "ZZZ > 150")))

meta.group$Depth_q_new <- factor(meta.group$Depth_q_new, c("SRF","DCM","ZZZ < 150","150m", "ZZZ > 150", "1000m","Deep"))

table(meta.group$Depth_q_new)
```

    ## 
    ##       SRF       DCM ZZZ < 150      150m ZZZ > 150     1000m      Deep 
    ##        60         3        17        37        12        18         3

``` r
table(meta.group$Depth_q)
```

    ## 
    ## DCM SRF ZZZ 
    ##   3  28 119

## Format water mass

``` r
meta.group <- meta.group %>% 
  mutate(Water_mass =
           factor(Water_mass,
                  levels = c("AASW",
                             "SASW",
                             "STSW",
                             "AASW-WW",
                             "WW",
                             "AASW-AAIW",
                             "SASW-AAIW",
                             "STSW-AAIW",
                             "AAIW-UCDW",
                             "WW-CDW",
                             "CDW",
                             "LCDW-AABW",
                             "AABW",
                             "DSW"))) %>% 
  mutate(Water_mass_simplified =
           case_when(Water_mass %in% c("AASW-WW","WW") ~ "WW influenced",
                     Water_mass %in% c("AASW-AAIW","SASW-AAIW","STSW-AAIW") ~ "SW-AAIW",
                     Water_mass %in% c("AAIW-UCDW","WW-CDW","CDW","LCDW-AABW") ~ "CDW influenced",
                     Water_mass == "SASW" ~ "SASW",
                     Water_mass == "STSW" ~ "STSW",
                     Water_mass == "AASW" ~ "AASW",
                     Water_mass == "DSW" ~ "DSW",
                     Water_mass == "AABW" ~ "AABW")) %>%
  mutate(Water_mass_simplified = factor(Water_mass_simplified, levels = c("AASW","SASW", "STSW","WW influenced","SW-AAIW","CDW influenced","DSW","AABW")))

head(meta.group)
```

    ##         Event_number Depth_m Depth_q Longhurst_Prov MertzGlacier Water_mass
    ## 75_15             75      15     ZZZ           SANT        FALSE       SASW
    ## 75_60             75      60     DCM           SANT        FALSE       SASW
    ## 75_150            75     150     ZZZ           SANT        FALSE       SASW
    ## 123_15           123      15     ZZZ           SANT        FALSE       SASW
    ## 123_60           123      60     DCM           SANT        FALSE       SASW
    ## 123_150          123     150     ZZZ           SANT        FALSE       SASW
    ##         Depth_q_new Water_mass_simplified
    ## 75_15           SRF                  SASW
    ## 75_60           DCM                  SASW
    ## 75_150         150m                  SASW
    ## 123_15          SRF                  SASW
    ## 123_60          DCM                  SASW
    ## 123_150        150m                  SASW

Because it is a qualitative variable with many levels, we need to create
a palette for water masses.

``` r
scales::show_col(colorRampPalette(pal_npg()(10))(14))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Create%20palette%20for%20water%20masses-1.png)<!-- -->

``` r
wm.pal <- c("#2FB0B7","#7C98A3","#E64B35","#049A87","#2E6587","#907483","#8BB3BB","#E19987","#9592AB","#A2A095","#D6100E","#A23B2C","#8D735A","#B09C85")
scales::show_col(wm.pal) # For all water masses
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Create%20palette%20for%20water%20masses-2.png)<!-- -->

``` r
wm.pal.sim <- pal_npg()(8)[c(2,3,1,4:8)]
scales::show_col(wm.pal.sim)
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Create%20palette%20for%20water%20masses-3.png)<!-- -->

``` r
meta.group <- meta.group %>% 
  mutate(Water_mass.col = Water_mass,
         Water_mass_simplified.col = Water_mass_simplified)
levels(meta.group$Water_mass.col) <- wm.pal
levels(meta.group$Water_mass_simplified.col) <- wm.pal.sim
```

# Missing values

## Explore missing values

``` r
dim(na.omit(meta.sum)) # A lot of lines are not complete
```

    ## [1] 10 50

``` r
gg_miss_var(meta.sum)
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Explore%20NA-1.png)<!-- -->

``` r
meta.sum %>% miss_var_summary()
```

    ## # A tibble: 50 × 3
    ##    variable            n_miss pct_miss
    ##    <chr>                <int>    <num>
    ##  1 Fe_Ligands.nM          109     72.7
    ##  2 log_Kcond_FeL          109     72.7
    ##  3 Humics_µg.SRFAeq.L      98     65.3
    ##  4 CSP_µg.BSA.eq.L         96     64  
    ##  5 TEP_µg.XG.eq.L          95     63.3
    ##  6 Synechos_cell.mL        93     62  
    ##  7 Picoeuk1_cell.mL        92     61.3
    ##  8 Picoeuk2_cell.mL        92     61.3
    ##  9 Nanoeuk_cell.mL         92     61.3
    ## 10 Cryptomonas_cell.mL     91     60.7
    ## # ℹ 40 more rows

``` r
vis_miss(meta.sum, sort_miss = TRUE)
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Explore%20NA-2.png)<!-- -->

``` r
res <- summary(aggr(meta.sum, sortVar = TRUE, cex.axis = .4))$combinations
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Explore%20NA-3.png)<!-- -->

    ## 
    ##  Variables sorted by number of missings: 
    ##               Variable      Count
    ##          Fe_Ligands.nM 0.72666667
    ##          log_Kcond_FeL 0.72666667
    ##     Humics_µg.SRFAeq.L 0.65333333
    ##        CSP_µg.BSA.eq.L 0.64000000
    ##         TEP_µg.XG.eq.L 0.63333333
    ##       Synechos_cell.mL 0.62000000
    ##       Picoeuk1_cell.mL 0.61333333
    ##       Picoeuk2_cell.mL 0.61333333
    ##        Nanoeuk_cell.mL 0.61333333
    ##    Cryptomonas_cell.mL 0.60666667
    ##                  Fv_Fm 0.55333333
    ##    Chlorophyll_chemtax 0.53333333
    ##        Crypto1_chemtax 0.53333333
    ##         Cyano2_chemtax 0.53333333
    ##          DiatA_chemtax 0.53333333
    ##          DiatB_chemtax 0.53333333
    ##          DinoA_chemtax 0.53333333
    ##         Hapto8_chemtax 0.53333333
    ##        Hapto67_chemtax 0.53333333
    ##          Pras3_chemtax 0.53333333
    ##         Pelago_chemtax 0.53333333
    ##            dCu_nmol.kg 0.52666667
    ##            dNi_nmol.kg 0.52666667
    ##            dCd_nmol.kg 0.42666667
    ##            dFe_nmol.kg 0.42000000
    ##            dZn_nmol.kg 0.42000000
    ##  Bacteria_HDNA.cell.mL 0.42000000
    ##  Bacteria_LDNA.cell.mL 0.42000000
    ##   Tot_bacteria_cell.mL 0.42000000
    ##              TPZT_µM.C 0.41333333
    ##                 POC_µM 0.37333333
    ##                   d13C 0.37333333
    ##               PON.µmol 0.37333333
    ##                   d15N 0.37333333
    ##     Fluorescence_mg.m3 0.28666667
    ##          PAR_µmol.m2.s 0.28666667
    ##                 Bsi_µM 0.14666667
    ##             NOx_µmol.L 0.10000000
    ##    Silicic_acid_µmol.L 0.10000000
    ##         Nitrite_µmol.L 0.10000000
    ##       Phosphate_µmol.L 0.10000000
    ##        Ammonium_µmol.L 0.10000000
    ##         Nitrate_µmol.L 0.10000000
    ##          Density_kg.m2 0.07333333
    ##            dO2_µmol.kg 0.07333333
    ##         Salinity_PSS78 0.04000000
    ##      Temperature_ITS90 0.03333333
    ##                Depth_m 0.00000000
    ##               Latitude 0.00000000
    ##              Longitude 0.00000000

``` r
head(res[rev(order(res[,2])),])
```

    ##                                                                                           Combinations
    ## 36 0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:1:1:1:1:1:1:1:1:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0
    ## 1  0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0
    ## 51 0:0:0:0:0:0:0:0:0:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1
    ## 57 0:0:0:0:0:0:0:1:1:0:0:0:0:0:0:0:0:0:0:1:1:1:1:0:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1
    ## 52 0:0:0:0:0:0:0:1:1:0:0:0:0:0:0:0:0:0:0:1:1:1:1:0:1:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0
    ## 26 0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:1:0:0:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1:1
    ##    Count  Percent
    ## 36    10 6.666667
    ## 1     10 6.666667
    ## 51     8 5.333333
    ## 57     6 4.000000
    ## 52     5 3.333333
    ## 26     5 3.333333

Biological variables are often missing simultaneously. dCd, dFe, dZn,
dCu, dNi are usually missing together. Fe_Ligands, log_Kcond_FeL, TPZT,
Humics include a lot of missing values

``` r
meta.sum.cat <- apply(meta.sum, FUN=function(x){if(is.na(x)) "m" else "o"}, MARGIN = c(1,2))

res.mca <- MCA(meta.sum.cat, graph = F)
plot(res.mca, invis = "ind", title = "MCA graph of the categories", cex = 0.5)
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Multiple%20Correspondance%20Analysis-1.png)<!-- -->

## Imputation

Preprocess the matrix by removing variables that have over 50% NA, since
it would not really make sense to impute more than half the values.

``` r
var <- meta.sum %>%
  miss_var_summary() %>%
  filter(pct_miss < 50) %>%
  pull(variable) # Variables with less than 50% NA

meta.sum.preprocess <- meta.sum[,var] # %NA < 50%
meta.sum.preprocess %>% miss_var_summary() # All variables have below 43% NA
```

    ## # A tibble: 27 × 3
    ##    variable              n_miss pct_miss
    ##    <chr>                  <int>    <num>
    ##  1 dCd_nmol.kg               64     42.7
    ##  2 dFe_nmol.kg               63     42  
    ##  3 dZn_nmol.kg               63     42  
    ##  4 Bacteria_HDNA.cell.mL     63     42  
    ##  5 Bacteria_LDNA.cell.mL     63     42  
    ##  6 Tot_bacteria_cell.mL      63     42  
    ##  7 TPZT_µM.C                 62     41.3
    ##  8 POC_µM                    56     37.3
    ##  9 d13C                      56     37.3
    ## 10 PON.µmol                  56     37.3
    ## # ℹ 17 more rows

``` r
set.seed(3)
meta.sum.bag <- preProcess(meta.sum.preprocess, method = c("bagImpute")) # Centered and scaled later on
meta.sum.bagImp <- predict(meta.sum.bag, meta.sum.preprocess)
```

Relevant post regarding how to compare the raw matrix and matrix with
imputation of NAs:
<http://juliejosse.com/wp-content/uploads/2018/06/DataAnalysisMissingR.html>.

> Bagging (Bootstrap aggregating) was originally proposed by Leo
> Breiman. It is one of the earliest ensemble methods (L, n.d.). When
> used in missing value imputation, it will use the remaining variables
> as predictors to train a bagging tree and then use the tree to predict
> the missing values. Although theoretically, the method is powerful,
> the computation is much more intense than KNN. In practice, there is a
> trade-off between computation time and the effect. If a median or mean
> meet the modeling needs, even bagging tree may improve the accuracy a
> little, but the upgrade is so marginal that it does not deserve the
> extra time. The bagging tree itself is a model for regression and
> classification. (see
> <https://scientistcafe.com/ids/missing-values.html>).

One way to compare imputed matrix and raw matrix would be to impute
missing values with the mean of the column, using the $k$ nearest
neighbours, and a bagging tree in order to compare the output and how it
deviates from the initial matrix.

``` r
meta.sum.knn <- preProcess(meta.sum.preprocess, method = c("knnImpute")) # knn automatically scales and centers the data
meta.sum.knnImp <- predict(meta.sum.knn, meta.sum.preprocess)
```

``` r
NA2mean <- function(x){replace(x, is.na(x), mean(x, na.rm = TRUE))}
meta.sum.mean <- meta.sum.preprocess %>%
  mutate(across(.cols = everything(), .fns = NA2mean))
```

## Inspect imputed data sets

``` r
formatLong <- function(table, method){
  if (method != "knn"){
    table %>%
      scale() %>%
      as.data.frame %>% 
      rownames_to_column("Sample") %>% 
      select(-c("Depth_m","Longitude","Latitude")) %>%
      pivot_longer(cols = !c("Sample"), names_to = "Variable", values_to = "Value") %>%
      mutate(Method = as.factor(method))
  } else {
    table %>%
      as.data.frame %>% 
      rownames_to_column("Sample") %>% 
      select(-c("Depth_m","Longitude","Latitude")) %>%
      pivot_longer(cols = !c("Sample"), names_to = "Variable", values_to = "Value") %>%
      mutate(Method = as.factor(method))
  }
}

comp.raw <- formatLong(meta.sum[,var], method = "Raw") 
comp.bag <- formatLong(meta.sum.bagImp, method = "Bag")
comp.knn <- formatLong(meta.sum.knnImp, method = "knn")
comp.mean <- formatLong(meta.sum.mean, method = "Mean")

comp <- comp.raw %>% 
  bind_rows(comp.bag, comp.knn, comp.mean)
```

``` r
plotMethod <- function(method){
  ggplot(data = comp %>% filter(Method == "Raw" | Method == method)) +
    geom_boxplot(aes(x = Variable, y = Value, fill = Method),
                 outlier.shape = "diamond",
                 notch = TRUE) +
    scale_fill_manual(values = c("#44015499","#35B77999")) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

(plot.mean <- plotMethod("Mean"))
```

    ## Warning: Removed 896 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?

![](02-ACE_CTD_Metadata_files/figure-gfm/Compare%20distributions-1.png)<!-- -->

``` r
(plot.bag <- plotMethod("Bag"))
```

    ## Warning: Removed 896 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?

![](02-ACE_CTD_Metadata_files/figure-gfm/Compare%20distributions-2.png)<!-- -->

``` r
(plot.knn <- plotMethod("knn"))
```

    ## Warning: Removed 896 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?

![](02-ACE_CTD_Metadata_files/figure-gfm/Compare%20distributions-3.png)<!-- -->

``` r
# plot_grid(plot.mean, plot.knn, plot.bag, nrow = 3)

(plot.comp <- ggplot(data = comp %>% mutate(Variable = factor(Variable, levels=(miss_var_summary(meta.sum[,var]) %>% pull(variable))))) + # Order by number of NAs
  geom_boxplot(aes(x = Variable, y = Value, fill = Method),
               outlier.shape = "diamond",
               notch = TRUE) +
  scale_fill_npg(alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.background = element_rect(fill='transparent', color=NA)))
```

    ## Warning: Removed 896 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?
    ## Notch went outside hinges
    ## ℹ Do you want `notch = FALSE`?

![](02-ACE_CTD_Metadata_files/figure-gfm/Compare%20distributions-4.png)<!-- -->

We see the imputation generally is bad when NA are replaced by the mean
of the column. The comparison between the k nearest neighbours and
bagging tree is more tricky. It seems that there are less outliers
introduced by the bagging tree method. As it seems that the notch
usually has a better match between the raw data with the bagging tree
method, this is the one we will further use.

Extract some stats to compare the imputations, especially number of
outliers as it is important to not introduce outliers during imputation,
the diff summaries…

``` r
# to check the diff between quartiles, median and mean in the two tables
tmp <- do.call(cbind, lapply(meta.sum.preprocess %>% scale() %>% as.data.frame(), summary))
```

    ## Warning in base::cbind(...): number of rows of result is not a multiple of
    ## vector length (arg 25)

``` r
tmp2 <- do.call(cbind, lapply(meta.sum.bagImp %>% scale() %>% as.data.frame(), summary))
tmp[1:6,]-tmp2
```

    ##           dCd_nmol.kg   dFe_nmol.kg   dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## Min.     5.737254e-03 -6.629558e-02 -1.177908e-01          8.396407e-02
    ## 1st Qu.  9.713653e-02 -1.390101e-01 -9.850417e-02         -2.398674e-02
    ## Median   2.663943e-01 -3.435732e-03  2.070891e-01          4.570108e-02
    ## Mean     1.267135e-15  2.563467e-16  7.177018e-16         -1.354983e-16
    ## 3rd Qu. -2.192620e-01  1.881463e-01 -7.849654e-02          3.784534e-01
    ## Max.    -3.377618e-01 -9.191893e-01 -3.982263e-01         -8.148506e-01
    ##         Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL     TPZT_µM.C        POC_µM
    ## Min.             2.078216e-02         5.941975e-02  3.441615e-01 -9.312663e-02
    ## 1st Qu.         -6.952591e-02         9.828463e-02 -1.228736e-01 -6.242232e-02
    ## Median          -2.481811e-02        -1.985100e-02  4.330339e-02  1.388920e-01
    ## Mean            -1.918516e-16         9.716748e-16  9.259596e-16  1.570493e-16
    ## 3rd Qu.          7.991186e-02         2.712425e-02  1.653747e-01 -1.216111e-01
    ## Max.            -1.052482e+00        -9.089018e-01 -1.178685e+00 -6.589537e-01
    ##                  d13C      PON.µmol          d15N Fluorescence_mg.m3
    ## Min.     4.732323e-01 -7.605107e-02  6.962981e-01       1.486723e-01
    ## 1st Qu. -1.131770e-01 -9.176385e-02 -9.656553e-03      -2.080589e-01
    ## Median  -3.342022e-02  9.870551e-02 -2.702054e-01       6.191080e-02
    ## Mean    -6.041818e-15  7.702507e-16  3.893025e-16       1.098291e-16
    ## 3rd Qu.  3.891516e-01 -6.532424e-02  7.403967e-02      -7.877803e-02
    ## Max.    -4.589928e-01 -6.356844e-01 -3.388965e-01      -8.832689e-01
    ##         PAR_µmol.m2.s        Bsi_µM    NOx_µmol.L Silicic_acid_µmol.L
    ## Min.     3.743393e-02  8.116577e-02  1.407041e-01        3.547616e-02
    ## 1st Qu.  3.743393e-02  7.888773e-02 -6.482475e-02       -8.548377e-02
    ## Median   1.036942e-02  8.523750e-02 -2.473458e-02        3.832020e-02
    ## Mean     1.868979e-16  1.956768e-16 -1.167790e-16       -1.805140e-17
    ## 3rd Qu. -2.677625e-01 -1.283780e-01 -1.225677e-02        2.083346e-02
    ## Max.    -9.207940e-01  5.851467e-01 -5.340992e-02       -1.057960e-01
    ##         Nitrite_µmol.L Phosphate_µmol.L Ammonium_µmol.L Nitrate_µmol.L
    ## Min.      7.416629e-02     1.184163e-01    8.016922e-03   1.488703e-01
    ## 1st Qu.  -1.940539e-02    -5.548997e-02    8.016922e-03  -3.775706e-02
    ## Median   -5.890860e-02     3.092709e-03   -2.992340e-02  -1.907422e-01
    ## Mean      1.754975e-16    -9.865360e-16   -1.400526e-16   1.559247e-16
    ## 3rd Qu.   1.298892e-01     6.181416e-02    1.018596e-02   1.242489e-02
    ## Max.     -3.305671e-02    -7.806096e-02   -3.142740e-01  -3.067303e-02
    ##         Density_kg.m2   dO2_µmol.kg Salinity_PSS78 Temperature_ITS90 Depth_m
    ## Min.    -1.722583e-02  7.079725e-02   1.363233e-02      2.870599e-02       0
    ## 1st Qu. -2.322334e-02 -1.045638e-01  -3.665762e-03      8.386944e-03       0
    ## Median   1.154692e-02  5.784513e-03  -7.710750e-03     -7.972038e-03       0
    ## Mean     1.741397e-14 -1.812566e-17  -4.967631e-15     -5.806339e-19       0
    ## 3rd Qu. -2.039281e-02  1.737397e-02  -8.941245e-03     -8.904068e-02       0
    ## Max.    -2.009899e-01 -2.703834e-02  -3.664203e-02      3.003191e-02       0
    ##         Latitude Longitude
    ## Min.           0         0
    ## 1st Qu.        0         0
    ## Median         0         0
    ## Mean           0         0
    ## 3rd Qu.        0         0
    ## Max.           0         0

``` r
tmp <- do.call(cbind, lapply(meta.sum.preprocess %>% scale() %>% as.data.frame(), summary))
```

    ## Warning in base::cbind(...): number of rows of result is not a multiple of
    ## vector length (arg 25)

``` r
tmp2 <- do.call(cbind, lapply(meta.sum.knnImp %>% as.data.frame(), summary))
tmp[1:6,]-tmp2
```

    ##          dCd_nmol.kg  dFe_nmol.kg   dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## Min.    6.661338e-16  0.000000000  6.661338e-16          2.220446e-16
    ## 1st Qu. 1.629046e-01 -0.009749753 -4.475651e-02         -4.152774e-03
    ## Median  5.536308e-01  0.142010979  5.915958e-01          1.994454e-01
    ## Mean    1.630843e-01  0.256731537  2.139976e-01          1.302703e-01
    ## 3rd Qu. 3.237095e-02  0.489903769  1.536640e-01          2.511306e-01
    ## Max.    0.000000e+00  0.000000000  1.110223e-15         -4.440892e-16
    ##         Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL     TPZT_µM.C       POC_µM
    ## Min.             4.440892e-16        -2.220446e-16  1.776357e-15 2.220446e-16
    ## 1st Qu.          4.042377e-02         1.034505e-01 -5.929440e-02 1.346849e-01
    ## Median           6.917551e-02         1.080784e-01  2.102566e-01 3.079729e-01
    ## Mean             1.692983e-01         1.544522e-01  1.236143e-01 2.097335e-01
    ## 3rd Qu.          2.757383e-01         2.186964e-01  2.044634e-01 1.994018e-01
    ## Max.            -4.440892e-16         0.000000e+00  8.881784e-16 8.881784e-16
    ##                  d13C   PON.µmol          d15N Fluorescence_mg.m3 PAR_µmol.m2.s
    ## Min.    -5.773160e-15 0.00000000 -8.881784e-16       1.110223e-16  2.775558e-16
    ## 1st Qu. -2.189695e-01 0.09636218 -2.513063e-01      -4.099383e-01  2.775558e-16
    ## Median   8.496692e-02 0.22237427 -2.020242e-01       4.114298e-02 -7.365470e-03
    ## Mean     5.460500e-02 0.20546541 -1.262638e-01      -3.692395e-02 -5.885303e-03
    ## 3rd Qu.  5.643325e-01 0.24090546 -1.370762e-02      -1.620005e-02 -3.054629e-01
    ## Max.    -5.773160e-15 0.00000000  4.440892e-16       0.000000e+00  0.000000e+00
    ##                Bsi_µM   NOx_µmol.L Silicic_acid_µmol.L Nitrite_µmol.L
    ## Min.     2.220446e-16 0.000000e+00        6.661338e-16     0.00000000
    ## 1st Qu. -6.791036e-03 8.534131e-02       -1.050021e-01    -0.08792836
    ## Median  -3.553976e-02 1.461096e-02       -5.044980e-02    -0.08792836
    ## Mean    -1.545431e-01 4.415964e-02        8.455974e-03    -0.04167153
    ## 3rd Qu. -1.316329e-01 1.494303e-02        3.525996e-02     0.00000000
    ## Max.     0.000000e+00 4.440892e-16        4.440892e-16     0.00000000
    ##         Phosphate_µmol.L Ammonium_µmol.L Nitrate_µmol.L Density_kg.m2
    ## Min.        1.332268e-15   -1.110223e-16   4.440892e-16  1.507683e-13
    ## 1st Qu.     8.803663e-02   -1.110223e-16   8.669153e-02 -1.062138e-02
    ## Median      2.515332e-02   -3.139635e-02   1.002062e-02  4.609217e-02
    ## Mean        5.137722e-02    8.947572e-03   4.451645e-02  3.655738e-02
    ## 3rd Qu.     8.908469e-02   -7.849087e-03   1.878866e-02  1.768576e-02
    ## Max.        1.998401e-15    0.000000e+00   4.440892e-16  1.483258e-13
    ##           dO2_µmol.kg Salinity_PSS78 Temperature_ITS90       Depth_m
    ## Min.     0.000000e+00   2.531308e-14     -2.220446e-16 -3.885781e-16
    ## 1st Qu. -1.415930e-01  -8.037361e-03     -2.049873e-02 -3.885781e-16
    ## Median   5.918555e-05   1.370005e-02     -3.702244e-02 -3.330669e-16
    ## Mean    -1.624531e-02   5.297484e-03     -2.284202e-02 -3.745152e-16
    ## 3rd Qu.  1.358731e-02   2.137207e-02     -7.583638e-02 -3.330669e-16
    ## Max.     0.000000e+00   2.642331e-14      1.332268e-15  0.000000e+00
    ##              Latitude     Longitude
    ## Min.    -2.220446e-15 -4.440892e-16
    ## 1st Qu. -2.109424e-15  1.110223e-16
    ## Median  -1.679212e-15  2.498002e-16
    ## Mean    -1.634248e-15  3.397282e-16
    ## 3rd Qu. -1.332268e-15  3.330669e-16
    ## Max.    -8.881784e-16  4.440892e-16

``` r
computeOutliers <- function(table){
  iqr <- do.call(cbind, lapply(table %>% select(-c("Depth_m", "Latitude", "Longitude")), IQR, na.rm = TRUE))
  lowerq <- do.call(cbind, lapply(table, quantile, na.rm = TRUE))[2,] # lower quantile
  upperq <- do.call(cbind, lapply(table, quantile, na.rm = TRUE))[4,] # upper quantile
  outliers.stats <- rbind(iqr, lowerq, upperq)
  rownames(outliers.stats) <- c("iqr", "lowerq", "upperq")
  
  outliers.stats <- outliers.stats %>% 
    t() %>% 
    as.data.frame() %>%
    # rownames_to_column(var = "Variable") %>% 
    mutate(across(.cols = everything(), .fns = as.numeric),
           Mild_threshold_upper = iqr*1.5 + upperq,
           Mild_threshold_lower = lowerq - iqr*1.5) %>%
    t() %>%
    as.data.frame()
  
  n_outliers <- c()
  for (i in 1:length(1:ncol(table))){
    ind <- which(table[,i] > outliers.stats[4,i] | table[,i] < outliers.stats[5,i])
    n_outliers <- c(n_outliers, length(ind))
  }
  
  outliers.stats <- rbind(outliers.stats, n_outliers)
  rownames(outliers.stats) <- c(rownames(outliers.stats)[-6],"Number_outliers")
  return(outliers.stats)
}
```

``` r
(outliers.stats.raw <-  computeOutliers(table = meta.sum.preprocess))
```

    ## Warning in base::rbind(...): number of columns of result is not a multiple of
    ## vector length (arg 2)

    ## Warning in rbind(deparse.level, ...): number of columns of result, 24, is not a
    ## multiple of vector length 27 of arg 2

    ##                      dCd_nmol.kg dFe_nmol.kg dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## iqr                     0.449750     0.27320     4.50450                225500
    ## lowerq                  0.323250     0.03665     0.77800                113000
    ## upperq                  0.773000     0.30985     5.28250                338500
    ## Mild_threshold_upper    1.447625     0.71965    12.03925                676750
    ## Mild_threshold_lower   -0.351375    -0.37315    -5.97875               -225250
    ## Number_outliers         0.000000     1.00000     0.00000                     3
    ##                      Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL TPZT_µM.C
    ## iqr                                 202550               330000     3.390
    ## lowerq                               79450               221000     5.150
    ## upperq                              282000               551000     8.540
    ## Mild_threshold_upper                585825              1046000    13.625
    ## Mild_threshold_lower               -224375              -274000     0.065
    ## Number_outliers                          5                    6     1.000
    ##                        POC_µM      d13C PON.µmol    d15N Fluorescence_mg.m3
    ## iqr                   3.82250   3.58250   0.6850  3.5150           0.093650
    ## lowerq                1.87250 -28.94000   0.3300 -0.5450           0.015000
    ## upperq                5.69500 -25.35750   1.0150  2.9700           0.108650
    ## Mild_threshold_upper 11.42875 -19.98375   2.0425  8.2425           0.249125
    ## Mild_threshold_lower -3.86125 -34.31375  -0.6975 -5.8175          -0.125475
    ## Number_outliers       5.00000   0.00000   1.0000  3.0000           6.000000
    ##                      PAR_µmol.m2.s    Bsi_µM NOx_µmol.L Silicic_acid_µmol.L
    ## iqr                        1.62690  1.276750      10.28             60.4850
    ## lowerq                     0.00000  0.174500      20.88              9.8550
    ## upperq                     1.62690  1.451250      31.16             70.3400
    ## Mild_threshold_upper       4.06725  3.366375      46.58            161.0675
    ## Mild_threshold_lower      -2.44035 -1.740625       5.46            -80.8725
    ## Number_outliers           25.00000 12.000000       0.00              0.0000
    ##                      Nitrite_µmol.L Phosphate_µmol.L Ammonium_µmol.L
    ## iqr                            0.22           0.6550          0.4650
    ## lowerq                         0.01           1.4150          0.0000
    ## upperq                         0.23           2.0700          0.4650
    ## Mild_threshold_upper           0.56           3.0525          1.1625
    ## Mild_threshold_lower          -0.32           0.4325         -0.6975
    ## Number_outliers                0.00           0.0000         11.0000
    ##                      Nitrate_µmol.L Density_kg.m2 dO2_µmol.kg Salinity_PSS78
    ## iqr                          10.410        1.6084     84.2650       0.771525
    ## lowerq                       20.745     1026.8637    225.5260      33.785425
    ## upperq                       31.155     1028.4720    309.7910      34.556950
    ## Mild_threshold_upper         46.770     1030.8846    436.1885      35.714237
    ## Mild_threshold_lower          5.130     1024.4511     99.1285      32.628138
    ## Number_outliers               0.000       23.0000      0.0000       3.000000
    ##                      Temperature_ITS90
    ## iqr                             3.6940
    ## lowerq                         -0.0874
    ## upperq                          3.6066
    ## Mild_threshold_upper            9.1476
    ## Mild_threshold_lower           -5.6284
    ## Number_outliers                14.0000

``` r
(outliers.stats.bag <- computeOutliers(table = meta.sum.bagImp))
```

    ## Warning in base::rbind(...): number of columns of result is not a multiple of
    ## vector length (arg 2)
    ## Warning in base::rbind(...): number of columns of result, 24, is not a multiple
    ## of vector length 27 of arg 2

    ##                      dCd_nmol.kg dFe_nmol.kg dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## iqr                    0.4924882  0.18122400    4.101088             128242.99
    ## lowerq                 0.2712618  0.04319188    0.675514             110756.42
    ## upperq                 0.7637500  0.22441588    4.776602             238999.41
    ## Mild_threshold_upper   1.5024823  0.49625188   10.928235             431363.90
    ## Mild_threshold_lower  -0.4674705 -0.22864412   -5.476119             -81608.06
    ## Number_outliers        0.0000000  4.00000000    0.000000                 12.00
    ##                      Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL TPZT_µM.C
    ## iqr                               142525.0             295374.3    2.0650
    ## lowerq                             79225.0             173375.7    5.4300
    ## upperq                            221750.0             468750.0    7.4950
    ## Mild_threshold_upper              435537.5             911811.4   10.5925
    ## Mild_threshold_lower             -134562.5            -269685.7    2.3325
    ## Number_outliers                        9.0                 10.0   13.0000
    ##                         POC_µM       d13C   PON.µmol      d15N
    ## iqr                   3.602082   1.971291  0.6047954  2.749711
    ## lowerq                1.590418 -28.477833  0.3052046  0.417500
    ## upperq                5.192500 -26.506543  0.9100000  3.167211
    ## Mild_threshold_upper 10.595623 -23.549606  1.8171931  7.291777
    ## Mild_threshold_lower -3.812704 -31.434770 -0.6019885 -3.707066
    ## Number_outliers       9.000000   9.000000  5.0000000 10.000000
    ##                      Fluorescence_mg.m3 PAR_µmol.m2.s    Bsi_µM NOx_µmol.L
    ## iqr                          0.07044931      7.054077  1.886718     9.4050
    ## lowerq                       0.04316753      0.000000  0.193000    21.6425
    ## upperq                       0.11361684      7.054077  2.079718    31.0475
    ## Mild_threshold_upper         0.21929080     17.635193  4.909795    45.1550
    ## Mild_threshold_lower        -0.06250643    -10.581116 -2.637077     7.5350
    ## Number_outliers              8.00000000     14.000000 17.000000     3.0000
    ##                      Silicic_acid_µmol.L Nitrite_µmol.L Phosphate_µmol.L
    ## iqr                              54.8100      0.1967084           0.5700
    ## lowerq                           13.2975      0.0200000           1.4525
    ## upperq                           68.1075      0.2167084           2.0225
    ## Mild_threshold_upper            150.3225      0.5117711           2.8775
    ## Mild_threshold_lower            -68.9175     -0.2750626           0.5975
    ## Number_outliers                   0.0000      0.0000000           3.0000
    ##                      Ammonium_µmol.L Nitrate_µmol.L Density_kg.m2 dO2_µmol.kg
    ## iqr                             0.44        9.57750      1.555075    75.08825
    ## lowerq                          0.00       21.43500   1026.863525   233.89925
    ## upperq                          0.44       31.01250   1028.418600   308.98750
    ## Mild_threshold_upper            1.10       45.37875   1030.751212   421.61988
    ## Mild_threshold_lower           -0.66        7.06875   1024.530912   121.26687
    ## Number_outliers                12.00        3.00000     23.000000     0.00000
    ##                      Salinity_PSS78 Temperature_ITS90
    ## iqr                        0.765225          4.050750
    ## lowerq                    33.780025         -0.012625
    ## upperq                    34.545250          4.038125
    ## Mild_threshold_upper      35.693087         10.114250
    ## Mild_threshold_lower      32.632188         -6.088750
    ## Number_outliers            3.000000         13.000000

``` r
(outliers.stats.knn <- computeOutliers(table = meta.sum.knnImp))
```

    ## Warning in base::rbind(...): number of columns of result is not a multiple of
    ## vector length (arg 2)
    ## Warning in base::rbind(...): number of columns of result, 24, is not a multiple
    ## of vector length 27 of arg 2

    ##                      dCd_nmol.kg dFe_nmol.kg dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## iqr                    1.7044620   1.0534850   1.7216336             1.0101904
    ## lowerq                -0.8899388  -0.8474908  -1.0261807            -0.6856167
    ## upperq                 0.8145232   0.2059942   0.6954529             0.3245737
    ## Mild_threshold_upper   3.3712162   1.7862217   3.2779033             1.8398593
    ## Mild_threshold_lower  -3.4466317  -2.4277182  -3.6086311            -2.2009023
    ## Number_outliers        0.0000000   3.0000000   0.0000000             3.0000000
    ##                      Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL  TPZT_µM.C
    ## iqr                              0.8491668           0.81597758  0.8914606
    ## lowerq                          -0.6979037          -0.71625752 -0.6086969
    ## upperq                           0.1512631           0.09972006  0.2827636
    ## Mild_threshold_upper             1.4250134           1.32368644  1.6199545
    ## Mild_threshold_lower            -1.9716540          -1.94022390 -1.9458878
    ## Number_outliers                  9.0000000          10.00000000  5.0000000
    ##                          POC_µM      d13C   PON.µmol       d15N
    ## iqr                   1.0202845  0.731096  1.1248432  0.9095984
    ## lowerq               -0.9082710 -0.490932 -0.8963576 -0.2926236
    ## upperq                0.1120135  0.240164  0.2284856  0.6169747
    ## Mild_threshold_upper  1.6424403  1.336808  1.9157504  1.9813723
    ## Mild_threshold_lower -2.4386978 -1.587576 -2.5836223 -1.6570212
    ## Number_outliers      11.0000000 17.000000  5.0000000  8.0000000
    ##                      Fluorescence_mg.m3 PAR_µmol.m2.s     Bsi_µM NOx_µmol.L
    ## iqr                           0.5695216   0.379826108  0.7028723  1.4358592
    ## lowerq                       -0.3762674  -0.382610733 -0.5234227 -0.5733032
    ## upperq                        0.1932542  -0.002784625  0.1794496  0.8625560
    ## Mild_threshold_upper          1.0475367   0.566954537  1.2337580  3.0163448
    ## Mild_threshold_lower         -1.2305499  -0.952349895 -1.5777311 -2.7270920
    ## Number_outliers               9.0000000  14.000000000 21.0000000  0.0000000
    ##                      Silicic_acid_µmol.L Nitrite_µmol.L Phosphate_µmol.L
    ## iqr                            1.7046330      1.8464956        1.3719042
    ## lowerq                        -0.8896978     -1.0323441       -0.5650570
    ## upperq                         0.8149353      0.8141515        0.8068472
    ## Mild_threshold_upper           3.3718848      3.5838949        2.8647036
    ## Mild_threshold_lower          -3.4466473     -3.8020875       -2.6229133
    ## Number_outliers                0.0000000      0.0000000        0.0000000
    ##                      Ammonium_µmol.L Nitrate_µmol.L Density_kg.m2 dO2_µmol.kg
    ## iqr                       0.49449249      1.4404638    0.50388698   1.2697541
    ## lowerq                   -0.43985897     -0.5736877   -0.53846581  -0.5916666
    ## upperq                    0.05463352      0.8667761   -0.03457883   0.6780875
    ## Mild_threshold_upper      0.79637226      3.0274717    0.72125165   2.5827186
    ## Mild_threshold_lower     -1.18159770     -2.7343833   -1.29429628  -2.4962977
    ## Number_outliers          11.00000000      0.0000000   24.00000000   0.0000000
    ##                      Salinity_PSS78 Temperature_ITS90
    ## iqr                       1.3799144         1.0680063
    ## lowerq                   -0.5451160        -0.6887507
    ## upperq                    0.8347984         0.3792556
    ## Mild_threshold_upper      2.9046701         1.9812649
    ## Mild_threshold_lower     -2.6149877        -2.2907601
    ## Number_outliers           3.0000000        13.0000000

``` r
(outliers.stats.mean <- computeOutliers(table = meta.sum.mean))
```

    ## Warning in base::rbind(...): number of columns of result is not a multiple of
    ## vector length (arg 2)
    ## Warning in base::rbind(...): number of columns of result, 24, is not a multiple
    ## of vector length 27 of arg 2

    ##                      dCd_nmol.kg dFe_nmol.kg dZn_nmol.kg Bacteria_HDNA.cell.mL
    ## iqr                      0.09650  0.09781523      1.0810              84412.87
    ## lowerq                   0.53100  0.08962500      2.7945             151500.00
    ## upperq                   0.62750  0.18744023      3.8755             235912.87
    ## Mild_threshold_upper     0.77225  0.33416307      5.4970             362532.18
    ## Mild_threshold_lower     0.38625 -0.05709784      1.1730              24880.69
    ## Number_outliers         54.00000 21.00000000     45.0000                 21.00
    ##                      Bacteria_LDNA.cell.mL Tot_bacteria_cell.mL TPZT_µM.C
    ## iqr                               94248.39            142911.95    0.7750
    ## lowerq                           108000.00            295250.00    6.6700
    ## upperq                           202248.39            438161.95    7.4450
    ## Mild_threshold_upper             343620.98            652529.89    8.6075
    ## Mild_threshold_lower             -33372.59             80882.07    5.5075
    ## Number_outliers                      15.00                25.00   48.0000
    ##                         POC_µM      d13C PON.µmol     d15N Fluorescence_mg.m3
    ## iqr                   1.355372   1.20250     0.24  1.32250          0.0372250
    ## lowerq                3.242500 -28.27250     0.53  0.41750          0.0683500
    ## upperq                4.597872 -27.07000     0.77  1.74000          0.1055750
    ## Mild_threshold_upper  6.630931 -25.26625     1.13  3.72375          0.1614125
    ## Mild_threshold_lower  1.209441 -30.07625     0.17 -1.56625          0.0125125
    ## Number_outliers      37.000000  32.00000    31.00 32.00000         36.0000000
    ##                      PAR_µmol.m2.s    Bsi_µM NOx_µmol.L Silicic_acid_µmol.L
    ## iqr                       8.370666  1.152633     9.4050             54.8100
    ## lowerq                    0.000000  0.193000    21.6425             13.2975
    ## upperq                    8.370666  1.345633    31.0475             68.1075
    ## Mild_threshold_upper     20.926666  3.074582    45.1550            150.3225
    ## Mild_threshold_lower    -12.556000 -1.535949     7.5350            -68.9175
    ## Number_outliers          14.000000 16.000000     3.0000              0.0000
    ##                      Nitrite_µmol.L Phosphate_µmol.L Ammonium_µmol.L
    ## iqr                            0.18           0.5700       0.4350741
    ## lowerq                         0.02           1.4525       0.0000000
    ## upperq                         0.20           2.0225       0.4350741
    ## Mild_threshold_upper           0.47           2.8775       1.0876852
    ## Mild_threshold_lower          -0.25           0.5975      -0.6526111
    ## Number_outliers                0.00           3.0000      13.0000000
    ##                      Nitrate_µmol.L Density_kg.m2 dO2_µmol.kg Salinity_PSS78
    ## iqr                         9.57750      1.623654     74.0750       0.755425
    ## lowerq                     21.43500   1026.899450    233.8992      33.789825
    ## upperq                     31.01250   1028.523104    307.9742      34.545250
    ## Mild_threshold_upper       45.37875   1030.958586    419.0867      35.678387
    ## Mild_threshold_lower        7.06875   1024.463969    122.7868      32.656687
    ## Number_outliers             3.00000     23.000000      0.0000       3.000000
    ##                      Temperature_ITS90
    ## iqr                           3.595050
    ## lowerq                       -0.012625
    ## upperq                        3.582425
    ## Mild_threshold_upper          8.975000
    ## Mild_threshold_lower         -5.405200
    ## Number_outliers              14.000000

``` r
diff.mean <- as.vector(as.data.frame(t(outliers.stats.mean))$Number_outliers - as.data.frame(t(outliers.stats.raw))$Number_outliers)
diff <- data.frame(as.vector(diff.mean))
rownames(diff) <- rownames(as.data.frame(t(outliers.stats.mean)))

diff$diff.knn <- as.vector(as.data.frame(t(outliers.stats.knn))$Number_outliers - as.data.frame(t(outliers.stats.raw))$Number_outliers)
diff$diff.bag <- as.vector(as.data.frame(t(outliers.stats.bag))$Number_outliers - as.data.frame(t(outliers.stats.raw))$Number_outliers)

#diff <- data.frame(Mean = diff.mean, knn = diff.knn, Bag = diff.bag)

colnames(diff) <- c("Mean","knn","Bag")
diff <- diff %>% 
  rownames_to_column(var = "Variable") %>% 
  pivot_longer(cols = !c("Variable"), names_to = "Method", values_to = "Diff_number_outliers")
diff
```

    ## # A tibble: 72 × 3
    ##    Variable              Method Diff_number_outliers
    ##    <chr>                 <chr>                 <dbl>
    ##  1 dCd_nmol.kg           Mean                     54
    ##  2 dCd_nmol.kg           knn                       0
    ##  3 dCd_nmol.kg           Bag                       0
    ##  4 dFe_nmol.kg           Mean                     20
    ##  5 dFe_nmol.kg           knn                       2
    ##  6 dFe_nmol.kg           Bag                       3
    ##  7 dZn_nmol.kg           Mean                     45
    ##  8 dZn_nmol.kg           knn                       0
    ##  9 dZn_nmol.kg           Bag                       0
    ## 10 Bacteria_HDNA.cell.mL Mean                     18
    ## # ℹ 62 more rows

``` r
ggplot(data = diff %>% mutate(Variable = factor(Variable, levels=(miss_var_summary(meta.sum[,var]) %>% pull(variable)))),
       aes(x = Variable, y = Diff_number_outliers, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#00A087CC", "#4DBBD5CC", "#E64B35CC")) +
  labs(y = "Diff number of outliers vs. raw") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.title.x = element_blank())
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Diff%20number%20outliers-1.png)<!-- -->

## Inspect effect of imputation in the TS space

``` r
plot.TS.all <- readRDS("./R_Data/TS_plot_all.rds")
```

``` r
plotTS <- function(metadata, col = meta.group$Water_mass){
  plot.TS.all +
    geom_point(data = metadata,
               aes(x = Salinity_PSS78, y = Temperature_ITS90,
                   col = col),
               size = 2.5) +
    scale_colour_manual(values = wm.pal) +
    labs(color = "Water mass",
         x = "Practical Salinity (psu)",
         y = "Conservative Temperature (°C)") +
    theme_bw()
}
  
plotTS(meta.sum)
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-1.png)<!-- -->

``` r
set.seed(13)
meta.sum.bag <- preProcess(meta.sum.preprocess, method = c("bagImpute")) # Centered and scaled later on
meta.sum.bagImp <- predict(meta.sum.bag, meta.sum.preprocess)
plotTS(meta.sum.bagImp)
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-2.png)<!-- -->

``` r
set.seed(37)
meta.sum.bag <- preProcess(meta.sum.preprocess, method = c("bagImpute")) # Centered and scaled later on
meta.sum.bagImp <- predict(meta.sum.bag, meta.sum.preprocess)
plotTS(meta.sum.bagImp)
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-3.png)<!-- -->

``` r
set.seed(1998)
meta.sum.bag <- preProcess(meta.sum.preprocess, method = c("bagImpute")) # Centered and scaled later on
meta.sum.bagImp <- predict(meta.sum.bag, meta.sum.preprocess)
plotTS(meta.sum.bagImp)
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-4.png)<!-- -->

``` r
ggplot() +
  geom_point(data = meta.sum.knnImp,
             aes(x = Salinity_PSS78, y = Temperature_ITS90,
                 col = meta.group$Water_mass),
             size = 3) +
  scale_colour_manual(values = wm.pal) +
  labs(color = "Water mass",
       x = "Practical Salinity (scaled)",
       y = "Temperature (scaled)") +
  theme_bw()  +
  labs()
```

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-5.png)<!-- -->

``` r
(plot.TS.wm.simplified <- plotTS(meta.sum, col = meta.group$Water_mass_simplified) +
  scale_colour_manual(values = wm.pal.sim) +
  #labs(title = "Temperature-Salinity diagram defining water masses for ACE CTD samples") +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) + 
    guides(col = guide_legend(title.position = "top", title.hjust = 0.5)))
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).
    ## Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-6.png)<!-- -->

``` r
(plot.TS.wm.simplified <- plotTS(meta.sum, col = meta.group$Water_mass_simplified) +
  scale_colour_manual(values = wm.pal.sim) +
  #labs(title = "Temperature-Salinity diagram defining water masses for ACE CTD samples") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)))
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).
    ## Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](02-ACE_CTD_Metadata_files/figure-gfm/TS%20plots-7.png)<!-- -->

**The bagging tree imputation sometimes introduces “outliers” in the T-S
plot (not necessarily in the temperature or salinity range but in the
association of the 2 variables). Therefore, we will use the *knn*
imputed dataset for further analyses. **

# Subset to metagenome samples

``` r
meta.metaG <- meta[meta$ACE_seq_name %in% metagenomes$ace_seq_name,]
dim(meta.metaG)
```

    ## [1] 206  93

``` r
meta.metaG.sum <- meta.metaG %>% 
  select(-c("CTD_cast_number","TM_station_number","DNA_vol_µl","DNA_Conc_ng.µl","DNA_qty_ng","RNA_vol_µl","RNA_Conc_ng.µl","RNA_qty_ng","Size_fraction","Event_start_date","Event_end_date")) %>% 
  group_by(Event_number, Depth_m,Depth_q, Longhurst_Prov) %>% 
  summarise_if(.predicate = is.numeric, .funs = mean) %>%
  mutate(Sample_name = paste(Event_number, Depth_m, sep = "_")) %>%
  filter(Sample_name %in% rownames(meta.sum)) %>% 
  column_to_rownames(var = "Sample_name")
dim(meta.metaG.sum)
```

    ## [1] 119  78

``` r
meta.metaG.bagImp <- meta.sum.bagImp[rownames(meta.sum.bagImp) %in% rownames(meta.metaG.sum),] # subset the bag imputed table
meta.metaG.knnImp <- meta.sum.knnImp[rownames(meta.sum.knnImp) %in% rownames(meta.metaG.sum),] # subset the bag imputed table
```

After grouping and summarizing, the subset of the table is done on the
imputed table so that the imputation was computed with all possible
samples and not on a reduced dataset.

``` r
meta.sum.knnImp.all <- meta.sum.knnImp

###### Run if only samples with metagenomes wanted ########
# meta.sum.knnImp <- meta.metaG.knnImp # don't run if running analyses on all samples
###########################################################

dim(meta.group) # needs to be subset accordingly if further using the subset of metagenomic samples
```

    ## [1] 150  10

``` r
dim(meta.sum.knnImp)
```

    ## [1] 150  27

``` r
meta.group.all <- meta.group

###### Run if only samples with metagenomes wanted ########
# meta.group <- meta.group[rownames(meta.group) %in% rownames(meta.sum.knnImp),]
###########################################################

dim(meta.group)
```

    ## [1] 150  10

# PCA

## All variables and all samples

``` r
plotNormalHistogram(meta.sum.knnImp) # centered and scaled while pre processing
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Transformation%20and%20NA-1.png)<!-- -->

``` r
meta.pca <- rda(meta.sum.knnImp, scale=FALSE) # knn is already centered and scaled
# Note that NaN values don't fail, why?
summary(meta.pca)$cont
```

    ## $importance
    ## Importance of components:
    ##                          PC1    PC2     PC3     PC4     PC5     PC6     PC7
    ## Eigenvalue            9.6018 4.0804 2.36453 1.67873 1.10472 0.81622 0.67532
    ## Proportion Explained  0.4029 0.1712 0.09921 0.07043 0.04635 0.03425 0.02833
    ## Cumulative Proportion 0.4029 0.5741 0.67327 0.74370 0.79005 0.82430 0.85263
    ##                           PC8     PC9    PC10    PC11    PC12    PC13     PC14
    ## Eigenvalue            0.64805 0.56439 0.44209 0.38404 0.31922 0.24213 0.194348
    ## Proportion Explained  0.02719 0.02368 0.01855 0.01611 0.01339 0.01016 0.008154
    ## Cumulative Proportion 0.87982 0.90350 0.92205 0.93817 0.95156 0.96172 0.969872
    ##                           PC15     PC16     PC17     PC18     PC19     PC20
    ## Eigenvalue            0.189256 0.160328 0.103304 0.080151 0.066589 0.037152
    ## Proportion Explained  0.007941 0.006727 0.004334 0.003363 0.002794 0.001559
    ## Cumulative Proportion 0.977812 0.984539 0.988873 0.992236 0.995030 0.996589
    ##                           PC21      PC22      PC23      PC24      PC25
    ## Eigenvalue            0.024429 0.0233605 0.0201995 0.0129527 3.540e-04
    ## Proportion Explained  0.001025 0.0009801 0.0008475 0.0005435 1.485e-05
    ## Cumulative Proportion 0.997614 0.9985940 0.9994415 0.9999850 1.000e+00
    ##                            PC26      PC27
    ## Eigenvalue            2.136e-06 1.512e-06
    ## Proportion Explained  8.961e-08 6.343e-08
    ## Cumulative Proportion 1.000e+00 1.000e+00

``` r
screeplot(meta.pca, bstick = TRUE, main = "Brocken stick - ACE CTD Metadata")
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Brocken%20stick-1.png)<!-- -->

According to the brocken stick model, we should then have a look at the
two first axes when considering only samples where we have metagenomes,
and at the third axis when considering all available CTD samples. With
all CTD samples:  
\* PC1 represents 40% of the variance \* PC2 represents 17% \* PC3
represents 10%  
With CTD samples where we have a metagenome:  
\* PC1 represents 37% \* PC2 represents 19%

``` r
renameVariables <- function(x) {
  # x vector of variable names
  if (length(x == 26)) {
    # Without factors
    df <- x %>% t() %>% as.data.frame()
    colnames(df) <- x
    df <- df %>%
      rename(
        #"LV6" = "LV5_mean",
        #"LV8" = "LV7_mean",
        "dCd" = "dCd_nmol.kg",
        "dFe" = "dFe_nmol.kg",
        "dZn" = "dZn_nmol.kg",
        "HDNA Bacteria" = "Bacteria_HDNA.cell.mL",
        "LDNA Bacteria" = "Bacteria_LDNA.cell.mL",
        "TPZT" = "TPZT_µM.C",
        "POC" = "POC_µM",
        "PON" = "PON.µmol",
        "Fluorescence" = "Fluorescence_mg.m3",
        "PAR" = "PAR_µmol.m2.s",
        "Biological silica" = "Bsi_µM",
        "Silicic acid" = "Silicic_acid_µmol.L",
        "Nitrite" = "Nitrite_µmol.L",
        "Nitrate" = "Nitrate_µmol.L",
        "Phosphate" = "Phosphate_µmol.L",
        "Ammonium" = "Ammonium_µmol.L",
        "Dissolved O2" = "dO2_µmol.kg",
        "Salinity" = "Salinity_PSS78",
        "Temperature" = "Temperature_ITS90",
        "Depth" = "Depth_m",
        "Latitude" = "Latitude",
        "Longitude" = "Longitude",
        "NO2+NO3" = "NOx_µmol.L",
        "Density" = "Density_kg.m2",
        "Total bacteria" = "Tot_bacteria_cell.mL",
        #"Distance to land" = "Distance_land"
      )
    renamed_x <- colnames(df)
    renamed_x
  }
}
```

``` r
(plot.pca.env <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Depth_q_new,
      shape.vec = meta.group$Longhurst_Prov,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = FALSE,
      scaling = 1) +
  scale_color_npg() +
  #scale_colour_viridis(discrete = TRUE) +
  labs(col = "Depth", shape = "Longhurst province"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-1.png)<!-- -->

``` r
(plot.pca.env.prov <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Longhurst_Prov,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = TRUE,
      scaling = 1,
      spe.names = renameVariables) +
  scale_color_npg() +
  #scale_colour_viridis(discrete = TRUE) +
  labs(col = "Longhurst province"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-2.png)<!-- -->

``` r
(plot.pca.env.depth <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Depth_q_new,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = TRUE,
      scaling = 1,
      spe.names = renameVariables) +
  scale_color_npg() +
  #scale_colour_viridis(discrete = TRUE) +
  labs(col = "Depth"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-3.png)<!-- -->

``` r
(plot.pca.env.wm <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = TRUE,
      scaling = 1,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal) +
  labs(col = "Water mass"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-4.png)<!-- -->

``` r
(plot.pca.env.wm <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass_simplified,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = FALSE,
      scaling = 1,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass")) +
  theme(legend.position = c(0.70,0.875),
        plot.title = element_blank()) +
  guides(col = guide_legend(ncol = 2))
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-5.png)<!-- -->

``` r
(plot.pca.env.wm <- ggPCA(meta.pca, 
      ax1 = 1, ax2 = 3,
      color.vec = meta.group$Water_mass_simplified,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      draw.circle = TRUE,
      scaling = 1,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass"))
```

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](02-ACE_CTD_Metadata_files/figure-gfm/Plot%20PCA-6.png)<!-- -->

The clustering according to Longhurst province is not obvious, mainly we
see ANTA and APLR clustering together but showing a depth gradient.
Moreover, some samples belonging to the ANTA province from the surface
or below 150 m mix with samples in the SSTC and SANT provinces. SSTC and
SANT provinces samples mix, and some of them are closer to the ANTA and
APLR samples, because they are deep samples.

Similarly when looking at the distribution of samples according to
depth, the discrimination is not obvious, especially for samples at the
surface and up to 150 m.

Water masses, however, seem to play a decisive role in determining where
the sample falls in the PCA space.

``` r
meta.pca.good <- goodness(meta.pca, model = "CA", choices = 1:3)

ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass_simplified,
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      select.spe = meta.pca.good[,2]>0.40,
      draw.circle = TRUE,
      scaling = 1) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass") +
  theme(legend.box.background = element_blank())
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Select%20variables-1.png)<!-- -->

``` r
ggPCA(meta.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Depth_q_new,
      shape.vec = meta.group$Longhurst_Prov,
      label.sites = FALSE,
      select.spe = meta.pca.good[,2]>0.40,
      scaling = 2) +
  scale_color_npg() +
  #scale_colour_viridis(discrete = TRUE) +
  labs(col = "Longhurst", shape = "Depth")
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Select%20variables-2.png)<!-- -->

``` r
ggPCA(meta.pca, 
      ax1 = 1, ax2 = 3,
      color.vec = meta.group$Depth_q_new,
      shape.vec = meta.group$Longhurst_Prov,
      label.sites = TRUE,
      select.spe = meta.pca.good[,3]>0.60,
      scaling = 2) +
  scale_color_npg() +
  #scale_colour_viridis(discrete = TRUE) +
  labs(col = "Longhurst", shape = "Depth") # PC3 driven by depth / density
```

    ## Warning: ggrepel: 142 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](02-ACE_CTD_Metadata_files/figure-gfm/Select%20variables-3.png)<!-- -->

## ANTA and APLR

``` r
ant <- which(meta.group$Longhurst_Prov == "ANTA" | meta.group$Longhurst_Prov == "APLR")
meta.ant <- meta.sum.knnImp[ant,]
dim(meta.ant)
```

    ## [1] 103  27

``` r
#meta.ant$log.depth <- log(meta.ant$sampling.event..elevation..depth.below.sea.surface)

meta.ant.pca <- rda(meta.ant, scale = FALSE)
screeplot(meta.ant.pca, bstick = TRUE, main = "Brocken stick - ACE CTD Metadata ANTA-APLR")
```

![](02-ACE_CTD_Metadata_files/figure-gfm/Subset%20ANTA%20APLR-1.png)<!-- -->

``` r
(plot.pca.ant <- ggPCA(meta.ant.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass_simplified[ant], #log(meta.ant$Depth_m),
      shape.vec = meta.group$Longhurst_Prov[ant],
      label.sites = TRUE,
      label.spe = TRUE,
      plot.spe = TRUE,
      scaling = 1,
      draw.circle = TRUE,
      spe.names = renameVariables) +
  scale_color_manual(values = wm.pal.sim) +
  #scale_colour_viridis(discrete = TRUE, na.value = "grey") +
  labs(col = "Water mass", shape = "Longhurst"))
```

    ## Warning: ggrepel: 92 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](02-ACE_CTD_Metadata_files/figure-gfm/Subset%20ANTA%20APLR-2.png)<!-- -->

``` r
(plot.pca.ant <- ggPCA(meta.ant.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = log(meta.group$Depth_m[ant]), #log(meta.ant$Depth_m),
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      scaling = 1,
      draw.circle = TRUE,
      spe.names = renameVariables) +
  #scale_color_manual(values = wm.pal.sim) +
  scale_colour_viridis(discrete = FALSE) +
  labs(col = "log(Depth)"))
```

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](02-ACE_CTD_Metadata_files/figure-gfm/Subset%20ANTA%20APLR-3.png)<!-- -->

## Surface samples

``` r
surf <- meta.group$Water_mass %in% c("AASW","SASW","STSW","AASW-WW","WW","AASW-AAIW","SASW-AAIW","STSW-AAIW")
meta.srf <- meta.sum.knnImp[surf,]
dim(meta.srf)
```

    ## [1] 101  27

``` r
meta.srf.pca <- rda(meta.srf, scale = FALSE)
screeplot(meta.srf.pca, bstick = TRUE, main = "Brocken stick - ACE CTD Metadata SW")
```

![](02-ACE_CTD_Metadata_files/figure-gfm/PCA%20Subset%20surface-1.png)<!-- -->

``` r
(plot.pca.srf <- ggPCA(meta.srf.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass[surf], #log(meta.ant$Depth_m),
      shape.vec = meta.group$Longhurst_Prov[surf],
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      scaling = 1,
      draw.circle = TRUE,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass", shape = "Longhurst"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/PCA%20Subset%20surface-2.png)<!-- -->

``` r
(plot.pca.srf <- ggPCA(meta.srf.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass_simplified[surf], #log(meta.ant$Depth_m),
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      scaling = 1,
      draw.circle = TRUE,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/PCA%20Subset%20surface-3.png)<!-- -->

``` r
(plot.pca.srf <- ggPCA(meta.srf.pca, 
      ax1 = 1, ax2 = 2,
      color.vec = meta.group$Water_mass[surf], #log(meta.ant$Depth_m),
      label.sites = FALSE,
      label.spe = TRUE,
      plot.spe = TRUE,
      scaling = 1,
      draw.circle = TRUE,
      spe.names = renameVariables) +
  scale_colour_manual(values = wm.pal.sim) +
  labs(col = "Water mass"))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/PCA%20Subset%20surface-4.png)<!-- -->

## With physical variables

# Match Latent Variables from Landwehr *et al.* (2021) with CTD samples

``` r
landwehr <- read.table("./Metadata/landwehr_latent_variables_mean_3h.csv", header = T, sep = ",")
head(landwehr)
```

    ##               timest_ LV0_mean  LV1_mean  LV2_mean LV3_mean    LV4_mean
    ## 1 2016-12-20 12:00:00 2.702255 0.2606801 -2.680330 3.755368  0.08723822
    ## 2 2016-12-20 15:00:00 2.798098 0.1492510 -2.428727 3.512605  0.07391819
    ## 3 2016-12-20 18:00:00 2.776233 0.1065375 -2.116134 3.494059  0.02079603
    ## 4 2016-12-20 21:00:00 3.235376 0.1305450 -0.998531 2.925723  0.75475145
    ## 5 2016-12-21 00:00:00 3.146485 1.8816363 -1.173363 1.424038 -0.15521118
    ## 6 2016-12-21 03:00:00 3.396193 2.2232706 -1.479160 1.403239  0.43120518
    ##      LV5_mean   LV6_mean    LV7_mean    LV8_mean   LV9_mean    LV10_mean
    ## 1 0.024078945 0.07026055 -0.03560512 -0.76792929 -2.1595833 -0.137428508
    ## 2 0.044210783 0.42671156  0.05610340 -0.75362937 -1.1890112 -0.002352388
    ## 3 0.019839662 0.59956406  0.05485859 -0.73850095  0.2519914 -0.027763009
    ## 4 0.008491765 0.58154201 -0.08022960 -0.59504812  0.3998022 -0.153841816
    ## 5 0.104737352 2.21049100  0.84218696 -0.03629023  0.4606314 -0.591503076
    ## 6 0.022630512 2.12106845  0.07866026 -0.65453774  0.4478866 -0.258207830
    ##     LV11_mean LV12_mean LV13_mean
    ## 1  0.39264345 0.4985240 0.6383805
    ## 2 -0.00324204 0.2574843 0.6397556
    ## 3  0.06702592 0.9301827 0.6507352
    ## 4  0.61712482 1.8235005 0.6762153
    ## 5  0.57704085 1.4396670 1.0564681
    ## 6  0.92710069 1.9284040 0.7478357

``` r
landwehr <- landwehr %>% 
  mutate(Date_time = ymd_hms(timest_)) %>%
  select(-"timest_") %>% 
  mutate(Date = date(Date_time))
```

``` r
ctd_timest <- meta %>% 
  select(c("Event_number","Event_end_date")) %>% 
  mutate(Date_time = ymd_hms(Event_end_date)) %>%
  group_by(Event_number) %>% 
  summarise(Date_time = mean(Date_time)) %>% # Time stamp is the same for a same event
  mutate(Sampling_type = as.factor("CTD")) # %>%  
  #mutate(Date = date(Date_time)) # %>% 
  # mutate(Date_time_new = round_date(Date_time, unit = "hour")))

head(ctd_timest)
```

    ## # A tibble: 6 × 3
    ##   Event_number Date_time           Sampling_type
    ##          <int> <dttm>              <fct>        
    ## 1           75 2016-12-28 08:00:00 CTD          
    ## 2          123 2016-12-30 14:02:00 CTD          
    ## 3          264 2017-01-07 06:12:00 CTD          
    ## 4          317 2017-01-09 10:59:00 CTD          
    ## 5          318 2017-01-09 13:31:00 CTD          
    ## 6          369 2017-01-11 07:24:00 CTD

``` r
meta_timest <- meta.udw %>%
  filter(ACE_seq_name %in% metagenomes$ace_seq_name) %>%
  select(c("Event_number","Event_end_date")) %>% 
  mutate(Date_time = ymd_hms(Event_end_date)) %>%
  group_by(Event_number) %>% 
  summarise(Date_time = mean(Date_time)) %>% 
  mutate(Sampling_type = as.factor("UDW")) %>%
  bind_rows(ctd_timest)
```

``` r
matched_timest <- meta_timest %>%
  full_join(landwehr, by = character()) %>%
  mutate(diff = abs((Date_time.x%--%Date_time.y) / hours(1))) %>%
  group_by(Date_time.y) %>%
  filter(diff <= 1.5) %>%
  ungroup() %>%
  select(!diff)
```

    ## Warning: Using `by = character()` to perform a cross join was deprecated in dplyr 1.1.0.
    ## ℹ Please use `cross_join()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
all.equal(matched_timest$Event_number, meta_timest$Event_number)
```

    ## [1] TRUE

``` r
colnames(matched_timest)
```

    ##  [1] "Event_number"  "Date_time.x"   "Sampling_type" "LV0_mean"     
    ##  [5] "LV1_mean"      "LV2_mean"      "LV3_mean"      "LV4_mean"     
    ##  [9] "LV5_mean"      "LV6_mean"      "LV7_mean"      "LV8_mean"     
    ## [13] "LV9_mean"      "LV10_mean"     "LV11_mean"     "LV12_mean"    
    ## [17] "LV13_mean"     "Date_time.y"   "Date"

``` r
matched_timest <- matched_timest %>% 
  rename(Date_time_meta = Date_time.x,
         Date_time_lv = Date_time.y)

head(matched_timest)[,c(1:4,17:18)]
```

    ## # A tibble: 6 × 6
    ##   Event_number Date_time_meta      Sampling_type LV0_mean LV13_mean
    ##          <int> <dttm>              <fct>            <dbl>     <dbl>
    ## 1         1104 2017-01-31 07:03:00 UDW              -2.19   -2.65  
    ## 2         1349 2017-02-09 07:00:00 UDW              -2.24   -1.08  
    ## 3         1350 2017-02-09 10:00:00 UDW              -1.39    0.0532
    ## 4         1352 2017-02-09 16:00:00 UDW              -1.29   -1.01  
    ## 5         1354 2017-02-09 22:00:00 UDW              -1.67   -2.72  
    ## 6         1480 2017-02-11 15:00:00 UDW              -1.03   -0.111 
    ## # ℹ 1 more variable: Date_time_lv <dttm>

``` r
meta.sum.lv.ctd <- meta %>%
  select(c("Event_number","Event_end_date","Depth_m")) %>%
  mutate(Date_time = ymd_hms(Event_end_date)) %>%
  group_by(Event_number, Depth_m) %>%
  summarise(Date_time = mean(Date_time)) %>% # 54 Event numbers because there are duplicated Event numbers between CTD and UDW samples
  inner_join(matched_timest %>% filter(Sampling_type == "CTD")) %>% 
  ungroup() %>% 
  select(!c("Date_time","Date"))
```

    ## `summarise()` has grouped output by 'Event_number'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Event_number)`

``` r
meta.sum.lv <- meta.udw %>% 
  filter(ACE_seq_name %in% metagenomes$ace_seq_name) %>% 
  select(c("Event_number","Event_end_date","Depth_m")) %>% 
  mutate(Date_time = ymd_hms(Event_end_date)) %>%
  group_by(Event_number, Depth_m) %>%
  summarise(Date_time = mean(Date_time)) %>%
  inner_join(matched_timest %>% filter(Sampling_type == "UDW")) %>% 
  ungroup() %>% 
  select(!c("Date_time","Date")) %>%
  bind_rows(meta.sum.lv.ctd)
```

    ## `summarise()` has grouped output by 'Event_number'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Event_number)`

``` r
tail(meta.sum.lv[,c(1:4,10,12,19)])
```

    ## # A tibble: 6 × 7
    ##   Event_number Depth_m Date_time_meta      Sampling_type LV5_mean LV7_mean
    ##          <int>   <int> <dttm>              <fct>            <dbl>    <dbl>
    ## 1         3114     150 2017-03-16 07:25:41 CTD            -0.0503     1.48
    ## 2         3114    1000 2017-03-16 07:25:41 CTD            -0.0503     1.48
    ## 3         3117       5 2017-03-16 10:52:00 CTD             0.573      2.61
    ## 4         3117      15 2017-03-16 10:52:00 CTD             0.573      2.61
    ## 5         3117      30 2017-03-16 10:52:00 CTD             0.573      2.61
    ## 6         3117    1000 2017-03-16 10:52:00 CTD             0.573      2.61
    ## # ℹ 1 more variable: Date_time_lv <dttm>

``` r
saveRDS(object = meta.sum.lv, file = "./R_Data/matched_latent_variables.rds")
```

``` r
# Merge with metagenomes to discriminate size fractions

meta.sum.lv.metaG <- meta %>% 
  select(c("ACE_seq_name","Size_fraction","Event_number","Depth_m")) %>% 
  inner_join(metagenomes, by = c("ACE_seq_name" = "ace_seq_name")) %>% 
  select(-c("r1","r2","DNA.Genoscope.code")) %>% 
  rename(Sample = sample, Group = group) %>%
  inner_join(meta.sum.lv.ctd)
```

    ## Joining with `by = join_by(Event_number, Depth_m)`

``` r
meta.sum.lv.metaG <- meta.udw %>% 
  select(c("ACE_seq_name","Size_fraction","Event_number","Depth_m")) %>% 
  inner_join(metagenomes, by = c("ACE_seq_name" = "ace_seq_name")) %>% 
  select(-c("r1","r2","DNA.Genoscope.code")) %>% 
  rename(Sample = sample, Group = group) %>%
  inner_join(meta.sum.lv %>% filter(Sampling_type == "UDW")) %>% 
  bind_rows(meta.sum.lv.metaG) %>% 
  mutate(Metadata_sample_name = paste(Event_number, Depth_m, sep = "_")) %>% 
  filter(Sampling_type == "CTD" | Metadata_sample_name != "2970_4") # 2970_4 is commmon for CTD and UDW 
```

    ## Joining with `by = join_by(Event_number, Depth_m)`

``` r
saveRDS(object = meta.sum.lv.metaG, file = "./R_Data/matched_latent_variables_metagenomes.rds")
```

# Hierarchical clustering

``` r
meta.sum.ts <- meta.sum.knnImp %>% select(c("Temperature_ITS90", "Salinity_PSS78")) %>% na.omit()

dist.ts <- vegdist(meta.sum.ts, method = "euclidean")
hclust.ts <- hclust(dist.ts, method = "average")
# plot(hclust.ts, hang = -1)

# ggdendrogram(hclust.ts)

dendro <- as.dendrogram(hclust.ts)
dendro.data <- dendro_data(dendro)
dendro.labels <- dendro.data$labels %>% 
  mutate(Event_number = as.numeric(sapply(str_split(label, "_"), getElement, 1)),
         Depth = as.numeric(sapply(str_split(label, "_"), getElement, 2))) %>% 
  arrange(Event_number, Depth)

(plot.dendro <- ggplot(dendro.data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro.labels,
            aes(x = x, y = y, label = label, col = meta.group$Water_mass_simplified),
            hjust = 1, angle = 90, size = 3) +
  ylim(-0.5, 3) +
  #scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(values = wm.pal.sim) +
  #labs(color = "TS k-mean cluster") +
  labs(color = "Water mass") +
  theme_dendro() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank()))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/hclust%20TS-1.png)<!-- -->

``` r
dist.all <- vegdist(meta.sum.knnImp, method = "euclidean") # %>% scale() if bag impute
hclust.all <- hclust(dist.all, method = "average")

ggdendrogram(hclust.all)
```

![](02-ACE_CTD_Metadata_files/figure-gfm/hclust%20all%20variables-1.png)<!-- -->

``` r
dendro <- as.dendrogram(hclust.all)
dendro.data <- dendro_data(dendro)

dendro.labels <- dendro.data$labels %>% 
  mutate(Event_number = as.numeric(sapply(str_split(label, "_"), getElement, 1)),
         Depth = as.numeric(sapply(str_split(label, "_"), getElement, 2))) %>% 
  arrange(Event_number, Depth)

(plot.dendro <- ggplot(dendro.data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro.labels,
            aes(x = x, y = y, label = label, col = meta.group$Water_mass_simplified),
            hjust = 1, angle = 90, size = 3) +
  #ylim(-0.5, 3) +
  #scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(values = wm.pal.sim) +
  labs(color = "Water mass") +
  theme_dendro() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank()))
```

![](02-ACE_CTD_Metadata_files/figure-gfm/hclust%20all%20variables-2.png)<!-- -->

``` r
permutest(betadisper(dist(meta.sum.knnImp), group = meta.group$Water_mass_simplified))
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq Mean Sq     F N.Perm Pr(>F)    
    ## Groups      7  85.897 12.2710 13.08    999  0.001 ***
    ## Residuals 142 133.222  0.9382                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# PERMANOVA should not conducted if this test is significant, however PERMANOVA is not too sensitive to dispersion effects

# adonis2(meta.sum.knnImp[!is.na(meta.group$Water_mass),] ~ meta.group$Water_mass, method = "euc", permutations = how(nperm=999))
```

# Merged CTD and UDW metadata files for metagenomic samples

The aim here is to have a metadata file with 218 rows and describing all
metagenomic samples for downstream biological analyses. As CTD and UDW
sampling are not comparable with regards to the metadata variables,
hence we will only merge tables for the grouping variables.

``` r
# Manually check on a map for each event in 
# meta.udw[meta.udw$ACE_seq_name %in% metagenomes$ace_seq_name, ] %>% pull(Event_number)
# 1104 1104 1349 1349 1350 1350 1352 1352 1354 1354 1480 1878 2970 2970

#1104: TRUE.Surround?
#Other: FALSE

# For now just assign a "FALSE" to sample_34 which is in common with CTD data (cf. below)
```

``` r
metaG.CTD <- meta[meta$ACE_seq_name %in% metagenomes$ace_seq_name, ] # 206 samples x 93 variables

metaG.CTD.id <- metaG.CTD %>% 
  mutate(Metadata_sample_name = paste(Event_number, Depth_m, sep = "_"),
         Sampling_type = as.factor("CTD")) %>% 
  inner_join(meta.group, by = c("Event_number" = "Event_number",
                                "Depth_m" = "Depth_m",
                                "Depth_q" = "Depth_q",
                                "Longhurst_Prov" = "Longhurst_Prov",
                                "MertzGlacier" = "MertzGlacier")) %>%
  select(c("ACE_seq_name","Sampling_type","Metadata_sample_name","Event_number","Depth_m","Depth_q","Depth_q_new","Size_fraction","Event_start_date","Event_end_date","Latitude","Longitude","Longhurst_Prov","MertzGlacier","Water_mass.y", "Water_mass_simplified", "Water_mass.col", "Water_mass_simplified.col")) %>% 
  inner_join(metagenomes, by = c("ACE_seq_name" = "ace_seq_name")) %>% 
  select(-c("r1","r2","DNA.Genoscope.code")) %>% 
  rename(Sample = sample, Group = group, Water_mass = Water_mass.y)


metaG.UDW.id <- meta.udw[meta.udw$ACE_seq_name %in% metagenomes$ace_seq_name, ] %>% # 14 samples x 87 variables (2 samples in common with CTD)
  mutate(Metadata_sample_name = paste(Event_number, Depth_m, sep = "_"),
         Sampling_type = as.factor("UDW"),
         Depth_q_new = Depth_q, # As UDW are only SRF samples
         Water_mass_simplified = Water_mass) %>% # As UDW are only SW
  select(c("ACE_seq_name","Sampling_type","Metadata_sample_name","Event_number","Depth_m","Depth_q","Depth_q_new","Size_fraction","Event_start_date","Event_end_date","Latitude","Longitude","Longhurst_Prov","Water_mass","Water_mass_simplified")) %>% 
  droplevels() %>% 
  inner_join(metagenomes, by = c("ACE_seq_name" = "ace_seq_name")) %>% 
  select(-c("r1","r2","DNA.Genoscope.code")) %>% 
  rename(Sample = sample, Group = group) %>% 
  mutate(Water_mass.col = Water_mass,
         Water_mass_simplified.col = Water_mass_simplified)
levels(metaG.UDW.id$Water_mass.col) <- c("#2FB0B7","#7C98A3")
levels(metaG.UDW.id$Water_mass_simplified.col) <- c("#4DBBD5FF","#00A087FF")

metaG.merged.id <- metaG.CTD.id %>% 
  full_join(metaG.UDW.id)
```

    ## Joining with `by = join_by(ACE_seq_name, Sampling_type, Metadata_sample_name,
    ## Event_number, Depth_m, Depth_q, Depth_q_new, Size_fraction, Event_start_date,
    ## Event_end_date, Latitude, Longitude, Longhurst_Prov, Water_mass,
    ## Water_mass_simplified, Water_mass.col, Water_mass_simplified.col, Sample,
    ## Group)`

``` r
metaG.CTD.id[which(metaG.CTD.id$ACE_seq_name %in% metaG.UDW.id$ACE_seq_name),] # Samples in common between UDW and CTD
```

    ##         ACE_seq_name Sampling_type Metadata_sample_name Event_number Depth_m
    ## 190 ACE3_2970_74A_4m           CTD               2970_4         2970       4
    ## 191 ACE3_2970_74B_4m           CTD               2970_4         2970       4
    ##     Depth_q Depth_q_new Size_fraction        Event_start_date
    ## 190     SRF         SRF         >3 µm 2017-03-14T22:21:00Z+00
    ## 191     SRF         SRF      0.2-3 µm 2017-03-14T22:21:00Z+00
    ##              Event_end_date  Latitude Longitude Longhurst_Prov MertzGlacier
    ## 190 2017-03-14T22:38:00Z+00 -48.99636  9.002573           SANT        FALSE
    ## 191 2017-03-14T22:38:00Z+00 -48.99636  9.002573           SANT        FALSE
    ##     Water_mass Water_mass_simplified Water_mass.col Water_mass_simplified.col
    ## 190       SASW                  SASW        #7C98A3                 #00A087FF
    ## 191       SASW                  SASW        #7C98A3                 #00A087FF
    ##         Sample Group
    ## 190  sample_61   G61
    ## 191 sample_171  G171

``` r
metaG.merged.id <- metaG.merged.id %>% 
  filter(Sampling_type == "CTD" | Metadata_sample_name != "2970_4") %>% # To not duplicate 2970_4
  mutate(Event_start_date = ymd_hms(Event_start_date),
         Event_end_date = ymd_hms(Event_end_date))

# Assign MertzGlacier variable for sample 34 that is in common between CTD and UDW 
metaG.merged.id[which(metaG.merged.id$Sample == "sample_34"),]$MertzGlacier = as.factor("FALSE")
metaG.merged.id[which(metaG.merged.id$Sampling_type == "UDW" & metaG.merged.id$Event_number != 1104),]$MertzGlacier = as.factor("FALSE")
metaG.merged.id[which(metaG.merged.id$Sampling_type == "UDW" & metaG.merged.id$Event_number == 1104),]$MertzGlacier = as.factor("TRUE.Surround")
```

    ## Warning in `[<-.factor`(`*tmp*`, iseq, value = structure(c(1L, 1L), class =
    ## "factor", levels = "TRUE.Surround")): invalid factor level, NA generated

``` r
metaG.merged.id <- metaG.merged.id %>% 
  arrange(Sample)
```

``` r
# Compute distance to land for CTD and UDW samples
metaG.merged.id <- dist2land(metaG.merged.id, lon = "Longitude", lat = "Latitude") %>% 
  rename(Distance_land = ldist)
```

    ## Using DecimalDegree as land shapes.

    ## Calculating distances...

    ## Returning great circle spherical distances from land as kilometers.

``` r
write.csv(metaG.merged.id, file = "./Metadata/metadata_metaG_merged_ID.csv")
saveRDS(metaG.merged.id, file = "./R_Data/metadata_metaG_grouping_variables.rds")
```

# CTD and UDW metadata table for metagenomic samples

The aim here is to build a metadata table that matches CTD metagenomic
samples, including the size fractions. We will build a similar table for
UDW samples, but those tables do not need to be matched since the
environmental variables are not comparable.

This was also done in the script where it was needed but we write it
here so as to only load one table at the beginning of the downstream
scripts, creating a table that includes size fractions (as this will be
needed in analyses involving metagenomes).

``` r
meta.metaG.knnImp <- meta.sum.knnImp %>% 
  rownames_to_column(var = "Metadata_sample_name") %>% 
  inner_join(metaG.merged.id, by = "Metadata_sample_name") %>% 
  select(-c(Depth_m.y, Latitude.y, Longitude.y)) %>% # keep centered-scaled duplicated variables
  rename(Depth_m = Depth_m.x,
         Longitude = Longitude.x,
         Latitude = Latitude.x)

meta.metaG <- meta.sum.preprocess %>% 
  rownames_to_column(var = "Metadata_sample_name") %>% 
  inner_join(metaG.merged.id, by = "Metadata_sample_name") %>% 
  select(-c(Depth_m.y, Latitude.y, Longitude.y)) %>% # keep centered-scaled duplicated variables
  rename(Depth_m = Depth_m.x,
         Longitude = Longitude.x,
         Latitude = Latitude.x)

meta.metaG.bagImp <- meta.sum.bagImp %>% 
  rownames_to_column(var = "Metadata_sample_name") %>% 
  inner_join(metaG.merged.id, by = "Metadata_sample_name") %>% 
  select(-c(Depth_m.y, Latitude.y, Longitude.y)) %>% # keep centered-scaled duplicated variables
  rename(Depth_m = Depth_m.x,
         Longitude = Longitude.x,
         Latitude = Latitude.x)
```

``` r
saveRDS(meta.metaG.knnImp, file = "./R_Data/metadata_metaG_knnImp.rds")
saveRDS(meta.metaG.bagImp, file = "./R_Data/metadata_metaG_bagImp.rds")
saveRDS(meta.metaG, file = "./R_Data/metadata_metaG_preprocessed.rds")
```

``` r
meta.metaG %>% 
  group_by(Water_mass_simplified) %>% 
  summarize_if(.predicate = is.numeric, .funs = list(mean = mean, sd = sd), na.rm = TRUE)
```

    ## # A tibble: 8 × 59
    ##   Water_mass_simplified dCd_nmol.kg_mean dFe_nmol.kg_mean dZn_nmol.kg_mean
    ##   <fct>                            <dbl>            <dbl>            <dbl>
    ## 1 AASW                             0.374           0.0762            1.95 
    ## 2 SASW                             0.304           0.0232            0.467
    ## 3 STSW                             0.129           0.0868            0.430
    ## 4 WW influenced                    0.724           0.148             3.29 
    ## 5 SW-AAIW                          0.368           0.0713            1.23 
    ## 6 CDW influenced                   0.814           0.276             5.49 
    ## 7 DSW                              0.729           0.192             4.99 
    ## 8 AABW                           NaN             NaN               NaN    
    ## # ℹ 55 more variables: Bacteria_HDNA.cell.mL_mean <dbl>,
    ## #   Bacteria_LDNA.cell.mL_mean <dbl>, Tot_bacteria_cell.mL_mean <dbl>,
    ## #   TPZT_µM.C_mean <dbl>, POC_µM_mean <dbl>, d13C_mean <dbl>,
    ## #   PON.µmol_mean <dbl>, d15N_mean <dbl>, Fluorescence_mg.m3_mean <dbl>,
    ## #   PAR_µmol.m2.s_mean <dbl>, Bsi_µM_mean <dbl>, NOx_µmol.L_mean <dbl>,
    ## #   Silicic_acid_µmol.L_mean <dbl>, Nitrite_µmol.L_mean <dbl>,
    ## #   Phosphate_µmol.L_mean <dbl>, Ammonium_µmol.L_mean <dbl>, …

``` r
# Number of observations
meta.metaG %>% 
  group_by(Water_mass_simplified) %>% 
  summarize(Nb_observations = n()) # this includes NAs because there are a different number for each variable
```

    ## # A tibble: 8 × 2
    ##   Water_mass_simplified Nb_observations
    ##   <fct>                           <int>
    ## 1 AASW                               94
    ## 2 SASW                               25
    ## 3 STSW                               22
    ## 4 WW influenced                      15
    ## 5 SW-AAIW                            11
    ## 6 CDW influenced                     32
    ## 7 DSW                                 7
    ## 8 AABW                                1
