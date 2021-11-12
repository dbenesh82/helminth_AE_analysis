Testing terms in Table 1
================

-   [Terms](#terms)
    -   [Base model](#base-model)
    -   [Development/survival metric](#developmentsurvival-metric)
-   [Environmental variables](#environmental-variables)
    -   [Mean temperature](#mean-temperature)
    -   [Mean temp x exact location](#mean-temp-x-exact-location)
    -   [Temperature variation](#temperature-variation)
    -   [Temp variation x exact
        location](#temp-variation-x-exact-location)
    -   [Latitude](#latitude)
    -   [Latitude x global
        distribution](#latitude-x-global-distribution)
    -   [Habitat - aquatic vs
        freshwater](#habitat---aquatic-vs-freshwater)
    -   [Habitat x latitude](#habitat-x-latitude)
-   [Host-parasite characteristics](#host-parasite-characteristics)
    -   [Stage in/out of host](#stage-inout-of-host)
    -   [Plant vs animal parasite](#plant-vs-animal-parasite)
    -   [Endotherm in life cycle](#endotherm-in-life-cycle)
-   [R<sup>2</sup> Table](#r2-table)

This document fits the series of models analyzing activation energies.
The results are summarized in Table 1.

Number of AE estimates:

    ## [1] 129

Number of species:

    ## [1] 87

# Terms

## Base model

## Development/survival metric

    ##                                         post.mean    l-95% CI   u-95% CI
    ## (Intercept)                           0.670061196  0.34128602 0.97289896
    ## simplified_metricmin                 -0.005723259 -0.10235035 0.09667503
    ## simplified_metricother/not specified  0.126466037 -0.01477946 0.27688497
    ##                                      eff.samp        pMCMC
    ## (Intercept)                          4800.000 0.0004166667
    ## simplified_metricmin                 4828.533 0.9066666667
    ## simplified_metricother/not specified 4800.000 0.0858333333

Not significant.

# Environmental variables

## Mean temperature

    ##                     post.mean    l-95% CI     u-95% CI eff.samp        pMCMC
    ## (Intercept)        0.69447069  0.40167267 0.9698179551     4800 0.0002083333
    ## mean_ann_temp_cen -0.00680309 -0.01352976 0.0002836502     4800 0.0579166667

Not significant.

## Mean temp x exact location

    ##                                                            post.mean
    ## (Intercept)                                              0.716473302
    ## mean_ann_temp_cen                                       -0.009932947
    ## sampling_location_availableyes, exact                   -0.033400115
    ## mean_ann_temp_cen:sampling_location_availableyes, exact  0.004165892
    ##                                                             l-95% CI
    ## (Intercept)                                              0.428653790
    ## mean_ann_temp_cen                                       -0.022198497
    ## sampling_location_availableyes, exact                   -0.138531161
    ## mean_ann_temp_cen:sampling_location_availableyes, exact -0.009558486
    ##                                                            u-95% CI eff.samp
    ## (Intercept)                                             1.010427590 4800.000
    ## mean_ann_temp_cen                                       0.001535114 4569.941
    ## sampling_location_availableyes, exact                   0.068393078 4800.000
    ## mean_ann_temp_cen:sampling_location_availableyes, exact 0.018584262 4471.145
    ##                                                                pMCMC
    ## (Intercept)                                             0.0002083333
    ## mean_ann_temp_cen                                       0.1000000000
    ## sampling_location_availableyes, exact                   0.5287500000
    ## mean_ann_temp_cen:sampling_location_availableyes, exact 0.5595833333

Not significant.

## Temperature variation

    ##                  post.mean     l-95% CI    u-95% CI eff.samp        pMCMC
    ## (Intercept)    0.679039696  0.378053470 0.981488727     4800 0.0002083333
    ## range_temp_cen 0.004503831 -0.000101902 0.009334732     4800 0.0650000000

Marg significant.

## Temp variation x exact location

    ##                                                         post.mean     l-95% CI
    ## (Intercept)                                           0.696670015  0.377079861
    ## range_temp_cen                                        0.003955204 -0.004828882
    ## sampling_location_availableyes, exact                -0.021026177 -0.116238633
    ## range_temp_cen:sampling_location_availableyes, exact  0.000845825 -0.009851303
    ##                                                        u-95% CI eff.samp
    ## (Intercept)                                          1.01503700 4559.678
    ## range_temp_cen                                       0.01314565 4800.000
    ## sampling_location_availableyes, exact                0.07685411 4879.564
    ## range_temp_cen:sampling_location_availableyes, exact 0.01099544 4800.000
    ##                                                             pMCMC
    ## (Intercept)                                          0.0002083333
    ## range_temp_cen                                       0.3870833333
    ## sampling_location_availableyes, exact                0.6762500000
    ## range_temp_cen:sampling_location_availableyes, exact 0.8683333333

Not significant.

## Latitude

    ##               post.mean     l-95% CI    u-95% CI eff.samp        pMCMC
    ## (Intercept) 0.692944206  0.394223913 0.977233619     4800 0.0004166667
    ## lat_abs_cen 0.000570807 -0.002757024 0.004227457     4800 0.7570833333

Not significant.

## Latitude x global distribution

    ##                                  post.mean     l-95% CI    u-95% CI eff.samp
    ## (Intercept)                   0.7062117641  0.391473531 1.001763075  4800.00
    ## lat_abs_cen                   0.0003268419 -0.004916177 0.005241343  4800.00
    ## globalnot global             -0.0214697833 -0.124997400 0.083472711  4800.00
    ## lat_abs_cen:globalnot global  0.0007422928 -0.005941174 0.007725447  4576.79
    ##                                     pMCMC
    ## (Intercept)                  0.0002083333
    ## lat_abs_cen                  0.9020833333
    ## globalnot global             0.6895833333
    ## lat_abs_cen:globalnot global 0.8333333333

Not significant.

## Habitat - aquatic vs freshwater

    ##                       post.mean    l-95% CI  u-95% CI eff.samp        pMCMC
    ## (Intercept)           0.6439462  0.31593338 0.9883916 4800.000 0.0004166667
    ## habitat_d1Terrestrial 0.1236610 -0.07751604 0.3261097 4536.443 0.2262500000

Not significant.

## Habitat x latitude

    ##                                      post.mean     l-95% CI    u-95% CI
    ## (Intercept)                       0.6445647368  0.295261910 0.969485680
    ## lat_abs_cen                       0.0001005484 -0.005796860 0.006323868
    ## habitat_d1Terrestrial             0.1388375952 -0.073756610 0.348182868
    ## lat_abs_cen:habitat_d1Terrestrial 0.0011096790 -0.006475623 0.008564729
    ##                                   eff.samp        pMCMC
    ## (Intercept)                       4800.000 0.0008333333
    ## lat_abs_cen                       4800.000 0.9608333333
    ## habitat_d1Terrestrial             4800.000 0.1912500000
    ## lat_abs_cen:habitat_d1Terrestrial 5263.483 0.7641666667

Not significant.

# Host-parasite characteristics

## Stage in/out of host

    ##                                post.mean   l-95% CI  u-95% CI eff.samp
    ## (Intercept)                   0.68783334  0.3554960 0.9950765 5271.424
    ## stage_in_out_hostoutside host 0.01192063 -0.1020254 0.1336555 4735.636
    ##                                      pMCMC
    ## (Intercept)                   0.0002083333
    ## stage_in_out_hostoutside host 0.8358333333

Not significant.

## Plant vs animal parasite

    ##                           post.mean   l-95% CI  u-95% CI eff.samp        pMCMC
    ## (Intercept)               0.6998694  0.3940683 1.0063042 4800.000 0.0002083333
    ## plant_anim_parasiteplant -0.1083868 -0.4215243 0.2037003 4481.542 0.4704166667

Not significant.

## Endotherm in life cycle

    ##                                        post.mean   l-95% CI  u-95% CI eff.samp
    ## (Intercept)                           0.69009485  0.3833816 0.9995382 4800.000
    ## endo_ecto2no endotherm in life cycle  0.05071574 -0.1012746 0.1994260 4800.000
    ## endo_ecto2plant                      -0.07875530 -0.3965075 0.2533927 4304.802
    ##                                             pMCMC
    ## (Intercept)                          0.0002083333
    ## endo_ecto2no endotherm in life cycle 0.5041666667
    ## endo_ecto2plant                      0.6045833333

Not significant

# R<sup>2</sup> Table

| model                    |  df |    DIC | R2m             | R2c                |
|:-------------------------|----:|-------:|:----------------|:-------------------|
| base                     |   3 | -40.97 | 0 \[0-0\]       | 0.78 \[0.58-0.89\] |
| metric                   |   5 | -42.07 | 0.01 \[0-0.05\] | 0.79 \[0.6-0.89\]  |
| mean temp                |   4 | -40.71 | 0 \[0-0.07\]    | 0.76 \[0.55-0.89\] |
| temp x location          |   6 | -38.91 | 0.02 \[0-0.08\] | 0.77 \[0.58-0.9\]  |
| temp range               |   4 | -44.86 | 0 \[0-0.04\]    | 0.8 \[0.61-0.9\]   |
| temp range x location    |   6 | -42.43 | 0.01 \[0-0.06\] | 0.79 \[0.62-0.91\] |
| latitude                 |   4 | -39.15 | 0 \[0-0.02\]    | 0.77 \[0.57-0.89\] |
| latitude x global        |   6 | -36.46 | 0 \[0-0.04\]    | 0.77 \[0.58-0.9\]  |
| habitat                  |   4 | -44.90 | 0 \[0-0.11\]    | 0.83 \[0.64-0.92\] |
| habitat x latitude       |   6 | -41.45 | 0.01 \[0-0.12\] | 0.85 \[0.63-0.92\] |
| stage in/out of host     |   4 | -39.78 | 0 \[0-0.03\]    | 0.79 \[0.59-0.9\]  |
| plant vs animal parasite |   4 | -39.49 | 0 \[0-0.15\]    | 0.78 \[0.57-0.9\]  |
| endotherm in life cycle? |   5 | -40.78 | 0 \[0-0.15\]    | 0.76 \[0.6-0.92\]  |
