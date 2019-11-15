Determinants of helminth activation energy (Sharpe-Schoolfield)
================
Dan Benesh
November 15, 2019

The purpose of this notebook is to explore and model the activation
energies for helminth development. I specifically use activation
energies fit with the **Sharpe-Schoolfield** model.

# Get data organized

First, I get the data organized for modelling. From the table of curve
fits, I take the **Sharpe-Schoolfield** activation energies and add them
to the table including the predictors.

How many AE estimates are there?

    ## [1] 142

How many unique species are in the data?

    ## [1] 88

So there are a fair number of species with multiple curves.

### Look at the response - activation energy

What is the distribution of activation energies?

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

The median activation energy was 1.04.

How confident are we in each activation energy estimate? For this, we
look at the distribution of standard error estimates for the activation
energies. Most are small, but there is an outlier.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

All else equal, we should have more confidence in curves fit with a
larger number of temperatures. To reflect this, standard error should
decrease with sample size (i.e. number of temperatures). But this is not
obvious when we plot standard error as a function of temperatures.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

When we zoom in to exclude the outlier, there is still not a clear
relationship between standard error and the number of temperatures.
However, the largest standard errors are observed with the lowest sample
sizes (left-hand part of plot), which is reassuring.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

With less precision (higher standard errors), we would expect a broader
distribution of activation energy estimates. Let’s look at the
distribution of activation energies split by whether standard error was
high or low (above or below the median SE). Surprisingly, the estimates
with high standard errors do not have a wider distribution, opposite to
expectations.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

The median of the distribution also seems to vary with standard error;
high standard error is associated with higher activation energy.

    ## # A tibble: 2 x 2
    ##   se_cat              median_AE
    ##   <chr>                   <dbl>
    ## 1 high standard error     1.21 
    ## 2 low standard error      0.786

Thus, when we weight the data points by the standard error of the
activation energies, it should shift the overall mean down.

## Same experiment, different metrics

In some cases, multiple TPCs were fit to data from the same experiment.
For example, one curve fit for minimum developmental time and one for
maximum developmental time. These are not independent because the same
individual parasites contribute to both metrics. In these cases AE
measurements should be highly correlated.

Let’s make a plot for the species with multiple activation energy
measurements. When we do this, we see that AE estimates from the same
experiment tend to cluster as expected.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

This indicates how these values are not independent. Since this is an
issue that only affects a few species, I’m not inclined to devise an
`experiment` random effect that accounts for this pseudoreplication
statistically. Rather, I would probably only take a single AE per
experiment. But which one?

We can zoom in on just these cases where there are multiple measurements
from the same experiment. Unlike for
[Arrhenius](06_fit_phylo_models_Arr_Nov4.md) activation energies, there
was not an obvious bias in how minimum vs max development times impacted
activation energy estimates.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

This is related to the concerns Peter has raised about *apparent*
developmental rates, in that development is not measured for all viable
worms (ones that die before maturing are excluded). We can examine
whether activation energy varies systematically with developmental
metric. It does, in that min devo times tend to yield higher AE
estimates.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Overall, most activation energies were based on mean developmental
rates. Thus, in the few cases where multiple TPCs were fit to data from
the same experiment, but different metrics, I would only retain the
activation energy for mean devo time.

    ## 
    ##           min 50% developed      midpoint          mean           max 
    ##            39            22             9            52             7 
    ## not specified 
    ##            13

In some cases, the source of repeated measures on a single species was
recorded. For example, if the temperature dependence of multiple
parasite populations or life stages was studied. We can also plot these
variables in a similar to above. I do not notice an obvious trend.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

After removing AEs from the same experiment, how many AE estimates per
species?

    ##    binomial               n        
    ##  Length:87          Min.   :1.000  
    ##  Class :character   1st Qu.:1.000  
    ##  Mode  :character   Median :1.000  
    ##                     Mean   :1.483  
    ##                     3rd Qu.:2.000  
    ##                     Max.   :8.000

Usually just one, but up to 8. For visualizing, I’ll randomly take one
AE per species.

This reduces the data from n = 129 to n = 87.

Make a plot of the distribution of activation energy across the tree.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

The bars on the right are the estimated activation energies. Maybe some
clades at the top (nematodes) have higher values than those at the
bottom (flatworms), but it is not conspicuous. I don’t think that is
super surprising, because I assume AE is a parameter that is subject to
a fair amount of measurement error.

# Phylogenetic mixed models

I’ll fit phylogenetic mixed models (see
[here](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2009.01915.x)
for an overview of these models) using the `MCMCglmm` library. These are
not simple models, but the Bayesian `MCMCglmm` package does allow us to
fit both within- and between-species effects in the same model. It also
allows us to weight the activation energies by their standard errors.

Before diving into the models, though, let’s outline a model-building
strategy. First, I fit a ‘base’ model with the tangential variables.
Here, I would include phylogeny and within-species variation. From the
outset, I weight each data point by the standard error of activation
energy. This weight was incorporated into the model as suggested
[here](https://ourcodingclub.github.io/2018/01/22/mcmcglmm.html) and in
the `MCMCglmm` course notes. Second, I add biologically interesting
variables to this base model. Here, I would categorize interesting
variables into two groups: (1) characteristics of the external
environment (like latitude, aquatic vs terrestrial, etc.) and (2)
characteristics of the parasite (e.g. in or out of host, stage in the
life cycle, etc.). With this as a plan, let’s get modelling.

### Within-species effect

First, I create a null model that includes a parasite species random
effect, but not a phylogenetic effect. Throughout, data points are
weighted by their standard error.

After fitting an MCMC model, a first quality check is looking at chain
mixing. Essentially, we want the chain to bounce back and forth randomly
- it shouldn’t get stuck at particular parameter estimates.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

The within-species effect, `binomial`, mixes well and is non-zero,
suggesting that multiple AE measurements on the same species tend to be
similar. The variance component for `std-error` is fixed at 1, such that
the weights assigned do not change as the posterior distribution is
sampled by the Markov chain.

The model intercept (the average activation energy) is about 1.05.

    ## 
    ##  Iterations = 3001:22991
    ##  Thinning interval  = 10
    ##  Sample size  = 2000 
    ## 
    ##  DIC: 147.1013 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial    0.1332  0.05591   0.2133     1296
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1257  0.07568   0.1836     1472
    ## 
    ##  Location effects: E_SS ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## (Intercept)     1.038    0.925    1.134     1849 <5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the proportion of variation attributable to the within-species
effect.

    ## [1] 0.5239126

It is high. On the one hand this makes sense, as many species only had a
single measurement (and hence no residual value after accounting for a
‘species effect’). Here’s an attempt to visualize the variation in
activation energy estimates for the species that had multiple
measurements. The species are ordered by their median activation energy.
One can see that some species tend to have very similar measurements,
while for other species, there is a fair amount of spread in the
estimated activation energies.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

### Phylogenetic effect

Now, let’s add phylogenetic effects to this mixed model. To do this, we
create the phylogenetic covariance matrix and then include this as a
random effect.

For this model, the chain needs to be run for longer to get reasonable
mixing, partly because the two random effects are related. If the model
estimates higher values for the within-species variance component, it
estimates less for the phylogenetic component, and vice versa.

After fitting, we can again plot the chains. They look good. The
within-species effect (`binomial`) is smaller than the phylogenetic
effect (`tree_tips`) and it pushes up against zero.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

These results suggest that phylogenetically related species tend to have
similar AEs, and when a species is measured multiple times, the
estimated AEs are similar, though they are not far more similar than
what they are expected to be based on phylogeny (i.e. the within-species
effect does not explain much additional variation beyond phylogeny).

The phylogenetic effect is strong. It explained 66% of the variation in
AE.

By contrast, the ‘species’ effect explained 5% of the variation.

Here’s a plot splitting AE by parasite order. It does look like some
clades tend to have higher or lower AEs.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

#### Developmental metric

As seen above, the chosen metric of development can affect activation
energy estimates. To check this effect, I added a fixed factor to the
model that distinguishes between activation estimates based on mean,
min, or unspecific developmental times.

The differences among these groups were rather minimal (non-significant
p-values).

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 124.7292 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02294 0.002137  0.06035    922.2
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2632  0.05889   0.4835     1014
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1106  0.06948   0.1503     1200
    ## 
    ##  Location effects: E_SS ~ trait_detail3 
    ## 
    ##                            post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## (Intercept)                  0.98718  0.52002  1.49566     1433 <8e-04 ***
    ## trait_detail3min             0.02496 -0.14888  0.21720     1043  0.822    
    ## trait_detail3not specified   0.01646 -0.23592  0.26216     1200  0.887    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Thus, it is not clear how valuable this term is in the model. Therefore,
I have not included it in further analyses.

### Environmental variables

Now, we move onto environmental variables of interest, though it is
important to keep phylogenetic effects in mind during this modeling
exercies. Essentially, any explanatory variable that is phylogenetically
structured will “compete” with the phylogenetic random effect to explain
variation in AE.

I have tried to find a series of logical steps to testing environmental
correlates of helminth AE. First and foremost, we expect AE to be
related to environmental temperature, so I start with temperature
variables.

#### Mean temperature

There are three temperature measurements in the data: mean, max, and min
for the location. I would expect the three temp variables to be tightly
correlated, so let’s start by just adding the mean to the model. Note
that I have centered the continuous variables around their median, such
that the model parameters are estimated at e.g. the data’s median annual
temperature.

The effect of mean locale temperature is negative but not significant.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 122.9785 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02558 0.002244  0.06203      768
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2401  0.06138   0.4614     1182
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1097  0.07497   0.1541     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen 
    ## 
    ##                   post.mean  l-95% CI  u-95% CI eff.samp   pMCMC   
    ## (Intercept)        1.013962  0.581701  1.525805     1200 0.00167 **
    ## mean_ann_temp_cen -0.010675 -0.024824  0.001682     1344 0.10000 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

For species living at lower temperatures, development speeds up more as
temperature increases. It also does not look like this depends on
whether an exact study location was given.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

#### Temperature variation

A mean temperature does not indicate the temperature range a parasite is
exposed to. Thus, as the next modelling step, we’ll include location min
and max temps. None of these terms are significant and the model DIC
decreases, suggesting adding temperature extremes does not improve the
model.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 125.1544 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02757 0.002494  0.06896    960.9
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2443  0.04023    0.471    986.4
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1089  0.06813    0.154     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + max_month_temp_cen + min_month_temp_cen 
    ## 
    ##                     post.mean   l-95% CI   u-95% CI eff.samp  pMCMC    
    ## (Intercept)         0.9771662  0.4261862  1.4740155     1200 <8e-04 ***
    ## mean_ann_temp_cen  -0.0008651 -0.1154674  0.1131404     1113  0.988    
    ## max_month_temp_cen -0.0019614 -0.0709731  0.0640072     1180  0.940    
    ## min_month_temp_cen -0.0064563 -0.0635448  0.0479945     1200  0.833    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Interestingly, the parameter estimates for max and min temps are
positive, which is opposite to the trend seen with mean temp. But when
we plot AE as a function of max and min, we see a negative relationship
just like for mean temp: low temps, high AE.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

This is suggestive that, once controlling for mean temp, activation
energies increase with higher max and min temps. However, as neither min
nor max seems to explain additional variation in activation energy, we
won’t consider these temperature extremes further.

#### Seasonality

Variation in temperature is one component of seasonality, but not the
only one. Season length, photoperiod, rainfall, etc. may all vary with
seasons and perhaps impact the thermal-dependence of parasite
development. The dataset contains two variables that might capture
seasonality: latitude and distribution. Distribution zones were defined
as ‘tropics’, ‘temperate’, ‘polar’, etc., which means this variable
overlaps substantially with latitude. Consequently, I examine these two
variables in parallel (i.e. I won’t put both in the same model).

##### Latitude

Let’s look at latitude first, as it is a single continuous variable
while distribution is multiple categories. In some cases, a study gave a
detailed location, so our data contains exact latitudes. In other cases,
no site was given and the author’s location was taken as a rough proxy
of the study’s location. Most records were from the Northern hemisphere,
but a few were from the Southern hemisphere. Thus, I took the absolute
value of latitude, essentially the distance from the equator.

There is not a clear relationship between AE and latitude.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 124.4942 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02458 0.001812  0.06168    931.1
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2265  0.04645   0.4268     1023
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1118  0.07219   0.1478     1162
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + lat_abs 
    ## 
    ##                    post.mean   l-95% CI   u-95% CI eff.samp  pMCMC    
    ## (Intercept)        1.368e+00  6.377e-01  2.058e+00     1200 <8e-04 ***
    ## mean_ann_temp_cen -2.224e-02 -4.403e-02  7.486e-05     1200 0.0583 .  
    ## lat_abs           -7.360e-03 -1.799e-02  4.885e-03     1200 0.2317    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Surprisingly, the parameter estimate for latitude is negative, but when
we plot it, we see that activation energies tend to increase with
latitude.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

This suggests that, like for `min` and `max` temperature, correlations
among variables (here between mean temp and latitude) might cause
counterintuitive or spurious results. When we remove the effect of mean
temperature on activation energy (i.e. we take the residuals of the
model with just mean temp), and then plot the residual variation as a
function of latitude, we see essentially no relationship between
latitude and activation energy.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

Let’s look at the next variable, which is related to latitude.

##### Climate zones

Our latitude variable is not a perfect representation of a parasite’s
‘ancestral’ habitat. Many species in the data are associated with
people or livestock and have thus been moved all around the world. This
is reflected in the distribution zone variable: it contains a category
for parasites with a broad distribution.

However, the distribution categories are still tightly linked to
latitude.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

Nonetheless, let’s add this variable to the model instead of latitude.

Several contrasts among these categories are marginally significant and
the model DIC decreases.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 119.2207 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.03176 0.002352  0.07043    839.9
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2288  0.04246   0.4361    904.2
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1028   0.0634   0.1415     1072
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone 
    ## 
    ##                      post.mean  l-95% CI  u-95% CI eff.samp   pMCMC   
    ## (Intercept)           0.952918  0.411748  1.485703     1337 0.00333 **
    ## mean_ann_temp_cen    -0.021485 -0.037456 -0.005617     1200 0.00333 **
    ## dist_zonetemperate   -0.061617 -0.399094  0.247721     1200 0.71667   
    ## dist_zonesubtropical  0.260994 -0.208490  0.715162     1080 0.27667   
    ## dist_zonetropical     0.441468 -0.084978  0.959347     1091 0.10333   
    ## dist_zoneglobal       0.143062 -0.157606  0.467471     1103 0.36500   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Consistent with the pattern for temperature and latitude, polar species
tend to have higher AE. But weirdly, a non-linear relationship with
latitude seems to emerge, with AE peaking in polar habitats, decreasing
in temperate habitats, and then increasing again toward the equator.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

Latitude and distribution overlap and presumably explain the same
variation in AE, but which explains it better? We can compare the models
with latitude vs distribution zone using information criteria. Here is
the progression of model fits thus far (low DIC, better model):

    ## [1] "DIC for model with phylogeny: 122.4"

    ## [1] "DIC for model with just mean temp: 123"

    ## [1] "DIC for model with latitude: 124.5"

    ## [1] "DIC for model with distribution zone: 119.2"

According to this measure, the model with distribution is better than
the one with latitude, which is not surprising given the seemingly
non-linear relationship with latitude. Another explanation for this, is
that including distribution zone strengthens the estimated effect of
mean annual temperature; it explains residual variation around the
temp-AE relationship.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

Thus, I’ll retain distribution in the model as a proxy for seasonality,
but I’ll drop latitude.

#### Habitat - aquatic vs freshwater

As a final environmental variable, we can add parasite habitat: aquatic
vs terrestrial. Water dampens temperature swings, such that seasonality
is less pronounced in aquatic habitats compared to terrestrial ones.
Thus, this variable might modify the effect of seasonality proxies like
distribution zone. (Note: for a couple TPCs, I filled in habitat and
distribution data.)

Adding habitat to the model improves DIC, but the parameter is not
significant.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 113.7393 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02333 0.001804  0.05849    920.3
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2885  0.08083   0.5578     1147
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.09787  0.06448   0.1365     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone + habitat_d1 
    ## 
    ##                       post.mean l-95% CI u-95% CI eff.samp   pMCMC   
    ## (Intercept)             0.83209  0.18187  1.38634     1200 0.00833 **
    ## mean_ann_temp_cen      -0.02340 -0.03891 -0.00785     1300 0.00667 **
    ## dist_zonetemperate     -0.05825 -0.39143  0.26210     1200 0.74000   
    ## dist_zonesubtropical    0.25197 -0.20630  0.67250     1200 0.25500   
    ## dist_zonetropical       0.46725 -0.04093  0.95169     1200 0.07500 . 
    ## dist_zoneglobal         0.10686 -0.20042  0.44292     1200 0.55500   
    ## habitat_d1terrestrial   0.31766 -0.05434  0.68286     1366 0.07833 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Terrestrial species have higher AE on average than aquatic species
(combines both freshwater and marine species).

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

We might expect the aquatic - terrestrial dichotomy to be particularly
important in more seasonal locations. That is, there may be an
interaction between our seasonality proxy (distribution) and habitat.
Unfortunately, this cuts the data quite thin, as there are not many
freshwater and terrestrial species in each distribution category.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

When we add this interaction to the model, it is not an improvement.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 116.4335 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02567 0.001941  0.06624    794.3
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2932  0.06779   0.5605    977.9
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.09748  0.06235   0.1325     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp + dist_zone * habitat_d1 
    ## 
    ##                                          post.mean l-95% CI u-95% CI
    ## (Intercept)                                1.15984  0.49935  1.80042
    ## mean_ann_temp                             -0.02563 -0.04368 -0.01037
    ## dist_zonetemperate                        -0.18385 -0.70363  0.30771
    ## dist_zonesubtropical                       0.39014 -0.05150  0.90122
    ## dist_zonetropical                          0.53643 -0.21679  1.27108
    ## dist_zoneglobal                           -0.07868 -0.55400  0.47658
    ## habitat_d1terrestrial                      0.11079 -0.48694  0.73893
    ## dist_zonetemperate:habitat_d1terrestrial   0.23352 -0.40760  0.89546
    ## dist_zonetropical:habitat_d1terrestrial   -0.01408 -0.90697  0.87495
    ## dist_zoneglobal:habitat_d1terrestrial      0.31366 -0.28045  0.96350
    ##                                          eff.samp  pMCMC    
    ## (Intercept)                                1013.9  0.005 ** 
    ## mean_ann_temp                              1200.0 <8e-04 ***
    ## dist_zonetemperate                          836.2  0.490    
    ## dist_zonesubtropical                       1200.0  0.115    
    ## dist_zonetropical                          1216.5  0.178    
    ## dist_zoneglobal                            1200.0  0.763    
    ## habitat_d1terrestrial                      1200.0  0.773    
    ## dist_zonetemperate:habitat_d1terrestrial    850.1  0.495    
    ## dist_zonetropical:habitat_d1terrestrial    1346.4  0.967    
    ## dist_zoneglobal:habitat_d1terrestrial      1200.0  0.312    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

To summarize, AE tends to be higher with low mean locale temperatures
and terrestrial habitats.

### Host-parasite characteristics

Now, we turn from the environment to variables describing the
host-parasite interaction. I consider a few variables: (1) host type,
(2) stage in or out of host, and (3) target host. Compared to the
environmental variables, these variables are not as obviously
confounded.

#### Plant vs animal parasite

A major dichotomy is between between plant and animal parasites. This
might be relevant if plant and animal hosts have different activation
energies or if they differ in how much they ‘protect’ parasites from
temperature variation.

Distinguishing plants and animals is not significant, though the model
DIC goes down. Plant parasites tend to have lower AE than animal
parasites. Also, I imagine this distinction is strongly conflated with
phylogeny.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 114.0268 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial    0.0242 0.001716  0.06102     1017
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.2923  0.06862    0.553     1002
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.09754  0.06626   0.1359     1101
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone + habitat_d1 + plant_anim 
    ## 
    ##                        post.mean   l-95% CI   u-95% CI eff.samp   pMCMC   
    ## (Intercept)            0.8562047  0.3088427  1.4702493     1200 0.01333 * 
    ## mean_ann_temp_cen     -0.0237974 -0.0402798 -0.0091996     1200 0.00833 **
    ## dist_zonetemperate    -0.0576272 -0.3870943  0.2625728     1200 0.72833   
    ## dist_zonesubtropical   0.2788682 -0.1724334  0.7000874     1200 0.21333   
    ## dist_zonetropical      0.4664205 -0.0696012  0.9668356     1200 0.08167 . 
    ## dist_zoneglobal        0.0954996 -0.2184935  0.4475779     1200 0.57167   
    ## habitat_d1terrestrial  0.3612250  0.0005105  0.7926231     1040 0.06333 . 
    ## plant_animplant       -0.2177108 -0.7919777  0.3646671     1200 0.47000   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the plot comparing plant and animal parasites.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-85-1.png)<!-- -->

#### Endotherm in life cycle

Besides the plant - animal distinction, parasites can be distinguished
by whether or not they have an endotherm (a bird or mammal) in the life
cycle. Since endotherm body temperatures are commonly higher than
ambient temperatures, we might expect AE to differ between species with
or without such a host in the life cycle.

This variable basically retains the previous plant-animal dichotomy, but
it additionally splits animal parasites into two groups.

    ##                             
    ##                              animal plant
    ##   endotherm in life cycle        84     0
    ##   no endotherm in life cycle     21    24

We can thus combine these into a single variable to add to the model.

The term contrasting animal parasites with and without endotherms in the
cycle was significant, and it suggested that species without endotherms
have higher AEs.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 108.1567 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02495 0.002042  0.07062    729.4
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.3134  0.05352   0.5917    944.6
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.09161  0.06027   0.1273     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 
    ## 
    ##                                       post.mean   l-95% CI   u-95% CI
    ## (Intercept)                           0.7436482  0.1186256  1.4316211
    ## mean_ann_temp_cen                    -0.0256220 -0.0405775 -0.0087099
    ## dist_zonetemperate                   -0.0547950 -0.3544871  0.2575825
    ## dist_zonesubtropical                  0.3187916 -0.1042593  0.7577990
    ## dist_zonetropical                     0.4891294 -0.0104361  0.9511119
    ## dist_zoneglobal                       0.1522220 -0.1787731  0.4872490
    ## habitat_d1terrestrial                 0.3703497  0.0006022  0.7731769
    ## endo_ecto2no endotherm in life cycle  0.2571359 -0.0276460  0.5292279
    ## endo_ecto2plant                      -0.0796691 -0.7005917  0.5192805
    ##                                      eff.samp   pMCMC   
    ## (Intercept)                              1200 0.04000 * 
    ## mean_ann_temp_cen                        1200 0.00167 **
    ## dist_zonetemperate                       1200 0.75333   
    ## dist_zonesubtropical                     1200 0.15833   
    ## dist_zonetropical                        1200 0.05833 . 
    ## dist_zoneglobal                          1200 0.35500   
    ## habitat_d1terrestrial                    1200 0.05333 . 
    ## endo_ecto2no endotherm in life cycle     1200 0.07500 . 
    ## endo_ecto2plant                          1200 0.76500   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

However, this trend is not visible in the raw data, suggesting it is
contingent on the covariates in the data.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-91-1.png)<!-- -->

When we plot these groups as a function of mean annual temperature, we
then see that, at a given temperature, species exclusively infecting
ectotherms have higher AE. This pattern is not strong, but it is
consistent with the idea that environmental stability lowers AE
(endotherms provide a stable temperature for helminths).

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

#### Stage in or out of host?

Stages inside a host may be more shielded from temperature variation
than stages in the environment. Let’s add this distinction to the model.

This term is not significant and the model DIC is worse.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 108.8492 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02589  0.00182  0.06927      711
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.3115  0.05222   0.6094    777.8
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.09092  0.05921   0.1327     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 + devo_type 
    ## 
    ##                                      post.mean  l-95% CI  u-95% CI
    ## (Intercept)                           0.725810  0.076112  1.347027
    ## mean_ann_temp_cen                    -0.025441 -0.041897 -0.010049
    ## dist_zonetemperate                   -0.055142 -0.376022  0.253646
    ## dist_zonesubtropical                  0.311335 -0.132067  0.753723
    ## dist_zonetropical                     0.491417  0.029672  1.014888
    ## dist_zoneglobal                       0.147803 -0.169794  0.475007
    ## habitat_d1terrestrial                 0.372614 -0.019080  0.767690
    ## endo_ecto2no endotherm in life cycle  0.269254 -0.007155  0.538280
    ## endo_ecto2plant                      -0.072333 -0.656251  0.624007
    ## devo_typeoutside host                 0.072716 -0.130338  0.283947
    ##                                      eff.samp   pMCMC   
    ## (Intercept)                              1326 0.03500 * 
    ## mean_ann_temp_cen                        1200 0.00167 **
    ## dist_zonetemperate                       1200 0.74667   
    ## dist_zonesubtropical                     1026 0.18333   
    ## dist_zonetropical                        1048 0.04500 * 
    ## dist_zoneglobal                          1222 0.37667   
    ## habitat_d1terrestrial                    1200 0.05667 . 
    ## endo_ecto2no endotherm in life cycle     1200 0.05500 . 
    ## endo_ecto2plant                          1200 0.81500   
    ## devo_typeoutside host                    1200 0.49833   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Though the difference is not significant, worms outside the host have
slightly higher AEs, which is also consistent with the idea that
environmental variation favors higher AE. The plot also shows that
habitat and in/outside host are conflated, as proportionally more of the
terrestrial values are outside the host and more of the aquatic values
are inside the host.

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->

#### Target host

Stages free in the environment are generally propagules that target the
next host in the life cycle. Are the parasites targeting a vertebrate or
an invertebrate? It is not obvious to me that this should matter. But
nonetheless, we can check it.

Distinguishing between invertebrates and vertebrates as next host
slightly improves the model, but the contrasts are not significant.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: 106.8255 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.02373 0.002749  0.06257    666.7
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.3314  0.03582   0.6324    819.1
    ## 
    ##                ~idh(stderr_ss):units
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## stderr_ss.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.08924  0.05613   0.1272     1200
    ## 
    ##  Location effects: E_SS ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 + devo_type3 
    ## 
    ##                                            post.mean  l-95% CI  u-95% CI
    ## (Intercept)                                 0.709785  0.082003  1.338254
    ## mean_ann_temp_cen                          -0.022714 -0.038049 -0.006558
    ## dist_zonetemperate                         -0.023062 -0.345818  0.267455
    ## dist_zonesubtropical                        0.242872 -0.262111  0.674594
    ## dist_zonetropical                           0.462716 -0.048367  0.993707
    ## dist_zoneglobal                             0.166620 -0.164122  0.484147
    ## habitat_d1terrestrial                       0.482846  0.093425  0.880331
    ## endo_ecto2no endotherm in life cycle        0.131881 -0.183084  0.414702
    ## endo_ecto2plant                            -0.108034 -0.799525  0.507372
    ## devo_type3in plant                         -0.143994 -0.516228  0.271594
    ## devo_type3outside host: targets invert      0.175411 -0.131047  0.483747
    ## devo_type3outside host: targets vertebrate -0.301663 -0.727548  0.124546
    ##                                            eff.samp   pMCMC   
    ## (Intercept)                                    1200 0.04000 * 
    ## mean_ann_temp_cen                              1200 0.00833 **
    ## dist_zonetemperate                             1200 0.90500   
    ## dist_zonesubtropical                           1200 0.27500   
    ## dist_zonetropical                              1200 0.07333 . 
    ## dist_zoneglobal                                1200 0.31667   
    ## habitat_d1terrestrial                          1200 0.01500 * 
    ## endo_ecto2no endotherm in life cycle           1422 0.38833   
    ## endo_ecto2plant                                1200 0.69667   
    ## devo_type3in plant                             1200 0.43167   
    ## devo_type3outside host: targets invert         1200 0.28333   
    ## devo_type3outside host: targets vertebrate     1232 0.15333   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![](06_fit_phylo_models_SS_Nov4_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->

# Conclusions

In this analysis of activation energies in helminths, I account for
measurement variance (standard error of activation energy estimates) and
I account for two sources of pseudoreplication: (1) multiple
measurements on a single parasite species and (2) shared ancestry. A
mixed model including both factors suggested that shared ancestry has a
clear effect on AEs: related species tend to exhibit similar activation
energies.

Next, I evaluated various environmental characteristics and found that
activation energies tend to be higher at low temperatures and high
latitudes. Finally, I assessed characteristics of the host-parasite
interaction, and I found that animal parasites with endotherms in the
life cycle tend to have lower AE, consistent with the notion of
temperature stability affecting temperature sensitivities.

Finally, let’s quantitatively summarize the results. For each of the
fitted models, I compiled the df, the DIC, and the R<sup>2</sup>. I
calculated R<sup>2</sup> according to [Nakagawa and
Schielzeth 2013](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x%4010.1111/%28ISSN%292041-210X.STATSTOO).
They distinguish between marginal and conditional R<sup>2</sup>.
Marginal R<sup>2</sup> is the proportion of variation explained by the
fixed effects while conditional R<sup>2</sup> is the variation explained
by both fixed and random effects. When the term is followed by an
asterisk, I retained it in the model.

| model                                 | df | df\_used |    DIC |  R2m |  R2c |
| :------------------------------------ | -: | -------: | -----: | ---: | ---: |
| just within-species\*                 |  3 |       NA | 147.10 | 0.00 | 0.56 |
| \+ phylogeny\*                        |  4 |        1 | 122.37 | 0.00 | 0.72 |
| \+ mean annual temp\*                 |  5 |        1 | 122.98 | 0.00 | 0.71 |
| \+ min/max temps                      |  7 |        2 | 125.15 | 0.01 | 0.70 |
| \+ latitude                           |  6 |      \-1 | 124.49 | 0.01 | 0.74 |
| \+ distribution zone\*                |  9 |        3 | 119.22 | 0.05 | 0.77 |
| \+ habitat (aquatic vs terrestrial)\* | 10 |        1 | 113.74 | 0.12 | 0.82 |
| \+ distribution x habitat interacton  | 13 |        3 | 116.43 | 0.14 | 0.80 |
| \+ plant vs animal parasite\*         | 11 |      \-2 | 114.03 | 0.17 | 0.86 |
| \+ endotherm in life cycle?\*         | 12 |        1 | 108.16 | 0.17 | 0.83 |
| \+ stage in/out of host\*             | 13 |        1 | 108.85 | 0.14 | 0.85 |
| \+ target host: invert vs vert        | 15 |        2 | 106.83 | 0.17 | 0.85 |

The most complex model, one that includes phylogeny, seemingly important
environmental variables, and characteristics of the host-parasite
system, explained over 89% of the variation in activation energy. Much
of this is due to the random effects (phylogeny and repeated measures
account for \>70% of the variation), but around 15% of the variation can
be explained by fixed effects.
