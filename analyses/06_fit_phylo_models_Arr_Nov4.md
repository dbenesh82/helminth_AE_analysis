Determinants of helminth activation energy (Arrhenius)
================
Dan Benesh
November 15, 2019

The purpose of this notebook is to explore and model the activation
energies for helminth development. I specifically use activation
energies fit with the **Arrhenius** model.

# Get data organized

First, I get the data organized for modelling. From the table of curve
fits, I take the **Arrhenius** activation energies and add them to the
table including the predictors.

How many AE estimates are there?

    ## [1] 142

How many unique species are in the data?

    ## [1] 88

So there are a fair number of species with multiple curves.

### Look at the response - activation energy

What is the distribution of activation energies?

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

The median activation energy was 0.65.

How confident are we in each activation energy estimate? For this, we
look at the distribution of standard error estimates for the activation
energies. Most are small, but there is an outlier.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

All else equal, we should have more confidence in curves fit with a
larger number of temperatures. To reflect this, standard error should
decrease with sample size (i.e. number of temperatures). But this is not
obvious when we plot standard error as a function of temperatures.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

When we zoom in to exclude the outlier, there is still not a clear
relationship between standard error and the number of temperatures.
However, the largest standard errors are observed with the lowest sample
sizes (left-hand part of plot), which is reassuring.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

With less precision (higher standard errors), we would expect a broader
distribution of activation energy estimates. Let’s look at the
distribution of activation energies split by whether standard error was
high or low (above or below the median SE). Surprisingly, the estimates
with lower standard errors have a broader distribution, opposite to
expectations.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

The median of the distribution also seems to vary with standard error;
high standard error is associated with higher activation energy.

    ## # A tibble: 2 x 2
    ##   se_cat              median_AE
    ##   <chr>                   <dbl>
    ## 1 high standard error     0.721
    ## 2 low standard error      0.574

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

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

This indicates how these values are not independent. Since this is an
issue that only affects a few species, I’m not inclined to devise an
`experiment` random effect that accounts for this pseudoreplication
statistically. Rather, I would probably only take a single AE per
experiment. But which one?

We can zoom in on just these cases where there are multiple measurements
from the same experiment. If the metric was minimum development time,
activation energy is higher; if it was maximum development time,
activation energy is lower. This demonstrates that the developmental
metric impacts activation energy estimates.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

This also highlights the concerns Peter has raised about *apparent*
developmental rates, in that development is not measured for all viable
worms (ones that die before maturing are excluded). We can examine
whether activation energy varies systematically with developmental
metric. It does, but not conspicuously or in the same way as within an
experiment. For example, the ‘max’ group is higher than the ‘min’ group
overall, while this was reversed within experiments.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

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

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

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

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

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

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

The within-species effect, `binomial`, mixes well and is non-zero,
suggesting that multiple AE measurements on the same species tend to be
similar. The variance component for `std-error` is fixed at 1, such that
the weights assigned do not change as the posterior distribution is
sampled by the Markov chain.

The model intercept (the average activation energy) is about 0.66

    ## 
    ##  Iterations = 3001:22991
    ##  Thinning interval  = 10
    ##  Sample size  = 2000 
    ## 
    ##  DIC: -22.66692 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.03944  0.01802   0.0626     1451
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.03212  0.01901  0.04745     1672
    ## 
    ##  Location effects: E_Arr ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## (Intercept)    0.6509   0.6008   0.7079     2000 <5e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the proportion of variation attributable to the within-species
effect.

    ## [1] 0.5546039

It is high. On the one hand this makes sense, as many species only had a
single measurement (and hence no residual value after accounting for a
‘species effect’). Here’s an attempt to visualize the variation in
activation energy estimates for the species that had multiple
measurements. The species are ordered by their median activation energy.
One can see that some species tend to have very similar measurements,
while for other species, there is a fair amount of spread in the
estimated activation energies.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

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
effect (`tree_tips`).

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

These results suggest that phylogenetically related species tend to have
similar AEs, and when a species is measured multiple times, the
estimated AEs are similar, though they are not far more similar than
what they are expected to be based on phylogeny (i.e. the within-species
effect does not explain much additional variation beyond phylogeny).

The phylogenetic effect is strong. It explained 61% of the variation in
AE.

By contrast, the ‘species’ effect explained 10% of the variation.

Here’s a plot splitting AE by parasite order. It does look like some
clades tend to have higher or lower AEs.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

#### Developmental metric

As seen above, the chosen metric of development can affect activation
energy estimates, at least within experiments. This was less clear
across experiments. To check this effect, I added a fixed factor to the
model that distinguishes between activation estimates based on mean,
min, or unspecific developmental times.

These differences were rather minimal (non-significant p-values).
Moreover, the contrast between mean and min metrics was negative (min
with lower AE than mean), which was opposite to that seen within
experiments.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -44.79256 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01004 0.001954  0.02066     1078
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.07123  0.01829   0.1315     1101
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02779  0.01797  0.03919     1200
    ## 
    ##  Location effects: E_Arr ~ trait_detail3 
    ## 
    ##                            post.mean  l-95% CI  u-95% CI eff.samp  pMCMC
    ## (Intercept)                 0.660748  0.401555  0.941275     1200 <8e-04
    ## trait_detail3min           -0.002775 -0.096603  0.097183     1200  0.907
    ## trait_detail3not specified  0.116393 -0.015031  0.277114     1538  0.110
    ##                               
    ## (Intercept)                ***
    ## trait_detail3min              
    ## trait_detail3not specified    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Thus, it is not clear how valuable this term is in the model. It does
not “correct” activation energies in the direction we would expect.
Therefore, I have not included it in further analyses.

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

The effect of mean locale temperature is negative and significant.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -47.63916 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01233 0.002543  0.02566     1098
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.05237  0.01411   0.1072    951.3
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02761  0.01685  0.03819     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen 
    ## 
    ##                   post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
    ## (Intercept)        0.688645  0.483050  0.924100     1200 <8e-04 ***
    ## mean_ann_temp_cen -0.008826 -0.015781 -0.002149     1200  0.005 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

For species living at lower temperatures, development speeds up more as
temperature increases. It also does not look like this depends on
whether an exact study location was given.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

#### Temperature variation

A mean temperature does not indicate the temperature range a parasite is
exposed to. Thus, as the next modelling step, we’ll include location min
and max temps. None of these terms are significant, though the model DIC
decreases, suggesting together they explain some variation.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -51.97266 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01385 0.002804  0.02786     1030
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.05681  0.01028   0.1123     1260
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02572  0.01587  0.03705     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + max_month_temp_cen + min_month_temp_cen 
    ## 
    ##                    post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
    ## (Intercept)         0.708671  0.468021  0.991658     1230 <8e-04 ***
    ## mean_ann_temp_cen  -0.031247 -0.088675  0.032841     1200  0.325    
    ## max_month_temp_cen  0.017194 -0.016606  0.054760     1200  0.340    
    ## min_month_temp_cen  0.008869 -0.020086  0.038679     1328  0.550    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Interestingly, the parameter estimates for max and min temps are
positive, which is opposite to the trend seen with mean temp. But when
we plot AE as a function of max and min, we see a negative relationship
just like for mean temp: low temps, high AE.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

This is suggestive that, once controlling for mean temp, activation
energies increase with higher max and min temps. Let’s see if we can
uncover this pattern in the data. We’ll take the residuals from the
previous model, which just included mean temperature, and plot them
again max and min temp.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->
![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

We uncover the positive relationships, but they are clearly weak. Thus,
we are in a situation when one metric (DIC) suggests adding `max` and
`min` temps improve the model, while other metrics, like p-values and
plots do not indicate strong relationships. In cases like this, we need
to decide if there is a good biological reason to include extreme temps,
even if their effects are ambiguous. On the one hand, we could expect
that temperature variability affects activation energy. On the other
hand, a lower min temperature might not matter much if parasites do not
develop below a certain threshold that is unlinked to the minimum
temperature. Also, if temperature variability is what drives activation
energy, then we would expect the parameter estimates for min and max to
go in opposite directions, i.e. larger temperature extremes are
associated with higher (or lower) AE values. Since (i) these results are
ambiguous, (ii) they are confounded with mean temperature, and (iii) we
test several other proxies of seasonality/temperature variability, my
preference is to drop min and max temp from the model.

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

There is a marginal relationship between AE and latitude.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -46.99349 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial    0.0112 0.002069  0.02196     1094
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.04359 0.007695  0.08582    924.5
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02749  0.01788  0.03872     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + lat_abs 
    ## 
    ##                   post.mean  l-95% CI  u-95% CI eff.samp   pMCMC    
    ## (Intercept)        1.029303  0.655740  1.359615     1200 < 8e-04 ***
    ## mean_ann_temp_cen -0.020925 -0.034327 -0.010114     1200 0.00333 ** 
    ## lat_abs           -0.007384 -0.014038 -0.001374     1200 0.02167 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Surprisingly, the parameter estimate for latitude is negative, but when
we plot it, we see that activation energies tend to increase with
latitude.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

This suggests that, like for `min` and `max` temperature, correlations
among variables (here between mean temp and latitude) might cause
counterintuitive or spurious results. When we remove the effect of mean
temperature on activation energy (i.e. we take the residuals of the
model with just mean temp), and then plot the residual variation as a
function of latitude, we see essentially no relationship between
latitude and activation energy.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

Thus, even though the term was significant, the effect of latitude is
hard to interpret and not an obvious improvement to the model. Let’s
look at the next variable, which is related to latitude.

##### Climate zones

Our latitude variable is not a perfect representation of a parasite’s
‘ancestral’ habitat. Many species in the data are associated with
people or livestock and have thus been moved all around the world. This
is reflected in the distribution zone variable: it contains a category
for parasites with a broad distribution.

However, the distribution categories are still tightly linked to
latitude.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

Nonetheless, let’s add this variable to the model instead of latitude.

Most of the contrasts among these categories are not significant.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -48.28212 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial    0.0139 0.004014  0.02744     1092
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.04449 0.008599  0.09394    904.9
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02661  0.01689  0.03798     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone 
    ## 
    ##                      post.mean  l-95% CI  u-95% CI eff.samp   pMCMC    
    ## (Intercept)           0.673514  0.415604  0.921991   1200.0 < 8e-04 ***
    ## mean_ann_temp_cen    -0.013554 -0.022004 -0.004936   1200.0 0.00167 ** 
    ## dist_zonetemperate   -0.058038 -0.236731  0.127131    997.9 0.50167    
    ## dist_zonesubtropical  0.112318 -0.131427  0.358740   1200.0 0.35000    
    ## dist_zonetropical     0.131854 -0.163470  0.401075   1200.0 0.37500    
    ## dist_zoneglobal       0.047345 -0.118749  0.235358   1200.0 0.60167    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Consistent with the pattern for temperature and latitude, polar species
tend to have higher AE.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

Latitude and distribution overlap and presumably explain the same
variation in AE, but which explains it better? We can compare the models
with latitude vs distribution zone using information criteria. Here is
the progression of model fits thus far (low DIC, better model):

    ## [1] "DIC for model with phylogeny: -44"

    ## [1] "DIC for model with just mean temp: -47.6"

    ## [1] "DIC for model with latitude: -47"

    ## [1] "DIC for model with distribution zone: -48.3"

According to this measure, the model with distribution is better than
the one with latitude. One explanation for this, is that including
distribution zone strengthens the estimated effect of mean annual
temperature; it explains residual variation around the temp-AE
relationship.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

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
    ##  DIC: -49.90825 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01287 0.002076  0.02478    951.6
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips    0.0569 0.005641   0.1177    878.5
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02602  0.01617  0.03624     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone + habitat_d1 
    ## 
    ##                       post.mean  l-95% CI  u-95% CI eff.samp   pMCMC    
    ## (Intercept)            0.644400  0.355889  0.947751     1200 0.00167 ** 
    ## mean_ann_temp_cen     -0.013992 -0.022700 -0.005731     1200 < 8e-04 ***
    ## dist_zonetemperate    -0.052732 -0.233651  0.122878     1091 0.55833    
    ## dist_zonesubtropical   0.112776 -0.163163  0.343743     1075 0.37500    
    ## dist_zonetropical      0.143146 -0.122848  0.408720     1200 0.33500    
    ## dist_zoneglobal        0.036191 -0.144671  0.205379     1107 0.69167    
    ## habitat_d1terrestrial  0.088498 -0.110212  0.294509     1099 0.39333    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Terrestrial species have higher AE on average than aquatic species
(combines both freshwater and marine species). This plot also
demonstrates how weights (standard errors) reduced the estimated effect
of habitat. The raw difference in AE between habitats was about 0.11,
but the model estimated it to be lower, 0.07.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

We might expect the aquatic - terrestrial dichotomy to be particularly
important in more seasonal locations. That is, there may be an
interaction between our seasonality proxy (distribution) and habitat.
Unfortunately, this cuts the data quite thin, as there are not many
freshwater and terrestrial species in each distribution category.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

When we add this interaction to the model, it is not an improvement.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -48.51064 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01396 0.002892  0.02823     1047
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.05603 0.007263   0.1172    753.7
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02581  0.01522  0.03612     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp + dist_zone * habitat_d1 
    ## 
    ##                                          post.mean  l-95% CI  u-95% CI
    ## (Intercept)                               0.825902  0.448981  1.114126
    ## mean_ann_temp                            -0.014988 -0.024652 -0.007075
    ## dist_zonetemperate                       -0.125854 -0.396295  0.153194
    ## dist_zonesubtropical                      0.181740 -0.113907  0.461167
    ## dist_zonetropical                         0.207935 -0.176665  0.625625
    ## dist_zoneglobal                          -0.029729 -0.287784  0.241044
    ## habitat_d1terrestrial                    -0.010996 -0.376541  0.317489
    ## dist_zonetemperate:habitat_d1terrestrial  0.146181 -0.216783  0.513023
    ## dist_zonetropical:habitat_d1terrestrial  -0.061943 -0.528063  0.428981
    ## dist_zoneglobal:habitat_d1terrestrial     0.130850 -0.211289  0.500786
    ##                                          eff.samp   pMCMC    
    ## (Intercept)                                 975.8 < 8e-04 ***
    ## mean_ann_temp                              1200.0 0.00167 ** 
    ## dist_zonetemperate                         1021.0 0.37833    
    ## dist_zonesubtropical                       1047.4 0.21833    
    ## dist_zonetropical                          1055.0 0.30333    
    ## dist_zoneglobal                            1200.0 0.79667    
    ## habitat_d1terrestrial                      1043.3 0.94833    
    ## dist_zonetemperate:habitat_d1terrestrial   1350.5 0.42000    
    ## dist_zonetropical:habitat_d1terrestrial    1262.7 0.80333    
    ## dist_zoneglobal:habitat_d1terrestrial      1771.2 0.46000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

To summarize, AE tends to be higher with low mean locale temperatures.
This happens at high latitudes, but this effect is not independent of
temperature.

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
    ##  DIC: -50.41971 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01373 0.002072  0.02787    633.3
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.04821 0.006711   0.1059    927.3
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02574   0.0169  0.03636     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone + habitat_d1 + plant_anim 
    ## 
    ##                       post.mean  l-95% CI  u-95% CI eff.samp   pMCMC    
    ## (Intercept)            0.657965  0.411077  0.936679     1200 0.00167 ** 
    ## mean_ann_temp_cen     -0.014114 -0.023138 -0.006064     1200 < 8e-04 ***
    ## dist_zonetemperate    -0.046449 -0.205130  0.133007     1200 0.55167    
    ## dist_zonesubtropical   0.144467 -0.083823  0.388000     1207 0.22333    
    ## dist_zonetropical      0.145581 -0.126935  0.413505     1200 0.31500    
    ## dist_zoneglobal        0.030420 -0.129098  0.223659     1200 0.71000    
    ## habitat_d1terrestrial  0.120660 -0.086550  0.310457     1200 0.19333    
    ## plant_animplant       -0.232284 -0.502326  0.035679     1200 0.10000 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the plot comparing plant and animal parasites.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-86-1.png)<!-- -->

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

Again, model DIC decreased, but the term contrasting animal parasites
with and without endotherms in the cycle was not significant (0.1). It
suggests species without endotherms have higher AEs.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -50.81563 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01454 0.002997  0.02858     1212
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.04752 0.005304   0.1041    757.3
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.0257  0.01605  0.03721     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 
    ## 
    ##                                      post.mean  l-95% CI  u-95% CI
    ## (Intercept)                           0.623763  0.336232  0.850788
    ## mean_ann_temp_cen                    -0.015189 -0.023348 -0.006195
    ## dist_zonetemperate                   -0.045680 -0.216683  0.113758
    ## dist_zonesubtropical                  0.163947 -0.073725  0.418963
    ## dist_zonetropical                     0.155336 -0.109802  0.435872
    ## dist_zoneglobal                       0.049105 -0.128906  0.212301
    ## habitat_d1terrestrial                 0.122355 -0.067562  0.339597
    ## endo_ecto2no endotherm in life cycle  0.083763 -0.063926  0.224321
    ## endo_ecto2plant                      -0.185874 -0.466608  0.092841
    ##                                      eff.samp  pMCMC    
    ## (Intercept)                            1200.0 <8e-04 ***
    ## mean_ann_temp_cen                      1200.0 <8e-04 ***
    ## dist_zonetemperate                     1200.0  0.622    
    ## dist_zonesubtropical                   1200.0  0.198    
    ## dist_zonetropical                      1200.0  0.278    
    ## dist_zoneglobal                        1200.0  0.565    
    ## habitat_d1terrestrial                  1200.0  0.232    
    ## endo_ecto2no endotherm in life cycle   1200.0  0.247    
    ## endo_ecto2plant                         964.9  0.198    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

However, this trend is not visible in the raw data, suggesting it is
contingent on the covariates in the data.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

When we plot these groups as a function of mean annual temperature, we
then see that, at a given temperature, species exclusively infecting
ectotherms have higher AE. This pattern is not strong, but it is
consistent with the idea that environmental stability lowers AE
(endotherms provide a stable temperature for helminths).

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->

#### Stage in or out of host?

Stages inside a host may be more shielded from temperature variation
than stages in the environment. Let’s add this distinction to the model.

This term is not significant and the model DIC is worse.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -49.82448 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01503 0.003221  0.02988    726.7
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.04626  0.00441   0.1048    710.5
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02573  0.01717   0.0371     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 + devo_type 
    ## 
    ##                                      post.mean  l-95% CI  u-95% CI
    ## (Intercept)                           0.609292  0.337560  0.888309
    ## mean_ann_temp_cen                    -0.014994 -0.024497 -0.006557
    ## dist_zonetemperate                   -0.042584 -0.212505  0.139692
    ## dist_zonesubtropical                  0.160288 -0.091369  0.417431
    ## dist_zonetropical                     0.161607 -0.116727  0.460418
    ## dist_zoneglobal                       0.052063 -0.141719  0.230135
    ## habitat_d1terrestrial                 0.118062 -0.072271  0.323096
    ## endo_ecto2no endotherm in life cycle  0.088953 -0.070787  0.222501
    ## endo_ecto2plant                      -0.181504 -0.487741  0.080050
    ## devo_typeoutside host                 0.025646 -0.071423  0.139446
    ##                                      eff.samp  pMCMC    
    ## (Intercept)                              1200 <8e-04 ***
    ## mean_ann_temp_cen                        1414  0.005 ** 
    ## dist_zonetemperate                       1200  0.618    
    ## dist_zonesubtropical                     1200  0.223    
    ## dist_zonetropical                        1200  0.242    
    ## dist_zoneglobal                          1200  0.560    
    ## habitat_d1terrestrial                    1003  0.233    
    ## endo_ecto2no endotherm in life cycle     1212  0.227    
    ## endo_ecto2plant                          1200  0.190    
    ## devo_typeoutside host                    1200  0.632    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Though the difference is not significant, worms outside the host have
slightly higher AEs, which is also consistent with the idea that
environmental variation favors higher AE. The plot also shows that
habitat and in/outside host are conflated, as proportionally more of the
terrestrial values are outside the host and more of the aquatic values
are inside the host.

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-97-1.png)<!-- -->

#### Target host

Stages free in the environment are generally propagules that target the
next host in the life cycle. Are the parasites targeting a vertebrate or
an invertebrate? It is not obvious to me that this should matter. But
nonetheless, we can check it.

Distinguishing between invertebrates and vertebrates as next host does
not improve the model.

    ## 
    ##  Iterations = 3001:62951
    ##  Thinning interval  = 50
    ##  Sample size  = 1200 
    ## 
    ##  DIC: -49.58007 
    ## 
    ##  G-structure:  ~binomial
    ## 
    ##          post.mean l-95% CI u-95% CI eff.samp
    ## binomial   0.01367 0.002217  0.02686     1006
    ## 
    ##                ~tree_tips
    ## 
    ##           post.mean l-95% CI u-95% CI eff.samp
    ## tree_tips   0.03903 0.002796  0.09519      775
    ## 
    ##                ~idh(stderr_arr):units
    ## 
    ##                  post.mean l-95% CI u-95% CI eff.samp
    ## stderr_arr.units         1        1        1        0
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units   0.02595  0.01672  0.03672     1200
    ## 
    ##  Location effects: E_Arr ~ mean_ann_temp_cen + dist_zone + habitat_d1 + endo_ecto2 + devo_type3 
    ## 
    ##                                             post.mean   l-95% CI
    ## (Intercept)                                 5.964e-01  3.426e-01
    ## mean_ann_temp_cen                          -1.265e-02 -2.078e-02
    ## dist_zonetemperate                         -2.625e-02 -1.992e-01
    ## dist_zonesubtropical                        1.172e-01 -1.451e-01
    ## dist_zonetropical                           1.436e-01 -1.475e-01
    ## dist_zoneglobal                             6.098e-02 -9.396e-02
    ## habitat_d1terrestrial                       1.940e-01 -6.174e-05
    ## endo_ecto2no endotherm in life cycle        9.115e-03 -1.362e-01
    ## endo_ecto2plant                            -2.169e-01 -4.965e-01
    ## devo_type3in plant                         -8.866e-02 -2.797e-01
    ## devo_type3outside host: targets invert      1.139e-01 -5.200e-02
    ## devo_type3outside host: targets vertebrate -1.748e-01 -3.577e-01
    ##                                              u-95% CI eff.samp  pMCMC    
    ## (Intercept)                                 8.855e-01     1200 <8e-04 ***
    ## mean_ann_temp_cen                          -3.467e-03     1200 0.0100 ** 
    ## dist_zonetemperate                          1.337e-01     1200 0.7767    
    ## dist_zonesubtropical                        3.387e-01     1200 0.3533    
    ## dist_zonetropical                           4.266e-01     1200 0.3283    
    ## dist_zoneglobal                             2.356e-01     1200 0.4883    
    ## habitat_d1terrestrial                       3.916e-01     1056 0.0483 *  
    ## endo_ecto2no endotherm in life cycle        1.625e-01     1294 0.9267    
    ## endo_ecto2plant                             6.357e-02     1428 0.1167    
    ## devo_type3in plant                          1.238e-01     1347 0.3883    
    ## devo_type3outside host: targets invert      2.734e-01     1049 0.1750    
    ## devo_type3outside host: targets vertebrate  3.278e-02     1312 0.0867 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![](06_fit_phylo_models_Arr_Nov4_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->

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

| model                                 | df | df\_used |     DIC |  R2m |  R2c |
| :------------------------------------ | -: | -------: | ------: | ---: | ---: |
| just within-species\*                 |  3 |       NA | \-22.67 | 0.00 | 0.53 |
| \+ phylogeny\*                        |  4 |        1 | \-43.99 | 0.00 | 0.70 |
| \+ mean annual temp\*                 |  5 |        1 | \-47.64 | 0.02 | 0.73 |
| \+ min/max temps                      |  7 |        2 | \-51.97 | 0.05 | 0.74 |
| \+ latitude                           |  6 |      \-1 | \-46.99 | 0.05 | 0.73 |
| \+ distribution zone\*                |  9 |        3 | \-48.28 | 0.08 | 0.72 |
| \+ habitat (aquatic vs terrestrial)\* | 10 |        1 | \-49.91 | 0.12 | 0.79 |
| \+ distribution x habitat interacton  | 13 |        3 | \-48.51 | 0.12 | 0.78 |
| \+ plant vs animal parasite\*         | 11 |      \-2 | \-50.42 | 0.18 | 0.75 |
| \+ endotherm in life cycle?\*         | 12 |        1 | \-50.82 | 0.23 | 0.74 |
| \+ stage in/out of host\*             | 13 |        1 | \-49.82 | 0.19 | 0.77 |
| \+ target host: invert vs vert        | 15 |        2 | \-49.58 | 0.24 | 0.74 |

The most complex model, one that includes phylogeny, seemingly important
environmental variables, and characteristics of the host-parasite
system, explained over 77% of the variation in activation energy. Much
of this is due to the random effects (phylogeny and repeated measures
account for \>50% of the variation), but around 20% of the variation can
be explained by fixed effects.
