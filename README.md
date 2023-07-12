
<!-- README.md is generated from README.Rmd. Please edit that file -->

# peathamstr

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

‘peathamstr’ extents the
[‘hamstr’](https://github.com/EarthSystemDiagnostics/hamstr) package
(Hierarchical Accumulation Modelling with Stan and R) (Dolman 2023) with
a modified version of Clymo’s Bog Growth Model (Clymo 1984).

The extension uses ‘hamstr’ to compute an age-depth model, but also
estimates layer-specific peat addition rates to the catotelm (PAR) and a
global exponential decomposition rate. In addition, peat mass fluxes net
carbon uptake (NCU) and net carbon release (NCR) as defined in Yu (2011)
are computed.

This allows to estimate catotelm PAR and decomposition rates also for
peat cores where the assumption of a constant PAR of Clymo’s Bog Growth
Model (Clymo 1984) is violated, e.g. when PAR declined over time (see
for example Yu et al. (2003)).

The package shamelessly recycles much of the code from the
[‘hamstr’](https://github.com/EarthSystemDiagnostics/hamstr) package and
much credit goes to its developers. When using ‘peathamstr’, do not
forget to also cite the ‘hamstr’ package!

## Installation

You can install the development version of ‘peathamstr’ like so:

``` r
remotes::install_github("henningte/peathamstr")
```

## Example

This example shows how to estimate the decomposition rate, net carbon
uptake (NCU) and net carbon release (NCR) for a peat core with
non-constant peat addition rate (PAR).

``` r
library(peathamstr)
library(pangaear) # to download example data
#> Warning: package 'pangaear' was built under R version 4.2.3
#> Registered S3 method overwritten by 'httr':
#>   method           from  
#>   print.cache_info hoardr
library(ggplot2)
```

Download and reshape data (Yu 2018b, 2018a) from the Pangaea database:

``` r
d1 <- pangaear::pg_data(doi = '10.1594/PANGAEA.890357')
#> Downloading 1 datasets from 10.1594/PANGAEA.890357
#> Processing 1 files
d2 <- pangaear::pg_data(doi = '10.1594/PANGAEA.890392')
#> Downloading 1 datasets from 10.1594/PANGAEA.890392
#> Processing 1 files

d <- 
  dplyr::full_join(
    d1[[1]]$data |>
      dplyr::select(1, 3, 6) |>
      setNames(nm = c("depth_midpoint", "bulk_density", "C")) |>
      dplyr::filter(!is.na(bulk_density)) |>
      dplyr::mutate(
        C = C/bulk_density,
        depth_midpoint = round(depth_midpoint * 100, 1), # assumed
        thickness = 1, # assumed
        mass = bulk_density * thickness * 10000/1000,
        cumulative_mass = c(0.00001, cumsum(mass)[-1])
      ),
    d2[[1]]$data |>
      dplyr::select(2, 4, 5) |>
      setNames(nm = c("depth_midpoint", "age", "age_sd")) |>
      dplyr::filter(!is.na(age)) |>
      dplyr::mutate(
        depth_midpoint = round(depth_midpoint * 100, 1),
        age = age * 1000,
        age_sd = age_sd * 1000
      ) |>
      hamstr::calibrate_14C_age(age.14C = "age", age.14C.se = "age_sd") |>
      dplyr::select(-age, -age_sd) |>
      dplyr::rename(
        age = "age.14C.cal",
        age_sd = "age.14C.cal.se"
      ),
    by = "depth_midpoint"
  ) |>
  dplyr::arrange(depth_midpoint) |>
  dplyr::filter(depth_midpoint <= max(depth_midpoint[!is.na(age)]) & depth_midpoint >= 44) |>
  dplyr::mutate(
    depth_upper = depth_midpoint - 0.5,
    depth_lower = depth_midpoint + 0.5
  )

# subsample the data to make the example faster to compute
set.seed(345345)
index <- c(which(!is.na(d$age)), which(is.na(d$age)) |> sample(size = 150, replace = FALSE)) |> sort()

d <- 
  d |>
  dplyr::slice(index) |>
  dplyr::mutate(
    depth_upper = c(depth_upper[[1]], depth_lower[-length(depth_lower)]) 
  )
```

Show the cumulative mass-age curve (here using only point age values and
cumulative masses for dated layers). It is clearly visible that the core
violates the assumptions of the original Bog Growth Model (see also Yu
et al. (2003) for a detailed description of the peatland).

``` r
d |>
  dplyr::filter(!is.na(age)) |>
  ggplot(aes(y = cumsum(mass * C), x = age)) +
  geom_path() +
  geom_point() +
  labs(y = expression("Cumulative carbon mass (kg m"^{-2}*")"), x = "age (yr BP)")
```

<img src="man/figures/README-example-3-1.png" width="60%" />

Estimate the posterior distribution (computes the age-depth model and
fits the modified Clymo model with PAR varying over time).

``` r
fit_1 <- 
  peat_hamstr(
    depth = 
      d |> 
      dplyr::filter(!is.na(age)) |> 
      dplyr::pull(depth_midpoint),
    obs_age =
      d |> 
      dplyr::filter(!is.na(age)) |> 
      dplyr::pull(age),
    obs_err = 
      d |> 
      dplyr::filter(!is.na(age)) |> 
      dplyr::pull(age_sd),
    cumulative_mass = 
      d |> 
      dplyr::filter(!is.na(cumulative_mass)) |> 
      dplyr::slice(-1) |> 
      dplyr::pull(cumulative_mass),
    cumulative_mass0 = d$cumulative_mass[[1]],
    depth2 =
      d |> 
      dplyr::filter(!is.na(cumulative_mass)) |> 
      dplyr::slice(-1) |> 
      dplyr::pull(depth_midpoint),
    depth2_upper =
      d |> 
      dplyr::filter(!is.na(cumulative_mass)) |> 
      dplyr::slice(-1) |> 
      dplyr::pull(depth_upper),
    depth2_lower =
      d |> 
      dplyr::filter(!is.na(cumulative_mass)) |> 
      dplyr::slice(-1) |> 
      dplyr::pull(depth_lower),
    min_age = 0,
    # the seed argument for the sampler is set here so that
    # this example always returns the same numerical result
    stan_sampler_args = 
      list(
        seed = 34564, 
        chains = 4, 
        iter = 4000, 
        cores = 4, 
        control = list(max_treedepth = 12)
      )
  )
#> Warning in validityMethod(object): The following variables have undefined
#> values: nmu_rep[126],The following variables have undefined values:
#> nmu_rep[127],The following variables have undefined values: nmu_rep[128],The
#> following variables have undefined values: nmu_rep[129],The following
#> variables have undefined values: nmu_rep[130],The following variables have
#> undefined values: nmu_rep[131],The following variables have undefined values:
#> nmu_rep[132],The following variables have undefined values: nmu_rep[133],The
#> following variables have undefined values: nmu_rep[134],The following
#> variables have undefined values: nmu_rep[135],The following variables have
#> undefined values: nmu_rep[136],The following variables have undefined values:
#> nmu_rep[137],The following variables have undefined values: nmu_rep[138],The
#> following variables have undefined values: nmu_rep[139],The following
#> variables have undefined values: nmu_rep[140],The following variables have
#> undefined values: nmu_rep[141],The following variables have undefined values:
#> nmu_rep[142],The following variables have undefined values: nmu_rep[143],The
#> following variables have undefined values: nmu_rep[144],The following
#> variables have undefined values: nmu_rep[145],The following variables have
#> undefined values: nmu_rep[146],The following variables have undefined values:
#> nmu_rep[147],The following variables have undefined values: nmu_rep[148],The
#> following variables have undefined values: nmu_rep[149],The following
#> variables have undefined values: nmu_rep[150],The following variables have
#> undefined values: nmu_rep[151],The following variables have undefined values:
#> nmu_rep[152],The following variables have undefined values: nmu_rep[153],The
#> following variables have undefined values: nmu_rep[154],The following
#> variables have undefined values: nmu_rep[155],The following variables have
#> undefined values: nmu_rep[156],The following variables have undefined values:
#> nmu_rep[157],The following variables have undefined values: nmu_rep[158],The
#> following variables have undefined values: nmu_rep[159],The following
#> variables have undefined values: nmu_rep[160],The following variables have
#> undefined values: nmu_rep[161],The following variables have undefined values:
#> nmu_rep[162],The following variables have undefined values: nmu_rep[163],The
#> following variables have undefined values: nmu_rep[164],The following
#> variables have undefined values: nmu_rep[165],The following variables have
#> undefined values: nmu_rep[166],The following variables have undefined values:
#> nmu_rep[167],The following variables have undefined values: nmu_rep[168],The
#> following variables have undefined values: nmu_rep[169],The following
#> variables have undefined values: nmu_rep[170],The following variables have
#> undefined values: nmu_rep[171],The following variables have undefined values:
#> nmu_rep[172],The following variables have undefined values: nmu_rep[173],The
#> following variables have undefined values: nmu_rep[174],The following
#> variables have undefined values: nmu_rep[175],The following variables have
#> undefined values: nmu_rep[176],The following variables have undefined values:
#> nmu_rep[177],The following variables have undefined values: nmu_rep[178],The
#> following variables have undefined values: nmu_rep[179],The following
#> variables have undefined values: nmu_rep[180],The following variables have
#> undefined values: nmu_rep[181],The following variables have undefined values:
#> nmu_rep[182],The following variables have undefined values: nmu_rep[183],The
#> following variables have undefined values: nmu_rep[184],The following
#> variables have undefined values: nmu_rep[185],The following variables have
#> undefined values: nmu_rep[186],The following variables have undefined values:
#> nmu_rep[187],The following variables have undefined values: nmu_rep[188],The
#> following variables have undefined values: nmu_rep[189],The following
#> variables have undefined values: nmu_rep[190],The following variables have
#> undefined values: nmu_rep[191],The following variables have undefined values:
#> nmu_rep[192],The following variables have undefined values: nmu_rep[193],The
#> following variables have undefined values: nmu_rep[194],The following
#> variables have undefined values: nmu_rep[195],The following variables have
#> undefined values: nmu_rep[196],The following variables have undefined values:
#> nmu_rep[197],The following variables have undefined values: nmu_rep[198],The
#> following variables have undefined values: nmu_rep[199],The following
#> variables have undefined values: nmu_rep[200],The following variables have
#> undefined values: nmu_rep[201],The following variables have undefined values:
#> nmu_rep[202],The following variables have undefined values: nmu_rep[203],The
#> following variables have undefined values: nmu_rep[204],The following
#> variables have undefined values: nmu_rep[205],The following variables have
#> undefined values: nmu_rep[206],The following variables have undefined values:
#> nmu_rep[207],The following variables have undefined values: nmu_rep[208],The
#> following variables have undefined values: nmu_rep[209],The following
#> variables have undefined values: nmu_rep[210],The following variables have
#> undefined values: nmu_rep[211],The following variables have undefined values:
#> nmu_rep[212],The following variables have undefined values: nmu_rep[213],The
#> following variables have undefined values: nmu_rep[214],The following
#> variables have undefined values: nmu_rep[215],The following variables have
#> undefined values: nmu_rep[216],The following variables have undefined values:
#> nmu_rep[217],The following variables have undefined values: nmu_rep[218],The
#> following variables have undefined values: nmu_rep[219],The following
#> variables have undefined values: nmu_rep[220],The following variables have
#> undefined values: nmu_rep[221],The following variables have undefined values:
#> nmu_rep[222],The following variables have undefined values: nmu_rep[223],The
#> following variables have undefined values: nmu_rep[224],The following
#> variables have undefined values: nmu_rep[225],The following variables have
#> undefined values: nmu_rep[226],The following variables have undefined values:
#> nmu_rep[227],The following variables have undefined values: nmu_rep[228],The
#> following variables have undefined values: nmu_rep[229],The following
#> variables have undefined values: nmu_rep[230],The following variables have
#> undefined values: nmu_rep[231],The following variables have undefined values:
#> nmu_rep[232],The following variables have undefined values: nmu_rep[233],The
#> following variables have undefined values: nmu_rep[234],The following
#> variables have undefined values: nmu_rep[235],The following variables have
#> undefined values: nmu_rep[236],The following variables have undefined values:
#> nmu_rep[237],The following variables have undefined values: nmu_rep[238],The
#> following variables have undefined values: nmu_rep[239],The following
#> variables have undefined values: nmu_rep[240],The following variables have
#> undefined values: nmu_rep[241],The following variables have undefined values:
#> nmu_rep[242],The following variables have undefined values: nmu_rep[243],The
#> following variables have undefined values: nmu_rep[244],The following
#> variables have undefined values: nmu_rep[245],The following variables have
#> undefined values: nmu_rep[246],The following variables have undefined values:
#> nmu_rep[247],The following variables have undefined values: nmu_rep[248],The
#> following variables have undefined values: nmu_rep[249],The following
#> variables have undefined values: nmu_rep[250],The following variables have
#> undefined values: nmu_rep[251],The following variables have undefined values:
#> nmu_rep[252],The following variables have undefined values: nmu_rep[253],The
#> following variables have undefined values: nmu_rep[254],The following
#> variables have undefined values: nmu_rep[255],The following variables have
#> undefined values: nmu_rep[256],The following variables have undefined values:
#> nmr_rep[134],The following variables have undefined values: nmr_rep[135],The
#> following variables have undefined values: nmr_rep[136],The following va
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
```

Show the age-depth model (this is the same function as in the original
‘hamstr’ package):

``` r
plot(fit_1, type = "default")
```

<img src="man/figures/README-example-5-1.png" width="100%" />

Show estimated and measured cumulative masses versus depth:

``` r
# plot(fit_1, type = "cumulative_mass_profile")
plot(
  fit_1, 
  type = "cumulative_carbon_mass_profile", 
  carbon_content = 
    d |> 
    dplyr::filter(!is.na(cumulative_mass)) |> 
    dplyr::slice(-1) |> 
    dplyr::pull(C)
)
#> Warning in get_posterior_cumulative_carbon_mass(object = object, depth =
#> depth, : * Cumulative carbon contents start at 0 for the first measured layer.
#> Warning in get_posterior_cumulative_carbon_mass(object = object, depth =
#> depth, : * Measured layers do not cover completely modeled layers (are not
#> contiguous). This results in `NA` carbon contents. `NA` carbon contents are
#> filled in order to compute cumulative carbon masses.
```

<img src="man/figures/README-example-6-1.png" width="60%" style="display: block; margin: auto;" />

Plot a histogram of the estimate exponential decomposition rate:

``` r
as.data.frame(fit_1$fit, pars = "clymo_alpha") |>
  ggplot(aes(x = clymo_alpha)) +
  geom_histogram(bins = 30) +
  labs(y = "Count", x = expression(alpha~"("*yr^{-1}*")"))
```

<img src="man/figures/README-example-7-1.png" width="60%" style="display: block; margin: auto;" />

Plot net carbon uptake (NCU) and net carbon release (NCR). Gaps in
default plots are due to incomplete coverage of modeled depth layers by
measured depth layers:

``` r
p1 <- 
  plot(
  fit_1, 
  type = "carbon_fluxes", 
  carbon_content = 
    d |> 
    dplyr::filter(!is.na(cumulative_mass)) |> 
    dplyr::slice(-1) |> 
    dplyr::pull(C)
) + 
  facet_wrap(~ variable, ncol = 1L)

p1
```

<img src="man/figures/README-example-8-1.png" width="60%" style="display: block; margin: auto;" />

To add apparent carbon accumulation rates (aCAR), you can add:

``` r
d_acar <- 
  predict(
        object = fit_1, 
        type = "apparent_carbon_accumulation_rates", 
        depth = "data", 
        carbon_content = d |> 
          dplyr::filter(!is.na(cumulative_mass)) |> 
          dplyr::slice(-1) |> 
          dplyr::pull(C)
      ) |>
  dplyr::group_by(depth_lower) |>
  dplyr::summarise(
    mean = mean(acar, na.rm = TRUE),
    `2.5%` = quantile(acar, probs = 0.025, na.rm = TRUE),
    `97.5%` = quantile(acar, probs = 0.975, na.rm = TRUE),
    age = mean(age_lower, na.rm = TRUE),
    .groups = "drop"
  )

p1 +
  geom_ribbon(
    data = d_acar, 
    aes(x = age, ymin = `2.5%`, ymax = `97.5%`), 
    color = NA, fill = "grey", alpha = 0.3
  ) +
  geom_path(data = d_acar, aes(x = age, y = mean))
#> Warning: Removed 1 row(s) containing missing values (geom_path).
```

<img src="man/figures/README-example-9-1.png" width="60%" style="display: block; margin: auto;" />

Plot for net mass uptake (NMU) and net mass release (NMR) can also be
created:

``` r
plot(fit_1, type = "mass_fluxes") + 
  facet_wrap(~ variable, ncol = 1L)
```

<img src="man/figures/README-example-10-1.png" width="60%" style="display: block; margin: auto;" />

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Clymo.1984" class="csl-entry">

Clymo, R S. 1984. “The Limits to Peat Bog Growth.” *Philosophical
Transactions of the Royal Society of London. Series B, Biological
Sciences* 303 (1117): 51.

</div>

<div id="ref-Dolman.2023" class="csl-entry">

Dolman, Andrew. 2023. “<span class="nocase">hamstr</span>: Hierarchical
Accumulation Modelling with Stan and R.”

</div>

<div id="ref-Yu.2011" class="csl-entry">

Yu, Zicheng. 2011. “Holocene Carbon Flux Histories of the World’s
Peatlands: Global Carbon-Cycle Implications.” *The Holocene* 21 (5):
761–74. <https://doi.org/10.1177/0959683610386982>.

</div>

<div id="ref-Yu.2018a" class="csl-entry">

———. 2018a. “Age Determination of Upper\_Pinto Peat Core.” PANGAEA.
<https://doi.org/10.1594/PANGAEA.890392>.

</div>

<div id="ref-Yu.2018" class="csl-entry">

———. 2018b. “Geochemistry of Upper\_Pinto Peat Core.” PANGAEA.
<https://doi.org/10.1594/PANGAEA.890357>.

</div>

<div id="ref-Yu.2003" class="csl-entry">

Yu, Zicheng, Dale H Vitt, Ian D Campbell, and Michael J Apps. 2003.
“Understanding Holocene Peat Accumulation Pattern of Continental Fens in
Western Canada.” *Can. J. Bot.* 81: 16.

</div>

</div>
