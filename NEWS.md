# peathamstr 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

## New models

* Add a new model version which estimates a separate decomposition rate for the acrotelm and catotelm and the time a newly formed layer spends in the acrotelm.

## Changes
* Adapt the functions to this new model.

## Bug fixes
* Remove the distinction between C mass fluxes and mass fluxes because this makes no sense: If we want C fluxes, we need to compute the model with C masses.
