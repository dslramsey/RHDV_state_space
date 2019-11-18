Emerging RHDV2 suppresses the impact of endemic and novel strains of
RHDV on wild rabbit populations
================

## Notes

This repository contains data and code from:

Ramsey, D.S.L., Cox, T., Strive, T., Forsyth, D.M., Stuart, I., Hall,
R., Elsworth, P., and Campbell, S. (2019). “Emerging RHDV2 suppresses
the impact of endemic and novel strains of RHDV on wild rabbit
populations” *Journal of Applied
Ecology*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3517266.svg)](https://doi.org/10.5281/zenodo.3517266)

## Getting started

**File descriptions**:

*R/.* – R code used to process the spotlight and seroprevalence data and
fit Stan models.

*Data/.* – Contains csv files for the rabbit spotlight and serology
datasets.

*Stan/.* – Stan code for the rabbit N-mixture state-space model as well
as the MAR(1) and age-specific trend models fitted to the seroprevalence
data for the three strains.

## Prerequisites

The R scripts require packages *tidyverse*, *lubridate*, *gridExtra*,
*rstan* and *MCMCvis*.
