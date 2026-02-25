# `tergm`: Fit, Simulate and Diagnose Models for Network Evolution Based on Exponential-Family Random Graph Models

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/tergm?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/tergm)](https://cran.r-project.org/package=tergm)
[![Coverage status](https://codecov.io/gh/statnet/tergm/branch/master/graph/badge.svg)](https://codecov.io/github/statnet/tergm?branch=master)
[![R build status](https://github.com/statnet/tergm/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/tergm/actions)

An integrated set of extensions to the 'ergm' package to analyze and simulate dynamic temporal networks based on exponential-family random graph models (ERGM). 'tergm' is a part of the Statnet suite of packages for network analysis.

TERGMs are a broad, flexible class of models for representing temporal edge dependence structure and dynamics.  They extend the generalized linear model framework by relaxing the assumption that observations (dyads) are independent.  Dependence is explicitly represented by model terms (e.g. mutuality, degree distributions, types of triads, etc) that can be fit along with dyad-independent terms (e.g., homophily).  Statistical inference supports significance testing for the model and parameter estimates.  `tergm` can be used for both estimation from and simulation of dynamic networks -- both workflows rely on Statnet's central MCMC algorithm, which can be tuned with a detailed set of control parameters.

Temporal network data come in different forms, continuous and discrete, completely observed and sampled.  The tergm package is designed to work with binary or continuously valued edges for:

* Network panel data 
* One cross-sectional network with continuous edge duration information
* One cross-sectional, egocentrically sampled network with continuous edge duration

Tools are included for implementing a standard statistical workflow:

* Exploratory data analysis -- via the related Statnet packages [`tsna`](https://github.com/statnet/tsna/) and [`ndtv`](https://github.com/statnet/ndtv/)
* Model estimation (data-type specific)
* Model diagnostics for assessing convergence and "Goodness of Fit"
* Simulating dynamic networks (from fitted or hypothetical models)

## Docs and examples

A brief introduction to the TERGM statistical framework, along with many examples of how `tergm` tools might be used in a network data analysis workflow can be found in our self-guided workshop materials: [Temporal Exponential Random Graph Models (TERGMs) for dynamic networks](https://statnet.org/workshop-tergm/).

Note that with the release of `tergm` 4.0, the API has been modified to allow for the specification of a much wider class of models.  A crosswalk of the differences can be found in the [`tergm` Conversion vignette](https://cran.r-project.org/web/packages/tergm/vignettes/tergm4_conversion.html) 

## Citation

Citation information (with bibtex entries) can be found [here](https://cran.r-project.org/web/packages/tergm/citation.html)

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/tergm`
* A private repository `statnet/tergm-private`

The intention is that all developments in `statnet/tergm-private` will eventually make their way into `statnet/tergm` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.github.io/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

A set of binaries is built after every commit to the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. They can be downloaded through the following links:

* [MacOS binary (a `.tgz` file in a `.zip` file)](https://nightly.link/statnet/tergm/workflows/R-CMD-check.yaml/master/macOS-rrelease-binaries.zip)
* [Windows binary (a `.zip` file in a `.zip` file)](https://nightly.link/statnet/tergm/workflows/R-CMD-check.yaml/master/Windows-rrelease-binaries.zip)

You will need to extract the MacOS `.tgz` or the Windows `.zip` file from the outer `.zip` file before installing. These binaries are usually built under the latest version of R and their operating system and may not work under other versions.

You may also want to install the corresponding latest binaries for packages on which `tergm` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`network`](https://github.com/statnet/network), [`networkDynamic`](https://github.com/statnet/networkDynamic), [`ergm`](https://github.com/statnet/ergm), and [`ergm.multi`](https://github.com/statnet/ergm.multi).
