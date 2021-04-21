# `tergm`: Fit, Simulate and Diagnose Models for Network Evolution Based on Exponential-Family Random Graph Models

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/tergm?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/tergm)](https://cran.r-project.org/package=tergm)
[![Coverage status](https://codecov.io/gh/statnet/tergm/branch/master/graph/badge.svg)](https://codecov.io/github/statnet/tergm?branch=master)
[![R build status](https://github.com/statnet/tergm/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/tergm/actions)

An integrated set of extensions to the 'ergm' package to analyze and simulate network evolution based on exponential-family random graph models (ERGM). 'tergm' is a part of the 'statnet' suite of packages for network analysis.

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

You may also want to install the corresponding latest binaries for packages on which `tergm` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`network`](https://github.com/statnet/network), [`networkDynamic`](https://github.com/statnet/networkDynamic), and [`ergm`](https://github.com/statnet/ergm).
