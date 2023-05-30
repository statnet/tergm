# BasketballAnalyzeR

<details>

* Version: 0.5.0
* GitHub: https://github.com/sndmrc/BasketballAnalyzeR
* Source code: https://github.com/cran/BasketballAnalyzeR
* Date/Publication: 2020-06-26 09:00:11 UTC
* Number of recursive dependencies: 74

Run `revdep_details(, "BasketballAnalyzeR")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘circlize’ ‘hexbin’ ‘scales’ ‘sna’
      All declared Imports should be used.
    ```

# Blaunet

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2022-09-27 08:10:08 UTC
* Number of recursive dependencies: 75

Run `revdep_details(, "Blaunet")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘digest’ ‘ergm’ ‘foreign’ ‘gWidgets2’ ‘gWidgets2tcltk’ ‘haven’
      ‘plot3D’ ‘plot3Drgl’ ‘rgl’ ‘sna’ ‘statnet.common’
      All declared Imports should be used.
    ```

# EpiModel

<details>

* Version: 2.3.2
* GitHub: https://github.com/EpiModel/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2023-02-16 23:30:02 UTC
* Number of recursive dependencies: 120

Run `revdep_details(, "EpiModel")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘EpiModel-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: as.data.frame.netdx
    > ### Title: Extract Timed Edgelists for netdx Objects
    > ### Aliases: as.data.frame.netdx
    > ### Keywords: extract
    > 
    > ### ** Examples
    > 
    ...
    Maximizing the pseudolikelihood.
    Finished MPLE.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in simulate.formula_lhs(object = TARGET_STATS ~ edges, nsim = if (dynamic ==  : 
      No applicable method for LHS of type ‘NULL’.
    Calls: netdx ... eval -> eval -> <Anonymous> -> simulate.formula_lhs
    Execution halted
    ```

*   checking tests ...
    ```
      Running ‘test-all.R’
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
      ── Error ('test-utils.R:38:3'): color_tea ──────────────────────────────────────
      Error in `UseMethod("get.vertex.attribute")`: no applicable method for 'get.vertex.attribute' applied to an object of class "NULL"
      Backtrace:
          ▆
       1. └─EpiModel::netsim(est, param, init, control) at test-utils.R:38:2
       2.   └─EpiModel::crosscheck.net(x, param, init, control)
       3.     ├─get_vertex_attribute(nw, "status") %in% c("s", "i", "r")
       4.     └─EpiModel::get_vertex_attribute(nw, "status")
       5.       └─network::get.vertex.attribute(...)
      
      [ FAIL 37 | WARN 1 | SKIP 119 | PASS 376 ]
      Error: Test failures
      In addition: Warning message:
      replacing previous import 'ergm::snctrl' by 'ergm.multi::snctrl' when loading 'tergm' 
      Execution halted
    ```

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘attributes-and-summary-statistics.Rmd’ using rmarkdown
    
    Quitting from lines 303-338 [tracker-list] (attributes-and-summary-statistics.Rmd)
    Error: processing vignette 'attributes-and-summary-statistics.Rmd' failed with diagnostics:
    no applicable method for 'get.vertex.attribute' applied to an object of class "NULL"
    --- failed re-building ‘attributes-and-summary-statistics.Rmd’
    
    --- re-building ‘Intro.Rmd’ using rmarkdown
    --- finished re-building ‘Intro.Rmd’
    ...
    --- failed re-building ‘model-parameters.Rmd’
    
    --- re-building ‘network-objects.Rmd’ using rmarkdown
    --- finished re-building ‘network-objects.Rmd’
    
    SUMMARY: processing the following files failed:
      ‘attributes-and-summary-statistics.Rmd’ ‘model-parameters.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# ergm

<details>

* Version: 4.5.0
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2023-05-28 09:50:07 UTC
* Number of recursive dependencies: 88

Run `revdep_details(, "ergm")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.0Mb
      sub-directories of 1Mb or more:
        doc    1.7Mb
        libs   2.6Mb
        R      1.0Mb
    ```

# latentnet

<details>

* Version: 2.10.6
* GitHub: https://github.com/statnet/latentnet
* Source code: https://github.com/cran/latentnet
* Date/Publication: 2022-05-11 12:30:05 UTC
* Number of recursive dependencies: 115

Run `revdep_details(, "latentnet")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘ergm.userterms’
    ```

*   checking Rd files ... NOTE
    ```
    checkRd: (-1) simulate.ergmm.Rd:35: Escaped LaTeX specials: \$
    checkRd: (-1) simulate.ergmm.Rd:36: Escaped LaTeX specials: \$
    ```

*   checking Rd cross-references ... NOTE
    ```
    Unknown package ‘ergm.userterms’ in Rd xrefs
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# networkDynamic

<details>

* Version: 0.11.3
* GitHub: NA
* Source code: https://github.com/cran/networkDynamic
* Date/Publication: 2023-02-16 08:20:02 UTC
* Number of recursive dependencies: 39

Run `revdep_details(, "networkDynamic")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running ‘activate_tests.R’
      Running ‘age.at_tests.R’
      Running ‘classTests.R’
      Running ‘converter_tests.R’
      Running ‘get_tests.R’
      Running ‘import_tests.R’
      Running ‘network_tests.R’
      Running ‘pid_tests.R’
      Running ‘query_tests.R’
      Running ‘reconcile.activityTests.R’
    ...
      > expect_equal(network.size(net),5)
      > expect_true(is.networkDynamic(net))
      > expect_equal(unlist(get.vertex.activity(net,as.spellList=TRUE)[4:5,1:2]),c(1,1,2,2),check.names=FALSE)
      Error: unlist(get.vertex.activity(net, as.spellList = TRUE)[4:5, 1:2]) (`actual`) not equal to c(1, 1, 2, 2) (`expected`).
      
      `names(actual)` is a character vector ('onset1', 'onset2', 'terminus1', 'terminus2')
      `names(expected)` is absent
      In addition: Warning message:
      Unused arguments (check.names = FALSE) 
      Execution halted
    ```

# pkggraph

<details>

* Version: 0.2.3
* GitHub: https://github.com/talegari/pkggraph
* Source code: https://github.com/cran/pkggraph
* Date/Publication: 2018-11-15 09:50:03 UTC
* Number of recursive dependencies: 71

Run `revdep_details(, "pkggraph")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

# statnet

<details>

* Version: 2019.6
* GitHub: https://github.com/statnet/statnet
* Source code: https://github.com/cran/statnet
* Date/Publication: 2019-06-14 08:00:06 UTC
* Number of recursive dependencies: 98

Run `revdep_details(, "statnet")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘networksis’
    ```

*   checking startup messages can be suppressed ... NOTE
    ```
    unable to reach CRAN
    
    It looks like this package (or a package it requires) has a startup
    message which cannot be suppressed: see ?packageStartupMessage.
    ```

# statnet.common

<details>

* Version: 4.9.0
* GitHub: https://github.com/statnet/statnet.common
* Source code: https://github.com/cran/statnet.common
* Date/Publication: 2023-05-24 16:10:02 UTC
* Number of recursive dependencies: 20

Run `revdep_details(, "statnet.common")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    Found the following possibly unsafe calls:
    File ‘statnet.common/R/control.utilities.R’:
      unlockBinding("snctrl", environment(snctrl))
      unlockBinding("snctrl", environment(update_my_snctrl))
    ```

