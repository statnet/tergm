# EpiModel

<details>

* Version: 2.0.5
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2021-05-15 18:50:03 UTC
* Number of recursive dependencies: 103

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly fixed

*   checking examples ... ERROR
    ```
    Running examples in ‘EpiModel-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: as.data.frame.netdx
    > ### Title: Extract Timed Edgelists netdx Objects
    > ### Aliases: as.data.frame.netdx
    > ### Keywords: extract
    > 
    > ### ** Examples
    > 
    ...
    Stopping at the initial estimate.
    Warning: Using x$coef to access the coefficient vector of an ergm is deprecated. Use coef(x) instead.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in is.durational(formation) : 
      could not find function "is.durational"
    Calls: netdx -> simulate -> simulate.network -> is.lasttoggle
    Execution halted
    ```

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
        4.       └─EpiModel:::netsim_loop(x, param, init, control, s)
        5.         ├─base::withCallingHandlers(...)
        6.         ├─base::do.call(...)
        7.         └─(function (x, param, init, control, s) ...
        8.           └─tergmLite::init_tergmLite(dat)
        9.             └─tergmLite:::stergm_prep(...)
       10.               ├─ergm::ergm_proposal(...)
       11.               └─ergm:::ergm_proposal.formula(...)
       12.                 ├─ergm::ergm_proposal(conlist, arguments, nw, ..., term.options = term.options)
       13.                 └─ergm:::ergm_proposal.ergm_conlist(...)
       14.                   └─ergm:::select_ergm_proposal(...)
      
      [ FAIL 6 | WARN 3 | SKIP 82 | PASS 379 ]
      Error: Test failures
      Execution halted
    ```

# ergm

<details>

* Version: 4.0.1
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2021-06-21 07:20:02 UTC
* Number of recursive dependencies: 81

Run `revdep_details(, "ergm")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 10.2Mb
      sub-directories of 1Mb or more:
        R      2.0Mb
        doc    2.2Mb
        help   1.5Mb
        libs   3.9Mb
    ```

# tergmLite

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/tergmLite
* Date/Publication: 2020-07-22 16:50:03 UTC
* Number of recursive dependencies: 70

Run `revdep_details(, "tergmLite")` for more info

</details>

## Newly broken

*   checking whether package ‘tergmLite’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: replacing previous import ‘ergm::snctrl’ by ‘tergm::snctrl’ when loading ‘tergmLite’
    See ‘/srv/scratch/z3528859/github/statnet/tergm/revdep/checks/tergmLite/new/tergmLite.Rcheck/00install.out’ for details.
    ```

*   checking R code for possible problems ... NOTE
    ```
    stergm_getMCMCsample: no visible global function definition for
      ‘stergm_MCMC_slave’
    Undefined global functions or variables:
      stergm_MCMC_slave
    ```

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Error: The combination of class (f), model constraints and hints ('.attributes' and 'sparse'), reference measure ("Bernoulli"), proposal weighting (default), and conjunctions and disjunctions is not implemented. Check your arguments for typos. 
      Backtrace:
          █
       1. └─EpiModel::initialize.net(est, param, init, control, s = 1) test-updateModelTermInputs.R:493:4
       2.   └─tergmLite::init_tergmLite(dat)
       3.     └─tergmLite:::stergm_prep(...)
       4.       ├─ergm::ergm_proposal(...)
       5.       └─ergm:::ergm_proposal.formula(...)
       6.         ├─ergm::ergm_proposal(conlist, arguments, nw, ..., term.options = term.options)
       7.         └─ergm:::ergm_proposal.ergm_conlist(...)
       8.           └─ergm:::select_ergm_proposal(...)
      
      [ FAIL 23 | WARN 2 | SKIP 0 | PASS 0 ]
      Error: Test failures
      Execution halted
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

