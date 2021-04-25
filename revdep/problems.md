# EpiModel

<details>

* Version: 2.0.3
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2020-11-09 21:40:13 UTC
* Number of recursive dependencies: 102

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly broken

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
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Warning in simulate_formula.network(object, nsim = nsim, seed = seed, ...,  :
      For dynamic simulation in ‘tergm’ you must pass ‘dynamic=TRUE’.  Attempting ‘ergm’ simulation instead...
    Error in get(name, pos = environment(new)) : 
      argument "set.control.stergm" is missing, with no default
    Calls: netdx ... .handle.auto.constraints -> nonsimp_update.formula -> assign -> get
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
       17.               ├─tergm:::simulate_formula.network(...)
       18.               │ └─base::eval.parent(mc)
       19.               │   └─base::eval(expr, p)
       20.               │     └─base::eval(expr, p)
       21.               └─ergm:::.simulate_formula.network(...)
       22.                 └─ergm:::.handle.auto.constraints(...)
       23.                   └─statnet.common::nonsimp_update.formula(...)
       24.                     ├─base::assign(name, get(name, pos = environment(new)), pos = e)
       25.                     └─base::get(name, pos = environment(new))
      
      [ FAIL 33 | WARN 34 | SKIP 78 | PASS 231 ]
      Error: Test failures
      In addition: Warning message:
      replacing previous import 'vctrs::data_frame' by 'tibble::data_frame' when loading 'dplyr' 
      Execution halted
    ```

# tergmLite

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/tergmLite
* Date/Publication: 2020-07-22 16:50:03 UTC
* Number of recursive dependencies: 67

Run `revdep_details(, "tergmLite")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
        9.     └─ergm:::simulate.formula_lhs_network(...)
       10.       ├─ergm::simulate_formula(...)
       11.       ├─tergm:::simulate_formula.network(...)
       12.       │ └─base::eval.parent(mc)
       13.       │   └─base::eval(expr, p)
       14.       │     └─base::eval(expr, p)
       15.       └─ergm:::.simulate_formula.network(...)
       16.         └─ergm:::.handle.auto.constraints(...)
       17.           └─statnet.common::nonsimp_update.formula(...)
       18.             ├─base::assign(name, get(name, pos = environment(new)), pos = e)
       19.             └─base::get(name, pos = environment(new))
      
      [ FAIL 24 | WARN 25 | SKIP 0 | PASS 0 ]
      Error: Test failures
      Execution halted
    ```

*   checking whether package ‘tergmLite’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: replacing previous import ‘ergm::snctrl’ by ‘tergm::snctrl’ when loading ‘tergmLite’
    See ‘/homes/morrism/GitHub/StatnetOrganization/tergm/revdep/checks/tergmLite/new/tergmLite.Rcheck/00install.out’ for details.
    ```

*   checking R code for possible problems ... NOTE
    ```
    stergm_getMCMCsample: no visible global function definition for
      ‘stergm_MCMC_slave’
    Undefined global functions or variables:
      stergm_MCMC_slave
    ```

