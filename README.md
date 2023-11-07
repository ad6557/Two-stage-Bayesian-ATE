# Two-stage-Bayesian-ATE
- Download the code
- open and run `run.R`

`risk.designfunction.R`: first stage
`risk.analysisfunction.R`: second stage
    `outcome~.-psvector` in regression: regress on treatment adjusted by covariates
    `outcome~treatment` in regression: regress on treatment only
`run.R`: run simulations.
