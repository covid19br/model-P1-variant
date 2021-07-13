README
================

This repository provides the script to reproduce the fitting analysis
and figures from Coutinho et al. 2021

# Notes:

-   This code runs best on R version 4.0.3 or above.
-   The code is also optimized for \[linux\].
-   When you first run the code from `Run_analysis.R`, R will install
    all required packages.
-   Running `all.params.exploration()` and `fit_minima()` takes a large
    number of hours to complete. For convenience, we also provide the
    tables resulting from these function as .csv files.

Run the file `Run_analysis.R` to reproduce the paper\`s main analysis
and sensitivity analysis.

## The following values correspond to the parameters presented in Table 1:

-   **REL.BETA.RATE2**: Relative transmission rate for the new variant
-   **PROB.REINFEC.V2**: Relative force of reinfection of P.1
-   **prevalence**: Prevalence of previous infection (2020-11-01) (%)
-   **INIT.V2.FRAC**: Initial fraction of the new variant (2020-11-01)
-   **r**: Intrinsic growth rate
-   **IHR.V2.PROP**: odds ratio of P.1 parameters relative to wild
    variant

# References:

Coutinho, R. M., Marquitti, F. M. D., Ferreira, L. S., Borges, M. E., da
Silva, R. L. P., Canton, O., … & Prado, P. I. (2021). Model-based
estimation of transmissibility and reinfection of SARS-CoV-2 P. 1
variant. medRxiv.
