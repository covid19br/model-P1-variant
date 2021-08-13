#' The following variable names corresponds to the following parameters shown in the paper:
#' REL.BETA.RATE2: Relative transmission rate for the new variant
#' PROB.REINFEC.V2: Relative force of reinfection of P.1
#' prevalence: Prevalence of previous infection (2020-11-01) (%)
#' INIT.V2.FRAC: Initial fraction of the new variant (2020-11-01) 
#' r: Intrinsic growth rate
#' IHR.V2.PROP: odds ratio of P.1 parameters relative to wild variant

### Set up
source('functions/ajuste_funs.R')
source('functions/fit_minima.R')
source('functions/mle.helpers.R')
source('functions/CI_estimator.R')
source('functions/get.reinfections.R')

LOCATION <- "results" # Set Output folder

###################################################################
#####################  Main Analysis    ###########################
###################################################################

####################################
# Run extensive parameter fitting
####################################

## Set upper and lower limits for parameter fitting

lower <- c(INIT.V2.FRAC = 1e-6,
           REL.BETA.RATE2 = 1.4,
           PROB.REINFEC.V2 = 1e-3,
           prevalence = 40,
           r = 0)

upper <- c(INIT.V2.FRAC = 0.01,
           REL.BETA.RATE2 = 3.4,
           PROB.REINFEC.V2 = 0.9,
           prevalence= 85,
           r = 0.08)

# This process takes about ~ 12 h in the cluster with 24 threads and 32 gigabytes of RAM
# main_result <- all.params.exploration(prep.guess(lower), prep.guess(upper),
#                                       conv.guess=prep.guess.inv, LHS = F, n=1e6)

## Save results

# write.csv(main_result, file = paste0(LOCATION, "/param_exploration.csv"), row.names = FALSE)

main_result <-read.csv(paste0(LOCATION, "/param_exploration.csv.xz"))

## calculate and plot profiles.

profile <- param_profiles(main_result)
profile_main_plot <- do.call(grid.arrange, c(profile$plots, list(ncol = 2, nrow = 3)))
profile_main_plot
#ggsave(profile_main_plot, file = "results/profile_main_plot.pdf", width = 8, height = 8)

###################
# Refine best fits
###################

# Obs: This method may take a couple hours

# best.fits <- fit_minima(paste0(LOCATION, "/param_exploration.csv.xz"),
#                         helper.mle.fun.2,
#                         paste0(LOCATION, "/best_fits_main.csv"),
#                         n = 100,
#                         optimizer = "optimx",
#                         method = "bobyqa")

best.fits <- read.csv(paste0(LOCATION, "/best.fits.csv"))

## Plot profile from best fits

profile_bestfits <- param_profiles(best.fits %>% select(-loglike))
profile_bestfits$plots[[1]] <- profile_bestfits$plots[[1]] + xlim(0.00033,0.00037) + ylim(0,0.004)#+ xlim(0.00033,0.00037) + ylim(0,0.0015)
profile_bestfits$plots[[2]] <- profile_bestfits$plots[[2]] + xlim(2.604,2.6059) + ylim(0,0.0001)
profile_bestfits$plots[[3]] <- profile_bestfits$plots[[3]] + xlim(0.0278,0.0280) + ylim(0,0.00015)
profile_bestfits$plots[[4]] <- profile_bestfits$plots[[4]] + xlim(72,74.7) + ylim(0,0.1)
profile_bestfits$plots[[5]] <- profile_bestfits$plots[[5]] + xlim(0.0307,0.0309) + ylim(0,0.0009)

profile_bestfits_plot <- do.call(grid.arrange, c(profile_bestfits$plots, list(ncol = 2, nrow = 3)))
profile_bestfits_plot
#ggsave(profile_bestfits_plot, file = "results/profile_bestfits_plot.pdf", width = 8, height = 8)


# Create best fit object
best <- best.fits[1,]
best$loglike <- NULL
best.fit <- mle2(helper.mle.fun.2,
                 start = as.list(prep.guess(best)),
                 optimizer = "optimx", method = "bobyqa") # Depending on your computer it can take a while

## Table of parameter estimates and 95% CI (Table 1: main fitting)

# Load output (To avoid take all the time on computing the mle2)
# load(paste0(LOCATION, "/bestfit.RData")) # You can alternatively load the best.fit object

best.parms <- prep.coef(best.fit)               # Get estimates
best.conf <- prep.conf(best.fit, method="quad") # Get parameters confidence interval
print(cbind(estimate = best.parms, best.conf))  # Print result table

# Save output
#write.csv(cbind(estimate = best.parms, best.conf),
#          file = paste0(LOCATION, "/best_coefs.csv"))

# Calculate solution confidence interval
## It take a while for greater n, n=1000 will take about several minutes
sol.CI <- solution_CI(best.fit, n=1000, conv.logit=c(1,3,4))

## Plot results (Figure 2) 
plot.3.CI(sol.CI)

# Estimate number of reinfections

sol.r <- get.reinfections(best)
reinfection_proportion <- sol.r["reinf"] / sum(sol.r["C2"])
print(reinfection_proportion)

###################################################################
##############  Sensitivity Analyses    ###########################
###################################################################


###########################
### SA 1 - Free IHR 
###########################

# Sensitivity to different pathogenicity of the variant P.1 is explored
# by repeating the fit assuming IHR as a free parameter (SA1).

####################################
# Run extensive parameter fitting
####################################

## Set upper and lower limits for parameter fitting
## IHR.V2.PROP is the odds ratio of P.1 parameters relative to wild type

lower <- c(INIT.V2.FRAC = 1e-6,
           REL.BETA.RATE2 = 1.4,
           PROB.REINFEC.V2 = 1e-3,
           prevalence = 40,
           r = 0,
           IHR.V2.PROP = 0.8)

upper <- c(INIT.V2.FRAC = 0.01,
           REL.BETA.RATE2 = 3.4,
           PROB.REINFEC.V2 = 0.9,
           prevalence= 85,
           r = 0.08,
           IHR.V2.PROP = 2.0)

# This process takes about ~ 12 h in the cluster with 24 threads and 32 gigabytes of RAM
# SA1_result <- all.params.exploration(prep.guess(lower), prep.guess(upper),
#                              conv.guess=prep.guess.inv, LHS = F, n=5e5)

## Save results
# write.csv(SA1_result, file = paste0(LOCATION, "/param_exploration_ihr_free_test.csv"),
#          row.names = FALSE)

SA1_result <- read.csv(paste0(LOCATION, "/param_exploration_ihr_free_test.csv"))

## calculate and plot profiles

profile_sa1 <- param_profiles(SA1_result)
do.call(grid.arrange, c(profile_sa1$plots, list(ncol = 2, nrow = 3)))

###################
# Refine best fits
###################

# Obs: This method may take a couple hours

# best.fits.sa1 <- fit_minima(paste0(LOCATION, "/param_exploration_ihr_free_test.csv"),
#                         helper.mle.fun.6,
#                         paste0(LOCATION, "/best.fits_ihr_free_test.csv"),
#                         n = 100,
#                         optimizer = "optimx",
#                         method = "bobyqa")

best.fits.sa1 <- read.csv(paste0(LOCATION, "/best.fits_ihr_free.csv"))

# Create best fit object
best.sa1 <- best.fits.sa1[1,]
best.sa1$loglike <- NULL
best.fit.sa1 <- mle2(helper.mle.fun.6,
                 start = as.list(prep.guess(best.sa1)),
                 optimizer = "optimx", method = "bobyqa") # Depending on your computer it may take a while

## Table of parameter estimates and 95% CI (Table 1: SA1)

best.parms <- prep.coef(best.fit.sa1)               # Get estimates
best.conf <- prep.conf(best.fit.sa1, method="quad") # Get parameters confidence interval
print(bind(estimate = best.parms, best.conf))       # Print result table

# Save output
#write.csv(cbind(estimate = best.parms, best.conf),
#          file = paste0(LOCATION, "/best_coefs_ihr_free.csv"))

###################################
### SA 2 - Health system collapse
###################################

# Sensitivity to period prior to the health system collapse is explored
# by filtering data after 10-jan-2021 (SA2)

# This filters the data after January 10
covid.zoo <- window(covid.zoo, end="2021-01-10")

###################
# Refine best fits
###################

# Obs: This method may take a couple hours
# We used best fit of the main analysis and used as initial guess for bobyqa method

best.fits.sa2 <- fit_minima(paste0(LOCATION, "/param_exploration.csv"),
                        helper.mle.fun.2,
                        paste0(LOCATION, "/best.fits_up_to_mid-jan.csv"),
                        n = 100,
                        optimizer = "optimx",
                        method = "bobyqa")

# Create best fit object
best.sa2 <- best.fits.sa2[1,]
best.sa2$loglike <- NULL
best.fit.sa2 <- mle2(helper.mle.fun.2,
                 start = as.list(prep.guess(best.sa2)),
                 optimizer = "optimx", method = "bobyqa") # Depending on your computer it may take a while

# table of parameter estimates and 95% CI (Table 1: SA2)

best.parms <- prep.coef(best.fit.sa2)               # Get estimates
best.conf <- prep.conf(best.fit.sa2, method="quad") # Get parameters confidence interval
print(cbind(estimate = best.parms, best.conf))      # Print result table
