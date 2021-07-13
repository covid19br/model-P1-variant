######################################
##### TWO-STRAIN SEIRHD MODEL ########
######################################
# This file contain the definition of the differential equations and the
# numerical integrator.

if (!require(deSolve)) install.packages(deSolve); library(deSolve)
if (!require(zeallot)) install.packages(zeallot); library(zeallot)
if (!require(dplyr)) install.packages(dplyr); library(dplyr)

# Reads the epidemiological parameters for COVID-19 and Brazilian demographics
source('functions/load_manaus_data.R')

# Parameters of the numerical integration
SIM.DURATION.DAYS = 120 # Days of simulation
TIME.VECTOR <- seq(0, SIM.DURATION.DAYS) 

func.factory <- function(
    BETA_RATE1 = BETA.RATE1,
    REL_BETA_RATE2 = REL.BETA.RATE2, CONT_REDUC_FRAC = CONT.REDUC.FRAC, CONTACT_M =
    CONTACT.M, DEATH_FRAC = DEATH.FRAC, EXPOSURE_PERIOD_DAYS =
    EXPOSURE.PERIOD.DAYS, FIRST_DOSE_REL_EFFIC = FIRST.DOSE.REL.EFFIC,
    POP_DISTR = POP.DISTR, POP_ESTADO_REL_FRAC = POP.ESTADO.REL.FRAC,
    POP_TOTAL_NUM = POP.TOTAL.NUM, REL_INFEC_PRESYMP = REL.INFEC.PRESYMP,
    ASYMPTOMATIC_FRAC = ASYMPTOMATIC.FRAC,SEVERE_CONT_REDUC_FRAC = SEVERE.CONT.REDUC.FRAC,
    SEVERE_PERIOD_DAYS = SEVERE.PERIOD.DAYS, SEVERITY_FRAC = SEVERITY.FRAC,
    SICKNESS_PERIOD_DAYS = SICKNESS.PERIOD.DAYS, PROB_REINFEC_V2 = PROB.REINFEC.V2,
    CONT_REDUC_FRAC_V2 = CONT.REDUC.FRAC.V2, DEATH_FRAC_V2 = DEATH.FRAC.V2,
    EXPOSURE_PERIOD_DAYS_V2 = EXPOSURE.PERIOD.DAYS.V2,
    REL_INFEC_PRESYMP_V2 = REL.INFEC.PRESYMP.V2, ASYMPTOMATIC_FRAC_V2 = ASYMPTOMATIC.FRAC.V2,
    SEVERE_CONT_REDUC_FRAC_V2 = SEVERE.CONT.REDUC.FRAC.V2,
    SEVERE_PERIOD_DAYS_V2 = SEVERE.PERIOD.DAYS.V2, SEVERITY_FRAC_V2 = SEVERITY.FRAC.V2,
    SICKNESS_PERIOD_DAYS_V2 = SICKNESS.PERIOD.DAYS.V2, REPORT_PROB = REPORT.PROB) {
    # code that generates this beauty:
    # source('./parms_vaccine.R')
    # objs = setdiff(ls(), lsf.str())
    # paste0(lapply(objs, function(x) paste0(gsub('\\.', '_', x), " = ", x)), collapse=", ")
    #
    # WARNING: several of those are not actual parameters! They are only used
    # in parms_* to calculate other parameters

    diffEqs = function(t, POP, parms) {
        c(POP.S, POP.E1, POP.A1, POP.I1, POP.H1, POP.C1, POP.R1, POP.D1,
                 POP.E2, POP.A2, POP.I2, POP.H2, POP.C2, POP.R2, POP.D2) %<-%
        split(POP, ceiling(seq_along(POP)/3))
    
        debug = FALSE
        if(debug && (any(POP < -1e-10) || ! all(is.finite(POP)))){
            cat(paste("Negative population!",
                      paste("time t = ", t),
                      paste(c(POP.S, POP.E1), collapse = ', '),
                      paste(c(POP.A1, POP.I1, POP.H1), collapse = ', '),
                      paste(c(POP.R1, POP.D1), collapse = ', '),
                      paste(c(POP.E2, POP.A2, POP.I2, POP.H2), collapse = ', '),
                      paste(c(POP.R2, POP.D2), collapse = ', '),
                      sep="\n"))
            stop()
        }
    
        # Some repeated factors:
        # Total infectious
        lambda1 = (BETA_RATE1 / POP_TOTAL_NUM) * CONTACT_M %*%
            (POP.A1 + REL_INFEC_PRESYMP*POP.E1  + (1-CONT_REDUC_FRAC)*POP.I1 + (1-SEVERE_CONT_REDUC_FRAC)*POP.H1)
        lambda2 = (REL_BETA_RATE2 * BETA_RATE1 / POP_TOTAL_NUM) * CONTACT_M %*%
            (POP.A2 + REL_INFEC_PRESYMP_V2*POP.E2  + (1-CONT_REDUC_FRAC_V2)*POP.I2 + (1-SEVERE_CONT_REDUC_FRAC_V2)*POP.H2)

        # Susceptible
        dS    = - POP.S * (lambda1 + lambda2) # Getting infected

        ## "wild" strain
        # Pre-symptomatic
        dE1   =   POP.S * lambda1 - # Getting infected
                  POP.E1 / EXPOSURE_PERIOD_DAYS # Becoming asymptomatic or mild case
        # Assymptomatic
        dA1   =   ASYMPTOMATIC_FRAC * (1 - SEVERITY_FRAC) * POP.E1 / EXPOSURE_PERIOD_DAYS - # Becoming asymptomatic
                  POP.A1 / SICKNESS_PERIOD_DAYS # Recovering
        # Mild Infectious
        dI1   =   (1 - ASYMPTOMATIC_FRAC) * (1 - SEVERITY_FRAC) * POP.E1 / EXPOSURE_PERIOD_DAYS - # Becoming mild case
                  POP.I1 / SICKNESS_PERIOD_DAYS # Recovering
        # Severe Case/Hospitalization
        dH1   =   SEVERITY_FRAC * POP.E1 / EXPOSURE_PERIOD_DAYS - # Becoming severe case
                  POP.H1 / SEVERE_PERIOD_DAYS
        # Cases counted
        dC1   =   c(0.0, 0.0, 0.0)
        # Recovered
        dR1   =   POP.A1 / SICKNESS_PERIOD_DAYS + # Asymptomatic recovering
                  POP.I1 / SICKNESS_PERIOD_DAYS + # Mild case recovering
                  (1-DEATH_FRAC)*POP.H1 / SEVERE_PERIOD_DAYS - # Severe case recovering
                  PROB_REINFEC_V2 * POP.R1 * lambda2 # Reinfection by the new strain
        # Dead
        dD1   =   DEATH_FRAC *POP.H1 / SEVERE_PERIOD_DAYS # Severe case dying

        ## new strain
        # Pre-symptomatic
        dE2   =   (POP.S + PROB_REINFEC_V2 * POP.R1) * lambda2 - # Getting infected
                  POP.E2 / EXPOSURE_PERIOD_DAYS_V2 # Becoming asymptomatic or mild case
        # Assymptomatic
        dA2   =   ASYMPTOMATIC_FRAC_V2 * (1 - SEVERITY_FRAC_V2) * POP.E2 / EXPOSURE_PERIOD_DAYS_V2 - # Becoming asymptomatic
                  POP.A2 / SICKNESS_PERIOD_DAYS_V2 # Recovering
        # Mild Infectious
        dI2   =   (1 - ASYMPTOMATIC_FRAC_V2) * (1 - SEVERITY_FRAC_V2) * POP.E2 / EXPOSURE_PERIOD_DAYS_V2 - # Becoming mild case
                  POP.I2 / SICKNESS_PERIOD_DAYS_V2 # Recovering
        # Severe Case/Hospitalization
        dH2   =   SEVERITY_FRAC_V2 * POP.E2 / EXPOSURE_PERIOD_DAYS_V2 - # Becoming severe case
                  POP.H2 / SEVERE_PERIOD_DAYS_V2
        dC2   =   c(0.0, 0.0, 0.0)
        # Recovered
        dR2   =   POP.A2 / SICKNESS_PERIOD_DAYS_V2 + # Asymptomatic recovering
                  POP.I2 / SICKNESS_PERIOD_DAYS_V2 + # Mild case recovering
                  (1-DEATH_FRAC_V2) * POP.H2 / SEVERE_PERIOD_DAYS_V2 # Severe case recovering
        # Dead
        dD2   =   DEATH_FRAC_V2 * POP.H2 / SEVERE_PERIOD_DAYS_V2 # Severe case dying


        if(debug && ! all(is.finite(c(dS, dE1, dA1, dI1, dH1, dR1, dD1,
                                          dE2, dA2, dI2, dH2, dR2, dD2)))){
            cat(paste(" *** Infinite Derivatives/NaN! ***",
                      paste("tempo t = ", t),
                      paste(c(dS, dE1), collapse = ', '),
                      paste(c(dA1, dI1, dH1), collapse = ', '),
                      paste(c(dR1, dD1), collapse = ', '),
                      sep="\n"))
            cat(paste(" *** Populations: ***",
                      paste(c(POP.S, POP.E1), collapse = ', '),
                      paste(c(POP.A1, POP.I1, POP.H1), collapse = ', '),
                      paste(c(POP.R1, POP.D1), collapse = ', '),
                      sep="\n"))
    
            stop()
        }
    
        return(list(c(dS, dE1, dA1, dI1, dH1, dC1, dR1, dD1,
                          dE2, dA2, dI2, dH2, dC2, dR2, dD2)))
    }

    return(diffEqs)
}


#diffEqs <- func.factory(BETA_RATE1 = BETA.RATE1, REL_BETA_RATE2 = REL.BETA.RATE2,
#                       PROB_REINFEC_V2 = 0)
#SOLUTION <- ode(y = POP0, times = TIME.VECTOR, func = diffEqs, parms = c())
#source('./functions/C_calculator.R')
#SOLUTION <- C_calculator(SOLUTION, EXPOSURE.PERIOD.DAYS.V2, SEVERITY.FRAC.V2)
