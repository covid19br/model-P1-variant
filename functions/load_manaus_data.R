################################################################################
##  location-specific data
################################################################################
if (!require(dplyr)) install.packages(dplyr); library(dplyr)
if (!require(zoo)) install.packages(zoo); library(zoo)

source('./functions/aggregate_age_matrix.R')
source('functions/logit.R')

if(!exists("init_condits",mode="function")) source('./functions/beta_init_condits2.R')
if(!exists("fitP.exp",mode="function"))  source('./functions/fitP.exp.R')

# Demographic data
# Source: https://demografiaufrn.net/laboratorios/lepp/
age.distr <- read.csv('DATA/Age_distribution_2020_Manaus.csv')
POP.TOTAL.NUM = sum(age.distr) # Total Population Size
POP.DISTR <- c(sum(age.distr[1:4]), sum(age.distr[5:12]), sum(age.distr[13:length(age.distr)]))

# Contact matrix    ( JJ, JA, JI, AJ, AA, AI, IJ, IA, II)
load('DATA/prem2020_bra.Rdata') # Load contact matrixes from Prem et al. 2020

full.mat.cont <- new_all
# alternate matrix considering different interventions
# full.mat.cont <- new_home + 0.7*new_others + 0.7*new_work + 0.2*new_school

# aggregated Contact Matrix
aggregate_indices <- list(c(1:4), c(5:12), c(13:16))
CONTACT.M <- aggregate_age_matrix(mat.cont = full.mat.cont,
                                  aggregate_indices = aggregate_indices,
                                  age_structure = age.distr)

# the names with a preceding \ indicate the parameter name in the model equations
# the names in () represent the respective parameter name in CoMo BR model
# references are in the CoMo BR model parameter table

EXPOSURE.PERIOD.DAYS   = 5.8  # (gamma) Average time between being infected and developing symptoms \gamma
SICKNESS.PERIOD.DAYS   = 9    # (nui) Average time between being infectious and recovering for asymptomatic and mild \nu_i
SEVERE.PERIOD.DAYS     = 8.4  # (nus) Average time between being infectious and recovering/dying for severe cases \nu_s

CONT.REDUC.FRAC        = 0.1  # Reduction on the expose of symptomatic (due to symptoms/quarantining). \tau 0 means the same level of exposure of asymptomatic and 1 means no expose whatsoever. 
SEVERE.CONT.REDUC.FRAC = 0.9  # Reduction on the expose of severe cases (due to hospitalization). \eta 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
REL.INFEC.PRESYMP      = 1    # (rho) relative infectiousness of pre-symptomatic individuals. \rho

#### Classification of cases related to the severity of cases of unvaccinated individuals
ASYMPTOMATIC.FRAC = c(0.67, 0.44, 0.31) # Fraction of asymptomatic cases in total cases (pclin) \alpha

# Fraction of severe cases/hospitalizations in symptomatic cases (IHR) \sigma
ihr <- read.csv('DATA/ihr.csv')
weighted_ihr <- ihr[,2]/100 * age.distr
SEVERITY.FRAC = c(sum(weighted_ihr[1:4]) / sum(age.distr[1:4]),
                  sum(weighted_ihr[5:12]) / sum(age.distr[5:12]),
                  sum(weighted_ihr[13:19]) / sum(age.distr[13:19]))

# Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) \mu
ihfr_estados <- read.csv('DATA/ihfr.csv')
DEATH.FRAC <- ihfr_estados %>%
    filter(sg_uf == "AM") %>%
    arrange(age_clas) %>%
    .$ihfr.covid

### NEW VARIANT parameters
REL.BETA.RATE2 = 1.3
REL.SEVERITY.FRAC.V2 = 1.0
PROB.REINFEC.V2 = 0.
EXPOSURE.PERIOD.DAYS.V2   = 5.8  # (gamma) Average time between being infected and developing symptoms \gamma
SICKNESS.PERIOD.DAYS.V2   = 9    # (nui) Average time between being infectious and recovering for asymptomatic and mild \nu_i
SEVERE.PERIOD.DAYS.V2     = 8.4  # (nus) Average time between being infectious and recovering/dying for severe cases \nu_s

CONT.REDUC.FRAC.V2        = 0.1  # Reduction on the expose of symptomatic (due to symptoms/quarantining). \tau 0 means the same level of exposure of asymptomatic and 1 means no expose whatsoever. 
SEVERE.CONT.REDUC.FRAC.V2 = 0.9  # Reduction on the expose of severe cases (due to hospitalization). \eta 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
REL.INFEC.PRESYMP.V2      = 1    # (rho) relative infectiousness of pre-symptomatic individuals. \rho
ASYMPTOMATIC.FRAC.V2 = ASYMPTOMATIC.FRAC
DEATH.FRAC.V2 = DEATH.FRAC
SEVERITY.FRAC.V2 = REL.SEVERITY.FRAC.V2*SEVERITY.FRAC

# prevalence
cases_dist_age <- read.csv('DATA/cases_dist_age_states.csv')
cases_dist_age <- cases_dist_age %>%
    filter(sg_uf == "AM") %>%
    arrange(age_clas) %>%
    .$N.covid

Recovered <- cases_dist_age / (ihr[,2] / 100)
POP.R1 <- c(sum(Recovered[1:4]),
            sum(Recovered[5:12]),
            sum(Recovered[13:19]))

# age distribution of recent hospitalized cases from Manaus
recent_cases_dist_age <- read.csv('DATA/recent_cases_dist_age_states.csv')
recent_cases_dist_age <- recent_cases_dist_age %>%
    filter(sg_uf == "AM") %>%
    arrange(age_clas) %>%
    .$N.covid

# Nowcasting from covid hospitalization in Manaus
now.covid <- (read.csv("DATA/now_covid_zoo_manaus.csv", header = T))
now.covid[,1] <- as.Date(now.covid[,1])

now.covid.zoo <- zoo(now.covid[,2], now.covid[,1])
now.covid.zoo <- now.covid.zoo[index(now.covid.zoo) < "2021-02-01",]

state.zoo.list <- list(now.covid.zoo = now.covid.zoo,
                       data.base = "2021_03_15") # data.base refers to SIVEP's database date used for the nowcasting projections

covid.zoo <- zoo(diff(c(state.zoo.list$now.covid.zoo)), index(state.zoo.list$now.covid.zoo)[-1])
index(covid.zoo) <- as.Date(index(covid.zoo))

new.hosp <- as.numeric(covid.zoo[as.Date("2020-11-01")]) *
    recent_cases_dist_age / sum(recent_cases_dist_age)

# probability of reporting hospitalized case
# (could be, for instance, ratio of Covid/SRAG reported)
REPORT.PROB = 1.

# exponential fit
#r <- fitP.exp(tail(covid.zoo, n = 3), only.coef = FALSE)$coefficients[2]

#new.hosp <- coredata(covid.zoo[max(as.Date(time(covid.zoo)))]) * recent_cases_dist_age / sum(recent_cases_dist_age)

# init.conds <- init_condits(r, new.hosp, PREVALENCE = POP.R1/POP.DISTR, POP.DISTR,
#                            CONTACT.M, EXPOSURE.PERIOD.DAYS, SICKNESS.PERIOD.DAYS,
#                            SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC, SEVERE.CONT.REDUC.FRAC,
#                            REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC, SEVERITY.FRAC, DEATH.FRAC,
#                            V2.FRAC = 0)

# POP0 <- init.conds$POP0
# Average contact rate and success between unvaccinated susceptibles and
# infectious individuals \beta
# calculated from hospitalized cases in the last 3 weeks
# BETA.RATE1 <- init.conds$BETA.RATE
# names(POP0) <- paste0(rep(c("POP.S", "POP.E1", "POP.A1", "POP.I1", "POP.H1",
#                             "POP.C1", "POP.R1", "POP.D1",
#                             "POP.E2", "POP.A2", "POP.I2", "POP.H2",
#                             "POP.C2", "POP.R2", "POP.D2"), each=3), rep(1:3, 15))


## Removed variables used in this code but not needed in the future
## in order to keep the workspace clear
# rm(ihfr_estados, weighted_ihr, age.distr, new_all, new_home,
#    new_others, new_school, new_work, u_all, u_home, u_work, u_others, u_school,
#    ihr, full.mat.cont, aggregate_indices, r_all, r_home, r_others,
#    r_school, r_work, cases_dist_age, recent_cases_dist_age, Recovered, POP.R1,
#    new.hosp)

# http://www.genomahcov.fiocruz.br/presenca-das-linhagens-por-estado/
d <- read.csv('DATA/freq_P1_manaus.csv')
d$t_start <- as.Date(d$t_start)
d$t_end <- as.Date(d$t_end)
d$freq = d$P.1 / (d$total)
# TODO: avoid infinite values using totals
d$freq.logit <- logit(d$freq)

TIME.VECTOR <- seq(0, 91) # 13 weeks
d0 <- as.Date("2020-11-01") # start date of simulation