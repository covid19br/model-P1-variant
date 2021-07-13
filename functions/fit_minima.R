source('functions/find_minima.R')

helper.mle.fun.2 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                               prevalence, r) {
        l <- as.list(match.call())
        l[[1]] <- NULL
        params <- unlist(l)
        params[c(1,3,4)] <- inv.logit(params[c(1,3,4)])
        params[4] <- 100*params[4] ## the model is fitted with prevalence in 0-100 scale
        ## params[c(2,5)] <- exp(params[c(2,5)])  ## log-link for strictly positive values, not working with solver
        solution <- model_solution(params)
        LL.cases <- with(solution$cases, -sum(dpois(covid.zoo, lambda = C1+C2, log = TRUE)))
        LL.freq <- with(solution$freq, -sum(dbinom(P.1, size = total, prob = pred.freq, log = TRUE)))
        return(LL.cases + LL.freq)
}

helper.mle.fun.5 <- function(INIT.V2.FRAC, REL.BETA.RATE2, PROB.REINFEC.V2,
                               prevalence, r, IHR.V2.PROP) {
        l <- as.list(match.call())
        l[[1]] <- NULL
        params <- unlist(l)
        params[c(1,3,4)] <- inv.logit(params[c(1,3,4)])
        params[4] <- 100*params[4] ## the model is fitted with prevalence in 0-100 scale
        ## params[c(2,5)] <- exp(params[c(2,5)])  ## log-link for strictly positive values, not working with solver
        solution <- model_solution(params)
        LL.cases <- with(solution$cases, -sum(dpois(covid.zoo, lambda = C1+C2, log = TRUE)))
        LL.freq <- with(solution$freq, -sum(dbinom(P.1, size = total, prob = pred.freq, log = TRUE)))
        return(LL.cases + LL.freq)
}


fit_minima <- function(datafile, helper, outfile, n = 100, ignore.cols = c(), ...) {
    dados <- read.csv(datafile)
    dados <- dados %>% dplyr::select(-all_of(ignore.cols))

    minima <- find_minima(dados)
    minima.10 <- minima %>% slice_min(neg.loglike, n = n) %>% dplyr::select(-c("neg.loglike"))

    minima.10 <- prep.guess(minima.10)
    print("mle will start")
    best.fits <- fit.mle.guesses(minima.10, helper, ...)
    save(best.fits, file = "bestfits_.Rdata")
    #load(file = "bestfits_.Rdata")
    print("mle ok")
    best.fits <- as.data.frame(lapply(best.fits,unlist)) %>%
        filter(! is.na(loglike)) %>%
        arrange(desc(loglike)) %>%
        prep.guess.inv()
    if (! missing(outfile))
        write.csv(best.fits, outfile, row.names=F)
    return(best.fits)
}

# example
# source('functions/ajuste_funs.R')
## CUIDADO: carregar dados DEPOIS
# source('manaus_data.R')
# best.fits <- fit_minima('results/param_exploration_ihr2.csv.xz',
# 'results/best.fits.csv', n = 100)

