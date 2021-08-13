if (!require(dplyr)) install.packages(dplyr); library(dplyr)
if (!require(ggplot2)) install.packages(ggplot2); library(ggplot2)
if (!require(gridExtra)) install.packages(gridExtra); library(gridExtra)

# long data format
#dados <- read.csv('param_exploration.csv')

local_minima <- function(x, boundaries = TRUE) {
    x1 <- c(FALSE, diff(sign(diff(x))) == 2, FALSE)
    if (boundaries) {
        if (x[1] < x[2])
            x1[1] <- TRUE
        if (x[length(x)] < x[length(x)-1])
            x1[length(x)] <- TRUE
    }
    x1
}

find_minima <- function(dados) {
    value.col <- "neg.loglike"
    margin.cols <- setdiff(colnames(dados), value.col)

    # this looks very wasteful but might work??
    minima <- list()
    for (col in margin.cols) {
        minima[[col]] <- dados %>%
            group_by(across(-all_of(c(value.col, col)))) %>%
            filter(local_minima(neg.loglike, boundaries=F)) %>%
            ungroup %>%
            as.data.frame
    }
    all_minima <- do.call(rbind, minima)
    all_minima <- all_minima %>%
        group_by_all() %>%
        filter(n() == length(margin.cols)) %>%
        dplyr::slice(1) %>%
        ungroup %>%
        as.data.frame
}

# prof <- param_profiles(dados)
# do.call(grid.arrange, c(prof$plots, list(ncol = 2, nrow =3)))
# calcula perfis
param_profiles <- function(dados, plot = TRUE, x_lab = c("Initial fraction of the new variant", "Relative transmission rate\nfor the new variant",
                                                                  "Relative force of reinfection of P.1", "Prevalence of the previous infection",
                                                                  "Intrinsic growth rate")) {
    value.col <- "neg.loglike"
    margin.cols <- setdiff(colnames(dados), value.col)
    perfil <- list()
    plots.perfis <- list()
    for (col in margin.cols) {
        perfil[[col]] <- dados %>%
            group_by_at(col) %>%
            summarise(neg.loglike = min(neg.loglike)) %>%
            as.data.frame
        perfil[[col]]$neg.loglike <-  perfil[[col]]$neg.loglike - min(perfil[[col]]$neg.loglike)
        if (plot)
            plots.perfis[[col]] <- ggplot(perfil[[col]]) + geom_point(aes_string(x=col, y=value.col)) +
            theme_minimal() + ylab("Negative log likelihood\n") + xlab(paste0("\n",x_lab[which(col == margin.cols)]))
    }
    if (plot) {
        return(list(profiles = perfil, plots = plots.perfis))
    } else {
        return(perfil)
    }
}

## rewrite as 5-d array
## nt all that useful
#dims <- list()
#for (col in margin.cols)
#    dims[[col]] <- dados %>% select(col) %>% unique %>% arrange %>% unlist %>% as.vector
#
#dim.lengths <- sapply(dims, length)
#A = array(dim = dim.lengths, dim.names = dims)
## big expensive loop follows
