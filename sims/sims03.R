#!/usr/bin/env Rscript

suppressMessages({
    library(asthmasims)
    library(parallel)
    library(extraDistr)
    library(splines2)
    library(data.table)
    library(optparse)
    library(qs)
    library(pbmcapply)
})


# ----- Command line arguments -----

option_list <- list(
    make_option(c("-c", "--cores"),
        type = "integer", default = 10,
        help = "number of cores to use [default %default]", metavar = "number"
    ),
    make_option(c("-n", "--nsim"),
        type = "integer", default = 10,
        help = "number of simulations to run under each configuration [default %default]", metavar = "number"
    ),
    make_option(c("-f", "--filename"),
        type = "character", default = "sims_scenario03_",
        help = "the output file name for the simulations [default %default]", metavar = "character"
    )
)
opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims <- opt$nsim
file_name <- opt$filename


# ----- PRNGs -----

RNGkind("L'Ecuyer-CMRG")
set.seed(357767)
mc.reset.stream()


# ----- Functions -----

generate_spline_basis <- function(x, iknots = 9, quants = FALSE, ...) {
    knot_seq <- seq(min(x), max(x), length.out = iknots + 2)
    knots <- floor(knot_seq[2:(iknots + 1)])
    return(iSpline(x, knots = knots, degree = 4, ...))
}


# ----- Configs -----

pt1 <- ddgamma(0:100, 6, 1 / 15.9)
pt1 <- pt1 / sum(pt1)
alp <- asthmasims:::ord_logit(pt1)
cfg <- CJ(
    sims = num_sims,
    n_seq = list(c(100, 150, 200)),
    alpha = list(alp),
    eta = list(
        c(0, 4.173, 0),
        c(0, 4.173, 0.902),
        c(0, 4.173, 1.827),
        c(0, 4.173, 2.870),
        c(0, 4.173, 4.173)
    ),
    sorted = FALSE
)


# ------ Model ------

mmod <- asthmasims:::compile_cmdstanr_mod("cumulative_logistic_smooth3")
Xsp <- generate_spline_basis(1:100, intercept = TRUE)
Xsp <- sweep(Xsp, 2, colMeans(Xsp))
mdat <- list(
    P = 3,
    beta_sd = 2.5 / sqrt(2),
    prior = 0,
    K = 101,
    P_Isp = ncol(Xsp),
    X_Isp = Xsp,
    y_weights = 0:100,
    alpha_int_sd = 100,
    theta_sd = rep(100, ncol(Xsp))
)
mod <- list(mmod, mdat)


# ------ Loop over configurations and save results ------
run_row <- seq_len(nrow(cfg))
for (z in run_row) {
    start_time <- Sys.time()

    res <- pbmclapply(seq_len(cfg[z][["sims"]]), function(j) {
        sim_asthma_trial(
            mod,
            n_seq = cfg[z][["n_seq"]][[1]],
            alpha = cfg[z][["alpha"]][[1]],
            eta = cfg[z][["eta"]][[1]],
            chains = 1,
            iter_warmup = 500,
            iter_sampling = 1000,
            adapt_delta = 0.98,
            show_messages = FALSE,
            refresh = 0
        )
    }, mc.cores = num_cores, ignore.interactive = TRUE)

    resl_alpha <- rbindlist(lapply(res, \(x) x[["alpha"]]), idcol = "trial")
    resl_contr <- rbindlist(lapply(res, \(x) x[["contr"]]), idcol = "trial")
    resl_trial <- rbindlist(lapply(res, \(x) x[["trial"]]), idcol = "trial")
    resl_yobs <- rbindlist(lapply(res, \(x) x[["yobs"]]), idcol = "trial")
    resl_alpha[, analysis := as.numeric(analysis)]
    resl_contr[, analysis := as.numeric(analysis)]
    resl_trial[, analysis := as.numeric(analysis)]
    resl_yobs[, analysis := as.numeric(analysis)]

    end_time <- Sys.time()

    qsave(
        list(
            cfg = cfg[z],
            alpha = resl_alpha,
            contr = resl_contr,
            trial = resl_trial,
            yobs = resl_yobs,
            runtime = end_time - start_time
        ),
        paste0(
            "out/",
            file_name,
            formatC(z, width = 2, flag = "0"),
            ".qs"
        )
    )
}
