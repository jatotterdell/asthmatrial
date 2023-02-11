#' @title update_model
#' @description Update the model used in clarity2sims
#' @param mod A compiled cmdstanr stanmodel
#' @param moddat Model data for use with mod
#' @param ... Additional arguments to sample
#' @return Posterior draws
#' @import posterior
update_model <- function(mod, moddat, ...) {
  drp <- utils::capture.output(fit <- mod$sample(
    data = moddat, ...
  ))
  draws <- posterior::as_draws_rvars(fit$draws(c("alpha", "beta", "mu")))
  diags <- unlist(fit$diagnostic_summary())
  return(list(draws = draws, diags = diags))
}


#' @title sim_asthma_trial
#' @description Simulate asthma trial
#' @details
#' Simulates asthma trial.
#' @param mod A list (model, model data)
#' @param n_seq Sequence of interim analysis sample sizes
#' @param p_assign Assignment probabilities to treatment arms
#' @param alpha True intercept parameter
#' @param eta True eta parameter
#' @param eff_eps Effectiveness threshold
#' @param ... Other arguments to cmdstanr::sample, e.g. adapt_delta, chains, etc.
#' @return A list of trial related data.tables
#' @export
sim_asthma_trial <- function(mod,
                             n_seq = 200,
                             p_assign = rep(1 / 3, 3),
                             alpha = stats::qlogis(cumsum(rep(1 / 101, 101)))[-101],
                             eta = c(0, 0, 0),
                             eff_eps = 0.975,
                             ...) {
  P <- length(p_assign)
  K <- length(alpha) + 1
  n_max <- max(n_seq)
  n_new <- diff(c(0, n_seq))
  n_int <- length(n_seq)


  # Storage
  labs <- list("analysis" = seq_len(n_int), "arm" = seq_len(P))
  n_obs <- matrix(0, n_int, P, dimnames = labs)
  y_obs <- array(0,
    dim = c(n_int, P, K),
    dimnames = list(
      "analysis" = seq_len(n_int),
      "arm" = seq_len(P),
      "level" = stringr::str_pad(seq_len(K) - 1, max(nchar(seq_len(K) - 1)), pad = "0")
    )
  )
  e_alpha <- matrix(0, n_int, K - 1,
    dimnames = list(labs[[1]], c("cut" = 1:(K - 1)))
  )
  v_alpha <- e_alpha
  e_beta <- matrix(0, n_int, P, dimnames = labs)
  v_beta <- e_beta
  e_mu <- e_beta
  v_mu <- e_beta
  e_ctr <- matrix(0, n_int, 3,
    dimnames = list(
      "analysis" = seq_len(n_int),
      "contrast" = c("3v1", "3v2", "2v1")
    )
  )
  diags <- matrix(
    0, n_int, 3,
    dimnames = list(
      "analysis" = seq_len(n_int),
      "diagnostic" = c("num_divergent", "num_max_treedepth", "ebfmi")
    )
  )
  pr_ctr_lt10 <- pr_ctr_lt0 <- hi_ctr <- lo_ctr <- v_ctr <- e_ctr

  # Generate data
  x <- permuted_block_rand(p_assign, n_max, 2 * P)[["trt"]]
  y <- generate_data_eta(alpha, eta[x])

  for (i in seq_len(n_int)) {
    idx <- seq_len(n_seq[i])
    xc <- x[idx]
    yc <- y[idx]
    n_obs[i, ] <- table(xc)
    y_obs[i, , ] <- table(xc, factor(yc - 1, levels = 1:K - 1))
    mod[[2]]$x <- xc
    mod[[2]]$y <- yc
    mod[[2]]$yc <- yc - 1
    mod[[2]]$N <- n_seq[i]

    fitdat <- update_model(mod[[1]], mod[[2]], ...)
    fit <- fitdat$draws
    diags[i, ] <- fitdat$diags
    # Define contrasts
    ctr <- c(
      fit$mu[3] - fit$mu[1],
      fit$mu[3] - fit$mu[2],
      fit$mu[2] - fit$mu[1]
    )
    # Summarise contrasts
    e_ctr[i, ] <- posterior::E(ctr)
    v_ctr[i, ] <- posterior::var(ctr)
    tmp <- posterior::quantile2(ctr, c(0.025, 0.975))
    lo_ctr[i, ] <- tmp[1, ]
    hi_ctr[i, ] <- tmp[2, ]
    pr_ctr_lt0[i, ] <- posterior::Pr(ctr < 0)
    pr_ctr_lt10[i, ] <- posterior::Pr(ctr < 10)
    # Other summaries
    e_alpha[i, ] <- posterior::E(fit$alpha)
    v_alpha[i, ] <- posterior::var(fit$alpha)
    e_beta[i, ] <- posterior::E(fit$beta)
    v_beta[i, ] <- posterior::var(fit$beta)
    e_mu[i, ] <- posterior::E(fit$mu)
    v_mu[i, ] <- posterior::var(fit$mu)
  }

  out_alpha <- list(
    e_alpha = e_alpha,
    v_alpha = v_alpha
  )
  out_arm <- list(
    n_obs = n_obs,
    e_beta = e_beta,
    v_beta = v_beta,
    e_mu = e_mu,
    v_mu = v_mu
  )
  out_ctr <- list(
    e_ctr = e_ctr,
    v_ctr = v_ctr,
    lo_ctr = lo_ctr,
    hi_ctr = hi_ctr,
    pr_ctr_lt0 = pr_ctr_lt0,
    pr_ctr_lt10 = pr_ctr_lt10
  )
  out_diags <- list(
    diags = diags
  )

  return(list(
    alpha = list_as_dt(out_alpha),
    trial = list_as_dt(out_arm),
    contr = list_as_dt(out_ctr),
    diags = list_as_dt(out_diags),
    yobs = dcast(
      array_as_dt(y_obs), ... ~ paste0("y", level),
      value.var = "y_obs"
    )
  ))
}
