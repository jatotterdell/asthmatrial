# Packages ----

suppressPackageStartupMessages({
    library(asthmasims)
    library(cmdstanr)
    library(posterior)
    library(extraDistr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggdist)
    library(splines2)
})


# Misc ----

theme_set(theme_bw(base_size = 10, base_family = "Palatino"))

labs <- c("Standard", "Maximum", "Investigational")

rvar_inv_ord_logit <- function(a) {
    cp <- exp(a) / (1 + exp(a))
    p <- diff(c(rvar(0), cp, rvar(1)))
    return(p)
}


# Model ----

mod <- asthmasims:::compile_cmdstanr_mod("cumulative_logistic_smooth3")


# DGP Parameters ----

pt1 <- dnorm(0:100, 100, 10)
pt1 <- pt1 / sum(pt1)
alp1 <- qlogis(cumsum(pt1))[-length(pt1)]
pt2 <- dnorm(0:100, 50, 10)
pt2 <- pt2 / sum(pt2)
alp2 <- qlogis(cumsum(pt2))[-length(pt2)]
pt3 <- dnorm(0:100, 70, 20)
pt3 <- pt3 / sum(pt3)
alp3 <- qlogis(cumsum(pt3))[-length(pt3)]



# Simulate observed data ----

set.seed(5136)
n <- 33
x <- rep(1:3, each = n)
Xdes <- contr.treatment(3)
y <- c(
    asthmasims:::generate_data_eta(alp1, rep(0, n)),
    asthmasims:::generate_data_eta(alp2, rep(0, n)),
    asthmasims:::generate_data_eta(alp3, rep(0, n))
)
y_weights <- y - 1
iknots <- 9
knots <- floor(seq(1, 100, length.out = iknots + 2))[2:(iknots + 1)]
Xsp <- iSpline(1:100, knots = knots, degree = 3)
Xsp <- sweep(Xsp, 2, colMeans(Xsp))
p <- aggregate(y, by = list(x), \(z) prop.table(table(factor(z, 0:100))))


# Frames for plotting ----

dtru <- tibble(
    x = rep(1:3, each = 101),
    xlab = factor(rep(labs, each = 101), levels = labs),
    y = rep(0:100, times = 3),
    p = c(pt1, pt2, pt3)
) |>
    group_by(x) |>
    mutate(cp = cumsum(p)) |>
    ungroup()

dtrum <- dtru |>
    group_by(xlab) |>
    summarise(m = sum(y * p))

dobs <- tibble(
    x = rep(1:3, each = 101),
    xlab = factor(rep(labs, each = 101), levels = labs),
    y = rep(0:100, times = 3),
    p = c(p[1, -1], p[2, -1], p[3, -1]),
) |>
    group_by(x) |>
    mutate(cp = cumsum(p)) |>
    ungroup()


# Update model ----

mdat <- list(
    N = n * 3,
    P = 3,
    x = x,
    y = y,
    beta_sd = 2.5 / sqrt(2),
    K = 101,
    P_Isp = ncol(Xsp),
    X_Isp = Xsp,
    y_weights = 0:100,
    alpha_int_sd = 100,
    theta_scale = 10,
    theta_sd = rep(100, ncol(Xsp)),
    prior = 0
)

snk <- capture.output(res <- mod$sample(
    data = mdat,
    chains = 3,
    parallel_chains = 3,
    adapt_delta = 0.98,
    refresh = 0
))

drws <- as_draws_rvars(res$draws(c("alpha", "beta", "mu")))
p1 <- rvar_inv_ord_logit(drws$alpha - drws$beta[1])
p2 <- rvar_inv_ord_logit(drws$alpha - drws$beta[2])
p3 <- rvar_inv_ord_logit(drws$alpha - drws$beta[3])

pd <- tibble(
    y = rep(0:100, times = 3),
    p = c(p1, p2, p3),
    x = rep(1:3, each = 101),
    xlab = factor(rep(labs, each = 101), levels = labs),
    ptru = c(pt1, pt2, pt3),
    cptru = c(cumsum(pt1), cumsum(pt2), cumsum(pt3))
) |>
    group_by(x, y) |>
    mutate(
        med = median(p),
        lo = quantile(p, 0.025),
        hi = quantile(p, 0.975)
    ) |>
    group_by(x) |>
    mutate(cp = cumsum(p)) |>
    ungroup()
pdE <- pd |>
    group_by(x) |>
    summarise(Ey = rvar_sum(y * p)) |>
    ungroup() |>
    mutate(EyDStandard = Ey - Ey[1], EyDMaximum = Ey - Ey[2])


# Figures ----

p1 <- ggplot(dtru, aes(x = y, y = p)) +
    geom_line(aes(colour = factor(xlab))) +
    geom_vline(aes(xintercept = dtrum$m[1]), linetype = 2) +
    geom_vline(aes(xintercept = dtrum$m[2]), linetype = 2, colour = "red") +
    scale_colour_manual("Intervention", values = c("black", "red", "red")) +
    labs(x = "Days with symptoms", y = "Proportion")

p2 <- ggplot(dobs, aes(y, cp, colour = xlab, group = x)) +
    geom_step() +
    labs(
        x = "Days with symptoms",
        y = "Cumulative proportion",
        colour = "Intervention"
    ) +
    scale_colour_viridis_d()

p3 <- ggplot(pd, aes(x = y, ydist = p)) +
    facet_wrap(~xlab, ncol = 1) +
    stat_interval(size = 1) +
    geom_point(aes(y = median(p)), size = 0.5) +
    geom_line(aes(y = ptru), colour = "red") +
    geom_point(data = dobs, aes(x = y, y = p), shape = 21, colour = "black") +
    labs(
        x = "Days with symptoms",
        y = "Proportion",
        colour = "Credible interval"
    ) +
    theme(legend.position = "bottom")

p4 <- pdE |>
    filter(x == 3) |>
    pivot_longer(
        EyDStandard:EyDMaximum,
        names_prefix = "EyD",
        names_to = "Comparison"
    ) |>
    mutate(
        Comparison = paste0("Investigational\nvs\n", Comparison)
    ) |>
    ggplot(data = _, aes(xdist = value, y = Comparison)) +
    stat_halfeye(
        aes(
            fill =
                after_stat(cut_cdf_qi(
                    cdf,
                    .width = c(.5, .8, .95, 0.99),
                    labels = scales::percent_format()
                ))
        ),
        adjust = 1, n = 1001, .width = c(0.5, 0.8, 0.95)
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_brewer(
        palette = "Reds",
        direction = -1,
        na.translate = FALSE
    ) +
    labs(
        x = "Posterior difference in average days with symptoms",
        fill = "Credible\ninterval"
    )
