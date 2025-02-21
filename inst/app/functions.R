require(deSolve)

object_usage_linter(interpret_glue = TRUE, skip_with = TRUE)

get_final_size_analytic <- function(r_init, i_init, v_init, n, r_0, a, eps, q) {
  if (sum(i_init) == 0)
    i_init <- n / sum(n)

  s_init <- n - r_init - i_init - v_init
  fn <- (1 - eps) * a * n
  f <- fn / sum(fn)
  cij <- diag(eps) + outer((1 - eps), f)

  r_0i <- r_0 / eigen(a * q * cij)$values[1] * a * q

  z_rhs <- function(z) {
    c(s_init * (1 - exp(-r_0i * cij %*% ((
      z + i_init
    ) / n))))
  }

  optfn <- function(x) {
    ifelse(all(x > 0), max(abs(x - z_rhs(x))), Inf)
  }

  opt_val <- Inf
  while (opt_val > 1) {
    opt <- optim((0.01 + 0.98 * runif(length(n))) * n, optfn)
    opt_val <- opt$value
  }
  opt$par + i_init + r_init
}


exposure_sir <- function(time, state, pars) {
  with(as.list(c(time, state, pars)), {
    s <- c(S1, S2)
    i <- c(I1, I2)
    r <- c(R1, R2)

    n <- c(N1, N2)

    # beta is a 2x2 transmission matrix
    beta <- (1 - epsilon) *
      outer(activities, activities) / sum(c(N1, N2) * activities) +
      epsilon * activities / c(N1, N2) * diag(2)

    beta <- beta / scaling_factor

    d_s <- -(beta %*% i) * s
    d_i <- (beta %*% i) * s - gam * i
    d_r <- gam * i

    return(list(c(d_s, d_i, d_r)))
  })
}


rescale_r_0 <- function(beta, gam, pop_p, n, r_0_value) {
  ngm_unscaled <- n * pop_p * beta * 1 / gam
  dom_eigen <- as.numeric(eigen(ngm_unscaled)$values[1])
  ngm_scaled <- ngm_unscaled / dom_eigen * r_0_value
  scaling_factor <- beta / (ngm_scaled / (n * pop_p) * gam)

  return(scaling_factor)
}


sim_exposure_sir <- function(r_init, i_init, v_init, tm,
    n, r_0, gam, a, eps, q) {
  n_tot <- sum(n)

  # beta is a 2x2 transmission matrix
  beta <- (1 - eps) * outer(a, a) / sum(n * a) +
    eps * a / (n) * diag(2)

  scaling_factor <- rescale_r_0(
    beta = beta,
    gam = gam,
    pop_p = n / n_tot,
    n = n_tot,
    r_0_value = r_0
  )

  pars <- list(
    n = n_tot,
    N1 = n[1],
    N2 = n[2],
    activities = a,
    epsilon = eps,
    gam = gam,
    scaling_factor = scaling_factor
  )

  # define the times vector

  # nolint times <- seq(0, tm, len=1000) #use this if you need plotting points
  times <- c(0, tm)

  # define the initial conditions
  s_init <- c(S1 = n[1] - r_init[1] - i_init[1] - v_init[1],
    S2 = n[2] - r_init[2] - i_init[2] - v_init[2])
  i_init <- c(I1 = i_init[1], I2 = i_init[2])
  r_init <- c(R1 = r_init[1], R2 = r_init[2])

  if (sum(i_init) == 0) {
    i_init <- c(I1 = n[1] / n_tot, I2 = n[2] / n_tot)
    s_init <- s_init - i_init
  }

  y0 <- c(s_init, i_init, r_init)

  # the ode() command solves the ODEs,
  # and stores the data in a dataframe called
  # "simulation"
  simulation <- as.data.frame(deSolve::ode(y0, times, exposure_sir, pars))

  simulation
}

#' @export
getFinalSize <- function(vactime, # nolint
                         vac_portion,
                         pop_size,
                         r_0,
                         recovery_rate,
                         contact_ratio,
                         contact_within_group,
                         susc_ratio) {
  # vactime: time after first case at which all vaccinations are delivered
  # vac_portion: fraction of each population vaccinated
  # pop_size: size of each population
  # r_0: overall basic reproduction number
  # recovery_rate: inverse of mean infectious period
  # # (same time units as vactime)
  # contact_ratio: ratio of 2nd group's : 1st group's overall contact rate
  # contact_within_group: fraction of each group's contacts that are
  # exclusively within group
  # susc_ratio: ratio of 2nd group's : 1st group's susceptibility to
  # infection per contact



  i_sim1 <- c(0, 0)
  r_sim1 <- c(0, 0)

  if (vactime > 0) {
    sim1 <- sim_exposure_sir(
      r_init = c(0, 0),
      i_init = c(0, 0),
      v_init = c(0, 0),
      tm = vactime,
      n = pop_size,
      r_0 = r_0,
      gam = recovery_rate,
      a = c(1, contact_ratio),
      eps = contact_within_group,
      q = c(1, susc_ratio)
    )

    i_sim1 <- as.numeric(sim1[nrow(sim1), c("I1", "I2")])
    r_sim1 <- as.numeric(sim1[nrow(sim1), c("R1", "R2")])
  }

  get_final_size_analytic(
    r_init = r_sim1,
    i_init = i_sim1,
    v_init = pop_size * vac_portion,
    n = pop_size,
    r_0 = r_0,
    a = c(1, contact_ratio),
    eps = contact_within_group,
    q = c(1, susc_ratio)
  )

}
