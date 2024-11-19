
#' @importFrom deSolve ode
#' @importFrom stats optim runif


getFinalSizeAnalytic <- function(Rinit, Iinit, Vinit, N, R0, a, eps, q) {

  if (sum(Iinit) == 0) Iinit <- N / sum(N)

  Sinit <- N - Rinit - Iinit - Vinit
  fn <- (1 - eps) * a * N
  f <- fn / sum(fn)

  cij <- diag(eps) + outer((1 - eps), f)

  R0i <- R0 / eigen(a * q * cij)$values[1] * a * q

  Zrhs <- function(Z) c(Sinit * (1 - exp(-R0i * cij %*% ((Z + Iinit) / N))))

  optfn <- function(x) ifelse(all(x > 0), max(abs(x - Zrhs(x))), Inf)

  optVal <- Inf
  while (optVal > 1) {
    opt <- optim((0.01 + 0.98 * runif(length(N))) * N, optfn)
    optVal <- opt$value
  }
  opt$par + Iinit + Rinit
}

#' @export
exposure.SIR <- function(Time, state, Pars) {
  with(as.list(c(Time, state, Pars)), {

    S <- c(S1, S2)
    I <- c(I1, I2)
    R <- c(R1, R2)

    N <- c(N1, N2)

    # beta is a 2x2 transmission matrix
    beta <- (1 - epsilon) *
      outer(activities, activities) / sum(c(N1, N2) * activities) +
      epsilon * activities / c(N1, N2) * diag(2)

    beta <- beta / scaling.factor

    dS <- -(beta %*% I) * S
    dI <- (beta %*% I) * S - gam * I
    dR <- gam * I

    return(list(c(dS, dI, dR)))
  })
}

#' @export
rescale.R0 <- function(beta, gam, pop.p, N, R0.value) {
  NGM.unscaled <- N * pop.p * beta * 1 / gam
  dom.eigen <- as.numeric(eigen(NGM.unscaled)$values[1])
  NGM.scaled <- NGM.unscaled / dom.eigen * R0.value
  scaling.factor <- beta / (NGM.scaled / (N * pop.p) * gam)

  return(scaling.factor)
}

#' @export
sim.exposure.SIR <- function(Rinit, Iinit, Vinit, tm, N, R0, gam, a, eps, q) {

  Ntot <- sum(N)

  # beta is a 2x2 transmission matrix
  beta <- (1 - eps) * outer(a, a) / sum(N * a) +
    eps * a / (N) * diag(2)

  scaling.factor <- rescale.R0(beta = beta, gam = gam,
    pop.p = N / Ntot, N = Ntot,
    R0.value = R0)

  pars <- list(N = Ntot, N1 = N[1], N2 = N[2],
    activities = a, epsilon = eps, gam = gam,
    scaling.factor = scaling.factor)

  # define the times vector
  # times <- seq(0, tm, len=1000) #use this if you need plotting points
  times <- c(0, tm)

  # define the initial conditions
  S.init <- c(S1 = N[1] - Rinit[1] - Iinit[1] - Vinit[1], S2 = N[2] - Rinit[2] - Iinit[2] - Vinit[2])
  I.init <- c(I1 = Iinit[1], I2 = Iinit[2])
  R.init <- c(R1 = Rinit[1], R2 = Rinit[2])

  if (sum(I.init) == 0) {
    I.init <- c(I1 = N[1] / Ntot, I2 = N[2] / Ntot)
    S.init <- S.init - I.init
  }

  y0 <- c(S.init, I.init, R.init)

  # the ode() command solves the ODEs, and stores the data in a dataframe called
  # "simulation"
  simulation <- as.data.frame(deSolve::ode(y0, times, exposure.SIR, pars))

  simulation
}

#'getFinalSize
#'
#'This function calculates the final sizes of the outbreaks.
#'
#'Run this function to run the model.
#'
#' @param vacTime time after first case at which all vaccinations are delivered
#' @param vacPortion fraction of each population vaccinated
#' @param popSize the size of each population
#' @param R0 overall basic reproduction number
#' @param recoveryRate inverse of mean infectious period (same time units as vacTime)
#' @param contactRatio ratio of 2nd group's : 1st group's overall contact rate
#' @param contactWithinGroup fraction of each group's contacts that are exclusively within group
#' @param suscRatio ratio of 2nd group's : 1st group's susceptibility to infection per contact
#' @export
getFinalSize <- function(vacTime, vacPortion, popSize, R0, recoveryRate,
                         contactRatio, contactWithinGroup, suscRatio) {


  Isim1 <- c(0, 0)
  Rsim1 <- c(0, 0)

  if (vacTime > 0) {
    sim1 <- sim.exposure.SIR(Rinit = c(0, 0), Iinit = c(0, 0), Vinit = c(0, 0), tm = vacTime,
      N = popSize, R0 = R0, gam = recoveryRate, a = c(1, contactRatio),
      eps = contactWithinGroup, q = c(1, suscRatio))

    Isim1 <- as.numeric(sim1[nrow(sim1), c("I1", "I2")])
    Rsim1 <- as.numeric(sim1[nrow(sim1), c("R1", "R2")])
  }
  getFinalSizeAnalytic(Rinit = Rsim1, Iinit = Isim1,
    Vinit = popSize * vacPortion, N = popSize, R0 = R0,
    a = c(1, contactRatio), eps = contactWithinGroup, q = c(1, suscRatio))

}

#' @export
finalSizeExample <- getFinalSize(vacTime = 100, vacPortion = c(0.1, 0.1),
  popSize = c(800000, 200000), R0 = 1.5, recoveryRate = 1 / 7,
  contactRatio = 1.7, contactWithinGroup = c(0.4, 0.4), suscRatio = 1)
