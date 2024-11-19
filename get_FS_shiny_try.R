## dependencies
require(deSolve)
require(shiny)
## dependent functions

#' @export
#' @import shiny
getFinalSizeAnalytic <- function(Rinit, Iinit, Vinit, N, R0, a, eps, q) {

  Sinit <- N - Rinit - Iinit - Vinit
  fn <- (1 - eps) * a * N
  f <- fn / sum(fn)

  cij <- diag(eps) + outer((1 - eps), f)

  R0i <- R0 / eigen(a * q * cij)$values[1] * a * q

  Zrhs <- function(Z) c(Sinit * (1 - exp(-R0i * cij %*% ((Z + Iinit) / N))))

  optfn <- function(x) ifelse(all(x / N > 0.001), max(abs(x - Zrhs(x))), Inf)

  optVal <- Inf
  counter <- 0
  while (optVal > 1 & counter < 5000) {
    opt <- optim((0.1 + 0.8 * runif(length(N))) * N, optfn)
    optVal <- opt$value
    counter <- counter + 1
  }
  if (counter >= 5000) {
    return("Error")
  } else {
    opt$par + Iinit + Rinit
  }
}

exposure.SIR <- function(Time, state, Pars) {
  with(as.list(c(Time, state, Pars)), {

    S <- c(S1, S2)
    I <- c(I1, I2)
    R <- c(R1, R2)

    N <- c(N1, N2)

    # beta is a 2x2 transmission matrix
    beta <- (1 - epsilon) * outer(activities, activities) / sum(c(N1, N2) * activities) +
      epsilon * activities / c(N1, N2) * diag(2)

    beta <- beta / scaling.factor

    dS <- -(beta %*% I) * S
    dI <- (beta %*% I) * S - gam * I
    dR <- gam * I

    return(list(c(dS, dI, dR)))
  })
}

rescale.R0 <- function(beta, gam, pop.p, N, R0.value) {
  NGM.unscaled <- N * pop.p * beta * 1 / gam
  dom.eigen <- as.numeric(eigen(NGM.unscaled)$values[1])
  NGM.scaled <- NGM.unscaled / dom.eigen * R0.value
  scaling.factor <- beta / (NGM.scaled / (N * pop.p) * gam)

  return(scaling.factor)
}

sim.exposure.SIR <- function(Rinit, Iinit, Vinit, tm, N, R0, gam, a, eps, q) {

  Ntot <- sum(N)

  # beta is a 2x2 transmission matrix
  beta <- (1 - eps) * outer(a, a) / sum(N * a) +
    eps * a / (N) * diag(2)

  scaling.factor <- rescale.R0(beta = beta, gam = gam, pop.p = N / Ntot, N = Ntot,
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
  simulation <- as.data.frame(ode(y0, times, exposure.SIR, pars))

  simulation
}

getFinalSize <- function(vacTime, vacPortion, popSize, R0, recoveryRate,
                         contactRatio, contactWithinGroup, suscRatio) {
  # vacTime: time after first case at which all vaccinations are delivered
  # vacPortion: fraction of each population vaccinated
  # popSize: size of each population
  # R0: overall basic reproduction number
  # recoveryRate: inverse of mean infectious period (same time units as vacTime)
  # contactRatio: ratio of 2nd group's : 1st group's overall contact rate
  # contactWithinGroup: fraction of each group's contacts that are exlusively within group
  # suscRatio: ratio of 2nd group's : 1st group's susceptibility to infection per contact

  Isim1 <- c(0, 0)
  Rsim1 <- c(0, 0)

  if (vacTime >= 0) {
    sim1 <- sim.exposure.SIR(Rinit = c(0, 0), Iinit = c(0, 0), Vinit = c(0, 0), tm = vacTime,
      N = popSize, R0 = R0, gam = recoveryRate, a = c(1, contactRatio),
      eps = contactWithinGroup, q = c(1, suscRatio))

    Isim1 <- as.numeric(sim1[nrow(sim1), c("I1", "I2")])
    Rsim1 <- as.numeric(sim1[nrow(sim1), c("R1", "R2")])
  }
  if (any(is.na(c(Isim1, Rsim1)))) {
    return("Error")
  } else {
    fsout <- try(getFinalSizeAnalytic(Rinit = Rsim1, Iinit = Isim1, Vinit = popSize * vacPortion, N = popSize, R0 = R0,
      a = c(1, contactRatio), eps = contactWithinGroup, q = c(1, suscRatio)))
    if (any(grepl("Error", fsout))) {
      return("Error")
    } else {
      return(fsout)
    }
  }
}


## helper function for layout
textInput2 <- function(inputId, label, value = "", ...) {
  tagList(tags$label(label, `for` = inputId), tags$input(id = inputId,
    type = "text", value = value, ...))
}


## shiny app
ui <- fixedPage(
  titlePanel("Model Title"),

  fluidRow(
    textInput2("vacPortion_a", label = "fraction vaccinated a = ", value = "0.1", class = "input-small"),

    textInput2("vacPortion_b", label = "fraction vaccinated b = ", value = "0.1", class = "input-small"),
  ),
  fluidRow(
    textInput2("vacTime", label = "time after first case at which all vaccinations are delivered = ", value = "100", class = "input-small"),
  ),
  fluidRow(
    textInput2("popSize_a", label = "a popsize = ", value = "800000", placeholder = "", class = "input-small"),
    textInput2("popSize_b", label = "b popsize = ", value = "200000", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("R0", label = "overall basic reproduction number = ", value = "1.5", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("recoveryRate", label = "inverse of mean infectious period = ", value = "0.14", placeholder = "", class = "input-small"),
  ),
  fluidRow(

    textInput2("contactRatio", label = "ratio of 2nd group's : 1st group's overall contact rate = ", value = "1.7", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("contactWithinGroup_a", label = "within a group mixing fraction = ", value = "0.4", placeholder = "", class = "input-small"),
    textInput2("contactWithinGroup_b", label = "within b group mixing fraction = ", value = "0.4", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("suscRatio", label = "ratio of 2nd group's : 1st group's susceptibility to infection per contact = ", value = "1", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("amountToSpend", label = "amount available to spend on vaccine rollout = ", value = "1e5", placeholder = "", class = "input-small"),
  ),
  fluidRow(
    textInput2("vaccineCostRatio", label = "ratio of 2nd group's : 1st group's per-individual cost to vaccinate = ", value = "1", placeholder = "", class = "input-small"),
  ),
  # "then R0 = ", textOutput("r0",inline=T),
  # ", total people infected = ", textOutput("infs",inline=T),
  # ", and vaccine doses given = ", textOutput("vaxs",inline=T),

  plotOutput("plot", click = "plot_click", )
)


server <- function(input, output, session) {
  # output$infs<-renderText(round(SIRV.model(t=input$time,V1=input$V1,I1=as.numeric(input$Iinit),b=input$beta,g=input$gamma,m=input$mu)$R[
  #   length(SIRV.model(t=input$time,V1=input$V1,I1=as.numeric(input$Iinit),b=input$beta,g=input$gamma,m=input$mu)$R)]
  #   * 1726765)) ## the big number is the population size
  # output$r0<-renderText(as.character(round(as.numeric(input$beta)/as.numeric(input$gamma),2)))
  # output$vaxs<-renderText(as.character(round(round(max(SIRV.model(t=input$time,V1=input$V1,I1=as.numeric(input$Iinit),b=input$beta,g=input$gamma,m=input$mu)$V)* 1726765 )-
  #                                              round(min(SIRV.model(t=input$time,V1=input$V1,I1=as.numeric(input$Iinit),b=input$beta,g=input$gamma,m=input$mu)$V)* 1726765))
  # ))
  output$plot <- renderPlot(
    {
      popSize <- c(as.numeric(input$popSize_a), as.numeric(input$popSize_b))
      R0 <- as.numeric(input$R0)
      recoveryRate <- as.numeric(input$recoveryRate)
      contactRatio <- as.numeric(input$contactRatio)
      contactWithinGroup <- c(as.numeric(input$contactWithinGroup_a), as.numeric(input$contactWithinGroup_b))
      suscRatio <- as.numeric(input$suscRatio)
      vacP <- c(as.numeric(input$vacPortion_a), as.numeric(input$vacPortion_b))
      vacTime <- as.numeric(input$vacTime)
      amountToSpend <- as.numeric(input$amountToSpend)
      vaccineCostRatio <- as.numeric(input$vaccineCostRatio)

      ## function using variables stored in environment
      getFS0 <- function(vacP) getFinalSize(vacTime = 0, vacPortion = vacP, popSize = popSize, R0 = R0,
        recoveryRate = recoveryRate, contactRatio = contactRatio,
        contactWithinGroup = contactWithinGroup, suscRatio = suscRatio)

      ## get baseline (needed to be in environment for following functions to run?)
      # fsNoVax <- getFS0(vacP = c(0,0))
      ## just embeded getFS0 in the getStats function below, instead
      getStats <- function(vacP) {
        fs <- getFS0(vacP)
        infPrev <- getFS0(vacP = c(0, 0)) - fs

        list(finalSize = fs, finalSizeProportion = fs / popSize, infectionsPrevented = infPrev)
      }

      getstatsout <- try(getStats(vacP))
      # getstatsout<-""

      if (any(grepl("Error", getstatsout))) {
        plot(1, 1, main = "Error!")
      } else {
        ## define functions (which use variables in environment)

        # plotStats <- function(numVaccines){
        #   v2max <- min(numVaccines,popSize[2])
        #   v2 <- seq(0,v2max,len=100)
        #   v1 <- numVaccines - v2
        #
        #   fs <- matrix(0,100,2)
        #   fsp <- fs
        #   ip <- fs
        #   equityIndex <- rep(0,100)
        #
        #   for(i in 1:100){
        #     stats <- getStats(c(v1[i],v2[i])/popSize)
        #     fs[i,] <- stats$finalSize
        #     fsp[i,] <- stats$finalSizeProportion
        #     ip[i,] <- stats$infectionsPrevented
        #
        #     equityIndex[i] <- min(fsp[i,])/max(fsp[i,])
        #   }
        #
        #
        #   plot(v2/numVaccines, ip[,1]+ip[,2], type='l', ylim = c(0,max(ip[,1]+ip[,2])),
        #        xlab='Proportion of vaccines to minority group',
        #        ylab='Infections prevented (red maj, green min)')
        #   lines(v2/numVaccines, ip[,1],col='red')
        #   lines(v2/numVaccines, ip[,2],col='green')
        #
        #   plot(v2/numVaccines, equityIndex, type='l', xlab='Proportion of vaccines to minority group')
        # }

        plotStatsCost <- function(amountToSpend, vaxCosts) {

          spend2max <- min(amountToSpend, popSize[2] * vaxCosts[2])
          spend2 <- seq(0, spend2max, len = 100)
          spend1 <- amountToSpend - spend2

          v1 <- spend1 / vaxCosts[1]
          v2 <- spend2 / vaxCosts[2]

          fs <- matrix(0, 100, 2)
          fsp <- fs
          ip <- fs
          ipPerCost <- fs
          equityIndex <- rep(0, 100)

          for (i in 1:100) {
            stats <- getStats(c(v1[i], v2[i]) / popSize)
            fs[i, ] <- stats$finalSize
            fsp[i, ] <- stats$finalSizeProportion
            ip[i, ] <- stats$infectionsPrevented

            equityIndex[i] <- min(fsp[i, ]) / max(fsp[i, ])
          }
          par(mfrow = c(1, 2))
          plot(spend2 / amountToSpend, ip[, 1] + ip[, 2], type = "l", ylim = c(0, max(ip[, 1] + ip[, 2])), lwd = 2,
            xlab = "Proportion of spending to minority group",
            ylab = "Infections prevented (red maj, green min)")
          lines(spend2 / amountToSpend, ip[, 1], col = "red", lwd = 2)
          lines(spend2 / amountToSpend, ip[, 2], col = "green", lwd = 2)

          plot(spend2 / amountToSpend, equityIndex, type = "l", xlab = "Proportion of spending to minority group",
            ylim = c(0, 1), lwd = 2, ylab = "maj prev (red); min prev (green); equity (black)")

          lines(spend2 / amountToSpend, fsp[, 1], col = "red", lty = 2, lwd = 2)
          lines(spend2 / amountToSpend, fsp[, 2], col = "green", lty = 2, lwd = 2)

          list(majFS = fs[, 1], minFS = fs[, 2], majSpend = spend1, minSpend = spend2)
          par(mfrow = c(1, 1))

        }

        ## now run

        plotStatsCost(amountToSpend,
          c(1,
            vaccineCostRatio))

      }
      # ticks<-signif(seq(0,signif(max(getFSout),3),length.out=10),2)
      # axis(side = 2, at = ticks,las=1)


    },
    res = 96)
}

#' Runs the vaccine model
#' @export
run_model <- function() {
  shinyApp(ui, server)
}
