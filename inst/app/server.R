library(shiny)



server <- function(input, output, session) {
  get_table <- function(input) {
    pop_size <- c(
      as.numeric(input$pop_size_a),
      as.numeric(input$pop_size_b)
    )
    r_0 <- as.numeric(input$r_0)
    recovery_rate <- as.numeric(input$recovery_rate)
    contact_ratio <- as.numeric(input$contact_ratio)
    contact_within_group <- c(
      as.numeric(input$contact_within_group_a),
      as.numeric(input$contact_within_group_b)
    )
    susc_ratio <- as.numeric(input$susc_ratio)
    vac_p <- c(
      as.numeric(input$vac_portion_a),
      as.numeric(input$vac_portion_b)
    )
    vac_time <- as.numeric(input$vactime)
    hosp_prob <- c(input$hosp_prob_a, input$hosp_prob_b)
    hosp_death <- c(input$hosp_deathProb_a, input$hosp_deathProb_b)
    nonhosp_death <- c(
      input$nonhosp_deathProb_a,
      input$nonhosp_deathProb_b
    )

    # input validation
    shiny::validate(
      shiny::need(
        all(vac_p >= 0 & vac_p <= 1),
        "Vaccination portion must be between 0 and 1"
      ),
      shiny::need(
        all(hosp_prob >= 0 & hosp_prob <= 1),
        "Hospitalization probability must be between 0 and 1"
      ),
      shiny::need(
        all(hosp_death >= 0 & hosp_death <= 1),
        "Hospitalization death probability must be between 0 and 1"
      ),
      shiny::need(
        all(nonhosp_death >= 0 & nonhosp_death <= 1),
        "Non-hospitalization death probability must be between 0 and 1"
      )
    )

    fs <- getFinalSize(
      vactime = vac_time,
      vac_portion = vac_p,
      pop_size = pop_size,
      r_0 = r_0,
      recovery_rate = recovery_rate,
      contact_ratio = contact_ratio,
      contact_within_group = contact_within_group,
      susc_ratio = susc_ratio
    )

    hosp <- hosp_prob * fs
    death <- hosp * hosp_death + (1 - hosp) * nonhosp_death

    num_vax <- vac_p * pop_size
    tbl_vax <- c(sum(num_vax), num_vax)
    tbl_vax_cov <- tbl_vax / c(sum(pop_size), pop_size)
    tbl_inf <- c(sum(fs), fs)
    tbl_inf_pref <- tbl_inf / c(sum(pop_size), pop_size)
    tbl_hosp <- c(sum(hosp), hosp)
    tbl_hosp_prev <- tbl_hosp / c(sum(pop_size), pop_size)
    tbl_death <- c(sum(death), death)
    tbl_death_prev <- tbl_death / c(sum(pop_size), pop_size)

    tbl <- rbind(
      vaccines = tbl_vax,
      vaccinationPercentage = 100 * tbl_vax_cov,
      infections = tbl_inf,
      infectionPercentage = 100 * tbl_inf_pref,
      hospitalizations = tbl_hosp,
      hospitalizationPercentage = 100 * tbl_hosp_prev,
      deaths = tbl_death,
      deathPercentage = 100 * tbl_death_prev
    )

    colnames(tbl) <- c("Total", "GroupA", "GroupB")
    as.data.frame(tbl)
  }

  output$table <- shiny::renderTable(get_table(input),
    rownames = TRUE, digits = 2)
}
