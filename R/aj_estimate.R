create_breaks_values <- function(probs_vec, stratified_by, by) {
  if (stratified_by != "probability_threshold") {
    breaks <- stats::quantile(probs_vec, probs = rev(seq(0, 1, by = by)))
  } else {
    breaks <- round(
      seq(0, 1, by = by),
      digits = nchar(format(by, scientific = FALSE))
    )
  }

  breaks
}

create_non_parametric_reals_estimates_by_assumptions <- function(
  data_for_tidycmprsk,
  fixed_time_horizons,
  censored_assumption = "adjusted",
  competing_assumption = "adjusted"
) {



  compute_aj_estimates <- function(
    data,
    fixed_time_horizon,
    event_label,
    event_filter
  ) {


    if (!is.na(event_filter)) {
      data <- data |>
        dplyr::mutate(
          is_excluded = (times <= fixed_time_horizon & reals == event_label)
        )
    } else {
      data <- data |>
        dplyr::mutate(is_excluded = FALSE)
    }

    subset_data <- data |>
      dplyr::filter(!is_excluded)


    ajfit <- survival::survfit(
      survival::Surv(
        subset_data$times,
        subset_data$reals
      ) ~
        subset_data$strata,
      istate = rep("real_negatives", length(subset_data$times))
    )

    data_for_exclusion_count <- data |>
      dplyr::filter(is_excluded) |>
      dplyr::group_by(strata) |>
      dplyr::summarise(
        reals_estimate = sum(is_excluded, na.rm = TRUE),
        estimate = sum(is_excluded, na.rm = TRUE) / dplyr::n(),
        .groups = "drop"
      ) |>
      tidyr::complete(
        strata, fill = list(reals_estimate = 0, estimate = 0))

    aj_survival <- summary(
      ajfit,
      times = fixed_time_horizon,
      extend = TRUE
    )$pstate |>
      as.data.frame() |>
      purrr::set_names(ajfit$states) |>
      tibble::add_column(time = fixed_time_horizon) |>
      tibble::add_column(strata = levels(data$strata))|>
      tidyr::pivot_longer(
        cols = !c(strata, time),
        values_to = "estimate",
        names_to = "outcome"
      )

    # print(aj_survival)


    aj_survival <- aj_survival |>
      dplyr::rename(fixed_time_horizon = time) |>
      dplyr::left_join(
        subset_data |>
          dplyr::group_by(strata) |>
          dplyr::summarise(total_included = n())
      ) |>
      dplyr::mutate(reals_estimate = total_included * estimate) |>
      dplyr::select(-total_included)


    if (!is.na(event_filter)) {

      aj_survival <- aj_survival |>
        dplyr::filter(
          !(outcome == event_filter)) |>
        dplyr::add_row(
          data_for_exclusion_count |>
            dplyr::mutate(
              outcome = event_filter,
              fixed_time_horizon = fixed_time_horizon
            )
        )

      if (event_filter != "real_censored") {

        # print("data_for_exclusion_count")
        # print(data_for_exclusion_count)

        aj_survival <- aj_survival |>
          dplyr::add_row(
            tibble::tibble(
              strata = levels(data$strata),
              outcome = "real_censored",
              fixed_time_horizon = fixed_time_horizon,
              estimate = 0,
              reals_estimate = 0
            )
          )
      }
    } else {
      aj_survival <- aj_survival |>
        dplyr::add_row(
          tibble::tibble(
            strata = levels(data$strata),
            outcome = "real_censored",
            fixed_time_horizon = fixed_time_horizon,
            estimate = 0,
            reals_estimate = 0
          )
        )
    }

    aj_survival |>
      dplyr::mutate(
        censored_assumption = censored_assumption,
        competing_assumption = competing_assumption
      )
  }

  if (censored_assumption == "excluded" & competing_assumption == "excluded") {
    return(
      data_for_tidycmprsk |>
        tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons) |>
        dplyr::mutate(outcome = reals) |>
        dplyr::mutate(
          outcome = dplyr::case_when(
            times < fixed_time_horizon & reals == "real_negatives" ~
              "real_censored",
            times > fixed_time_horizon & reals != "real_negatives" ~
              "real_negatives",
            TRUE ~ outcome
          )
        ) |>
        dplyr::mutate(outcome = as.factor(outcome)) |> # format
        dplyr::count(
          fixed_time_horizon,
          strata,
          outcome,
          name = "reals_estimate",
          .drop = FALSE
        ) |>
        dplyr::group_by(fixed_time_horizon, strata) |>
        dplyr::mutate(estimate = reals_estimate / sum(reals_estimate)) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          censored_assumption = censored_assumption,
          competing_assumption = competing_assumption
        )
    )
  }

  if (censored_assumption == "excluded" & competing_assumption == "adjusted") {
    return(purrr::map_dfr(
      fixed_time_horizons,
      ~ compute_aj_estimates(
        data_for_tidycmprsk,
        .x,
        "real_negatives",
        "real_censored"
      )
    ))
  }

  if (censored_assumption == "adjusted" & competing_assumption == "excluded") {
    return(purrr::map_dfr(
      fixed_time_horizons,
      ~ compute_aj_estimates(
        data_for_tidycmprsk,
        .x,
        "real_competing",
        "real_competing"
      )
    ))
  }

  if (censored_assumption == "adjusted" & competing_assumption == "adjusted") {
    return(purrr::map_dfr(
      fixed_time_horizons,
      ~ compute_aj_estimates(data_for_tidycmprsk, .x, NA, NA)
    ))
  }
}

c("ppcr", "probability_threshold")
create_estimates_with_assumptions <- function(
  data_for_tidycmprsk,
  fixed_time_points,
  assumptions
) {
  create_non_parametric_reals_estimates_by_assumptions(
    data_for_tidycmprsk = data_for_tidycmprsk,
    fixed_time_horizons = fixed_time_points,
    competing_assumption = assumptions$competing,
    censored_assumption = assumptions$censored
  ) |>
    dplyr::left_join(
      data_for_tidycmprsk |>
        select(strata, stratified_by) |>
        distinct()
    )
}

add_cutoff_strata <- function( data, by ) {

  data |>
    dplyr::mutate(
      strata_probability_threshold = cut(
        probs,
        breaks = create_breaks_values(
          probs, stratified_by = "probability_threshold", by = by),
        include.lowest = TRUE
      ),
      strata_ppcr = dplyr::ntile( dplyr::desc(probs), 1 / by),
      strata_ppcr = factor((strata_ppcr) / (1 / by))
    )

}

create_aj_data <- function(probs, reals, times, fixed_time_horizons, stratified_by, by) {

  assumption_sets <- list(
    list(competing = "excluded", censored = "excluded"),
    list(competing = "adjusted", censored = "adjusted"),
    list(competing = "excluded", censored = "adjusted"),
    list(competing = "adjusted", censored = "excluded")
  )

  probs_reals_times_data <- tibble::tibble(
    probs = probs,
    reals = reals,
    times = times
  ) |>
    dplyr::mutate(
      reals = factor(
        dplyr::recode(
          reals,
          `0` = "real_negatives",
          `2` = "real_competing",
          `1` = "real_positives"
        ),
        levels = c("real_negatives", "real_competing", "real_positives") # Setting reference level
      )
    ) |>
    add_cutoff_strata(by = by) |>  # TODO: add stratified by
    tidyr::pivot_longer(
      cols = tidyr::starts_with("strata"),
      names_to = "stratified_by",
      values_to = "strata",
      names_prefix = "strata_")

  strata_combinations <- stratified_by |>
    purrr::map_df(
      ~create_strata_combinations(.x, 0.01)
    )

  aj_dat <- assumption_sets |>
    purrr::map_dfr(
      ~ create_estimates_with_assumptions(
        data_for_tidycmprsk = probs_reals_times_data,
        fixed_time_points = fixed_time_horizons,
        assumptions = .x
      )
    ) |>
    dplyr::arrange(fixed_time_horizon, outcome)

  aj_dat

}

# assumption_sets <- list(
#   list(competing = "excluded", censored = "excluded"),
#   list(competing = "adjusted", censored = "adjusted"),
#   list(competing = "excluded", censored = "adjusted"),
#   list(competing = "adjusted", censored = "excluded")
# )
#
#
# aj_dat <- assumption_sets |>
#   purrr::map_dfr(
#     ~ create_estimates_with_assumptions(
#       data_for_tidycmprsk = data_for_cmprsk,
#       fixed_time_points = fixed_time_horizons,
#       assumptions = .x
#     )
#   ) |>
#   dplyr::arrange(fixed_time_horizon, outcome)
#
# library(dplyr)
#
# aj_dat |>
#   count(
#     outcome
#   )
#
#
# aj_dat |>
#   count(
#     censored_assumption,
#     competing_assumption,
#     outcome
#   ) |>
#   View()
#
#
# aj_dat |>
#   group_by(
#     fixed_time_horizon,
#     censored_assumption,
#     competing_assumption
#   ) |>
#   summarise(sum(reals_estimate)) |>
#   View()
#
#
# library(plotly)
#
# aj_dat$censored_assumption
#
# adjadj <- aj_dat |>
#   dplyr::filter(
#     censored_assumption == "adjusted",
#     competing_assumption == "adjusted"
#   ) |>
#   plot_ly(
#     x =~ strata,
#     y =~ reals_estimate,
#     color =~ outcome,
#     frame =~ fixed_time_horizon
#   )
#
# adjexc <- aj_dat |>
#   dplyr::filter(
#     censored_assumption == "adjusted",
#     competing_assumption == "excluded"
#   ) |>
#   plot_ly(
#     x =~ strata,
#     y =~ reals_estimate,
#     color =~ outcome,
#     frame =~ fixed_time_horizon
#   )
#
# subplot(
#   adjadj,
#   adjexc,shareY = TRUE
# )



