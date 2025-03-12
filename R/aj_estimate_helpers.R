
create_aj_data_combinations <- function(
    probs, fixed_time_horizons, stratified_by, by) {

  strata_combinations <- stratified_by |>
    purrr::map_df(
      ~create_strata_combinations(.x, by)
    )

  assumption_sets <- list(
    list(competing = "excluded", censored = "excluded"),
    list(competing = "adjusted", censored = "adjusted"),
    list(competing = "excluded", censored = "adjusted"),
    list(competing = "adjusted", censored = "excluded")
  )

  reals <- factor(
    x = paste("real", c("negatives", "positives", "competing", "censored"), sep = "_"),
    levels = paste("real", c("negatives", "positives", "competing", "censored"), sep = "_"))

  tidyr::expand_grid(
   reference_group = names(probs),
   fixed_time_horizon = fixed_time_horizons,
   reals = reals,
   censoring_assumption = c("excluded", "adjusted"),
   competing_assumption = c(
     "excluded",
     "adjusted_as_negative",
     "adjusted_as_censored"),
   strata_combinations
  )

}

#' Create strata combinations
#'
#' @description
#' Returns strata combinations that are necessary for estimation
#'
#' @param stratified_by A single string, one of `"probability_threshold"` or `"ppcr"`.
#' @param by A numeric value.
#'
#' @returns
#' A tibble of strata combinations.
#'
#' @export
create_strata_combinations <- function( stratified_by, by) {

  if ( stratified_by == "probability_threshold" ) {

    strata_combinations <- tibble::tibble(
      upper_bound = create_breaks_values(NA, "probability_threshold", by),
      lower_bound = dplyr::lag(upper_bound, default = 0),
      mid_point = upper_bound - by/2,
      include_lower_bound = (lower_bound == 0 & upper_bound != 0),
      include_upper_bound = upper_bound != 0,
      strata = glue::glue(
        "{ifelse(include_lower_bound==TRUE,'[','(')}{lower_bound},\\
        {upper_bound}{ifelse(include_upper_bound==TRUE,']',')')}"
      ),
      chosen_cutoff = upper_bound
    )

  } else if ( stratified_by == "ppcr" ) {

    strata_combinations <- tibble::tibble(
      strata = create_breaks_values(NA, "probability_threshold", by)[-1],
      lower_bound = strata - by,
      upper_bound = strata + by,
      mid_point = upper_bound - by/2,
      include_lower_bound = TRUE,
      include_upper_bound = FALSE,
      chosen_cutoff = strata
    )

  }

  strata_combinations |>
    dplyr::mutate(
      strata = as.factor(strata),
      stratified_by )

}

#' Update administrative censoring
#'
#' @description
#' This function update outcomes according to fixed time horizons:
#' If no event happened until fixed time horizon the outcome will be updated
#' to `real_negatives`. If the original outcome is already defined as `real_negatives`
#' but the tiemes is updated prior to the fixed time horizon, the outcome will be
#' updated as `real_censored`.
#'
#' @inheritParams extract_crude_estimate
#'
#' @returns
#' A modified version of `data_to_adjust`.
#'
#' @export

update_administrative_censoring <- function( data_to_adjust ) {

  data_to_adjust |>
    dplyr::mutate(
      reals = dplyr::case_when(
        times > fixed_time_horizon & reals == "real_positives" ~ "real_negatives",
        times < fixed_time_horizon & reals == "real_negatives" ~ "real_censored",
        TRUE ~ reals
      )
    )

}


extract_aj_estimate <- function(data_to_adjust, fixed_time_horizons) {

  ajfit <- survival::survfit(
    survival::Surv(
      times,
      reals
    ) ~
      strata,
    data = data_to_adjust,
    istate = rep("real_negatives", length(data_to_adjust$times))
  )

  summary_aj_fit <- summary(
    ajfit,
    times = fixed_time_horizons,
    extend = TRUE
  )

  counts_by_strata <- tibble::tibble(
    strata = levels(summary_aj_fit$strata),
    n = summary_aj_fit$n
  ) |>
    dplyr::mutate(
      strata = gsub("strata=", "", strata) # TODO: define in survift prefix
    )

  print("summary_aj_fit$time")
  print(summary_aj_fit$time)

  summary_aj_fit$pstate |>
    as.data.frame() |>
    purrr::set_names(ajfit$states) |>
    tibble::add_column(fixed_time_horizon = summary_aj_fit$time) |>
    tibble::add_column(strata = gsub("strata=", "", summary_aj_fit$strata)) |>
    # tibble::add_column(
    #   total_included = rep(summary_aj_fit$n, each = length(fixed_time_horizons))
    #   ) |>
    tidyr::pivot_longer(
      cols = !c(strata, fixed_time_horizon),
      values_to = "estimate",
      names_to = "reals"
    ) |>
    dplyr::left_join(
      counts_by_strata
    ) |>
    dplyr::mutate(
      reals_estimate = estimate * n
    ) |>
    dplyr::mutate( strata = as.factor(strata) )

}


#' Extract crude estimate
#'
#' @description
#' This function extract the crude estimates of the
#' state occupancy in the end of the followup
#'
#' @param data_to_adjust A data frame contains original estimated risks,
#' estimated risks stratas, the time to events, and the stratification technique
#'
#' @returns
#' A data frame contains the crude estimates, i.e the counts for each state.
#'
#' @export
extract_crude_estimate <- function(data_to_adjust) {

  data_to_adjust |>
    dplyr::group_by(strata, reals, fixed_time_horizon) |>
    dplyr::summarise(
      reals_estimate = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    tidyr::complete(strata, reals, fixed_time_horizon,
             fill = list(reals_estimate = 0))

}

#' Extract adjusted estimate by assumptions
#'
#' @description
#' A short description...
#'
#' @inheritParams extract_crude_estimate
#' @param fixed_time_horizons Fixed time horizons.
#' @param censoring_assumption One of `"excluded"` or `"adjusted"`.
#' @param competing_assumption One of `"excluded"` or `"adjusted"`.
#'
#' @returns
#' An adjusted estimate.
#'
#' @export
extract_aj_estimate_by_assumptions <- function(
    data_to_adjust,
    fixed_time_horizons,
    censoring_assumption = "excluded",
    competing_assumption = "excluded") {

  if ( censoring_assumption == "excluded" &
    competing_assumption == "excluded" ) {

    aj_estimate_data <- data_to_adjust |>
      tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
      update_administrative_censoring() |>
      extract_crude_estimate()

  } else if ( censoring_assumption == "excluded" &
              competing_assumption == "adjusted_as_negative") {

    aj_estimate_data_excluded <- data_to_adjust |>
      tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
      update_administrative_censoring() |>
      dplyr::filter( reals == "real_censored" ) |>
      extract_crude_estimate()

    aj_estimate_data_adjusted <- purrr::map(
      fixed_time_horizons,
      ~data_to_adjust |>
        tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
        dplyr::filter( reals != "real_censored" ) |>
        dplyr::filter( fixed_time_horizon == .x) |>
        extract_aj_estimate(fixed_time_horizons = .x)
    ) |>
      purrr::set_names(fixed_time_horizons) |>
      purrr::list_rbind(
        names_to = "fixed_time_horizon") |>
      dplyr::mutate( fixed_time_horizon = as.numeric(fixed_time_horizon) )

    aj_estimate_data <- aj_estimate_data_excluded |>
      dplyr::bind_rows(aj_estimate_data_adjusted)



  } else if ( censoring_assumption == "adjusted" &
              competing_assumption == "excluded" ) {

    aj_estimate_data_excluded <- data_to_adjust |>
      tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
      update_administrative_censoring() |>
      dplyr::filter( reals == "real_competing" ) |>
      extract_crude_estimate()

    aj_estimate_data_adjusted <- purrr::map(
      fixed_time_horizons,
      ~data_to_adjust |>
        tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
        dplyr::filter( reals != "real_competing" ) |>
        dplyr::filter( fixed_time_horizon == .x) |>
        extract_aj_estimate(fixed_time_horizons = .x) |>
        dplyr::filter( reals != "real_competing" )
    ) |>
      purrr::set_names(fixed_time_horizons) |>
      purrr::list_rbind(
        names_to = "fixed_time_horizon") |>
      dplyr::mutate( fixed_time_horizon = as.numeric(fixed_time_horizon) )

    # TODO: add rows for reals_censoring = 0 for each strata

    aj_estimate_data <- aj_estimate_data_excluded |>
      dplyr::bind_rows(aj_estimate_data_adjusted)

  } else if (censoring_assumption == "adjusted" &
             competing_assumption == "adjusted_as_negative") {


    aj_estimate_data <- extract_aj_estimate(
      data_to_adjust,
      fixed_time_horizons = fixed_time_horizons
    )

  }

  aj_estimate_data |>
    dplyr::mutate(
      censoring_assumption = censoring_assumption,
      competing_assumption = competing_assumption
    )

}
