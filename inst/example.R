# First Try

library(survival)
library(dcurves)
library(gtsummary); library(dplyr); library(tidyr)


cox_model <- coxph(
  Surv(ttcancer, cancer) ~
    age + famhistory + marker, data = df_surv)

fixed_time_horizons <- c(0.2, 0.5, 1)


breaks_pt <- create_breaks_values(
  data_for_cmprsk$probs,
  stratified_by = "probability_threshold",
  by = 0.2
)

breaks_ppcr <- create_breaks_values(
  data_for_cmprsk$probs,
  stratified_by = "ppcr",
  by = 0.1
)

data_for_cmprsk |>
  add_cutoff_strata() |>


data_for_cmprsk <- df_surv |>
  dplyr::mutate(
    reals = dplyr::case_when(
      cancer_cr == "censor" ~ 0,
      cancer_cr == "diagnosed with cancer" ~ 1,
      cancer_cr == "dead other causes" ~ 2
    )
  ) |>
  dplyr::rename("times" = "ttcancer") |>
  dplyr::select(
    reals,
    times
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
  dplyr::bind_cols(
    probs = pred_1_5
  )

data_for_cmprsk <- data_for_cmprsk |>
  add_cutoff_strata()|>
  dplyr::mutate(strata = strata_pt )


assumption_sets <- list(
  list(competing = "excluded", censored = "excluded"),
  list(competing = "adjusted", censored = "adjusted"),
  list(competing = "excluded", censored = "adjusted"),
  list(competing = "adjusted", censored = "excluded")
)

aj_dat <- assumption_sets |>
  purrr::map_dfr(
    ~ create_estimates_with_assumptions(
      data_for_tidycmprsk = data_for_cmprsk,
      fixed_time_points = fixed_time_horizons,
      assumptions = .x,
      stratified_by = "strata"
    )
  ) |>
  dplyr::arrange(fixed_time_horizon, outcome)

# 2nd try



library(rticompeting)
library(survival)
library(dcurves)
library(gtsummary); library(dplyr); library(tidyr)
library(reactable)
cox_model <- coxph(
  Surv(ttcancer, cancer) ~
    age + famhistory + marker, data = df_surv)

thin_model <- coxph(
  Surv(ttcancer, cancer) ~
    age + marker, data = df_surv)

df_surv <- df_surv |>
  dplyr::mutate(
    reals = dplyr::case_when(
      cancer_cr == "censor" ~ 0,
      cancer_cr == "diagnosed with cancer" ~ 1,
      cancer_cr == "dead other causes" ~ 2
    )
  )

pred_1_5 <- broom::augment(
  cox_model,
  newdata = df_surv %>% mutate(ttcancer =  1.5),
  type.predict = "expected"
) |>
  mutate(
    pr_failure18 = 1 - exp(-.fitted)
  ) |>
  pull(pr_failure18)

pred_thin <- broom::augment(
  thin_model,
  newdata = df_surv %>% mutate(ttcancer =  1.5),
  type.predict = "expected"
) |>
  mutate(
    pr_failure18 = 1 - exp(-.fitted)
  ) |>
  pull(pr_failure18)




fixed_time_horizons <- c(1, 3, 5)

aj_dat <- create_aj_data(
  probs = pred_1_5,
  reals = df_surv$reals,
  times = df_surv$ttcancer,
  fixed_time_horizons = fixed_time_horizons,
  by = 0.1
)

full_aj_dat <- list(
  "thin" = pred_thin,
  "full" = pred_1_5
) |>
  purrr::map_df(
    ~create_aj_data(
      probs = .x,
      reals = df_surv$reals,
      times = df_surv$ttcancer,
      fixed_time_horizons = fixed_time_horizons,
      by = 0.1
    ), .id = "reference_group" )


rticompeting::create_ojs_summary_report(
  probs = list("Thin" = pred_thin, "Full" = pred_1_5),
  reals = df_surv$reals,
  times = df_surv$ttcancer,
  fixed_time_horizons = fixed_time_horizons,
  by = 0.1,
  stratified_by = c("probability_threshold", "ppcr")
)

# 3rd try

probs_cmprsk <- list(
  "thin" = pred_thin,
  "full" = pred_1_5
)

fixed_time_horizons <- c(1, 3, 5)



# Example usage

# stratified_by <- "probability_threshold"

stratified_by <- c("probability_threshold", "ppcr")

aj_data_combinations <- create_aj_data_combinations(
  probs_cmprsk,
  fixed_time_horizons,
  stratified_by,
  0.1
)


View(aj_data_combinations)



aj_data_combinations <- create_aj_data_combinations(
  list(
    "thin" = pred_thin
  ),
  fixed_time_horizons,
  stratified_by,
  0.1
)

aj_data_combinations |>
  View()

data_to_adjust <- tibble::tibble(
  probs = probs_cmprsk$thin,
  reals = df_surv$reals,
  times = df_surv$ttcancer
) |>
  add_cutoff_strata(by = 0.1) |>  # TODO: add stratified by
  tidyr::pivot_longer(
    cols = tidyr::starts_with("strata"),
    names_to = "stratified_by",
    values_to = "strata",
    names_prefix = "strata_") |>
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
  )

# Censored Excluded, Competing Excluded

levels(aj_data_combinations$reals)


adjusted_data <- data_to_adjust |>
  extract_aj_estimate_by_assumptions(
    censoring_assumption = "excluded",
    competing_assumption = "excluded",
    fixed_time_horizons = fixed_time_horizons
  ) |>
  mutate(
    reals = factor(reals, paste("real", c("negatives", "positives", "competing",
                                          "censored"), sep = "_")))


# TODO: fill (0,0) with reals_estimate = 0

aj_data_combinations |>
  dplyr::filter(
    censoring_assumption == "excluded",
    competing_assumption == "excluded") |>
  dplyr::left_join(
    adjusted_data,
    by = c(
      "strata", "reals", "censoring_assumption", "fixed_time_horizon",
      "competing_assumption"
    )
  ) |>
  View()


# Censored Excluded, Competing Adjusted as Negatives

aj_data_combinations |>
  dplyr::filter(
    censoring_assumption == "excluded",
    competing_assumption == "adjusted_as_negative") |>
  dplyr::left_join(
    data_to_adjust |>
      extract_aj_estimate_by_assumptions(
        censoring_assumption = "excluded",
        competing_assumption = "adjusted_as_negative",
        fixed_time_horizons = fixed_time_horizons
      ),
    by = c(
      "strata", "reals", "censoring_assumption", "fixed_time_horizon",
      "competing_assumption"
    )
  ) |>
  View()





# Censored Adjusted, Competing Excluded TODO: add rows for censored data
# TODO: add assertion by number of combinations

data_to_adjust |>
  extract_aj_estimate_by_assumptions(
    censoring_assumption = "adjusted",
    competing_assumption = "excluded",
    fixed_time_horizons = fixed_time_horizons
  ) |>
  View()

# Censored Adjusted, Competing Adjusted as Negative


adjusted_data <- data_to_adjust |>
  extract_aj_estimate_by_assumptions(
    censoring_assumption = "adjusted",
    competing_assumption = "adjusted_as_negative",
    fixed_time_horizons = fixed_time_horizons
  ) |>
  mutate(
    reals = factor(reals, paste("real", c("negatives", "positives", "competing",
                                          "censored"), sep = "_")))


aj_data_combinations |>
  dplyr::filter(
    censoring_assumption == "adjusted",
    competing_assumption == "adjusted_as_negative") |>
  dplyr::left_join(
    adjusted_data,
    # by = "strata"
    # by = c(
    # "strata", "reals", "censoring_assumption", "fixed_time_horizon",
    # "competing_assumption"
    # )
  ) |>
  View()

# Censored Adjusted, Competing Adjusted as Censored


adjusted_data <- data_to_adjust |>
  extract_aj_estimate_by_assumptions(
    censoring_assumption = "adjusted",
    competing_assumption = "adjusted_as_censored",
    fixed_time_horizons = fixed_time_horizons
  ) |>
  mutate(
    reals = factor(reals, paste("real", c("negatives", "positives", "competing",
                                          "censored"), sep = "_")))


aj_data_combinations |>
  dplyr::filter(
    censoring_assumption == "adjusted",
    competing_assumption == "adjusted_as_censored") |>
  dplyr::left_join(
    adjusted_data,
    # by = "strata"
    # by = c(
    # "strata", "reals", "censoring_assumption", "fixed_time_horizon",
    # "competing_assumption"
    # )
  ) |>
  View()




# TODO: long format and add crude_aj_estimates

# add total_included

data_to_adjust |>
  View()

# TODO: change labels to censored

data_to_adjust |>
  View()

# TODO: should work for all_included = TRUE

# If competing_assumption = 'excluded', censoring_assumption = 'exclucded'

fixed_time_horizons

# DONE: add fixed_time_horizons: expand grid

data_to_adjust |>
  tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
  update_administrative_censoring() |>
  extract_crude_estimate() |>
  View()

# Competing Excluded, Censored Adjusted

data_to_adjust |>
  tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
  update_administrative_censoring() |>
  dplyr::filter( reals == "real_competing" ) |>
  extract_crude_estimate() |>
  View()


# Competing Adjusted, Censored Excluded

# Excluded obs

data_to_adjust |>
  tidyr::expand_grid(fixed_time_horizon = fixed_time_horizons)  |>
  update_administrative_censoring() |>
  dplyr::filter( reals == "real_censored" ) |>
  extract_crude_estimate() |>
  View()

# Adjusted obs

purrr::map(
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
  View()


extract_aj_estimate(
  data_to_adjust,
  c(1, 2, 3)
) |>
  View()




# Competing Adjusted, Censored Adjusted

# return none


?dplyr::if_else

extract_crude_estimate(
  data_to_adjust) |>
  View()
# dplyr::group_by(strata) |>
# dplyr::mutate(
#   total_excluded = sum(reals_estimate),
#   estimate = reals_estimate / total_excluded) |>
# dplyr::ungroup() |>
# View()



# add reals estimate

extract_aj_estimate(
  data_to_adjust, c(1, 2, 3)
) |>
  View()





summary_aj_fit$time

summary_aj_fit <- summary(
  ajfit,
  times = fixed_time_horizons,
  extend = TRUE
)

summary_aj_fit$strata

summary_aj_fit$n



ajfit$n
ajfit$time
ajfit$n.risk
ajfit$states
