# Methods --------

#' Predict method for objects of class `peathamstr_fit`
#'
#' Either returns
#' posterior ages at specified depths, posterior age accumulation rates
#' (yr cm\eqn{^{-1}}) at the depths the model was constructed on, or posterior
#' mass fluxes (kg yr\eqn{^{-1}} m\eqn{^{-2}}) at the depths the model was
#' constructed on.
#'
#' @param object An object of class `peathamstr_fit`.
#'
#' @param type A character value defining what posterior variables are predicted.
#' One of:
#' \describe{
#'   \item{`"age_models"`}{Predicts posterior ages at specified depths.}
#'   \item{`"acc_rates"`}{Predicts posterior age accumulation rates
#'      (yr cm\eqn{^{-1}}) at the depths the model was constructed on}
#'   \item{`"mass_fluxes"`}{Predicts posterior mass fluxes (kg yr\eqn{^{-1}}
#'      m\eqn{^{-2}}) at the depths the model was constructed on.}
#' }
#'
#' @param depth A character value defining depth values for which predictions are
#' computed. Must be one of:
#' \describe{
#'   \item{`"modelled"`}{Predicts for depths the model was constructed on.}
#'   \item{`"data"`}{Predicts for depths of the input data to the model.}
#' }
#' Alternatively can be a numeric vector to specify depths (this only works if
#' `type == "age_models"`).
#'
#' @param ... additional arguments to `hamstr` predict methods
#'
#' @return A `tibble` of `peathamstr_fit` age-depth model realizations.
#'
#' @export
predict.peathamstr_fit <- function(
    object,
    type = c("age_models", "acc_rates", "mass_fluxes", "cumulative_mass", "apparent_mass_accumulation_rates"),
    depth = c("modelled", "data"),
    ...)
{

  stopifnot(type %in% c("age_models", "acc_rates", "mass_fluxes", "cumulative_mass", "apparent_mass_accumulation_rates"))

  if(type %in% c("age_models", "acc_rates")) {
    NextMethod()
  } else {

    res <-
      switch(
        type,
        "mass_fluxes" = {
          stopifnot(depth %in% c("modelled"))

          res <-
            get_posterior_mass_fluxes(object)

          res_age <-
            object |>
            structure(class = setdiff(class(object), "peathamstr_fit")) |>
            predict(type = "age_models", depth = res$depth[res$iter == 1])

          res <-
            dplyr::left_join(
              res,
              res_age,
              by = c("iter", "depth")
            )

          res

        },
        "cumulative_mass" = {

          res <-
            get_posterior_cumulative_mass(
              object = object,
              depth = depth
            )

          res_age <-
            object |>
            structure(class = setdiff(class(object), "peathamstr_fit")) |>
            predict(type = "age_models", depth = res$depth[res$iter == 1])

          dplyr::left_join(
            res,
            res_age,
            by = c("iter", "depth")
          )

        },
        "apparent_mass_accumulation_rates" = {

          if(depth != "data") {
            stop('Making predictions for `type = "apparent_mass_accumulation_rates"` is currently only possible with `depth = "data"`.')
          }

          # first: extract posterior ages
          res_age_upper <-
            object |>
            structure(class = setdiff(class(object), "peathamstr_fit")) |>
            predict(type = "age_models", depth = object$data$depth2_upper)

          res_age_lower <-
            object |>
            structure(class = setdiff(class(object), "peathamstr_fit")) |>
            predict(type = "age_models", depth = object$data$depth2_lower)

          res_age <-
            res_age_upper |>
            dplyr::rename(
              age_upper = "age",
              depth_upper = "depth"
            ) |>
            dplyr::mutate(
              age_lower = res_age_lower$age,
              depth_lower = res_age_lower$depth,
              duration = age_lower - age_upper
            ) |>
            dplyr::left_join(
              tibble::tibble(
                depth_upper = object$data$depth2_upper,
                mass = {
                  .x <- c(object$data$cumulative_mass0, object$data$layer_mass)
                  .x[-1] - .x[-length(.x)]
                }
              ),
              by = "depth_upper"
            ) |>
            dplyr::mutate(
              amar = mass / duration
            )

          res_age |>
            dplyr::select(dplyr::all_of(c("iter", "depth_upper", "depth_lower", "age_upper", "age_lower", "duration", "mass", "amar")))

        })

    return(res)

  }

}


#' Summarizes objects of class `peathamstr_fit`
#'
#' @inheritParams predict.peathamstr_fit
#'
#' @param ... further parameters passed to individual summary functions.
#'
#' @return A tibble with summarized posterior draws.
#'
#' @export
summary.peathamstr_fit <- function(
    object,
    type = c("age_models", "acc_rates", "pars", "mass_fluxes", "cumulative_mass"),
    #probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975),
    #tau = 0, kern = c("U", "G", "BH"),
    ...
) {

  if(type %in% c("age_models", "acc_rates", "pars")) {
    NextMethod()
  } else {
    switch(
      type,
      "mass_fluxes" = summarise_hamstr_mass_fluxes(object, ...)
    )
  }

}



#' Computes posterior mass fluxes
#'
#' @keywords internal
#' @noRd
get_posterior_mass_fluxes <- function(object) {

  # add ids
  depths <-
    tibble::as_tibble(
      object$data[c("depth2_upper", "depth2_lower")]
    ) |>
    dplyr::mutate(
      depth = .data$depth2_upper,
      idx = seq_along(depth)
    )

  dplyr::bind_cols(
    as.data.frame(object$fit, pars = "clymo_par") |>
      tibble::as_tibble() |>
      dplyr::mutate(
        iter = seq_len(dplyr::n())
      ) |>
      tidyr::pivot_longer(
        cols = -tidyr::contains("iter"),
        names_to = "par",
        values_to = "nmu"
      ) |>
      dplyr::mutate(
        nmu = nmu * object$data$p3_clymo_par
      ),
    as.data.frame(object$fit, pars = "nmr_rep") |>
      tibble::as_tibble() |>
      dplyr::mutate(
        iter = seq_len(dplyr::n())
      ) |>
      tidyr::pivot_longer(
        cols = -tidyr::contains("iter"),
        names_to = "par",
        values_to = "nmr"
      ) |>
      dplyr::select(nmr)
  ) |>
    dplyr::mutate(
      nmb = nmu + nmr, #---note: net mass balance ---todo: check sign here: NCB should be smaller than NCU, hence in the formula in @Yu.2011, we have to assume NCR to have positive values
      idx = get_par_idx(.data$par),
      par = "nmu"
    ) |>
    dplyr::right_join(depths, by = "idx") |>
    dplyr::arrange("par", "iter", "idx", "depth") |>
    dplyr::select("iter", "idx", "depth", "depth2_upper", "depth2_lower", "nmu", "nmr", "nmb")

}


#' Computes posterior cumulative masses
#'
#' @keywords internal
#' @noRd
get_posterior_cumulative_mass <- function(object, depth = c("modelled", "data")) {

  switch(
    depth,
    "modelled" = {

      # add ids
      depths <-
        tibble::as_tibble(
          object$data[c("depth2_upper", "depth2_lower")]
        ) |>
        dplyr::mutate(
          depth = .data$depth2_upper,
          idx = seq_along(depth)
        )

      as.data.frame(object$fit, pars = "layer_cumulative_mass_rep") |>
        tibble::as_tibble() |>
        dplyr::mutate(
          iter = seq_len(dplyr::n())
        ) |>
        tidyr::pivot_longer(
          cols = -tidyr::contains("iter"),
          names_to = "par",
          values_to = "cumulative_mass"
        ) |>
        dplyr::mutate(
          idx = get_par_idx(.data$par),
          par = "cumulative_mass"
        ) |>
        dplyr::right_join(depths, by = "idx") |>
        dplyr::arrange("par", "iter", "idx", "depth") |>
        dplyr::select("iter", "idx", "depth", "depth2_upper", "depth2_lower", "cumulative_mass")

    },
    "data" = {
      tibble::tibble(
        depth = object$data$depth2_upper,
        cumulative_mass =
          as.data.frame(object$fit, pars = "Mod_layer_mass_cumulative") |>
          purrr::map(
            function(.x) tibble::tibble(iter = seq_along(.x), cumulative_mass = .x)
          )
      ) |>
        tidyr::unnest(cumulative_mass)
    })

}


summarise_hamstr_mass_fluxes <- function(
    object,
    probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)
) {

  res <-
    stats::predict(object, type = "mass_fluxes", depth = "modelled") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("nmu", "nmr", "nmb")),
      names_to = "variable"
    ) |>
    dplyr::group_by(.data$variable, .data$depth, .data$depth2_upper, .data$depth2_lower, .data$idx)

  dplyr::bind_cols(
    res |>
      summarise_q(var = .data$value, probs = probs) |>
      dplyr::ungroup(),
    res |>
      dplyr::summarise(
        age = mean(age, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::select(age)
  ) |>
    dplyr::arrange(.data$variable, .data$depth)

}
