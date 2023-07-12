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
#'   \item{`"carbon_fluxes"`}{Predicts posterior carbon mass fluxes
#'     (kg yr\eqn{^{-1}} m\eqn{^{-2}}) at the depths the model was constructed
#'     on. For this, argument `carbon_content` needs to be specified.}
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
#' @param carbon_content A numeric vector with the same length as there are
#' `depth2` values in the input data. This only needs to be specified when
#' `type == "carbon_fluxes"` and otherwise can be set to `NULL`. This will be
#' used to interpolate carbon contents for layers the model was constructed on.
#' The interpolated carbon contents are multiplied with the mass fluxes (as
#' predicted when setting `type` to `"mass_fluxes"`) to compute respective
#' carbon fluxes.
#'
#' @param ... additional arguments to `hamstr` predict methods
#'
#' @return A `tibble` of `peathamstr_fit` age-depth model realizations.
#'
#' @export
predict.peathamstr_fit <- function(
    object,
    type = c("age_models", "acc_rates", "mass_fluxes", "carbon_fluxes", "cumulative_mass", "cumulative_carbon_mass", "apparent_mass_accumulation_rates", "apparent_carbon_accumulation_rates"),
    depth = c("modelled", "data"),
    carbon_content = NULL,
    ...)
{

  stopifnot(type %in% c("age_models", "acc_rates", "mass_fluxes", "carbon_fluxes", "cumulative_mass", "cumulative_carbon_mass", "apparent_mass_accumulation_rates", "apparent_carbon_accumulation_rates"))
  if(type %in% c("carbon_fluxes", "cumulative_carbon_mass", "apparent_carbon_accumulation_rates")) {
    stopifnot(!is.null(carbon_content) & length(carbon_content) == object$data$N2)
  }

  if(type %in% c("age_models", "acc_rates")) {
    NextMethod()
  } else {

    res <-
    switch(
      type,
      "mass_fluxes" = ,
      "carbon_fluxes" = {
        stopifnot(depth %in% c("modelled"))

        # first: extract posterior ages
        res_age <-
          object |>
          structure(class = setdiff(class(object), "peathamstr_fit")) |>
          predict(type = "age_models", depth = depth)

        # second: add mass fluxes
        res_mass_fluxes <-
          get_posterior_mass_fluxes(object, res_age)

        res <-
          dplyr::left_join(
            res_mass_fluxes,
            res_age,
            by = c("iter", "depth")
          )

        if(type == "carbon_fluxes") {

          d_increments <-
            res_mass_fluxes |>
            dplyr::select(dplyr::all_of(c("c_depth_top", "c_depth_bottom"))) |>
            dplyr::filter(! duplicated(.data$c_depth_top)) |>
            dplyr::mutate(
              thickness = .data$c_depth_bottom - .data$c_depth_top
            )

          M <-
            do.call(
              "cbind",
              purrr::map(seq_along(object$data$depth2_upper), function(i) {
                translate_depths(
                  layers_from =
                    tibble::tibble(
                      depth_upper = object$data$depth2_upper,
                      depth_lower = object$data$depth2_lower,
                      thickness = .data$depth_lower - .data$depth_upper
                    ) |>
                    dplyr::slice(i) |>
                    as.matrix(),
                  layers_to = as.matrix(d_increments)
                )
              })
            )

          d_increments$C <- as.numeric(M %*% carbon_content/rowSums(M))

          res <-
            dplyr::left_join(
              res,
              d_increments |>
                dplyr::select(! dplyr::all_of(c("thickness", "c_depth_bottom"))),
              by = "c_depth_top"
            ) |>
            dplyr::mutate(
              dplyr::across(dplyr::all_of(c("nmu", "nmr", "nmb")), function(.x) .x * .data$C)
            ) |>
            dplyr::rename(
              ncu = "nmu",
              ncr = "nmr",
              ncb = "nmb"
            )

        }

        res

      },
      "cumulative_mass" = {

        # first: extract posterior ages
        res_age <-
          object |>
          structure(class = setdiff(class(object), "peathamstr_fit")) |>
          predict(type = "age_models", depth = depth)

        res <-
          get_posterior_cumulative_mass(
            object = object,
            depth = depth
          )

        dplyr::left_join(
          res,
          res_age,
          by = c("iter", "depth")
        )

      },
      "cumulative_carbon_mass" = {

        # first: extract posterior ages
        res_age <-
          object |>
          structure(class = setdiff(class(object), "peathamstr_fit")) |>
          predict(type = "age_models", depth = depth)

        res <-
          get_posterior_cumulative_carbon_mass(
            object = object,
            depth = depth,
            carbon_content = carbon_content
          )

        dplyr::left_join(
          res,
          res_age,
          by = c("iter", "depth")
        )

      },
      "apparent_mass_accumulation_rates" =,
      "apparent_carbon_accumulation_rates" = {

        if(depth != "data") {
          stop('Making predictions for `type = "apparent_carbon_accumulation_rates"` or `type = "apparent_mass_accumulation_rates"` is currently only possible with `depth = "data"`.')
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
                .x <- c(object$data$cumulative_mass0, object$data$cumulative_mass)
                .x[-1] - .x[-length(.x)]
              },
              C =
                if(type == "apparent_carbon_accumulation_rates") {
                  carbon_content
                } else {
                  NA_real_
                }
            ),
            by = "depth_upper"
          ) |>
          dplyr::mutate(
            amar = mass / duration,
            acar = mass * C / duration
          )

        if(type == "apparent_carbon_accumulation_rates") {
          res_age |>
            dplyr::select(dplyr::all_of(c("iter", "depth_upper", "depth_lower", "age_upper", "age_lower", "duration", "mass", "C", "acar")))
        } else {
          res_age |>
            dplyr::select(dplyr::all_of(c("iter", "depth_upper", "depth_lower", "age_upper", "age_lower", "duration", "mass", "amar")))
        }

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
    type = c("age_models", "acc_rates", "pars", "mass_fluxes", "carbon_fluxes", "cumulative_mass"),
    #probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975),
    #tau = 0, kern = c("U", "G", "BH"),
    ...
) {

  if(type %in% c("age_models", "acc_rates", "pars")) {
    NextMethod()
  } else {
    switch(
      type,
      "mass_fluxes" = summarise_hamstr_mass_fluxes(object, ...),
      "carbon_fluxes" = summarise_hamstr_carbon_fluxes(object, ...)
    )
  }

}



#' Computes posterior mass fluxes
#'
#' @keywords internal
#' @noRd
get_posterior_mass_fluxes <- function(object, res_age) {

  # add ids
  depths <-
    tibble::as_tibble(
      object$data[c("c", "c_depth_top", "c_depth_bottom")]
    ) |>
    dplyr::mutate(depth = .data$c_depth_top) |>
    dplyr::rename(idx = "c")

  dplyr::bind_cols(
    as.data.frame(object$fit, pars = "nmu_rep") |>
      tibble::as_tibble() |>
      dplyr::mutate(
        iter = seq_len(dplyr::n())
      ) |>
      tidyr::pivot_longer(
        cols = -tidyr::contains("iter"),
        names_to = "par",
        values_to = "nmu"
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
    dplyr::select("iter", "idx", "depth", "c_depth_top", "c_depth_bottom", "nmu", "nmr", "nmb")

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
          object$data[c("c", "c_depth_top", "c_depth_bottom")]
        ) |>
        dplyr::mutate(depth = .data$c_depth_top) |>
        dplyr::rename(idx = "c")

      as.data.frame(object$fit, pars = "c_cumulative_mass") |>
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
        dplyr::select("iter", "idx", "depth", "c_depth_top", "c_depth_bottom", "cumulative_mass")

    },
    "data" = {
      tibble::tibble(
        depth = object$data$depth2,
        cumulative_mass =
          as.data.frame(object$fit, pars = "Mod_cumulative_mass") |>
          purrr::map(
            function(.x) tibble::tibble(iter = seq_along(.x), cumulative_mass = .x)
          )
      ) |>
        tidyr::unnest(cumulative_mass)
    })

}

#' Computes posterior cumulative carbon masses
#'
#' Note: This will start the cumulative carbon mass at 0 for the first modeled
#' layer because the carbon content of previous layers is unknown.
#'
#' @keywords internal
#' @noRd
get_posterior_cumulative_carbon_mass <- function(object, depth = c("modelled", "data"), carbon_content) {

  warning("* Cumulative carbon contents start at 0 for the first measured layer.")
  depth_difference <- (object$data$depth2_lower - object$data$depth2_upper)[-length(object$data$depth2_upper)] - (object$data$depth2_lower[-1] - object$data$depth2_lower[-length(object$data$depth2_upper)])
  if(! isTRUE(all.equal(depth_difference, rep(0.0, length(depth_difference))))) {
    warning("* Measured layers do not cover completely modeled layers (are not contiguous). This results in `NA` carbon contents. `NA` carbon contents are filled in order to compute cumulative carbon masses.")
  }

  switch(
    depth,
    "modelled" = {

      # add ids
      depths <-
        tibble::as_tibble(
          object$data[c("c", "c_depth_top", "c_depth_bottom")]
        ) |>
        dplyr::mutate(depth = .data$c_depth_top) |>
        dplyr::rename(idx = "c") |>
        dplyr::mutate(
          thickness = .data$c_depth_bottom - .data$c_depth_top
        )

      # assign carbon contents to increments
      M <-
        do.call(
          "cbind",
          purrr::map(seq_along(object$data$depth2_upper), function(i) {
            translate_depths(
              layers_from =
                tibble::tibble(
                  depth_upper = object$data$depth2_upper,
                  depth_lower = object$data$depth2_lower,
                  thickness = .data$depth_lower - .data$depth_upper
                ) |>
                dplyr::slice(i) |>
                as.matrix(),
              layers_to = as.matrix(depths)
            )
          })
        )

      depths <-
        depths |>
        dplyr::mutate(
          C = as.numeric(M %*% carbon_content/rowSums(M)),
          C_filled = .data$C
        ) |>
        tidyr::fill("C_filled", .direction = "updown")

      # compute cumulative C mass
      res <-
        as.data.frame(object$fit, pars = "c_cumulative_mass")

      res <-
        res[, -1] - res[-ncol(res)]

      res <- res * matrix(depths$C_filled, nrow = nrow(res), ncol = ncol(res), byrow = TRUE)

      as.data.frame(t(apply(res, 1, cumsum))) |>
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
        dplyr::right_join(depths |> dplyr::select(-.data$thickness), by = "idx") |>
        dplyr::arrange("par", "iter", "idx", "depth") |>
        dplyr::select("iter", "idx", "depth", "c_depth_top", "c_depth_bottom", "cumulative_mass") |>
        dplyr::filter(c_depth_top >= min(object$data$depth2_upper)) |>
        dplyr::mutate(
          cumulative_mass = .data$cumulative_mass - .data$cumulative_mass[depth == min(depth)[[1]]]
        )

    },
    "data" = {

      res <-
        as.data.frame(object$fit, pars = "Mod_cumulative_mass")

      res <-
        res[, -1] - res[-ncol(res)]

      res <- res * matrix(carbon_content[-1], nrow = nrow(res), ncol = ncol(res), byrow = TRUE)

      tibble::tibble(
        depth = object$data$depth2[-1],
        cumulative_mass =
          as.data.frame(t(apply(res, 1, cumsum))) |>
          purrr::map(
            function(.x) tibble::tibble(iter = seq_along(.x), cumulative_mass = .x)
          )
      ) |>
        tidyr::unnest(cumulative_mass)
    })

}


#' Summarizes carbon fluxes from a `peathamstr_fit` object
#'
#' @inheritParams predict.peathamstr_fit
#'
#' @param probs A numeric vector with values in \\[0, 1\\] representing
#' quantiles for which to summarize the posterior distrbutions.
#'
#' @return A tibble containing the summarized posterior carbon fluxes.
#'
#' @noRd
#' @keywords internal
summarise_hamstr_carbon_fluxes <- function(
    object,
    carbon_content,
    probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)
) {

  res <-
    stats::predict(object, type = "carbon_fluxes", carbon_content = carbon_content, depth = "modelled") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("ncu", "ncr", "ncb")),
      names_to = "variable"
    ) |>
    dplyr::group_by(.data$variable, .data$depth, .data$c_depth_top, .data$c_depth_bottom, .data$idx)

  dplyr::bind_cols(
    res |>
      summarise_q(var = .data$value, probs = probs) |>
      dplyr::ungroup(),
    res |>
      dplyr::summarise(
        age = mean(age, na.rm = TRUE),
        C = mean(C, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::select(age, C)
  ) |>
    dplyr::arrange(.data$variable, .data$depth)

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
    dplyr::group_by(.data$variable, .data$depth, .data$c_depth_top, .data$c_depth_bottom, .data$idx)

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
