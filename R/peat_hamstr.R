#' Estimates an age-depth model and parameters from a modified Clymo Bog Growth Model
#'
#' @inheritParams hamstr::hamstr
#'
#' @param cumulative_mass A numeric vector with the cumulative mass of the peat
#' layers (kg m\eqn{^{-1}}).
#'
#' @param cumulative_mass0 A numeric value with the cumulative mass of all peat
#' layers above the first considered peat layer (kg m\eqn{^{-1}}).
#'
#' @param depth2 A numeric vector with the same length as `cumulative_mass`
#' defining the depth of the corresponding peat layers.
#'
#' @param depth2_upper,depth2_lower A numeric vector with the same length as
#' `cumulative_mass` defining the depth of the upper and lower boundaries of the
#' corresponding peat layers.
#'
#' @param p1_layer_mass_shape,p2_layer_mass_shape A numeric
#' value representing the shape and rate parameter of a Gamma distribution
#' modeling the shape parameter of the Gamma distribution assumed for modeled
#' peat cumulative masses.
#'
#' @param p1_clymo_par,p2_clymo_par A numeric
#' value representing the shape and rate parameter of the Gamma distribution
#' modeling peat mass addition rates (g m\eqn{^{-2}} yr\eqn{^{-1}}).
#'
#' @param p1_clymo_alpha,p2_clymo_alpha A numeric
#' value representing the shape and rate parameter of the Gamma distribution
#' modeling peat decomposition rate (in a simple exponential decomposition
#' model) (yr\eqn{^{-1}}).
#'
#' @export
peat_hamstr <-
  function(depth, obs_age, obs_err,
           min_age = 1950 - as.numeric(format(Sys.Date(), "%Y")),
           K = NULL,
           top_depth = NULL, bottom_depth = NULL,
           acc_mean_prior = NULL,
           acc_shape = 1.5,
           mem_mean = 0.5, mem_strength = 10,
           model_bioturbation = FALSE,
           n_ind = NULL,
           L_prior_mean = 10,
           L_prior_shape = 2,
           L_prior_sigma = NULL,
           model_displacement = FALSE,
           D_prior_scale = 10,
           model_hiatus = FALSE,
           H_top = NULL, H_bottom = NULL,
           layer_mass,
           cumulative_mass0,
           depth_clymo_par_constant,
           depth2_upper,
           depth2_lower,
           p1_layer_mass_shape = 5,
           p2_layer_mass_shape = 0.01,
           p3_layer_mass_shape = 100,
           p1_clymo_par = 3,
           p2_clymo_par = 0.02,
           p3_clymo_par = 100,
           p1_clymo_alpha_1 = 2,
           p2_clymo_alpha_1 = 50,
           p3_clymo_alpha_1 = 1/1000,
           p1_clymo_alpha_2 = 3,
           p2_clymo_alpha_2 = 3/0.005,
           p3_clymo_alpha_2 = 1/10000,
           p1_ac_age = 3,
           p2_ac_age = 3/80,
           p3_ac_age = 40,
           clymo_par_memory_p1 = 1,
           clymo_par_memory_p2 = 1,
           p1_age0 = 5,
           p2_age0 = 5/1,
           p3_age0 = 1,
           sample_posterior = TRUE,
           hamstr_control = list(),
           stan_sampler_args = list()
  ) {

    # # first: get hamstr_fit object
    # res <-
    #   hamstr::hamstr(
    #     depth = depth,
    #     obs_age = obs_age,
    #     obs_err = obs_err,
    #     min_age = min_age,
    #     K = K,
    #     top_depth = top_depth,
    #     bottom_depth = bottom_depth,
    #     acc_mean_prior = acc_mean_prior,
    #     acc_shape = acc_shape,
    #     mem_mean = mem_mean,
    #     mem_strength = mem_strength,
    #     model_bioturbation = model_bioturbation,
    #     n_ind = n_ind,
    #     L_prior_mean = L_prior_mean,
    #     L_prior_shape = L_prior_shape,
    #     L_prior_sigma = L_prior_sigma,
    #     model_displacement = model_displacement,
    #     D_prior_scale = D_prior_scale,
    #     model_hiatus = model_hiatus,
    #     H_top = H_top,
    #     H_bottom = H_bottom,
    #     sample_posterior = FALSE,
    #     hamstr_control = hamstr_control,
    #     stan_sampler_args = stan_sampler_args
    #   )

    # ## next: add additional parts

    # # data
    # res$data$cumulative_mass <- cumulative_mass
    # res$data$cumulative_mass0 <- cumulative_mass0
    # res$data$p1_layer_mass_shape <- p1_layer_mass_shape
    # res$data$p2_layer_mass_shape <- p2_layer_mass_shape
    # res$data$p1_clymo_par <- p1_clymo_par
    # res$data$p2_clymo_par <- p2_clymo_par
    # res$data$p1_clymo_alpha <- p1_clymo_alpha
    # res$data$p1_clymo_alpha <- p2_clymo_alpha
    # res$data$depth2 <- depth2
    # res$data$N2 <- length(res$data$depth2)
    # res$data$which_c2 <-
    #   sapply(
    #     res$data$depth2,
    #     function(d) {
    #       which.max((res$data$c_depth_bottom < d) * (res$data$c_depth_bottom - d))
    #     })

    # used_sampler_args <- do.call(hamstr:::get_stan_sampler_args, stan_sampler_args)


    # # set the seed here and not inside get_inits_hamstr, so that the chains are
    # # different

    # set.seed(used_sampler_args$seed)

    # inits <- replicate(used_sampler_args$chains, #---todo: avoid using internal functions
    #                    list(hamstr:::get_inits_hamstr(res$data)))

    # args <- #---todo
    #   list(
    #     object = stanmodels$peathamstr,
    #     data = res$data,
    #     init = inits
    #   )

    # args <- append(args, used_sampler_args)

    # if(sample_posterior) {
    #   res$fit <- do.call(rstan::sampling, args)
    # }


    ###

    stan_dat <- peat_make_stan_dat_hamstr()
    used_sampler_args <- do.call(hamstr:::get_stan_sampler_args, stan_sampler_args)
    set.seed(used_sampler_args$seed)
    inits <- replicate(used_sampler_args$chains, list(hamstr:::get_inits_hamstr(stan_dat)))
    args <- list(object = stanmodels$mm2, data = stan_dat,
                 init = inits)
    args <- append(args, used_sampler_args)
    if (sample_posterior) {
      fit <- do.call(rstan::sampling, args)
    }
    else if (sample_posterior == FALSE) {
      fit <- NA
    }
    stan_dat <- append(stan_dat, used_sampler_args)
    info <- list(version = utils::packageVersion("hamstr"), time = Sys.time())
    out <- list(fit = fit, data = stan_dat, info = info)
    class(out) <- append(c("peathamstr_fit", "hamstr_fit"), class(out))
    return(out)

  }


peat_make_stan_dat_hamstr <- function (...)
{
  l <- c(as.list(parent.frame()))
  default.args <- formals(peat_hamstr)
  default.args$cumulative_mass <- as.name("layer_mass")
  default.args$cumulative_mass0 <- as.name("cumulative_mass0")
  default.args$p1_layer_mass_shape <- as.name("p1_layer_mass_shape")
  default.args$p2_layer_mass_shape <- as.name("p2_layer_mass_shape")
  default.args$p3_layer_mass_shape <- as.name("p3_layer_mass_shape")
  default.args$p1_clymo_par <- as.name("p1_clymo_par")
  default.args$p2_clymo_par <- as.name("p2_clymo_par")
  default.args$p3_clymo_par <- as.name("p3_clymo_par")
  default.args$p1_clymo_alpha_1 <- as.name("p1_clymo_alpha_1")
  default.args$p2_clymo_alpha_1 <- as.name("p2_clymo_alpha_1")
  default.args$p3_clymo_alpha_1 <- as.name("p3_clymo_alpha_1")
  default.args$p1_clymo_alpha_2 <- as.name("p1_clymo_alpha_2")
  default.args$p2_clymo_alpha_2 <- as.name("p2_clymo_alpha_2")
  default.args$p3_clymo_alpha_2 <- as.name("p3_clymo_alpha_2")
  default.args$p1_ac_age <- as.name("p1_ac_age")
  default.args$p2_ac_age <- as.name("p2_ac_age")
  default.args$p3_ac_age <- as.name("p3_ac_age")
  default.args$depth2_upper <- as.name("depth2_upper")
  default.args$depth2_lower <- as.name("depth2_lower")
  default.args$clymo_par_memory_p1 <- as.name("clymo_par_memory_p1")
  default.args$clymo_par_memory_p2 <- as.name("clymo_par_memory_p2")
  default.args$clymo_par_memory_p1 <- as.name("clymo_par_memory_p1")
  default.args$p1_age0 <- as.name("p1_age0")
  default.args$p2_age0 <- as.name("p2_age0")
  default.args$p3_age0 <- as.name("p3_age0")
  default.arg.nms <- names(default.args)
  l <- l[lapply(l, is.null) == FALSE]
  default.args[names(l)] <- l
  l <- default.args
  hc.default.args <- formals(hamstr:::hamstr_control)
  hc.default.arg.nms <- names(hc.default.args)
  hc <- l$hamstr_control
  hc <- hc[lapply(hc, is.null) == FALSE]
  hc.default.args[names(hc)] <- hc
  hc <- hc.default.args
  l <- append(l, hc)
  l <- l[names(l) != "hamstr.control"]
  if (is.null(l$acc_mean_prior)) {
    d <- data.frame(depth = l$depth, obs_age = l$obs_age)
    acc_mean <- stats::coef(MASS::rlm(obs_age ~ depth, data = d))[2]
    acc_mean <- signif(acc_mean, 2)
    if (acc_mean <= 0) {
      warning("Estimated mean accumulation rate is negative - using value = 20")
      acc_mean <- 20
    }
    l$acc_mean_prior <- acc_mean
  }
  ord <- order(l$depth)
  l$depth <- l$depth[ord]
  l$obs_age <- l$obs_age[ord]
  l$obs_err <- l$obs_err[ord]
  if (l$model_bioturbation == TRUE) {
    if (is.null(l$L_prior_sigma) == FALSE)
      message("L_prior_shape is being overriden by L_prior_sigma.")
    if (length(c(l$L_prior_sigma, l$L_prior_shape)) == 0)
      stop("One of either L_prior_sigma or L_prior_shape must be specified.\n             Set either to 0 to impose a fixed mixing depth.")
    if ((length(l$n_ind) == 1 | length(l$n_ind) == length(l$obs_age)) ==
        FALSE)
      stop("n_ind must be either a single value or a vector the same length as obs_age")
    if (length(l$n_ind == 1))
      l$n_ind <- rep(l$n_ind, length(l$obs_age))
    l$n_ind <- l$n_ind[ord]
    if (is.null(l$L_prior_sigma) == FALSE) {
      if (l$L_prior_sigma == 0)
        l$L_prior_shape <- 0
      else l$L_prior_shape <- gamma_sigma_shape(mean = l$L_prior_mean,
                                                sigma = l$L_prior_sigma)$shape
    }
  }
  else if (l$model_bioturbation == FALSE) {
    l$n_ind <- numeric(0)
  }
  if (is.null(l$infl_sigma_sd)) {
    l$infl_sigma_sd <- 10 * mean(l$obs_err)
  }
  if (l$min_age > min(l$obs_age)) {
    warning("min_age is older than minimum obs_age")
  }
  if (is.null(l$top_depth))
    l$top_depth <- l$depth[1]
  if (is.null(l$bottom_depth))
    l$bottom_depth <- utils::tail(l$depth, 1)
  depth_range <- l$bottom_depth - l$top_depth
  if (l$top_depth > min(l$depth))
    stop("top_depth must be above or equal to the shallowest data point")
  if (l$bottom_depth < max(l$depth))
    stop("bottom_depth must be deeper or equal to the deepest data point")
  if (is.null(l$K)) {
    K_fine_1 <- l$bottom_depth - l$top_depth
    min.d.depth <- stats::median(diff(sort(unique(l$depth))))
    K_fine_2 <- round(16 * K_fine_1/min.d.depth)
    K_fine <- min(c(K_fine_1, K_fine_2))
    if (K_fine > 900)
      K_fine <- 900
    l$K <- hamstr:::default_K(K_fine)
  }
  l$N <- length(l$depth)
  stopifnot(l$N == length(l$obs_err), l$N == length(l$obs_age))
  alpha_idx <- hamstr:::alpha_indices(l$K)
  l$K_tot <- sum(alpha_idx$nK)
  l$K_fine <- utils::tail(alpha_idx$nK, 1)
  l$c <- 1:l$K_fine
  l$mem_alpha = l$mem_strength * l$mem_mean
  l$mem_beta = l$mem_strength * (1 - l$mem_mean)
  l$mem_mean = l$mem_mean
  l$mem_strength = l$mem_strength
  l$delta_c = depth_range/l$K_fine
  l$c_depth_bottom = l$delta_c * l$c + l$top_depth
  l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K_fine - 1)])
  l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d)))
  l <- append(l, alpha_idx)
  l$n_lvls <- length(l$K)
  l$scale_shape = as.numeric(l$scale_shape)
  l$model_bioturbation = as.numeric(l$model_bioturbation)
  l$model_displacement = as.numeric(l$model_displacement)
  l$smooth_s = as.numeric(l$smooth_s)
  l$model_hiatus = as.numeric(l$model_hiatus)
  if (is.null(l$H_top))
    l$H_top = l$top_depth
  if (is.null(l$H_bottom))
    l$H_bottom = l$bottom_depth
  l$smooth_i <- hamstr:::get_smooth_i(l, l$L_prior_mean)
  l$I <- nrow(l$smooth_i)

  ## clymo model

  # first: construct sections for mass data
  c2 <-
    tibble::tibble(
      sample_depth_upper = l$depth2_upper,
      sample_depth_lower = l$depth2_lower,
      index_has_mass_measurements = TRUE
    )

  # add missing sections to make the sequence contiguous
  c2 <-
    dplyr::bind_rows(
      c2,
      purrr::map_dfr(seq_len(nrow(c2))[-1], function(i) {
        tibble::tibble(
          sample_depth_upper = c2$sample_depth_lower[i - 1],
          sample_depth_lower = c2$sample_depth_upper[i],
          index_has_mass_measurements = FALSE
        )
      })
    ) |>
    dplyr::mutate(
      sample_thickness = sample_depth_lower - sample_depth_upper,
      sample_depth_mean = sample_depth_upper + sample_thickness * 0.5
    ) |>
    dplyr::filter(sample_thickness != 0.0) |>
    dplyr::arrange(sample_depth_upper)

  l$index_clymo_par_constant <- which(c2$sample_depth_lower >= l$depth_clymo_par_constant & c2$sample_depth_upper <= l$depth_clymo_par_constant)[[1]]
  l$N2 <- sum(c2$index_has_mass_measurements)
  l$N2_c <- nrow(c2)
  l$index_has_mass_measurements <- which(c2$index_has_mass_measurements)
  l$depth2_all <- c(c2$sample_depth_upper, tail(c2$sample_depth_lower, 1L))
  l$which_c2_all <- sapply(l$depth2_all, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d)))
  l$which_c2_lower <- sapply(c2$sample_depth_lower, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d)))

  return(l)
}


#' Plots objects of class `peathamstr_fit`
#'
#' @param x An object of class `peathamstr_fit`.
#'
#' @param type A character value. Must be one of the plot `type`s available for
#' `hamstr_fit` objects (see `hamstr::plot.hamstr_fit()`) or the additional types
#' provided by the 'peathamstr' package:
#' \describe{
#'   \item{`"cumulative_mass_profile"`}{Plots the modeled and measured
#'    cumulative masses versus depth.}
#'   \item{`"mass_fluxes"`}{Plots the modeled and measured
#'    net mass uptake (NMU) and net mass release (NMR) versus average estimated
#'    ages. Shaded areas are 50\% and 95\% posterior intervals.}
#' }
#'
#' @param summarise See `hamstr::plot.hamstr_fit()`.
#'
#' @param ... Additional paramters (currently ignored).
#'
#' @return An object of class `ggplot`.
#'
#' @export
plot.peathamstr_fit <- function(
    x,
    type =
      c("default",
        "age_models",
        "acc_rates",
        "hier_acc_rates",
        "acc_mean_prior_post",
        "mem_prior_post",
        "L_prior_post",
        "D_prior_post",
        "PDF_14C",
        "cumulative_mass_profile",
        "mass_fluxes"
      ),
    summarise,
    ...)
{

  if(type %in% c("default", "age_models", "acc_rates", "hier_acc_rates", "acc_mean_prior_post", "mem_prior_post", "L_prior_post", "D_prior_post", "PDF_14C")) {
    NextMethod()
  } else {
    switch(
      type,
      "cumulative_mass_profile" = {

          dplyr::bind_rows(
            tibble::tibble(
              iter = 1L,
              depth_upper = x$data$depth2_upper,
              depth_lower = x$data$depth2_lower,
              depth = x$data$depth2_upper,
              y = x$data$layer_mass,
              variable_type = "observed"
            ),
            predict(
              object = x,
              type = "cumulative_mass",
              depth = "modelled",
              carbon_content = NULL
            ) |>
              dplyr::rename(
                y = "cumulative_mass",
                depth_upper = "depth2_upper",
                depth_lower = "depth2_lower"
              ) |>
              dplyr::mutate(
                variable_type = "predicted"
              )
          ) |>
          dplyr::group_by(
            .data$variable_type, .data$depth, .data$depth_upper, .data$depth_lower
          ) |>
          summarise_q(var = .data$y, probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)) |>
          ggplot2::ggplot(ggplot2::aes(y = .data$mean, x = .data$depth)) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, group = .data$variable_type), color = "grey50", width = 0) +
          ggplot2::geom_point(ggplot2::aes(fill = .data$variable_type), shape = 21) +
          ggplot2::labs(y = expression("Cumulative mass (kg m"^{-2}*")"), x = "Depth of upper layer boundary (cm)") +
          ggplot2::guides(color = ggplot2::guide_legend(title = ""))

      },
      "mass_fluxes" = {

        summary(
          object = x,
          type = type,
          probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)
        ) |>
          dplyr::mutate(
            variable = toupper(variable)
          ) |>
          ggplot2::ggplot(ggplot2::aes(y = .data$mean, x = .data$age)) +
          geom_hline(yintercept = 0, color = "grey50") +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = .data$variable), color = NA, alpha = 0.3) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$`25%`, ymax = .data$`75%`, fill = .data$variable), color = NA, alpha = 0.3) +
          ggplot2::geom_path(ggplot2::aes(color = .data$variable)) +
          ggplot2::labs(y = expression("Mass flux (kg m"^{-2}~yr^{-1}*")"), x = "Average estimated age (yr)") +
          ggplot2::guides(fill = ggplot2::guide_legend(title = "Mass flux"), color = ggplot2::guide_legend(title = "Mass flux"))

      })
  }


}
