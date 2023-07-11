#' Extract Stan Parameter Index
#'
#' @param x vector of stan parameter names
#'
#' @return a vector
#' @keywords internal
#' @noRd
get_par_idx <- function(x) {

  a <- strsplit(x, "\\[")
  a <- sapply(a, function(x) x[2])
  b <- unlist(strsplit(a, "\\]"))

  as.numeric( b )
}


#' Translates between two sets of depth profiles with different layer boundaries
#'
#' @keywords internal
#' @noRd
translate_depths <- function(layers_from, layers_to) {

  J  = nrow(layers_from)
  N = nrow(layers_to)
  res = matrix(rep(0.0, N * J), nrow = N, ncol = J)
  from_contains_to_upper = numeric(J)
  from_contains_to_lower = numeric(J)

  for(n in 1:N) {

    index_start = 0L
    index_end = 0L

    for(j in 1:J) {
      from_contains_to_upper[j] = (layers_from[j, 2] >= layers_to[n, 1]) && (layers_from[j, 1] <= layers_to[n, 1])
      from_contains_to_lower[j] = (layers_from[j, 2] >= layers_to[n, 2]) && (layers_from[j, 1] <= layers_to[n, 2])
    }

    if(sum(from_contains_to_upper) == 0 && sum(from_contains_to_lower) == 0) {
      next
    }

    if(sum(from_contains_to_upper) > 0) {
      cur_index = 0
      while(!cur_index) {
        index_start = index_start + 1
        if(from_contains_to_upper[index_start]) {
          cur_index = 1
        }
      }
    } else {
      index_start = 1
    }

    if(sum(from_contains_to_lower) > 0) {
      cur_index = 0
      while(!cur_index) {
        index_end = index_end + 1
        if(from_contains_to_lower[index_end]) {
          cur_index = 1
        }
      }
    } else {
      index_end = J
    }

    if(index_start == index_end && sum(from_contains_to_lower) > 0) {
      res[n, index_end] = (layers_to[n, 2] - layers_to[n, 1]) / layers_from[index_end, 3]
    } else {
      res[n, index_start:index_end] = rep(1, index_end - index_start + 1)
      res[n, index_start] = (layers_from[index_start, 2] - layers_to[n, 1]) / layers_from[index_start, 3]
      if(sum(from_contains_to_lower) > 0) {
        res[n, index_end] = (layers_to[n, 2] - layers_from[index_end, 1]) / layers_from[index_end, 3]
      } else {
        res[n, index_end] = 1
      }
    }

  }

  res

}


#' Summarise to Quantiles and Moments
#'
#' @param dat a dataframe or tibble
#' @param var the variable to summarise (unquoted)
#' @param probs the quantiles at which to summarise
#'
#' @return a tibble containing the summarised posterior
#' @keywords internal
summarise_q <- function(dat,
                        var,
                        probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)
                        )
{
  dat |>
    dplyr::summarise(mean = mean({{ var }}, na.rm = TRUE),
                     sd = stats::sd({{ var }}, na.rm = TRUE),
                     x = stats::quantile({{ var }}, probs, na.rm = TRUE),
                     q = paste0(round(100*probs, 1), "%")) |>
    tidyr::pivot_wider(names_from = .data$q, values_from = .data$x) |>
    dplyr::as_tibble()
}
