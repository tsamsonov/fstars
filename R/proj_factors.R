#' Returns a stars object with projection factors calculated at each pixel as variables
#'
#' @param s stars object
#'
#' @return stars object
#' @export
get_factors <- function(s) {
  get_factors_stars(attr(s, "dimensions"))
}

#' Tests external proj functionality
#'
#' @return nothing
#' @export
test_rproj <- function() {
  test_proj()
}
