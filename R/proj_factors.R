#' Returns a stars object with projection factors calculated at each pixel as variables
#'
#' @param s stars object
#'
#' @return stars object
#' @export
get_factors <- function(s) {
  pf = get_factors_stars(attr(s, "dimensions"), st_crs(s)$proj4string)
  spf = s %>%
    mutate(meridional_scale = pf[['meridional_scale']],
           parallel_scale = pf[['parallel_scale']],
           areal_scale = pf[['areal_scale']],
           angular_distortion = pf[['angular_distortion']],
           meridian_parallel_angle = pf[['meridian_parallel_angle']],
           meridian_convergence = pf[['meridian_convergence']],
           parallel_convergence = pf[['parallel_convergence']],
           tissot_semimajor = pf[['tissot_semimajor']],
           tissot_semiminor = pf[['tissot_semiminor']],
           tissot_orientation = pf[['tissot_orientation']],
           dx_dlam = pf[['dx_dlam']],
           dx_dphi = pf[['dx_dphi']],
           dy_dlam = pf[['dy_dlam']],
           dy_dphi = pf[['dy_dphi']])
}

#' Filter stars object with selected kernel
#'
#' @param s a raster stars object
#' @param kernel string (name of the filter), matrix (with values)
#' @param size
#'
#' @return matrix
#' @export
#'
#' @examples
st_convolve <- function(s, kernel = 'mean', size = c(3, 3), fun = NULL) {
  krnl = matrix(rep(1/9, 9), nrow = 3)
  res = filter_matrix(s[[1]], krnl) %>%
    st_as_stars()
  attr(res, "dimensions")[[1]]$delta = attr(s, "dimensions")[[1]]$delta
  attr(res, "dimensions")[[2]]$delta = attr(s, "dimensions")[[2]]$delta
  set_names(attr(res, "dimensions"), names(attr(s, "dimensions")))
  return(res)
}

#' Tests external proj functionality
#'
#' @return nothing
#' @export
test_rproj <- function() {
  test_proj()
}
