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
           tissot_semimajor = pf[['tissot_semimajor']],
           tissot_semiminor = pf[['tissot_semiminor']],
           dx_dlam = pf[['dx_dlam']],
           dx_dphi = pf[['dx_dphi']],
           dy_dlam = pf[['dy_dlam']],
           dy_dphi = pf[['dy_dphi']])
}

#' Tests external proj functionality
#'
#' @return nothing
#' @export
test_rproj <- function() {
  test_proj()
}
