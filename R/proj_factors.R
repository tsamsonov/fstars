#' Extract value using a bilinear interpolation
#'
#' @param s stars object
#' @param x x coordinate to be interpolated
#' @param y y coordinate to be interpolated
#'
#' @return interpolated value
#' @export
#'
#' @examples
interpolate_xy <- function(s, x, y) {
  cpp_interpolate_xy(s[[1]], attr(s, "dimensions"), x, y)
}

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
st_convolve <- function(s, kernel = 'mean', size = 3, fun = NULL) {
  krnl = matrix(rep(size^-2, size^2), nrow = size)
  res = filter_matrix(s[[1]], krnl) %>%
    st_as_stars()
  attr(res, "dimensions")[[1]]$delta = attr(s, "dimensions")[[1]]$delta
  attr(res, "dimensions")[[2]]$delta = attr(s, "dimensions")[[2]]$delta
  attr(res, "dimensions")[[1]]$offset = attr(s, "dimensions")[[1]]$offset + attr(s, "dimensions")[[1]]$delta * size %/% 2
  attr(res, "dimensions")[[2]]$offset = attr(s, "dimensions")[[2]]$offset + attr(s, "dimensions")[[2]]$delta * size %/% 2
  attr(res, "dimensions")[[1]]$refsys = attr(s, "dimensions")[[1]]$refsys
  attr(res, "dimensions")[[2]]$refsys = attr(s, "dimensions")[[2]]$refsys
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
