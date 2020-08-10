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
  rcpp_interpolate_xy(s[[1]], attr(s, "dimensions"), x, y)
}

#' Returns a stars object with projection factors calculated at each pixel as variables
#'
#' @param s stars object
#'
#' @return stars object
#' @export
get_factors <- function(s) {
  pf = rcpp_get_factors(attr(s, "dimensions"), st_crs(s)$proj4string)
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

#' Calculate statistics using floating window
#'
#' @param s a raster stars object
#' @param stats a vector of stats names
#' @param size a size of the kernel. 3 by defaulr
#' @param adaptive should a kernel be adaptive to map projection distorions?
#'
#' @return a stars object
#' @export
#'
#' @examples
st_convolve <- function(s, stats = 'mean', size = 3, adaptive = FALSE) {
  krnl = matrix(rep(size^-2, size^2), nrow = size)

  res = rcpp_filter_matrix(s[[1]], attr(s, "dimensions"), st_crs(s)$proj4string, size, stats,
                           attr(attr(s, "dimensions"), "raster")[["curvilinear"]],
                           adaptive, FALSE) %>%
    st_as_stars()

  attr(res, "dimensions")[[1]]$delta = attr(s, "dimensions")[[1]]$delta
  attr(res, "dimensions")[[2]]$delta = attr(s, "dimensions")[[2]]$delta
  attr(res, "dimensions")[[1]]$offset = attr(s, "dimensions")[[1]]$offset
  attr(res, "dimensions")[[2]]$offset = attr(s, "dimensions")[[2]]$offset
  attr(res, "dimensions")[[1]]$refsys = attr(s, "dimensions")[[1]]$refsys
  attr(res, "dimensions")[[2]]$refsys = attr(s, "dimensions")[[2]]$refsys
  set_names(attr(res, "dimensions"), names(attr(s, "dimensions")))

  return(res)
}

#' Calculate surface derivatives
#'
#' @param s a raster stars object
#' @param stats a vector of derivatives names
#' @param method surface fitting method. Currentlyr EVANS, ZEVENBERGEN and FLORINSKY methods are supported. ZEVENBERGEN is used by default
#' @param adaptive should a kernel be adaptive to map projection distorions?
#'
#' @return a stars object
#' @export
#'
#' @examples
st_deriv <- function(s, stats = 'slope', method = 'ZEVENBERGEN', adaptive = FALSE) {

  res = rcpp_filter_matrix(s[[1]], attr(s, "dimensions"), st_crs(s)$proj4string, 3, stats,
                           attr(attr(s, "dimensions"), "raster")[["curvilinear"]],
                           adaptive && (method != "FLORINSKY"), TRUE, method) %>%
    st_as_stars()

  attr(res, "dimensions")[[1]]$delta = attr(s, "dimensions")[[1]]$delta
  attr(res, "dimensions")[[2]]$delta = attr(s, "dimensions")[[2]]$delta
  attr(res, "dimensions")[[1]]$offset = attr(s, "dimensions")[[1]]$offset
  attr(res, "dimensions")[[2]]$offset = attr(s, "dimensions")[[2]]$offset
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
