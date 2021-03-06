% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proj_factors.R
\name{st_deriv}
\alias{st_deriv}
\title{Calculate surface derivatives}
\usage{
st_deriv(s, stats = "slope", method = "ZEVENBERGEN", adaptive = FALSE)
}
\arguments{
\item{s}{a raster stars object}

\item{stats}{a vector of derivatives names}

\item{method}{surface fitting method. Currentlyr EVANS, ZEVENBERGEN and FLORINSKY methods are supported. ZEVENBERGEN is used by default}

\item{adaptive}{should a kernel be adaptive to map projection distorions?}
}
\value{
a stars object
}
\description{
Calculate surface derivatives
}
