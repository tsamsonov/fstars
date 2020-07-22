library(stars)
library(fstars)
library(sf)

tif = system.file("tif/L7_ETMs.tif", package = "stars")
x = read_stars(tif)
plot(x)

attr(x, "dimensions")[[1]][['refsys']][['wkt']]

get_factors(x)
test_rproj()
