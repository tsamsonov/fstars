library(stars)
library(fstars)

tif = system.file("tif/L7_ETMs.tif", package = "stars")
x = read_stars(tif)
plot(x)

attr(x, "dimensions")

get_factors(x)
