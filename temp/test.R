library(stars)
library(fstars)
library(tidyverse)
library(classInt)
library(mapview)

data(land, package = 'tmap')

# box = st_bbox(c(xmin = 30, xmax = 60, ymin = 50, ymax = 60), crs = st_crs(land))
# box = st_bbox(c(xmin = -12, xmax = 60, ymin = 30, ymax = 72), crs = st_crs(land))
# box = st_bbox(c(xmin = -70, xmax = --10, ymin = 75, ymax = 84), crs = st_crs(land))
box = st_bbox(c(xmin = -12, xmax = 60, ymin = 28, ymax = 75), crs = st_crs(land))

prj = '+proj=merc'

landp = land %>%
  st_crop(box) %>%
  st_warp(crs = prj)

fct = get_factors(landp)

pal = c("#003200", "#3C9600", "#006E00", "#556E19", "#00C800", "#8CBE8C",
        "#467864", "#B4E664", "#9BC832", "#EBFF64", "#F06432", "#9132E6",
        "#E664E6", "#9B82E6", "#B4FEF0", "#646464", "#C8C8C8", "#FF0000",
        "#FFFFFF", "#5ADCDC")

plot(fct['tissot_orientation'])

# ggplot() +
#   geom_stars(data = landp['cover']) +
#   scale_fill_manual(values = pal, guide = guide_legend(ncol = 3), name = NULL) +
#   coord_sf(crs = st_crs(landp)) +
#   theme(legend.position = 'bottom')

# vars = c(
  # "meridional_scale",
  # "parallel_scale",
  # "areal_scale",
  # "angular_distortion",
  # "meridian_parallel_angle",
  # "meridian_convergence"
  # "parallel_convergence",
  # "tissot_semimajor",
  # "tissot_semiminor",
  # "tissot_orientation"
  # "dx_dlam",
  # "dx_dphi",
  # "dy_dlam",
  # "dy_dphi"
# )


# for (var in vars) {
#   brks = unique(classIntervals(fct[[var]], n = 10, style = 'fisher')$brks)
#
#   cont = st_contour(
#     fct[var],
#     na.rm = TRUE,
#     contour_lines = TRUE,
#     breaks = brks
#   )
#
#   g = ggplot() +
#     geom_stars(data = cut(fct[var], breaks = brks)) +
#     geom_sf(data = cont) +
#     coord_sf(crs = st_crs(fct)) +
#     ggtitle(var)
#
#   print(g)
# }

f = st_convolve(landp['elevation'], size = 9)
f1 = st_convolve(landp['elevation'], size = 21, adaptive = TRUE)
f2 = st_convolve(landp['elevation'], size = 9, adaptive = TRUE)
f3 = st_convolve(landp['elevation'], size = 3, adaptive = TRUE)
# f3 = st_convolve(landp['elevation'], size = 7, adaptive = TRUE)
# f4 = st_convolve(landp['elevation'], size = 11, adaptive = TRUE)

plot(landp['elevation'])
plot(f3)
plot(f2)
plot(f1)
plot(f)

# df = f2 - f
#
# mapview(f)
# mapview(f2)
# mapview(df)

# plot(f3)
# plot(f4)

# ggplot() +
#   geom_stars(data = f)

# pt = st_point(c(45, 60), dim = "XY") %>%
#   st_sfc(crs = 4326) %>%
#   st_transform(prj)
#
# mapview(f) + mapview(pt)
#
# coord = st_coordinates(pt)
#
# interpolate_xy(f, coord[1], coord[2])
