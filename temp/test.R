library(stars)
library(fstars)
library(tidyverse)
library(classInt)
library(mapview)

data(land, package = 'tmap')

# land = read_stars('/Volumes/Data/Spatial/DEM/HYPSO/elevation/hypso5g.tif', proxy = TRUE) %>%
#   select(elevation = 1)

# box = st_bbox(c(xmin = -20, xmax = 60, ymin = 30, ymax = 70), crs = st_crs(land))
# box = st_bbox(c(xmin = 5, xmax = 13, ymin = 43.5, ymax = 48), crs = st_crs(land))
box = st_bbox(c(xmin = -70, xmax = -15, ymin = 60, ymax = 85), crs = st_crs(land))

# box = st_bbox(c(xmin = -160, xmax = 160, ymin = -60, ymax = 85), crs = st_crs(land))

cland = land %>%
  st_crop(box) %>%
  st_as_stars()

prj = 32622 #'+proj=eqc'
#
# landp = land %>%
#   st_warp(crs = prj)

clandp = cland %>%
  st_warp(crs = prj)

# fct = get_factors(cland)

# f0 = st_convolve(clandp['elevation'], size = 9)
# f1 = st_convolve(clandp['elevation'], size = 3, adaptive = TRUE)

f0f = st_deriv(cland['elevation'], 'slope', method = "FLORINSKY")
f0e = st_deriv(clandp['elevation'], 'slope', method = "EVANS")
f0z = st_deriv(clandp['elevation'], 'slope', method = "ZEVENBERGEN")
# f0ea = st_deriv(clandp['elevation'], 'slope', method = "EVANS", adaptive = TRUE)
# f0za = st_deriv(clandp['elevation'], 'slope', method = "ZEVENBERGEN", adaptive = TRUE)

# f1 = st_convolve(clandp['elevation'], 'mean', size = 7)
# f1a = st_convolve(clandp['elevation'], 'mean', size = 7, adaptive = TRUE)
# f1s = st_deriv(f1, 'aspect')
# f1as = st_deriv(cland, 'aspect', adaptive = TRUE)

plot(cland['elevation'])
plot(f0f)
plot(f0e)
plot(f0z)
# plot(f0ea)
# plot(f0za)

write_stars(cland['elevation'], 'temp/cland.tif')
write_stars(clandp['elevation'], 'temp/clandp.tif')
write_stars(f0f,  'temp/slope_flor.tif')
write_stars(f0za, 'temp/slope_zeven.tif')
write_stars(f0ea, 'temp/slope_evans.tif')

write_stars(f0z, 'temp/slope_zeven_fix.tif')
write_stars(f0e, 'temp/slope_evans_fix.tif')

# mapview(f0f)
# plot(f1)
# plot(f1a)
# plot(f1s)
# plot(f1as)

# write_stars(f1, 'temp/filtered.tif')
# write_stars(f1a, 'temp/filtered_a.tif')
#
# write_stars(f1s, 'temp/aspect.tif')
# write_stars(f1as, 'temp/aspect_a.tif')

# View(clandp[['elevation']])
# View(f1[[1]])
# View(f1a[[1]])

# f1s = st_convolve(f1, 'slope', adaptive = TRUE)



# f0h = st_convolve(f0, 'hill', adaptive = TRUE)
# write_stars(f1s, 'temp/slope.tif')


# plot(f1)
# plot(f0)
# plot(f1s)
# plot(f0h)
# plot(f1)
# plot(f1s)


#
# mapview(clandp['elevation']) + mapview(f1) + mapview(f1a) + mapview(f1a - f1)
# mapview(f1) + mapview(f1a)


# plot(f1)
# plot(f1s)

# pal = c("#003200", "#3C9600", "#006E00", "#556E19", "#00C800", "#8CBE8C",
#         "#467864", "#B4E664", "#9BC832", "#EBFF64", "#F06432", "#9132E6",
#         "#E664E6", "#9B82E6", "#B4FEF0", "#646464", "#C8C8C8", "#FF0000",
#         "#FFFFFF", "#5ADCDC")

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

# mapview(f3[1])


# plot(f2)
# plot(f1)
# plot(f)

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
