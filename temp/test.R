library(stars)
library(fstars)
library(tidyverse)
library(classInt)

data(land, package = 'tmap')
box = st_bbox(c(xmin = -12, xmax = 60, ymin = 30, ymax = 72), crs = st_crs(land))

landp = land %>%
  st_crop(box) %>%
  st_transform(crs = '+proj=eck3')

fct = get_factors(landp)

pal = c("#003200", "#3C9600", "#006E00", "#556E19", "#00C800", "#8CBE8C",
        "#467864", "#B4E664", "#9BC832", "#EBFF64", "#F06432", "#9132E6",
        "#E664E6", "#9B82E6", "#B4FEF0", "#646464", "#C8C8C8", "#FF0000",
        "#FFFFFF", "#5ADCDC")

ggplot() +
  geom_stars(data = landp['cover']) +
  scale_fill_manual(values = pal, guide = guide_legend(ncol = 3), name = NULL) +
  coord_sf(crs = st_crs(landp)) +
  theme(legend.position = 'bottom')

var = 'areal_scale'
brks = classIntervals(sample(fct[[var]], 10000), n = 10, style = 'fisher')$brks

cont = st_contour(
  fct[var],
  na.rm = TRUE,
  contour_lines = TRUE,
  breaks = brks
)

ggplot() +
  geom_stars(data = cut(fct[var], breaks = brks)) +
  geom_sf(data = cont)
  coord_sf(crs = st_crs(fct))

