
library(dplyr)
library(nhdplusTools)
library(sf)
library(tmap)

nwis_id <- 'USGS-08447410'
site <- sf::st_point(c(-101.4462223, 29.8029823)) %>%
  sf::st_sfc(crs = 4326)

start_comid <- nhdplusTools::discover_nhdplus_id(site)

# Get upstream mainstem and tributary flowlines for this point
flowline <- nhdplusTools::navigate_nldi(list(featureSource = "comid", 
                                             featureID = start_comid), 
                                        mode = "upstreamMain", 
                                        distance_km = 9999)


usgs_sites <- navigate_nldi(list(featureSource = "comid", 
                                 featureID = start_comid), 
                            mode = "upstreamMain", "nwissite", 
                            distance = 9999) #returns all upstream USGS gages

plot_nhdplus(as.integer(flowline$nhdplus_comid))


bbox <- sf::st_bbox(flowline)

bbox <- bbox + c(-.1, -.1, .1, .1)

dat <- plot_nhdplus()

outlet_fline <- dat$flowline %>%
  filter(Hydroseq == min(Hydroseq))

plot(st_transform(st_geometry(outlet_fline), 3857), add = TRUE, lwd = 3)

plot(st_transform(start_point, 3857), add = TRUE, pch = 19)
get_UM()
flowline <- navigate_nldi(cur_site, mode = "UM", 
                          data_source = 'flowlines') #returns all upstream flowlines

tm_shape(flowline) + 
  tm_lines() +
  tm_shape(usgs_sites) +
  tm_dots(col="red") 
