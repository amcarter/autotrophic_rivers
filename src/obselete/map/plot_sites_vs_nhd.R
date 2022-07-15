library(tidyverse)
library(nhdplusTools)
library(sf)
library(mapview)
mv <- mapview::mapview

site <- sf::st_point(c(-101.4462223, 29.8029823)) %>%
    sf::st_sfc(crs = 4326)

comid <- nhdplusTools::discover_nhdplus_id(site)

buf_dist = 10
site_buf <- sf::st_buffer(x = site,
                          dist = buf_dist)

site_box <- st_bbox(site_buf)

subset <- nhdplusTools::subset_nhdplus(bbox = site_box,
                                       output_file = tempfile(fileext = '.gpkg'),
                                       nhdplus_data = 'download',
                                       return_data = TRUE,
                                       overwrite = FALSE,
                                       out_prj = 4326) %>%
    suppressMessages()

## NHD HR stuff
huc12 <- get_huc12(site)
print(huc12$huc12)
huc4 <- substr(huc12$huc12, 1, 4)[1]
nhdplusTools::download_nhdplushr(nhd_hr_dir,
                                 hu_list = huc4) %>%
    invisible()

HRflowlines <- nhdplusTools::get_nhdplushr(file.path(nhd_hr_dir,
                                                     substr(huc4, 1, 2)),
                                           file.path(nhd_hr_dir,
                                                     paste0(huc4, '.gpkg')),
                                           layers = 'NHDFlowline',
                                           proj = 4326)$NHDFlowline

NHD_HR <- suppressWarnings(sf::st_crop(HRflowlines, site_box))

NHDPlus <- subset$NHDFlowline_Network

xx = mv(NHD_HR, color = 'darkslategray3') + mv(NHDPlus, color='deepskyblue4') + mv(site, color='red')
mapview_save_path <- file.path(mapview_save_dir,
                               paste0(dmn, '_', site_code, '.html'))
mapview::mapshot(xx,
                 url = mapview_save_path)
print(xx)
