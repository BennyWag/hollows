# figures -----------------------------------------------------------------

library(sf)
library(terra)
#library(tidyterra)
library(raster)
library(viridis)
library(shadowtext)
library(ggrepel)
library(ggspatial)
library(RColorBrewer)
library(tidyverse)
library(ggnewscale)
library(gtsummary)
library(gt)
library(patchwork)

# study area map - prep----------------------------------------------------

#load data = shapes

sites_spatial<-st_read('data/post_fire/spatial/sites_all.gpkg')%>%
  dplyr::mutate(site_ID = str_replace(site_ID, 'BBBR\\b', 'BBRR'),
                firesev = case_when(firesev == 0 ~ 2,
                                    TRUE ~ firesev))

sites_spatial_pre<-read.csv('data/pre_fire/plot_data.csv')%>%
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = crs(sites_spatial))%>%
  mutate(transect = case_when(transect == 'T3' ~ 'high-country',
                              transect == 'T2' ~ 'mid-hills',
                              transect == 'T1' ~ 'lowlands'),
         transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))


EG<-st_read('data/spatial/rfa/rfa25.shp')%>%st_transform(crs = crs(sites_spatial))

AU<-st_read('data/spatial/AU/Australia_proj.shp')%>%st_transform(crs = crs(sites_spatial))

places<-st_read('data/spatial/places/EG_places.shp')%>%st_transform(crs = crs(sites_spatial))%>%select(NAME, geometry)

st_x = function(x) st_coordinates(x)[,1]
st_y = function(x) st_coordinates(x)[,2]

places$x<-st_x(places)
places$y<-st_y(places)

places$x2<-places$x
places$y2<-places$y

places[1,5]<-148.85
places[5,5]<-148.95

places[3,5]<-148.05
places[3,6]<-(-37.82)

#rasters

DEM<-rast('data/spatial/DEM_VIC.tif')

#black_summer<-rast('data/spatial/fire_eg_comp.tif')


#crop to study area

EG_crop<-EG%>%filter(NAME %in% c('GIPPSLAND', 'EAST GIPPSLAND'))%>%st_crop(xmin = 147.5, xmax = 150,
                                                                           ymin = -36.5, ymax = -38)

#st_write(EG_crop, 'outputs/spatial/study_area.gpkg')

#fire_eg<-raster::crop(black_summer, EG_crop)
#fire_eg<-raster::mask(fire_eg, EG_crop)

DEM_EG<-terra::crop(DEM, EG_crop)
DEM_EG<-terra::mask(DEM_EG, vect(EG_crop))

DEM_df<-DEM_EG%>%as.data.frame(xy = T, na.rm = T)

# map ---------------------------------------------------------------------

sf_use_s2(FALSE)

map<-ggplot()+
  geom_raster(data = DEM_df, aes(x = x, y = y, fill = dem.9s))+
  geom_shadowtext(data = places, aes(x = x2, y = y2, label = NAME), 
                  size = 3, col = "white", fontface = "bold")+
  geom_sf(data=sites_spatial, size = 4, shape = 15,  aes(color = as.factor(firesev)))+
  scale_colour_brewer(palette = 'OrRd', name  = "Fire severity",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  new_scale_color()+
  geom_sf(data=sites_spatial_pre, size = 1.5, aes(color = as.factor(transect)))+
  scale_colour_brewer(palette = 'Set2', name  ='Pre-fire band',labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_fill_viridis(option = 'E', name = 'Elevation (m)')+
  labs(x = NULL, y = NULL, title = NULL)+
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.002, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  annotation_scale(location = 'br', style = 'ticks', text_col = 'black')+
  coord_sf(expand = F)+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,0,0),"pt"),
        axis.text.x = element_text(size = 12, face = 'bold'), 
        axis.text.y = element_text(size = 12, face = 'bold'),
        legend.title = element_text(size = 13, face = 'bold'),
        legend.text = element_text(size = 11, face = 'bold'))

#inset

#get extent as bounding box

ext<-extent(EG_crop)%>%as('SpatialPolygons')%>%st_as_sf()
st_crs(ext)<-crs(sites_spatial)

#create raster df

inset<-ggplot()+
  geom_sf(data = AU, show.legend = "point")+
  geom_sf(data = ext,lwd=1, color = 'red', fill = NA)+
  #annotate("text", x = 147, y = -34, label= "Study area", size = 5, col = 'black')+
  labs(x = '', y = '', title = '')+
  coord_sf(crs = 4326)+
  theme_void()

#combine

map+annotation_custom(ggplotGrob(inset), xmin = 149.4, xmax = 149.95, ymin = -37, ymax = -36.38)

ggsave('map_new.svg', path = 'figures',
       width = 25, height = 16, units = 'cm', dpi = 600)


# tables ------------------------------------------------------------------

t1<-tbl_regression(hollowmod_final2_dredge$best_model, intercept = T)
t2<-tbl_regression(hollowmod_final2_num_dredge$best_model, intercept = T)

hol_final_post <-
  tbl_merge(
    tbls = list(t1, t2),
    tab_spanner = c("**Hollow occurrence**", "**Hollow abundance**")
  )%>%as_gt()

gtsave(hol_final_post, 'tables/hollows_final_post.rtf')
