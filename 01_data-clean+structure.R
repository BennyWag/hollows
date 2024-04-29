# data clean + structure --------------------------------------------------

library(openxlsx)
library(sf)
library(raster)
library(terra)
library(fasterize)
library(RColorBrewer)
library(tidyverse)
library(patchwork)

sites<-read.csv('data/post_fire/transects_all.csv')%>%
  mutate('x' = 
           as.numeric(str_split_fixed(sites$phone_location, ',', n = 2)[,2]),
         'y' = 
           as.numeric(str_split_fixed(sites$phone_location, ',', n = 2)[,1]))%>%
  dplyr::mutate(site_ID = str_replace(site_ID, 'BBBR\\b', 'BBRR'))

sites_spatial<-st_read('data/post_fire/spatial/sites_all.gpkg')%>%
  dplyr::mutate(site_ID = str_replace(site_ID, 'BBBR\\b', 'BBRR'))

structure<-read.csv('data/post_fire/structure_final.csv')


# data prep ---------------------------------------------------------------


structure_fix<-structure%>%mutate(site_ID = str_trim(site_ID, 'right'))%>%
  dplyr::mutate(site_ID = str_replace(site_ID, 'BBR\\b', 'BBRR'),
                site_ID = case_when(site_ID == 'G3 F4' ~ 'G3F4',
                                    site_ID == 'G3 F0' ~ 'G3F0',
                                    site_ID == '\nT15P4' ~ 'T15P4',
                                    site_ID == '\nT1P1' ~ 'T1P1',
                                    site_ID == ' T1P3' ~ 'T1P3',
                                    site_ID == 'T10' ~ 'T1P3',
                                    site_ID == 'ARI BBR04' ~ 'ARI BBRR 04',
                                    site_ID == 'ARI BBRR 45a' ~ 'ARI BBRR 45',
                                    
                                    TRUE ~ site_ID))

#check site data

sites_struc<-as.data.frame(unique(structure_fix$site_ID))%>%rename(site_ID = 1)

sites_all_list<-sites%>%filter(type == 'measure')%>%dplyr::select(site_ID)

site_miss<-sites_struc%>%dplyr::filter(!site_ID %in% sites_all_list$site_ID)

#join in GGs and firesev

structure_join<-left_join(structure_fix%>%mutate(plot = as.factor(site_ID)), 
                          sites_spatial%>%mutate(plot = as.factor(site_ID))%>%
                            dplyr::select(plot, gg_pres, n_gg, firesev, 
                                          elevation, aspect, slope, mid_understory_density,
                                          layers, foliagecover), by = 'plot')%>%
  mutate(species = str_trim(species, 'right'))%>%
  select(-geom)

#tree species

unique(structure_join$species)

structure_join_fix<-structure_join%>%
  mutate(species = case_when(species == 'Eucalyptus croajingolensis6.' ~ 'Eucalyptus croajingolensis',
                             species == 'unknow' ~ 'unknown',
                             species == ' fastigata' ~ 'Eucalyptus fastigata',
                             species == 'Acacia dealbatal' ~ 'Acacia dealbata',
                             species == 'Eucalyptus  cypellocarpa' ~ 'Eucalyptus cypellocarpa',
                             species == 'Acacia melaxnomia' ~ 'Acacia melanoxylon',
                             species == 'Acacia melaxnomia' ~ 'Acacia melanoxylon',
                             species == 'Hydecarpia angustifolia' ~ 'Hedycarya angustifolia',
                             species == 'Tasmania lanceolata' ~ 'Tasmannia lanceolata',
                             species == 'Tasmania spp.' ~ 'Tasmannia spp',
                             species == 'Proathanthera lasianthas' ~ 'Prostanthera lasianthos',
                             species == 'banksia spp' ~ 'Banksia spp',
                             TRUE ~ species))

unique(structure_join_fix$species)

#measurements

structure_join_fix<-structure_join_fix%>%
  mutate(dbh = case_when(dbh>5000 ~ dbh/100,
                         dbh>300 ~ dbh/10,
                         TRUE ~ dbh),
         height = case_when(height > 70 ~ height/10,
                            height < 0 ~ 1,
                            TRUE ~ height),
         hollow_num = case_when(is.na(hollow_num) ~ 0,
                                TRUE ~ as.numeric(hollow_num)))

# add additional data -----------------------------------------------------

#spatial data

spatial<-rast('data/spatial/gg_all.tif')
climate<-rast('data/spatial/climate_model.tif')

crs(spatial)<-crs(climate)

#rasterize soil layer

gridsize <- res(climate)[1]

soil_shp<-st_read('data/spatial/soil_type_4326.gpkg')%>%
  mutate(ASC_O_num = as.numeric(as.factor(ASC_O)),
         ASC_SO_num = as.numeric(as.factor(ASC_SO)))

soil_rast1 <- raster(soil_shp, res=gridsize)
soil_rast2 <- raster(soil_shp, res=gridsize)


soil_rast_O<-fasterize(soil_shp, soil_rast1, field = 'ASC_O_num',
                 fun = 'max')

soil_rast_SO<-fasterize(soil_shp, soil_rast2, field = 'ASC_SO_num',
                       fun = 'max')

soil_stack<-c(rast(soil_rast_O), rast(soil_rast_SO))

soil_stack<-resample(soil_stack, climate, method = 'near')

rm(soil_shp)

writeRaster(soil_stack, 'outputs/spatial/soil_types.tif')

#combine

spatial<-c(spatial, climate, soil_stack)

names(spatial)<-c('AHMI', 'cond', 'et0', 'hn', 
                  'rain', 'vpd', 'FT', 'ph', 
                  'slope', 'TWI', 'fires', 'NDVI', 
                  'climsuit', 'ASC_O', 'ASC_SO')

sites_spatial_data<-terra::extract(spatial, vect(sites_spatial), 
                                   fun = NULL , method = 'simple',  xy = T, ID = F)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  mutate(fires = round(fires))

sites_spatial_data<-cbind(sites_spatial, sites_spatial_data)

#add soil nitrogen 

soilchem<-read.csv('data/post_fire/soils_tidy.csv')%>%
  mutate(Nitrate.Nitrogen = case_when(is.na(Nitrate.Nitrogen) ~ 0.9,
                                            TRUE ~ as.numeric(Nitrate.Nitrogen)))

sites_spatial_data<-left_join(sites_spatial_data, soilchem, by = 'site_ID')


#add leaf N (per tree and site average)

leafchem<-read.csv('data/post_fire/leafsample_structure.csv')

#plot averages

leafav<-leafchem%>%group_by(site_ID)%>%summarise(mean_leaf_n = mean(total_n),
                                              mean_leaf_c = mean(c))


sites_spatial_data<-left_join(sites_spatial_data, leafav, by = 'site_ID')%>%
  mutate(firesev = case_when(site_ID == 'T35P3' ~ 5,
                             site_ID == 'T1P1' ~ 3,
                             site_ID == 'T2P1' ~ 3,
                             site_ID == 'ARI BBRR 26' ~ 3,
                             site_ID == 'ARI BBRR 36' ~ 3,
                             site_ID == 'ARI BBRR 22' ~ 4,
                             site_ID == 'ARI BBRR 45' ~ 4,
                             site_ID == 'T1P2' ~ 6,
                             TRUE ~ firesev))

#add to single trees in structure sheet

structure_join_fix_leaf<-left_join(structure_join_fix, leafchem%>%select(c, total_n, bag_ID), by = 'bag_ID')

#combine and save

structure_join_fix_final<-left_join(structure_join_fix_leaf, 
                                    sites_spatial_data%>%select(site_ID, 35:77), 
                                    by = 'site_ID')%>%
  select(-geom)


write.csv(structure_join_fix_final, 'outputs/structure_fixed.csv', row.names = F)


# point centered quarter method calculations ------------------------------

struc<-structure_join_fix_final%>%mutate(firesev = case_when(firesev == 0 ~ 2,
                                                             TRUE ~ firesev),
                                         sampling_point = case_when(sampling_point>10 ~ 3,
                                                                    TRUE ~ as.numeric(sampling_point)))


#bitterlich calc

#Fi = pi*c^2*(dbh/100)^2

#distance = 21.8 m, dbh = 68.7 cm, radius = dbh/2 = 34.35

#angle (sin angle/2) = (dbh/2) / distance, (0.687/2) / 21.8 = 0.015

#expansion factor = 10000 * (sin angle/2)^2 = 10000*(0.015^2) = 2.25

#c = 50/sqrt(k) = 50/sqrt(2.25) = 33.33

#rep no of trees: pi*33.33^2*(68.7/100)^2 = Fi = 1647

# EFI = 10000/FI = 10000/1647 = 6.07


#number of trees per species accodring to PCQM documentation

struc_sample<-struc%>%filter(site_ID == 'ARI BBRR 03', dead_alive == 'A')%>%
  mutate(mean_dist = mean(distance),
         n_trees_m2 = 100/(mean_dist)^2)

struc_sample_summary<-struc_sample%>%group_by(site_ID, species, n_trees_m2)%>%
  summarise(n_trees = length(species),
            n_hollows = sum(hollow_num),
            total_trees = nrow(struc_sample))%>%
  mutate(n_quarter = n_trees/total_trees,
         n_ha = (n_quarter * n_trees_m2)*100,
         hollows_quarter = n_hollows/total_trees,
         hollows_ha = (hollows_quarter*n_trees_m2)*100)

#write.csv(struc_sample%>%select(1:5, 9:11, 18), 'outputs/PCQM_example-Lutz.csv', row.names = F)

#using distance to furthest tree to get expansion factor:
#http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Distance_based_plots

struc_sample_EFI<-struc_sample%>%select(1:5, 9:11, 18)%>%
  group_by(sampling_point)%>%
  mutate(dist = max(distance),
         EFI = 10000/dist^2)

#trial with bitterlich, see if similar no of hollows

struc_sample_bitt<-struc_sample%>%select(1:5, 9:11, 18)%>%
  mutate(angle = ((dbh/100)/2)/distance,
         k = 10000*(angle^2),
         c = 50/sqrt(k),
         FI = pi*c^2*(dbh/100)^2,
         EFI = 10000/FI,
         hollows_ha = hollow_num*EFI)

#trial with bitterlich but fixed c (furthest tree), see if similar no of hollows


struc_sample_bit_mix<-struc_sample%>%select(1:5, 9:11, 18)%>%
  group_by(sampling_point)%>%
  filter(distance == max(distance))%>%
  mutate(angle = ((dbh/100)/2)/distance,
         k = 10000*(angle^2),
         c = 50/sqrt(k))

struc_sample_bit_mix_cal<-left_join(struc_sample%>%select(1:5, 9:11, 18), 
                                    struc_sample_bit_mix%>%select(sampling_point, c), 
                                    by = 'sampling_point')%>%
  mutate(FI = pi*c^2*(dbh/100)^2,
         EFI = 10000/FI,
         hollows_ha = hollow_num*EFI)

struc_sample_bit_mix_cal%>%group_by(sampling_point)%>%
  summarise(meanh = mean(hollows_ha),
            medianh = median(hollows_ha),
            maxh = max(hollows_ha),
            minh = min(hollows_ha),
            sd = sd(hollows_ha))

#calc for all trees

struc_fact<-struc%>%select(1:5, 9:11, 18)%>%
  group_by(site_ID, sampling_point)%>%
  filter(distance == max(distance))%>%
  mutate(angle = ((dbh/100)/2)/distance,
         k = 10000*(angle^2),
         c_fact = 50/sqrt(k))


struc_ha<-left_join(struc,
                      struc_fact%>%select(site_ID, sampling_point, c_fact),
                      by = c('site_ID','sampling_point'))%>%
  mutate(FI = pi*c_fact^2*(dbh/100)^2,
         EFI = 10000/FI,
         hollows_ha = hollow_num*EFI,
         meancrown = (crown_NS+crown_EW)/2,
         crown_area = crown_EW*crown_NS,
         ba = (dbh/200)^2*3.142,
         volume = ba*(height/3))

write.csv(struc_ha, 'outputs/structure_fixed_perha.csv', row.names = F)


# by plot and distance ----------------------------------------------------

struc_PCQM<-struc%>%filter(dead_alive == 'A')%>%
  mutate(hollow_num = case_when(hollows == 'absent' ~ 0,
                                TRUE ~ hollow_num),
         hbt = case_when(hollow_num>0 ~ 1,
                         TRUE ~ 0),
         distance = case_when(distance>300 ~ distance/10,
                              TRUE ~ distance))%>%
  group_by(site_ID, sampling_point, firesev)%>%
  summarise(av_d_sq = mean(distance)^2,
            density = 10000/av_d_sq,
            av_hollow = sum(hollow_num)/4,
            hol_ha = av_hollow*av_d_sq,
            no_hbt = sum (hbt),
            av_hbt = sum(hbt)/4,
            hbt_ha = av_hbt*av_d_sq)

struc_PCQM_av<-struc_PCQM%>%na.omit()%>%
  group_by(site_ID, firesev)%>%
  summarise(av_dens = mean(density),
            av_hol = mean(hol_ha),
            av_hbt = mean(hbt_ha),
            prop_hbt = (sum(no_hbt)/40)*100)

ggplot(struc_PCQM_av, aes(x = as.factor(firesev), y = av_hbt))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  geom_point(aes(color = av_hbt), size = 2)+
  scale_color_gradient2(midpoint = mean(struc_PCQM_av$av_hbt), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(struc_PCQM_av, aes(x = as.factor(firesev), y = av_hol))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  geom_point(aes(color = av_hol), size = 2)+
  scale_color_gradient2(midpoint = mean(struc_PCQM_av$av_hol), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(struc_PCQM%>%na.omit(), aes(x = as.factor(firesev), y = hol_ha))+
  geom_boxplot(width = 0.5)+
  geom_point(aes(color = hol_ha), size = 2)+
  scale_color_gradient2(midpoint = mean(struc_PCQM$hol_ha), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

# significance tests --------------------------------------------------------------

hollows_tally<-struc_ha%>%filter(grepl("Euc",species))%>%
  mutate(hollows = case_when(hollows == 'present' ~ 1,
                             TRUE ~ 0))%>%
  group_by(plot, gg_pres, firesev)%>%
  summarise('trees measured' = length(species),
            'No. of species' = n_distinct(species),
            'trees/ha' = sum(EFI),
            'hollows_measured' = sum(hollow_num),
            'hollowspha' = sum(hollows_ha),
            'HBTs' = sum(EFI[hollows == 1]),
            'meanhol' = mean(hollows_ha),
            'meanhbt' = mean(EFI[hollows == 1]))%>%
  mutate(gg_pres = case_when(is.na(gg_pres) ~ 1,
                             TRUE ~ as.numeric(gg_pres)))


ggplot(hollows_tally, aes(x = as.factor(firesev), y = hollows_measured))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  geom_point(aes(color = hollows_measured), size = 2)+
  scale_color_gradient2(midpoint = mean(hollows_tally$hollows_measured), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(hollows_tally%>%filter(HBTs<3000), aes(x = as.factor(firesev), y = HBTs))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  geom_point(aes(color = HBTs), size = 2)+
  scale_color_gradient2(midpoint = mean(na.omit(hollows_tally$HBTs)), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

#combine 3+4 and 5+6 to medium and high sev?

ggplot(hollows_tally%>%filter(meanhol<400), aes(x = as.factor(firesev), y = meanhol))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  geom_point(aes(color = meanhol), size = 2)+
  scale_color_gradient2(midpoint = mean(na.omit(hollows_tally$meanhol)), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(hollows_tally, aes(x = as.factor(firesev), y = meanhbt, fill = as.factor(firesev)))+
  geom_boxplot(outlier.size = -1, width = 0.5, show.legend = F)+
  geom_point(aes(color = meanhbt), size = 2, show.legend = F)+
  scale_color_gradient2(midpoint = mean(na.omit(hollows_tally$meanhbt)), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )+
  labs(title=NULL,x = NULL, y = 'HBTs/ha', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

#by gg presence

ggplot(hollows_tally, aes(x = as.factor(firesev), y = hollows_measured))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  facet_wrap(~gg_pres)+
  geom_point(aes(color = hollows_measured), size = 2)+
  scale_color_gradient2(midpoint = mean(hollows_tally$hollows_measured), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(hollows_tally, aes(x = as.factor(firesev), y = meanhbt))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  facet_wrap(~gg_pres)+
  geom_point(aes(color = meanhbt), size = 2)+
  scale_color_gradient2(midpoint = mean(na.omit(hollows_tally$meanhbt)), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(hollows_tally%>%filter(meanhol<400), aes(x = as.factor(firesev), y = meanhol))+
  geom_boxplot(outlier.size = -1, width = 0.5)+
  #facet_wrap(~gg_pres)+
  geom_point(aes(color = meanhol), size = 2)+
  scale_color_gradient2(midpoint = mean(na.omit(hollows_tally$meanhol)), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )


na.omit(hollows_tally)%>%group_by(firesev)%>%summarise(mean = mean(hollows_measured),
                                                       max = max (hollows_measured),
                                                       min = min(hollows_measured),
                                                       sd = sd(hollows_measured))


holstats1<-na.omit(hollows_tally)%>%group_by(firesev)%>%summarise(mean = mean(meanhbt),
                                               max = max (meanhbt),
                                               min = min(meanhbt),
                                               sd = sd(meanhbt))

holstats2<-na.omit(hollows_tally%>%filter(meanhol<400),)%>%group_by(firesev)%>%summarise(mean = mean(meanhol),
                                                       max = max (meanhol),
                                                       min = min(meanhol),
                                                       sd = sd(meanhol))


summary(aov(hollows_measured~as.factor(firesev), data = hollows_tally))

TukeyHSD(aov(hollows_measured~as.factor(firesev), data = hollows_tally))

summary(aov(log(hollowspha)~as.factor(firesev), data = na.omit(hollows_tally)))

TukeyHSD(aov(log(hollowspha)~as.factor(firesev), data = na.omit(hollows_tally)))

summary(aov(log(HBTs)~as.factor(firesev), data = na.omit(hollows_tally)))

TukeyHSD(aov(log(HBTs)~as.factor(firesev), data = na.omit(hollows_tally)))

summary(aov(log(meanhol)~as.factor(firesev), data = na.omit(hollows_tally))) #used

TukeyHSD(aov(log(meanhol)~as.factor(firesev), data = na.omit(hollows_tally)))

summary(aov(log(meanhbt)~as.factor(firesev), data = na.omit(hollows_tally))) #used

TukeyHSD(aov(log(meanhbt)~as.factor(firesev), data = na.omit(hollows_tally)))

summary(aov(log(meanhol)~gg_pres, data = na.omit(hollows_tally))) #used

summary(aov(log(meanhbt)~gg_pres, data = na.omit(hollows_tally))) #used

summary(aov(hollows_measured~gg_pres, data = hollows_tally))

#  comparison between pre and post fire hollow occurrence

hol_pre<-read.csv('data/pre_fire/structure_extended_euc-only.csv', stringsAsFactors = T)%>%
  group_by(plot, GG, transect)%>%
  summarise('trees measured' = length(species),
            'No. of species' = n_distinct(species),
            'trees/ha' = sum(EFi),
            'hollows_measured' = sum(hollownumber),
            'hollowspha' = sum(hollows_ha),
            'HBTs' = sum(EFi[hollowpres == 1]),
            'meanhol' = mean(hollows_ha),
            'meanhbt' = mean(EFi[hollowpres == 1]),
            'time' = 'pre')%>%
  mutate(transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))%>%
  select(1, gg_pres = 2, 8, 9, 12)

hollows_tally_comp<-hollows_tally%>%select(1,2,hollowspha = 10, HBTs = 11)%>%
  mutate(time = 'post')%>%na.omit()

hol_comp<-rbind(hollows_tally_comp, hol_pre)%>%mutate(time = as.factor(time))%>%
  mutate(time = fct_relevel(time, 'pre', 'post'))

#AOVs

summary(aov(log(hollowspha)~time, data = hol_comp))

TukeyHSD(aov(log(hollowspha)~time, data = hol_comp))

summary(aov(log(HBTs)~time, data = hol_comp))

TukeyHSD(aov(log(HBTs)~time, data = hol_comp))

#Plots

hol_comp_prepost<-hol_comp%>%filter(!grepl('ARI', plot))%>%filter(!grepl('G3', plot))

ggplot(hol_comp, aes(x = time, y = hollowspha))+
  geom_boxplot(width = 0.5)+
  geom_point(aes(color = hollowspha), size = 2)+
  ylim(0,1200)+
  scale_color_gradient2(midpoint = mean(hol_comp$hollowspha), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

ggplot(hol_comp_prepost, aes(x = time, y = HBTs))+
  geom_boxplot(width = 0.5)+
  geom_point(aes(color = HBTs), size = 2)+
  ylim(0,150)+
  scale_color_gradient2(midpoint = mean(hol_comp_prepost$HBTs), 
                        low = "red", mid = 'green', high = "blue", space = "Lab" )

# descriptive stats/graphs - structure ------------------------------------

struc<-structure_join_fix_final%>%mutate(firesev = case_when(firesev == 0 ~ 2,
                                                             TRUE ~ firesev))

struc_pre<-read.csv('data/pre_fire/structure_extended_euc-only.csv', stringsAsFactors = T)%>%
  mutate(transect = fct_relevel(transect,'lowlands','mid-hills','high-country'))
  


#Fire Sev Classes:
#Canopy burnt (Class 6)- CB (> 20% canopy foliage consumed)
#High canopy scorch (5) - HCS (>80% of canopy foliage is scorched)
#Medium canopy scorch (4) - MCS (Canopy is a mosaic of both unburnt and scorched foliage, 20 - 80%)
#Low canopy scorch (3) - LCS (Canopy foliage is largely unaffected (<20% scorched), but the understorey has been burnt)
#Unburnt (2) - UB (Canopy and understorey foliage are largely (>90%) unburnt)

#dbh graphs

#1. Euc only

struc%>%filter(grepl("Euc",species))%>%
ggplot(aes(dbh, fill = as.factor(firesev)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~as.factor(firesev), labeller = as_labeller(c('2' = 'Unburnt', 
                                                           '3' = 'Low canopy scorch',
                                                           '4' = 'Medium canopy scorch',
                                                           '5' = 'High canopy scorch',
                                                           '6' = 'Canopy burnt')))+
  labs(title=NULL,x = 'DBH', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

#2. Euc HBTs only

holfig_post_1<-struc%>%filter(grepl("Euc",species), hollows == 'present')%>%
  ggplot(aes(dbh, fill = as.factor(firesev)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~as.factor(firesev), ncol = 2, labeller = as_labeller(c('2' = 'Unburnt', 
                                                           '3' = 'Low canopy scorch',
                                                           '4' = 'Medium canopy scorch',
                                                           '5' = 'High canopy scorch',
                                                           '6' = 'Canopy burnt')))+
  labs(title=NULL,x = 'DBH (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

#stats

#dbh

struc%>%filter(grepl("Euc",species), hollows == 'present', dead_alive == 'A')%>%
  group_by(firesev)%>%
  summarise(mean = mean(dbh),
            max = max (dbh),
            min = min(dbh),
            sd = sd(dbh))

summary(aov(dbh~as.factor(firesev), data = struc%>%filter(grepl("Euc",species), hollows == 'present')))
TukeyHSD(aov(dbh~as.factor(firesev), data = struc%>%filter(grepl("Euc",species), hollows == 'present')))

summary(aov(dbh~as.factor(gg_pres), data = struc%>%filter(grepl("Euc",species), hollows == 'present')))

#height

struc_ha%>%filter(grepl("Euc",species), hollows == 'present', dead_alive == 'A', height>0)%>%
  group_by(firesev)%>%
  summarise(mean = mean(height),
            max = max (height),
            min = min(height),
            sd = sd(height))

summary(aov(height~as.factor(firesev), data = struc_ha%>%filter(grepl("Euc",species), 
                                                                hollows == 'present',
                                                                dead_alive == 'A',
                                                                height>0)))

TukeyHSD(aov(dbh~as.factor(firesev), data = struc_ha%>%filter(grepl("Euc",species), 
                                                              hollows == 'present',
                                                              dead_alive == 'A',
                                                              height>0)))

summary(aov(height~as.factor(gg_pres), data = struc_ha%>%filter(grepl("Euc",species), 
                                                                hollows == 'present',
                                                                dead_alive == 'A',
                                                                height>0)))

#crown area

struc_ha%>%filter(grepl("Euc",species), hollows == 'present', dead_alive == 'A', crown_area>0)%>%
  group_by(firesev)%>%
  summarise(mean = mean(crown_area),
            max = max (crown_area),
            min = min(crown_area),
            sd = sd(crown_area))

summary(aov(crown_area~as.factor(firesev), data = struc_ha%>%filter(grepl("Euc",species), 
                                                                hollows == 'present',
                                                                dead_alive == 'A',
                                                                height>0)))

summary(aov(crown_area~as.factor(gg_pres), data = struc_ha%>%filter(grepl("Euc",species), 
                                                                    hollows == 'present',
                                                                    dead_alive == 'A',
                                                                    height>0)))

TukeyHSD(aov(crown_area~as.factor(firesev), data = struc_ha%>%filter(grepl("Euc",species), 
                                                              hollows == 'present',
                                                              dead_alive == 'A',
                                                              height>0)))


#2. Euc HBTs only - pre fire

#by transect

holfig_pre_1<-struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')%>%
ggplot(aes(dbh, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~transect, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                 'mid-hills' = 'Mid-hills', 
                                                 'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'DBH (cm)', y = 'Density', fill = 'Transect')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect",labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"),
        plot.margin = margin(0, 15, 0, 0))

#by gg

holfig_pre_2<-struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')%>%
  ggplot(aes(dbh,  fill = as.factor(GG)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~GG, labeller = as_labeller(c('0' = 'Not detected', 
                                                 '1' = 'Detected')))+
  labs(title=NULL,x = 'DBH (cm)', y = 'Density', fill = 'Transect')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect",labels=c("Not detected", "Detected"))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"),
        plot.margin = margin(0, 15, 0, 0))



struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')%>%
  group_by(transect)%>%
  summarise(mean = mean(dbh),
            max = max (dbh),
            min = min(dbh),
            sd = sd(dbh))


summary(aov(dbh~transect, data = struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')))
TukeyHSD(aov(dbh~transect, data = struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')))


summary(aov(dbh~GG, data = struc_pre))

summary(aov(dbh~GG, data = struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')))
TukeyHSD(aov(dbh~as.factor(GG), data = struc_pre%>%filter(grepl("Euc",species), hollows == 'Present')))


# structure - figures -----------------------------------------------------

holpatch<-(holfig_pre_1 / holfig_pre_2) | holfig_post_1 
  
holpatch + plot_annotation(tag_levels = 'A')

ggsave('hol_dbh.svg',path = 'figures/', width = 25, height = 15.6, units = 'cm', dpi = 600)


# soil data by firesev ----------------------------------------------------

sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  ggplot(aes(y = Total.Nitrogen, x = as.factor(firesev), fill = as.factor(firesev)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'Total Nitrogen (%)', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  ggplot(aes(y = Phosphorus.Colwell, x = as.factor(firesev), fill = as.factor(firesev)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'P levels', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  ggplot(aes(y = Potassium.Colwell, x = as.factor(firesev), fill = as.factor(firesev)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'K levels', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  ggplot(aes(y = pH.Level..CaCl2., x = as.factor(firesev), fill = as.factor(firesev)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'pHs', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

#pre fire soil data

sites_pre<-read.csv('data/pre_fire/plot_data.csv')%>%
  mutate(transect = case_when(transect == 'T3' ~ 'high-country',
                              transect == 'T2' ~ 'mid-hills',
                              transect == 'T1' ~ 'lowlands'))

summary(aov(soil_total.nitrogen~as.factor(transect), data = sites_pre))

sites_pre%>%
  mutate(transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))%>%
  ggplot(aes(y = soil_total.nitrogen, x = as.factor(transect), fill = as.factor(transect)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'Nitrogen (%)', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect",labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


sites_pre%>%
  mutate(transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))%>%
  ggplot(aes(y = as.numeric(soil_phosphorus.colwell), x = as.factor(transect), fill = as.factor(transect)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'P', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect",labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


sites_pre%>%
  mutate(transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))%>%
  ggplot(aes(y = as.numeric(soil_potassium.colwell), x = as.factor(transect), fill = as.factor(transect)))+
  geom_boxplot(width = 0.5, show.legend = F)+
  labs(title=NULL,x = NULL, y = 'K', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect",labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


#stats

sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  st_drop_geometry()%>%
  group_by(firesev)%>%
  summarise(mean = mean(Total.Nitrogen),
            max = max (Total.Nitrogen),
            min = min(Total.Nitrogen),
            sd = sd(Total.Nitrogen))

summary(aov(Total.Nitrogen~as.factor(firesev), data = sites_spatial_data%>%
              mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))))

TukeyHSD(aov(Total.Nitrogen~as.factor(firesev), data = sites_spatial_data%>%
               mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))))


sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  st_drop_geometry()%>%
  #group_by(firesev)%>%
  summarise(mean = mean(Total.Carbon),
            max = max (Total.Carbon),
            min = min(Total.Carbon),
            sd = sd(Total.Carbon))

sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  st_drop_geometry()%>%
  #group_by(firesev)%>%
  summarise(mean = mean(Phosphorus.Colwell),
            max = max (Phosphorus.Colwell),
            min = min(Phosphorus.Colwell),
            sd = sd(Phosphorus.Colwell))


sites_spatial_data%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,TRUE ~ firesev))%>%
  st_drop_geometry()%>%
  #group_by(firesev)%>%
  summarise(mean = mean(Potassium.Colwell),
            max = max (Potassium.Colwell),
            min = min(Potassium.Colwell),
            sd = sd(Potassium.Colwell))


# descriptive -------------------------------------------------------------

sites_spatial_data_new<-sites_spatial_data%>%filter(!(grepl("T",site_ID)))
sites_spatial_data_old<-sites_spatial_data%>%filter((grepl("T",site_ID)))


struc_post_euc<-struc%>%filter(grepl("Euc",species))
