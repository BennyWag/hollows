library(car)
library(MuMIn)
library(effects)
library(effectsize)
library(broom)
library(gt)
library(gtsummary)
library(sf)
library(lme4)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

# combined analyses - hollows/ha and structure ----------------------------

## functions

dredge_fix<-function(model, ylim = c(0, 1)){
  
  options(na.action = 'na.fail')
  dredger<-dredge(model, trace = 3,  extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  options(na.action = "na.omit")
  
  best<-get.models(dredger, 1)[[1]]
  
  list<-list(dredger, best)
  names(list)<-c('dredge', 'best_model')
  
  print(summary(best))
  plot(allEffects(mod = best), type = "response", ylim = ylim)
  
  return(list)
}

estimate_graph<-function(model){
  
  results <- tidy(model)
  
  fit_cis_95 <- confint(model, level = 0.95) %>% 
    data.frame() %>%
    rename("conf.low_95" = "X2.5..",
           "conf.high_95" = "X97.5..")
  fit_cis_90 <- confint(model, level = 0.90) %>% 
    data.frame() %>%
    rename("conf.low_90" = "X5..",
           "conf.high_90" = "X95..")
  
  results <- bind_cols(results, 
                       fit_cis_95, 
                       fit_cis_90) %>%
    rename(Variable = term,
           Coefficient = estimate,
           SE = std.error) %>%
    filter(Variable != "(Intercept)")
  
  results <- results %>% dplyr::select(-SE, 
                                       -statistic,
                                       -p.value)
  plot<-ggplot(results, 
               aes(x = Variable, y = Coefficient)) +
    geom_hline(yintercept = 0, 
               colour = gray(1/2), lty = 2) +
    geom_point(aes(x = Variable, 
                   y = Coefficient)) + 
    geom_linerange(aes(x = Variable, 
                       ymin = conf.low_90,
                       ymax = conf.high_90),
                   lwd = 1) +
    geom_linerange(aes(x = Variable, 
                       ymin = conf.low_95,
                       ymax = conf.high_95),
                   lwd = 1/2) + 
    coord_flip()+
    theme_bw()+
    theme(plot.title = element_text(size =18, face='bold'),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
          axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(
            size = 12, color = "black", face = "bold.italic"))
  
  plot
  
  list<-list(results,plot)
  names(list)<-c('table', 'plot')
  
  return(list)
  
}

## data setup/clean

post<-read.csv('outputs/structure_fixed.csv')%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,
                             TRUE ~ firesev),
         sampling_point = case_when(sampling_point>10 ~ 3,
                                    TRUE ~ as.numeric(sampling_point)),
         transect = case_when(elevation < 420 ~ 'lowlands',
                              elevation > 420 & elevation < 840 ~ 'mid-hills',
                              elevation >= 840 ~ 'high-country'),
         transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))%>%
  mutate(firesev = case_when(site_ID == 'T35P3' ~ 5,
                             site_ID == 'T1P1' ~ 3,
                             site_ID == 'T2P1' ~ 3,
                             site_ID == 'ARI BBRR 26' ~ 3,
                             site_ID == 'ARI BBRR 36' ~ 3,
                             site_ID == 'ARI BBRR 22' ~ 4,
                             site_ID == 'ARI BBRR 45' ~ 4,
                             site_ID == 'T1P2' ~ 6,
                             TRUE ~ firesev))

#fix data for model

post_model<-post%>%
  mutate(sev_simp = as.factor(case_when(firesev > 2 & firesev <=4 ~ 'low',
                                        firesev > 4 ~ 'high',
                                        TRUE ~ 'unburnt')),
         sev_simp = fct_relevel(sev_simp, 'unburnt', 'low', 'high'))%>%
  filter(grepl("Euc",species))%>%#filter(dead_alive == 'A')%>%
  mutate(hollows = case_when(hollows == 'present' ~ 1,
                             TRUE ~ 0),
         firesev_fact = as.factor(firesev))%>%mutate(layers_num = case_when(layers == '1' ~ 1,
                                                                            layers == '2' ~ 2,
                                                                            layers == '3' ~ 3),
                                                     layer_num = case_when(layer == 'canopy' ~ 1,
                                                                           layer == 'subcanopy' ~ 2,
                                                                           layer == 'midstorey' ~ 3),
                                                     ageclass_num = case_when(age_class == 'regrowth' ~ 1,
                                                                              age_class == 'mature' ~ 2,
                                                                              age_class == 'over-mature' ~ 3,
                                                                              age_class == 'late-mature' ~ 4,
                                                                              age_class == 'stag' ~ 5),
                                                     soilcode_1 = as.numeric(ASC_O),
                                                     soilcode_2 = as.numeric(ASC_SO),
                                                     hollow_num = case_when(hollow_num > 9 ~ 1,
                                                                            TRUE ~ as.numeric(hollow_num)),
                                                     hollow_num = case_when(hollows == 'absent' ~ 0,
                                                                            TRUE ~ hollow_num),
                                                     crown_area = crown_EW*crown_NS,
                                                     ba = (dbh/200)^2*3.142,
                                                     volume = ba*(height/3))

## site stats

post_sites<-st_read('data/post_fire/spatial/sites_all.gpkg')%>%
  dplyr::mutate(site_ID = str_replace(site_ID, 'BBBR\\b', 'BBRR'))%>%st_drop_geometry()%>%
  mutate(firesev = case_when(firesev == 0 ~ 2,
                             TRUE ~ firesev))%>%
  mutate(firesev = case_when(site_ID == 'T35P3' ~ 5,
                             site_ID == 'T1P1' ~ 3,
                             site_ID == 'T2P1' ~ 3,
                             site_ID == 'ARI BBRR 26' ~ 3,
                             site_ID == 'ARI BBRR 36' ~ 3,
                             site_ID == 'ARI BBRR 22' ~ 4,
                             site_ID == 'ARI BBRR 45' ~ 4,
                             site_ID == 'T1P2' ~ 6,
                             TRUE ~ firesev))%>%
  slice(1:40)%>%
  select(site_ID, firesev, gg_pres)%>%
  filter(site_ID %in% unique(pre$plot))

table(post_sites$firesev)
table(post_sites$gg_pres, post_sites$firesev)

# pre<-read.csv('data/pre_fire/structure_extended_euc-only_modeling-pres_abs.csv', stringsAsFactors = T)%>%
#   mutate(transect = ifelse(grepl('T3', plot, ignore.case = T), 'high-country',
#                            ifelse(grepl('T2', plot, ignore.case = T), 'mid-hills', 'lowlands')),
#          transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'),
#          Fi_re = pi*(((dbh)/4)^2),
#          EFi_re = 10000/Fi_re,
#          hollows_ha_re = hollownumber*EFi_re)

pre<-read.csv('data/pre_fire/tree_data.csv', stringsAsFactors = T)%>%
  mutate(plot = as.factor(str_remove(plot, '[.]')))%>%
  mutate(meancrown = (crownNS+crownEW)/2,
         crown_area = crownEW*crownNS,
         ba = (dbh/200)^2*3.142,
         Fi = pi*25^2*(dbh/100)^2,
         EFi = 10000/Fi,
         hollows_ha = hollownumber*EFi)

#add gg presabs

gg_pre<-read.csv('data/pre_fire/gg_pres_abs.csv', stringsAsFactors = T)%>%
  mutate(plot = as.factor(str_remove(Plot, '[.]')))%>%select(4,2,3)

pre_gg<-left_join(pre, gg_pre, by = 'plot')

#add plot-level data

pre_plots<-read.csv('data/pre_fire/plot_data.csv', stringsAsFactors = T)%>%
  select(3:9,14,15,18,19,21,22,25,32,33,39:41,46:48,52,64,65)

pre<-left_join(pre_gg, pre_plots, by = 'plot')%>%
  mutate(soil_total.nitrogen = ifelse(soil_total.nitrogen == '<1', 0.9, soil_total.nitrogen),
         soil_phosphorus.colwell = ifelse(soil_phosphorus.colwell == '<2', 2, soil_phosphorus.colwell))%>%
  mutate(obs = row_number(),
            hollowpres = as.factor(hollowpres),
            health_cat_fct = fct_relevel(health, 
                                         'Alive', 
                                         'Senescing', 
                                         'Dead'),
            health_cat = as.numeric(health_cat_fct),
            health_cat_factor = as.factor(as.numeric(health_cat)),
            layers_num = case_when(layers == '>3' ~ 4,
                                   layers == '1' ~ 1,
                                   layers == '2' ~ 2,
                                   layers == '3' ~ 3),
            layer_num = case_when(layer == 'Canopy' ~ 1,
                                  layer == 'Subcanopy' ~ 2,
                                  layer == 'Midstorey' ~ 3,
                                  layer == 'Understorey ' ~ 4),
            ageclass_num = case_when(ageclass == 'Regrowth' ~ 1,
                                     ageclass == 'Mature ' ~ 2,
                                     ageclass == 'Over-mature' ~ 3,
                                     ageclass == 'Late mature' ~ 4,
                                     ageclass == 'Stag' ~ 5),
         Fi_re = pi*(((dbh)/4)^2),
         EFi_re = 10000/Fi_re,
         hollows_ha_re = hollownumber*EFi_re)%>%
  mutate(transect = ifelse(grepl('T3', plot, ignore.case = T), 'high-country',
                           ifelse(grepl('T2', plot, ignore.case = T), 'mid-hills', 'lowlands')),
         transect = fct_relevel(transect, 'lowlands', 'mid-hills', 'high-country'))

pre_nopot<-pre%>%mutate(hollownumber_fix = case_when(hollowsfix == 'Potential ' ~ 0,
                                                 TRUE ~ hollownumber),
                        hollows_ha_fix = case_when(hollowsfix == 'Potential ' ~ 0,
                                          TRUE ~ hollows_ha),
                        hollowpres_fix = case_when(hollowpres == 2 ~ '0',
                                               TRUE ~ hollowpres))

write.csv(pre_nopot, 'outputs/pre_fixed.csv', row.names = F)

# get prop hollows and hollow/site/ha -------------------------------------

## post-fire

post_holha<-post%>%filter(dead_alive == 'A')%>%
  mutate(hollow_num = case_when(hollows == 'absent' ~ 0,
                                TRUE ~ hollow_num),
         hbt = case_when(hollow_num>0 ~ 1,
                         TRUE ~ 0),
         distance = case_when(distance>300 ~ distance/10,
                              TRUE ~ distance))%>%
  group_by(site_ID, gg_pres, sampling_point, firesev,transect)%>%
  summarise(av_d_sq = mean(distance)^2,
            density = 10000/av_d_sq,
            av_hollow = sum(hollow_num)/4,
            hol_ha = av_hollow*av_d_sq,
            no_hbt = sum(hbt),
            no_lhbt = sum(hbt[dbh>=100]),
            av_hbt = sum(hbt)/4,
            hbt_ha = av_hbt*av_d_sq)

ggplot(post_holha%>%na.omit(), aes(x = as.factor(firesev), y = hol_ha))+
  geom_boxplot(width = 0.5)+
  ylim(0, 350)

post_holha_av<-post_holha%>%na.omit()%>%#filter(hol_ha>0)%>%
  group_by(site_ID, gg_pres, firesev, transect)%>%
  summarise(av_dens = mean(density),
            av_hol = mean(hol_ha),
            av_hbt = mean(hbt_ha),
            prop_hbt = (sum(no_hbt)/40)*100,
            prop_lhbt = (sum(no_lhbt)/40)*100,
            prop_lhbt_hbt = (sum(no_lhbt)/sum(no_hbt))*100)%>%
  mutate(sev_simp = as.factor(case_when(firesev > 2 & firesev <=4 ~ 'low',
                              firesev > 4 ~ 'high',
                              TRUE ~ 'unburnt')),
         sev_simp = fct_relevel(sev_simp, 'unburnt', 'low', 'high'))

ggplot(post_holha_av%>%na.omit(), aes(x = as.factor(firesev), y = av_hol))+
  geom_boxplot(width = 0.5)

#not for each point but per transect

# post_holha_transect<-post%>%filter(dead_alive == 'A')%>%
#   mutate(hollow_num = case_when(hollows == 'absent' ~ 0,
#                                 TRUE ~ hollow_num),
#          hbt = case_when(hollow_num>0 ~ 1,
#                          TRUE ~ 0),
#          distance = case_when(distance>300 ~ distance/10,
#                               TRUE ~ distance))%>%
#   group_by(site_ID, firesev,transect)%>%
#   summarise(av_d_sq = mean(distance)^2,
#             density = 10000/av_d_sq,
#             av_hollow = sum(hollow_num)/40,
#             hol_ha = av_hollow*av_d_sq,
#             no_hbt = sum(hbt),
#             av_hbt = sum(hbt)/40,
#             hbt_ha = av_hbt*av_d_sq)

#check potential bias

plot(post$firesev, post$elevation)

table(post_holha_av$transect, post_holha_av$firesev)

## pre-fire

#bitterlich calc

# pre_holha<-pre%>%group_by(plot, GG, transect)%>%
#   summarise('trees measured' = length(species),
#             'No. of species' = n_distinct(species),
#             'trees/ha' = sum(EFi),
#             'hollows_measured' = sum(hollownumber),
#             'hollows' = sum(hollows_ha),
#             'HBTs' = sum(EFi[hollowpres == 1]))
# 
# #Craig's calc - same results
# 
# pre_holha_re<-pre%>%group_by(plot, GG, transect)%>%
#   summarise('trees measured' = length(species),
#             'No. of species' = n_distinct(species),
#             'trees/ha' = sum(EFi_re),
#             'hollows_measured' = sum(hollownumber),
#             'hollows' = sum(hollows_ha_re),
#             'HBTs' = sum(EFi_re[hollowpres == 1]))

#with potential HBTs converted to absences

pre_holha_nopot<-pre_nopot%>%filter(!health %in% c('Dead'))%>%
  group_by(plot, GG, transect)%>%
  summarise(n_trees = length(species),
            n_hbt = length(species[hollowpres_fix == 1]),
            n_lhbt = length(species[hollowpres_fix == 1 & dbh>=100]),
            trees_ha = sum(EFi),
            hollows_ha = sum(hollows_ha_fix),
            hbt_ha = sum(EFi[hollowpres_fix == 1]),
            prop_hbt = (n_hbt/n_trees)*100,
            prop_lhbt = (n_lhbt/n_trees)*100,
            prop_lhbt_hbt = (n_lhbt/n_hbt)*100)

ggplot(pre_holha_nopot, aes(x = as.factor(transect), y = hollows_ha))+
  geom_boxplot(width = 0.5)+
  ylim(0,350)

#match above by using averages?

# pre_holha_av<-pre%>%
#   group_by(plot, transect)%>%
#   summarise(#av_dens = mean(EFi),
#             #av_hol = mean(hollows_ha),
#             #av_hbt = mean(EFi[hollowpres == 1]),
#             n_trees = length(species),
#             n_hbt = length(species[hollowpres == 1]),
#             prop_hbt = (n_hbt/n_trees)*100)

## significance tests

#pre

summary(aov(prop_hbt~as.factor(transect), data = pre_holha_nopot))
TukeyHSD(aov(prop_hbt~as.factor(transect), data = pre_holha_nopot))

#post

summary(aov(prop_hbt~as.factor(transect), data = post_holha_av))
TukeyHSD(aov(prop_hbt~as.factor(transect), data = post_holha_av))

summary(aov(prop_hbt~as.factor(firesev), data = post_holha_av))
TukeyHSD(aov(prop_hbt~as.factor(firesev), data = post_holha_av))

# summary table on hollow abundance/occurrence/prop -----------------------

post_hbtprob<-post_holha_av%>%filter(!sev_simp %in% c('unburnt'))%>%
  group_by(transect, sev_simp)%>%
  summarise(mean_prophbt = mean(prop_hbt),
            sd_probhbt = sd(prop_hbt),
            mean_proplhbt = mean(prop_lhbt),
            sd_problhbt = sd(prop_lhbt),
            mean_proplhbt_hbt = mean(prop_lhbt_hbt),
            sd_problhbt_hbt = sd(prop_lhbt_hbt))

pre_hbtprob<-pre_holha_nopot%>%
  mutate(sev_simp = 'pre-fire')%>%
  group_by(transect, sev_simp)%>%
  summarise(mean_prophbt = mean(prop_hbt),
            sd_probhbt = sd(prop_hbt),
            mean_proplhbt = mean(prop_lhbt),
            sd_problhbt = sd(prop_lhbt),
            mean_proplhbt_hbt = mean(prop_lhbt_hbt),
            sd_problhbt_hbt = sd(prop_lhbt_hbt))

pre_post_hbtprob<-rbind(post_hbtprob, pre_hbtprob)

write.csv(pre_post_hbtprob, 'outputs/pre-post_hbtprob.csv', row.names = F)


# size characteristics HBTs pre/post --------------------------------------

## dbh

pre_nopot%>%filter(hollowpres_fix == 1)%>%
  group_by(transect)%>%
  summarise(mean = mean(dbh),
            max = max (dbh),
            min = min(dbh),
            sd = sd(dbh))

post_model%>%filter(hollows == 1, firesev > 2)%>%
  group_by(transect)%>%
  summarise(mean = mean(dbh),
            max = max (dbh),
            min = min(dbh),
            sd = sd(dbh))

#sign. test

summary(aov(log(dbh)~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1)))
TukeyHSD(aov(log(dbh)~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1)))

summary(aov(log(dbh)~as.factor(transect), data = post_model%>%filter(hollows == 1, firesev > 2)))
TukeyHSD(aov(dbh~as.factor(transect), data = post_model%>%filter(hollows == 1, firesev > 2)))

summary(aov(log(dbh)~as.factor(firesev), data = post_model%>%filter(hollows == 1)))
TukeyHSD(aov(dbh~as.factor(firesev), data = post_model%>%filter(hollows == 1)))

## height

pre_nopot%>%filter(hollowpres_fix == 1)%>%
  #group_by(transect)%>%
  summarise(mean = mean(height),
            max = max (height),
            min = min(height),
            sd = sd(height))

post_model%>%filter(hollows == 1, height > 0, firesev > 2)%>%
  #group_by(firesev)%>%
  summarise(mean = mean(height),
            max = max (height),
            min = min(height),
            sd = sd(height))

#sign. test

summary(aov(log(height)~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1)))
TukeyHSD(aov(height~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1)))

summary(aov(log(height)~as.factor(transect), data = post_model%>%filter(hollows == 1, height > 0, firesev > 2)))
TukeyHSD(aov(height~as.factor(transect), data = post_model%>%filter(hollows == 1, height > 0, firesev > 2)))

summary(aov(log(height)~as.factor(firesev), data = post_model%>%filter(hollows == 1, height > 0)))
TukeyHSD(aov(log(height)~as.factor(firesev), data = post_model%>%filter(hollows == 1, height > 0)))

## crown

pre_nopot%>%filter(hollowpres_fix == 1, crown_area>0)%>%
  #group_by(transect)%>%
  summarise(mean = mean(crown_area),
            max = max (crown_area),
            min = min(crown_area),
            sd = sd(crown_area))

post_model%>%filter(hollows == 1, crown_area > 0, firesev > 2)%>%
  #group_by(firesev)%>%
  summarise(mean = mean(crown_area),
            max = max (crown_area),
            min = min(crown_area),
            sd = sd(crown_area))


#sign. test

summary(aov(log(crown_area)~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1, crown_area>0)))
TukeyHSD(aov(crown_area~as.factor(transect), data = pre_nopot%>%filter(hollowpres_fix == 1, crown_area>0)))

summary(aov(log(crown_area)~as.factor(transect), data = post_model%>%filter(hollows == 1, crown_area > 0, firesev > 2)))
TukeyHSD(aov(crown_area~as.factor(transect), data = post_model%>%filter(hollows == 1, crown_area > 0, firesev > 2)))

summary(aov(log(crown_area)~as.factor(firesev), data = post_model%>%filter(hollows == 1, crown_area > 0)))
TukeyHSD(aov(log(crown_area)~as.factor(firesev), data = post_model%>%filter(hollows == 1, crown_area > 0)))

# dbh dist hbts  ---------------------------------------------

## post fire - fire severity 

post_model%>%filter(hollows == 1)%>%
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

## post fire - transect

post_hbt_dbh<-post_model%>%filter(hollows == 1, firesev > 2)%>%
  ggplot(aes(dbh, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'HBT diameter (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect")+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.025), labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.8, "cm", data = NULL))

## pre fire - transect

pre_hbt_dbh<-pre_nopot%>%filter(hollowpres_fix == 1)%>%
  ggplot(aes(dbh, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = NULL, y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.025),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_hbt_dbh / post_hbt_dbh + plot_annotation(tag_levels = 'A')

ggsave('hol_dbh_dist.svg',path = 'figures/', width = 25, height = 15.6, units = 'cm', dpi = 600)

### using only burnt sites

unique(post_model$site_ID)

post_model_pre_only<-post_model%>%filter(grepl("T",site_ID), hollows == 1, firesev > 2)

unique(post_model_pre_only$site_ID)

pre_model_post_only<-pre_nopot%>%filter(plot %in% unique(post_model_pre_only$site_ID),  hollowpres_fix == 1)

unique(pre_model_post_only$plot)

## post fire - transect (pre burnt only)

post_hbt_compare<-post_model_pre_only%>%
  ggplot(aes(dbh, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'HBT diameter (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect")+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.025), labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.8, "cm", data = NULL))

## pre fire - transect (pre burnt only)

pre_compare<-pre_model_post_only%>%
  ggplot(aes(dbh, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = NULL, y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.025),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_compare / post_hbt_compare + plot_annotation(tag_levels = 'A')

ggsave('hol_dbh_dist_compare.svg',path = 'figures/', width = 25, height = 15.6, units = 'cm', dpi = 600)


### histograms

## post fire - transect

post_hbt_dbh_hist<-post_model%>%filter(hollows == 1, firesev > 2)%>%
  ggplot(aes(dbh, fill = transect))+
  geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  #geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'HBT diameter (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect")+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 50))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.8, "cm", data = NULL))

## pre fire - transect

pre_hbt_dbh_hist<-pre_nopot%>%filter(hollowpres_fix == 1)%>%
  ggplot(aes(dbh, fill = transect))+
  geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  #geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = NULL, y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 50))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_hbt_dbh_hist / post_hbt_dbh_hist + plot_annotation(tag_levels = 'A')


# no of hollow density plots ----------------------------------------------

## post fire - transect

post_hbt_hollows<-post_model%>%filter(hollows == 1, firesev > 2)%>%
  ggplot(aes(hollow_num, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'No. of hollows per HBT', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect")+
  #scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.6), labels = scales::percent_format(scale = 100))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.8, "cm", data = NULL))

## pre fire - transect

pre_hbt_hollows<-pre_nopot%>%filter(hollowpres_fix == 1)%>%
  ggplot(aes(hollownumber_fix, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = NULL, y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  #scale_x_continuous(limits = c(0,300), breaks = seq(0,300,100))+
  scale_y_continuous(limits = c(0, 0.6),labels = scales::percent_format(scale = 100))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_hbt_hollows / post_hbt_hollows + plot_annotation(tag_levels = 'A')


# fire effect models ------------------------------------------------------

#fire sev effect

hollowmod_firesev<-glm(hollows~firesev_fact, 
                     family = 'binomial', 
                     data = post_model, 
                     na.action = 'na.fail')

summary(hollowmod_firesev)
plot(allEffects(mod = hollowmod_firesev), type = "response", ylim = c(0, 1))

est_holocc<-estimate_graph(hollowmod_firesev)


hollowmod_firesev_oc<-glm(hollow_num~firesev_fact,
                          family = 'poisson',
                          data = post_model,
                          na.action = 'na.fail')

summary(hollowmod_firesev_oc)
plot(allEffects(mod = hollowmod_firesev_oc), type = "response", ylim = c(0, 2))

est_holab<-estimate_graph(hollowmod_firesev_oc)


plot_est_holocc<-est_holocc$table%>%
  mutate(f = c(3,4,5,6),
         f_fact = as.factor(c('Low scorch', 'Medium scorch', 'High scorch', 'Canopy burnt')),
         f_fact = fct_relevel(f_fact, 'Canopy burnt', 'High scorch', 'Medium scorch', 'Low scorch'))%>%
  ggplot( 
    aes(x = f_fact, y = Coefficient)) +
  geom_hline(yintercept = 0, 
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = f_fact, 
                 y = Coefficient)) + 
  geom_linerange(aes(x = f_fact, 
                     ymin = conf.low_90,
                     ymax = conf.high_90),
                 lwd = 1) +
  geom_linerange(aes(x = f_fact, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                 lwd = 1/2) + 
  coord_flip()+
  labs(x = NULL, y = 'Estimates')+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 15, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


plot_est_holab<-est_holab$table%>%
  mutate(f = c(3,4,5,6),
         f_fact = as.factor(c('Low scorch', 'Medium scorch', 'High scorch', 'Canopy burnt')),
         f_fact = fct_relevel(f_fact, 'Canopy burnt', 'High scorch', 'Medium scorch', 'Low scorch'))%>%
  ggplot( 
    aes(x = f_fact, y = Coefficient)) +
  geom_hline(yintercept = 0, 
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = f_fact, 
                 y = Coefficient)) + 
  geom_linerange(aes(x = f_fact, 
                     ymin = conf.low_90,
                     ymax = conf.high_90),
                 lwd = 1) +
  geom_linerange(aes(x = f_fact, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                 lwd = 1/2) + 
  coord_flip()+
  labs(x = NULL, y = 'Estimates')+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))


plot_est_holocc + plot_est_holab + plot_annotation(tag_levels = 'A')

ggsave('firesev_est.svg',path = 'figures/', width = 25, height = 15.6, units = 'cm', dpi = 600)

#significance tests

#by firesev

summary(aov(hollows~as.factor(firesev_fact), data = post_model))
TukeyHSD(aov(hollows~as.factor(firesev_fact), data = post_model))

summary(aov(hollow_num~as.factor(firesev), data = post_model))
TukeyHSD(aov(hollow_num~as.factor(firesev), data = post_model))

#by simplified firesev:

summary(aov(hollows~sev_simp, data = post_model))
TukeyHSD(aov(hollows~sev_simp, data = post_model))

summary(aov(hollow_num~sev_simp, data = post_model))
TukeyHSD(aov(hollow_num~sev_simp, data = post_model))

#LHBTs only

summary(aov(hollows~as.factor(firesev), data = post_model%>%filter(hollows>0 & dbh>=100)))
TukeyHSD(aov(hollows~as.factor(firesev), data = post_model%>%filter(hollows>0 & dbh>=100)))

#trees with large hollows

summary(aov(hollows~as.factor(firesev), data = post_model%>%filter(hollows>0 & !hollow_size %in% c('<5cm'))))
TukeyHSD(aov(hollows~as.factor(firesev), data = post_model%>%filter(!hollow_size %in% c('<5cm'))))

# probability of being a HBT ----------------------------------------------

hollowmod_dbh_post<-glm(hollows~dbh, 
                       family = 'binomial', 
                       data = post_model, 
                       na.action = 'na.fail')

summary(hollowmod_dbh_post)
plot(allEffects(mod = hollowmod_dbh_post), type = "response", ylim = c(0, 1))

hollowmod_dbh_post_high<-glm(hollows~dbh, 
                        family = 'binomial', 
                        data = post_model%>%filter(firesev >4), 
                        na.action = 'na.fail')

plot(allEffects(mod = hollowmod_dbh_post_high), type = "response", ylim = c(0, 1))



hollowmod_dbh_pre<-glm(as.numeric(hollowpres_fix)~dbh, 
                        family = 'binomial', 
                        data = pre_nopot, 
                        na.action = 'na.fail')

summary(hollowmod_dbh_pre)
plot(allEffects(mod = hollowmod_dbh_pre), type = "response", ylim = c(0, 1))

#predict and visualise

hol_dbh_pred_pre<-expand.grid(seq(0,300,5), 'pre')%>%rename(dbh = 1, fire = 2)
hol_dbh_pred_post<-expand.grid(seq(0,300,5), 'post')%>%rename(dbh = 1, fire = 2)

hol_dbh_pred_pre$pred<-predict(hollowmod_dbh_pre, hol_dbh_pred_pre, type = 'response')
hol_dbh_pred_post$pred<-predict(hollowmod_dbh_post, hol_dbh_pred_post, type = 'response')

hol_dbh_pred<-rbind(hol_dbh_pred_pre, hol_dbh_pred_post)%>%mutate(transect = 'all')

ggplot(hol_dbh_pred, aes(x = dbh, y = pred*100, color = fire))+
  geom_line()+
  geom_point()

#plot diff

hol_dbh_pred_diff<-cbind(hol_dbh_pred_pre, hol_dbh_pred_post)%>%
  select(1, pre = 3, post = 6)%>%
  mutate(diff = post-pre)

ggplot(hol_dbh_pred_diff, aes(x = dbh, y = diff*100))+
  geom_line()+
  geom_point()


#models by elevation band

#function

holmod<-function(data,
                 response,
                 predictor,
                 fire,
                 transect){
  
  modeldat<-data%>%select(hollow = {{response}}, dbh = {{predictor}})
  
  mod<-glm(as.numeric(hollow)~dbh, 
           family = 'binomial',
           data = modeldat,
           na.action = 'na.fail')
  
  predmat<-expand.grid(seq(0,300,5), fire, transect)%>%rename(dbh = 1, fire = 2, transect = 3)
  
  predmat$pred<-predict(mod, predmat, type = 'response')
  
  ggplot(predmat, aes(x = dbh, y = pred*100, color = fire))+
    geom_line()+
    geom_point()
  
  return(predmat)
  
}



#post models

post_low<-holmod(data = post_model%>%filter(transect == 'lowlands'),
                 response = 'hollows',
                 predictor = 'dbh',
                 fire = 'post',
                 transect = 'lowlands')

post_mid<-holmod(data = post_model%>%filter(transect == 'mid-hills'),
                 response = 'hollows',
                 predictor = 'dbh',
                 fire = 'post',
                 transect = 'mid-hills')

post_high<-holmod(data = post_model%>%filter(transect == 'high-country'),
                 response = 'hollows',
                 predictor = 'dbh',
                 fire = 'post',
                 transect = 'high-country')

#pre models

pre_low<-holmod(data = pre_nopot%>%filter(transect == 'lowlands'),
       response = 'hollowpres_fix',
       predictor = 'dbh',
       fire = 'pre',
       transect = 'lowlands')

pre_mid<-holmod(data = pre_nopot%>%filter(transect == 'mid-hills'),
                response = 'hollowpres_fix',
                predictor = 'dbh',
                fire = 'pre',
                transect = 'mid-hills')

pre_high<-holmod(data = pre_nopot%>%filter(transect == 'high-country'),
                response = 'hollowpres_fix',
                predictor = 'dbh',
                fire = 'pre',
                transect = 'high-country')


a<-ggplot(rbind( pre_low,
                 pre_mid,
                 pre_high,
                 post_low, 
             post_mid, 
             post_high,
             hol_dbh_pred_pre%>%mutate(transect = 'all'),
             hol_dbh_pred_post%>%mutate(transect = 'all')), aes(x = dbh, y = pred*100, color = fire))+
  geom_line()+
  geom_point()+
  facet_grid(~transect, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                 'mid-hills' = 'Mid-elevation',
                                                 'high-country' = 'High elevation',
                                                 'all' = 'Combined')))+
  scale_color_manual(labels = c('Pre-fire', 'Post-fire'),
                     values = c('#FED976', '#800026'))+
  scale_x_continuous(breaks = seq(0,300,50))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  labs(x = NULL, y = 'Prob. of hollow occurence', color = NULL)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 12,  color = 'black'),
        axis.text.x = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"))
  
# difference

diff_fun<-function(post, pre){
  
  cbind(post,
        pre)%>%
    select(1, 3, pre = 4, post = 8)%>%
    mutate(diff = post-pre)%>%
    select(1,2,5)
}

diff_all<-rbind(diff_fun(pre_low, post_low),
                diff_fun(pre_mid, post_mid),
                diff_fun(pre_high, post_high),
                diff_fun(hol_dbh_pred_pre%>%mutate(transect = 'all')%>%select(1,2,4,3), 
                         hol_dbh_pred_post%>%mutate(transect = 'all')%>%select(1,2,4,3)))

b<-ggplot(diff_all, aes(x = dbh, y = diff*100))+
  geom_line()+
  geom_point()+
  facet_grid(~transect, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                 'mid-hills' = 'Mid-elevation',
                                                 'high-country' = 'High elevation',
                                                 'all' = 'Combined')))+
  scale_x_continuous(breaks = seq(0,300,50))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  labs(y = 'Change in prob.', x = 'Tree diameter (cm)')+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 12,  color = 'black'),
        axis.text.x = element_text(size = 12,  color = 'black'))

a/b + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

ggsave('dbh_holprob.svg',path = 'figures/', width = 30, height = 15.6, units = 'cm', dpi = 600)

# stats

diff_all%>%
  group_by(transect)%>%
  filter(diff == max(diff))
  


# sign. tests tree size ---------------------------------------------------

#significant difference/change in tree size/hbt size after fire?

#difference in hbt size by firesev

hbt_aov<-aov(log(dbh)~as.factor(firesev), data = post_model%>%filter(hollows>0)) #log transform these data for paper

#Anova(hbt_aov)

summary(hbt_aov)
hbt_aov$coefficients

TukeyHSD(hbt_aov)

#diff all trees

tree_aov<-aov(dbh~as.factor(firesev), data = post_model)

TukeyHSD(tree_aov)

#compare pre and post-fire HBT size

post_hbt<-post_model%>%filter(hollows>0)%>%dplyr::select(dbh, hollow_num, transect)%>%mutate(fire = 'post')# %>% sample_n(206)
pre_hbt<-pre_nopot%>%filter(hollowpres_fix>0)%>%dplyr::select(dbh, hollow_num = hollownum, transect)%>%mutate(fire = 'pre')
hbt_compare<-rbind(post_hbt, pre_hbt)

hbt_compare_aov<-aov(dbh~as.factor(fire), data = hbt_compare)

summary(hbt_compare_aov)

summary(glm(dbh~as.factor(fire), data = hbt_compare))

#compare pre + post tree size (all trees)

post_trees<-post_model%>%dplyr::select(dbh, hollow_num, transect)%>%mutate(fire = 'post')# %>% sample_n(206)
pre_trees<-pre_nopot%>%dplyr::select(dbh, hollow_num = hollownum, transect)%>%mutate(fire = 'pre')
tree_compare<-rbind(post_trees, pre_trees)

summary(glm(dbh~as.factor(fire), data = tree_compare))

# animal occurrence by elev (pre) and firesev (post) -----------------------

## GGs and n arboreal/hollow dep. species

#post

summary(glm(prop_hbt/100~as.factor(gg_pres), family = 'quasibinomial', data = post_holha_av))


summary(aov(prop_hbt~as.factor(gg_pres), data = post_holha_av))
TukeyHSD(aov(prop_hbt~as.factor(gg_pres), data = post_holha_av))


summary(aov(gg_pres~prop_hbt, data = post_holha_av))
summary(aov(gg_pres~prop_lhbt, data = post_holha_av))

summary(glm(gg_pres~prop_hbt, data = post_holha_av%>%filter(firesev > 2)))

summary(aov(gg_pres~firesev, data = post_holha_av))
summary(aov(gg_pres~sev_simp, data = post_holha_av))

summary(aov(gg_pres~transect, data = post_holha_av))

TukeyHSD(aov(gg_pres~transect, data = post_holha_av))

#pre

summary(aov(GG~prop_hbt, data = pre_holha_nopot))
summary(aov(GG~prop_lhbt, data = pre_holha_nopot))

summary(glm(GG~prop_hbt, data = pre_holha_nopot))


summary(aov(GG~transect, data = pre_holha_nopot))
TukeyHSD(aov(GG~transect, data = pre_holha_nopot))

## get n arboreal species

post_spot<-read.csv('data/post_fire/spot_fixed.csv')

post_count<-post_spot%>%group_by(Site.Name, species)%>%
  summarise(n_obs = sum(n_animals))%>%
  filter(!species %in% c("Dingo & Dog (feral)",
                         "Deer",
                         "Domestic Cat (feral)",
                         "Tawny Frogmouth",
                         "Australian Owlet-nightjar"))%>%
  group_by(Site.Name)%>%
  summarise(n_spec = length(species))%>%
  rename(site_ID = 1)

pre_spot<-read.csv('data/pre_fire/spotlighting.csv')
pre_spot$plot<-interaction(pre_spot$transectID,pre_spot$plotID, drop = T, sep = 'P')
pre_spot<-pre_spot%>%mutate(plot = as.factor(str_remove(plot, '[.]')),
                            plot = case_when(plot == 'T2 P3' ~ 'T2P3',
                                             TRUE ~ plot))

pre_count<-pre_spot%>%group_by(plot, species)%>%
  summarise(n_obs = sum(count))%>%
  filter(!species %in% c("Feral Cat",
                         "Tawny frogmouth ",
                         "Owlet nightjar ",
                         "Feral dog",
                         "Victoria smooth froglet "))%>%
  group_by(plot)%>%
  summarise(n_spec = length(species))

pre_count<-rbind(pre_count, data.frame(plot = c('T2P2', 'T25P5'), n_spec = c(0,0))) # add plots with 0 obs

unique(pre_spot$species)
unique(post_spot$species)

#merge

post_holha_av_spec<-left_join(post_holha_av, post_count, by = 'site_ID')%>%
  mutate(n_spec = case_when(is.na(n_spec) ~ 0,
                            TRUE ~ n_spec))

pre_holha_nopot_spec<-left_join(pre_holha_nopot, pre_count, by = 'plot')

#test post

summary(aov(n_spec~prop_hbt, data = post_holha_av_spec))
summary(aov(n_spec~prop_lhbt, data = post_holha_av_spec))

summary(aov(n_spec~firesev, data = post_holha_av_spec))
summary(aov(n_spec~sev_simp, data = post_holha_av_spec))

summary(aov(n_spec~transect, data = post_holha_av_spec))

TukeyHSD(aov(n_spec~transect, data = post_holha_av_spec))

#test pre

summary(aov(n_spec~prop_hbt, data = pre_holha_nopot_spec))
summary(aov(n_spec~prop_lhbt, data = pre_holha_nopot_spec))

summary(aov(n_spec~transect, data = pre_holha_nopot_spec))
TukeyHSD(aov(n_spec~transect, data = pre_holha_nopot_spec))

## differences in dbh/tree size in occ and unocc sites

#post

summary(aov(dbh~as.factor(gg_pres), data = post_model))
summary(aov(dbh~as.factor(gg_pres), data = post_model%>%filter(hollows>0)))
summary(aov(height~as.factor(gg_pres), data = post_model%>%filter(hollows>0)))


summary(aov(gg_pres~dbh, data = post_model))
summary(glm(gg_pres~dbh, data = post_model))
summary(glm(gg_pres~dbh, data = post_model%>%filter(hollows>0)))

#pre

summary(aov(dbh~as.factor(GG), data = pre_nopot))
summary(aov(dbh~as.factor(GG), data = pre_nopot%>%filter(hollowpres_fix>0)))
summary(aov(height~as.factor(GG), data = pre_nopot%>%filter(hollowpres_fix>0)))


summary(glm(GG~dbh, data = pre_nopot))
summary(glm(GG~dbh, data = pre_nopot%>%filter(hollowpres_fix>0)))

## plots

post_holha_av%>%
  ggplot(aes(x = as.factor(gg_pres), y = prop_hbt))+
  geom_boxplot()

pre_holha_nopot%>%
  ggplot(aes(x = as.factor(GG), y = prop_hbt))+
  geom_boxplot()

## post-fire av. dbh

pre_nopot%>%group_by(GG)%>%
  summarise(meandbh = mean(dbh),
            sddbh = sd(dbh))

post_model%>%filter(!site_ID %in% c('ARI BBRR 40'))%>%
  group_by(gg_pres)%>%
  summarise(meandbh = mean(dbh),
            sddbh = sd(dbh))

# multivariate hollow models ----------------------------------------------

## post-fire: model with best variable from each group

post_hollowpres_final<-glm(hollows~ageclass_num+layer_num+dbh+
                       AHMI+firesev_fact+TWI+
                       Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                     family = 'binomial', 
                     data = post_model,
                     na.action = 'na.fail')

post_hollowpres_final_dredge<-dredge_fix(post_hollowpres_final)

plot(allEffects(mod = post_hollowpres_final_dredge$best_model), type = "response", ylim = c(0, 1))

summary(post_hollowpres_final_dredge$best_model)

#using volume and crown

post_model_crown<-post_model%>%filter(crown_area>0)

post_hollowpres_final_crown<-glm(hollows~layer_num+volume+crown_area+ageclass_num+
                        AHMI+firesev_fact+TWI+fires+
                        Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                      family = 'binomial', 
                      data = post_model_crown)

post_hollowpres_final_crown_dredge<-dredge_fix(post_hollowpres_final_crown)

plot(allEffects(mod = post_hollowpres_final_crown_dredge$best_model), type = "response", ylim = c(0, 1))

summary(post_hollowpres_final_crown_dredge$best_model)

## pre-fire: model with bas variable from each group

pre_nopot_vol<-pre_nopot%>%
  mutate(volume = ba*(height/3))%>%
  filter(hollowpres %in% c('0', '1'))%>%
  rename(TWI = TWI250, Total.Nitrogen = soil_total.nitrogen)

pre_hollowpres_final<-glm(as.numeric(hollowpres_fix)~volume+layer_num+ageclass_num+AHMI+
                            TWI+Total.Nitrogen, 
                          family = 'binomial', 
                          data = pre_nopot_vol)

pre_hollowpres_final_dredge<-dredge_fix(pre_hollowpres_final)

plot(allEffects(mod = pre_hollowpres_final_dredge$best_model), type = "response", ylim = c(0, 1))

summary(pre_hollowpres_final_dredge$best_model)


## plotting partial dependence plots

postvars<-c('ageclass_num', 'crown_area', 'layer_num', 'pH.Level..CaCl2.', 'Total.Nitrogen')
prevars<-c('ageclass_num', 'volume', 'layer_num', 'TWI', 'Total.Nitrogen')

#functions for each variable

part_fun<-function(var, model){
  
  part<-model%>%
    pdp::partial(pred.var = c(var),  prob = T, rug = T)%>%
    pivot_longer(1)
  
  return(part)
  
}

part_plotfun<-function(dataset, varname, yname, ylim = c(0,1)){
  
  a<-dataset%>%
    ggplot()+
    geom_line(aes(y=yhat, x=value))+
    ylim(ylim)+
    labs(x = varname, y = yname)
  
  return(a)
  
}

#get estimates

post_part<-map(postvars, part_fun, post_hollowpres_final_crown_dredge$best_model)
names(post_part)<-postvars

pre_part<-map(prevars, part_fun, pre_hollowpres_final_dredge$best_model)
names(pre_part)<-prevars

#make generic plots

post_part_plot<-map2(post_part, postvars, part_plotfun, yname = 'Prob. of hollow occ.')
pre_part_plot<-map2(pre_part, prevars, part_plotfun, yname = 'Prob. of hollow occ.')

#final plots

postplot<-
  post_part_plot$ageclass_num+
  labs(x = 'Age class')+
  theme_bw() +
  theme(axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13), 
        axis.text.y = element_text(size = 12,  color = 'black'),
        axis.text.x = element_text(size = 12,  color = 'black'))+
  
  post_part_plot$layer_num+
  labs(x = 'Structural layer', y = NULL)+
  theme_bw() +
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
    
  post_part_plot$Total.Nitrogen+
  labs(x = 'Soil nitrogen', y = NULL)+
  theme_bw() + 
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
  
  post_part_plot$pH.Level..CaCl2.+
  labs(x = 'Soil pH')+
  theme_bw()+
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
  
  post_part_plot$crown_area+
  labs(x = bquote('Crown area ' (m^2)), y = NULL)+
  theme_bw()+
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))

preplot<-
  pre_part_plot$ageclass_num+
  labs(x = 'Age class')+
  theme_bw() +
  theme(axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13), 
        axis.text.y = element_text(size = 12,  color = 'black'),
        axis.text.x = element_text(size = 12,  color = 'black'))+
  
  pre_part_plot$layer_num+
  labs(x = 'Structural layer', y = NULL)+
  theme_bw() +
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
  
  pre_part_plot$Total.Nitrogen+
  labs(x = 'Soil nitrogen', y = NULL)+
  theme_bw() + 
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
  
  pre_part_plot$TWI+
  labs(x = 'TWI value')+
  theme_bw()+
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))+
  
  pre_part_plot$volume+
  labs(x = bquote('Tree volume ' (m^3)), y = NULL)+
  theme_bw()+
    theme(axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13), 
          axis.text.y = element_text(size = 12,  color = 'black'),
          axis.text.x = element_text(size = 12,  color = 'black'))

preplot / postplot + plot_annotation(tag_levels = 'A')

ggsave('part_dep_models.svg',path = 'figures/', width = 30, height = 20, units = 'cm', dpi = 600)


## model results tables

summary(pre_hollowpres_final_dredge$best_model)
summary(post_hollowpres_final_crown_dredge$best_model)

t1<-tbl_regression(pre_hollowpres_final_dredge$best_model, intercept = T)
t2<-tbl_regression(post_hollowpres_final_crown_dredge$best_model, intercept = T)

hol_table_pre_post <-
  tbl_merge(
    tbls = list(t1, t2),
    tab_spanner = c("**Pre-fire**", "**Post-fire**")
  )%>%as_gt()

gtsave(hol_table_pre_post, 'tables/hollows_prepost.rtf')

## get undstandardised odd ratio

OR_pre<-effectsize::standardize_parameters(pre_hollowpres_final_dredge$best_model,  exp = TRUE)%>%
  mutate(model = 'pre')

OR_post<-effectsize::standardize_parameters(post_hollowpres_final_crown_dredge$best_model,  exp = TRUE)%>%
  mutate(model = 'post')

OR_both <-rbind(OR_pre, OR_post)

write.csv(OR_both, 'outputs/hollows_prepost_OR.csv', row.names = F)

# Appendix Figures --------------------------------------------------------

## S3 - pre-fire crown and height dist of HBTS


pre_hbt_crown<-pre_nopot%>%filter(hollowpres_fix == 1)%>%
  ggplot(aes(crown_area, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL, x = bquote('HBT crown area ' (m^2)), y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_x_continuous(limits = c(0,500), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0, 0.01),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_hbt_height<-pre_nopot%>%filter(hollowpres_fix == 1)%>%
  ggplot(aes(height, fill = transect))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+ #aes(y=..scaled..)
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'HBT height (m)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_x_continuous(limits = c(0,70), breaks = seq(0,70,10))+
  scale_y_continuous(limits = c(0, 0.07),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

pre_hbt_crown / pre_hbt_height + plot_annotation(tag_levels = 'A')

ggsave('hol_crown_height_dist.svg',path = 'figures/appendix/', width = 25, height = 15.6, units = 'cm', dpi = 600)

## S4 - post-fire DBH, crown + height by firesev

post_hbt_dbh_firesevt<-post_model%>%filter(hollows == 1)%>%
  ggplot(aes(dbh, fill = as.factor(firesev)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~as.factor(firesev), ncol = 2, labeller = as_labeller(c('2' = 'Unburnt', 
                                                                     '3' = 'Low canopy scorch',
                                                                     '4' = 'Medium canopy scorch',
                                                                     '5' = 'High canopy scorch',
                                                                     '6' = 'Canopy burnt')))+
  labs(title=NULL, x = 'HBT DBH (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  scale_y_continuous(limits = c(0, 0.025),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))


post_hbt_crown_firesevt<-post_model%>%filter(hollows == 1, crown_area>0)%>%
  ggplot(aes(crown_area, fill = as.factor(firesev)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~as.factor(firesev), ncol = 2, labeller = as_labeller(c('2' = 'Unburnt', 
                                                                     '3' = 'Low canopy scorch',
                                                                     '4' = 'Medium canopy scorch',
                                                                     '5' = 'High canopy scorch',
                                                                     '6' = 'Canopy burnt')))+
  labs(title=NULL, x = bquote('HBT crown area ' (m^2)), y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  scale_y_continuous(limits = c(0, 0.045),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

post_hbt_height_firesevt<-post_model%>%filter(hollows == 1)%>%
  ggplot(aes(height, fill = as.factor(firesev)))+
  #geom_histogram(bins = 20, color = 'black', size = 0.8, show.legend = F)+
  geom_density(show.legend = F)+
  facet_wrap(~as.factor(firesev), ncol = 2, labeller = as_labeller(c('2' = 'Unburnt', 
                                                                     '3' = 'Low canopy scorch',
                                                                     '4' = 'Medium canopy scorch',
                                                                     '5' = 'High canopy scorch',
                                                                     '6' = 'Canopy burnt')))+
  labs(title=NULL, x = 'HBT height (m)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'OrRd', name  = "firesev",labels=c('2' = 'Unburnt', 
                                                                 '3' = 'Low canopy scorch',
                                                                 '4' = 'Medium canopy scorch',
                                                                 '5' = 'High canopy scorch',
                                                                 '6' = 'Canopy burnt'))+
  scale_y_continuous(limits = c(0, 0.09),labels = scales::percent_format(scale = 1000))+
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12,  color = 'black'),
        axis.text.y = element_text(size = 12,  color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"),
        strip.background =element_rect(fill="white"),
        panel.spacing = unit(0.8, "cm", data = NULL))

post_hbt_dbh_firesevt + post_hbt_crown_firesevt + post_hbt_height_firesevt + plot_annotation(tag_levels = 'A')

ggsave('hol_firesev_dist.svg',path = 'figures/appendix/', width = 60, height = 30, units = 'cm', dpi = 600)


# Appendix tables ---------------------------------------------------------

## S2 + S3: animals (and number) per plot

post_nspec<-post_spot%>%group_by(Site.Name, species)%>%
  summarise(n_obs = sum(n_animals))%>%
  filter(!species %in% c("Dingo & Dog (feral)",
                         "Deer",
                         "Domestic Cat (feral)",
                         "Tawny Frogmouth",
                         "Australian Owlet-nightjar"))%>%
  rename(site = 1)

pre_nspec<-pre_spot%>%group_by(plot, species)%>%
  summarise(n_obs = sum(count))%>%
  filter(!species %in% c("Feral Cat",
                         "Tawny frogmouth ",
                         "Owlet nightjar ",
                         "Feral dog",
                         "Victoria smooth froglet "))%>%
  rename(site = 1)


#add elev band and + firesev


post_elev_firesev<-post_model%>%select(site = site_ID, transect, elevation, firesev)%>%
  mutate(fireclass = case_when(firesev == 2 ~ 'Unburnt',
                               firesev == 3 ~ 'Low canopy scorch',
                               firesev == 4 ~ 'Medium canopy scorch',
                               firesev == 5 ~ 'High canopy scorch',
                               firesev == 6 ~ 'Canopy burnt'))%>%
  select(!firesev)

pre_elev_firesev<-pre_nopot%>%select(site = plot, transect, elevation)%>%
  mutate(fireclass = 'Unburnt')

#combine

post_nspec_transect<-left_join(post_nspec, unique(post_elev_firesev), by = 'site')

pre_nspec_transect<-left_join(pre_nspec, unique(pre_elev_firesev), by = 'site')

#merge 

#spot_nspec<-full_join(pre_nspec_transect, post_nspec_transect , by = 'site')

write.csv(post_nspec_transect, 'outputs/spot_nspec_post.csv', row.names = F)
write.csv(pre_nspec_transect, 'outputs/spot_nspec_pre.csv', row.names = F)


## S4: post-fire site details

pre_ele<-pre_nopot%>%select(site = plot, prele = elevation)%>%unique

species_mix_post<-post_model%>%
  mutate(species = case_when(species == 'Euccalyptus viminalis' ~ 'Eucalyptus viminalis',
                             TRUE ~ species))%>%
  filter(grepl("Euc",species))%>%
  group_by(site_ID, species)%>%
  summarise(n = n())%>%
  mutate(species = str_remove(species, 'Eucalyptus'))%>%
  group_by(site_ID)%>%
  summarise(mix = paste0(species, collapse = ", "))%>%
  rename(site = 1)


post_details<-post_model%>%select(site = site_ID, transect, firesev, elevation, slope, aspect, gg_pres)%>%unique()%>%
  mutate(firesev = case_when(firesev == 2 ~ 'Unburnt',
                               firesev == 3 ~ 'Low canopy scorch',
                               firesev == 4 ~ 'Medium canopy scorch',
                               firesev == 5 ~ 'High canopy scorch',
                               firesev == 6 ~ 'Canopy burnt'))

post_details_speces<-left_join(post_details, species_mix_post, by = 'site')
post_details_speces_preele<-left_join(post_details_speces, pre_ele, by = 'site')

write.csv(post_details_speces_preele, 'outputs/post_site-descr.csv', row.names = F)
