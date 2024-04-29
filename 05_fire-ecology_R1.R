# fire ecology review 1 -----------------------------------------------------

library(lme4)
library(performance)
library(sandwich)
library(lmtest)
library(tidyverse)
library(broom.mixed)

# fire effects model - account for nested design --------------------------

#cluster robust SE

post_model_noNA<-post_model%>%drop_na(sampling_point)

post_model_siteID<-post_model_noNA%>%mutate(site_no = as.numeric(as.factor(site_ID)))
  # %>%
  # group_by(firesev_fact)%>%
  # sample_n(100)


hollowmod_firesev<-glm(hollows~firesev_fact, 
                       family = 'binomial', 
                       data = post_model_siteID, 
                       na.action = 'na.fail')

summary(hollowmod_firesev)

tidy(hollowmod_firesev)

estimate_graph(hollowmod_firesev)

robust_se<-vcovCL(hollowmod_firesev, cluster = ~site_ID + sampling_point)

coef(hollowmod_firesev)

coeftest(hollowmod_firesev, vcov = robust_se)

hollowmod_firesev_oc<-glm(hollow_num~firesev_fact, 
                          family = 'poisson',
                          data = post_model, 
                          na.action = 'na.fail')

summary(hollowmod_firesev_oc)

est_holocc_rev<-estimate_graph(hollowmod_firesev)

est_holab_rev<-estimate_graph(hollowmod_firesev_oc)


plot_est_holocc_rev<-est_holocc_rev$table%>%
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


plot_est_holab_rev<-est_holab_rev$table%>%
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


plot_est_holocc_rev + plot_est_holab_rev + plot_annotation(tag_levels = 'A')

dir.create('figures/review')

ggsave('firesev_est_rev.svg',path = 'figures/review/', width = 25, height = 15.6, units = 'cm', dpi = 600)

#random effect

post_model_siteID<-post_model_noNA%>%mutate(site_no = as.numeric(as.factor(site_ID)))

hollowmod_firesev_lme <- glmer(hollows ~ firesev_fact + (1 | site_no:sampling_point), data = post_model_siteID, family = binomial)

hollowmod_firesev_lme_oc <- glmer(hollow_num ~ firesev_fact +  (1 | site_no:sampling_point), data = post_model_siteID, family = poisson)

summary(hollowmod_firesev_lme)

summary(hollowmod_firesev_lme_oc)


# proportion data ---------------------------------------------------------

#pre

summary(aov(prop_hbt~as.factor(transect), data = pre_holha_nopot))
TukeyHSD(aov(prop_hbt~as.factor(transect), data = pre_holha_nopot))

summary(glm(prop_hbt/100~as.factor(transect), family = 'quasibinomial', data = pre_holha_nopot))


#post

summary(aov(prop_hbt~as.factor(transect), data = post_holha_av))
TukeyHSD(aov(prop_hbt~as.factor(transect), data = post_holha_av))

summary(glm(prop_hbt/100~as.factor(transect), family = 'quasibinomial', data = post_holha_av))


summary(aov(prop_hbt~as.factor(firesev), data = post_holha_av))
TukeyHSD(aov(prop_hbt~as.factor(firesev), data = post_holha_av))

summary(glm(prop_hbt/100~as.factor(firesev), family = 'quasibinomial', data = post_holha_av))


# hollow models by dbh ----------------------------------------------------

hollowmod_dbh_post_mixed<-glmer(hollows~dbh + (1|site_no:sampling_point), 
                          family = 'binomial', 
                          data = post_model_siteID)

summary(hollowmod_dbh_post_mixed)

confint(hollowmod_dbh_post_mixed)

summary(hollowmod_dbh_post)


hollowmod_dbh_pre_mixed<-glmer(as.numeric(hollowpres_fix)~dbh + (1|plot), 
                                family = 'binomial', 
                                data = pre_nopot)


summary(hollowmod_dbh_pre_mixed)

confint(hollowmod_dbh_pre_mixed)

summary(hollowmod_dbh_pre)

#predict and visualise

hol_dbh_pred_pre_mixed<-expand.grid(seq(0,300,5), 'pre')%>%rename(dbh = 1, fire = 2)
hol_dbh_pred_post_mixed<-expand.grid(seq(0,300,5), 'post')%>%rename(dbh = 1, fire = 2)

hol_dbh_pred_pre_mixed$pred<-predict(hollowmod_dbh_pre_mixed, hol_dbh_pred_pre_mixed, type = 'response', ReForm = NA)
hol_dbh_pred_post_mixed$pred<-predict(hollowmod_dbh_post_mixed, hol_dbh_pred_post_mixed, type = 'response',  ReForm = NA)

hol_dbh_pred_mixed<-rbind(hol_dbh_pred_pre_mixed, hol_dbh_pred_post_mixed)%>%mutate(transect = 'all')

ggplot(hol_dbh_pred_mixed, aes(x = dbh, y = pred*100, color = fire))+
  geom_line()+
  geom_point()

#plot diff

hol_dbh_pred_diff_mixed<-cbind(hol_dbh_pred_pre_mixed, hol_dbh_pred_post_mixed)%>%
  select(1, pre = 3, post = 6)%>%
  mutate(diff = post-pre)

ggplot(hol_dbh_pred_diff_mixed, aes(x = dbh, y = diff*100))+
  geom_line()+
  geom_point()

#by elev. band

holmod_lme_pre<-function(data,
                 response,
                 predictor,
                 rf,
                 fire,
                 transect){
  
  modeldat<-data%>%select(hollow = {{response}}, dbh = {{predictor}}, site = {{rf}})
  
  mod<-glmer(as.numeric(hollow)~dbh + (1|site), 
           family = 'binomial',
           data = modeldat)
  
  predmat<-expand.grid(seq(0,300,5), fire, transect)%>%rename(dbh = 1, fire = 2, transect = 3)
  
  predmat$pred<-predict(mod, predmat, type = 'response', ReForm = NA)
  
  ggplot(predmat, aes(x = dbh, y = pred*100, color = fire))+
    geom_line()+
    geom_point()
  
  return(predmat)
  
}

holmod_lme_post<-function(data,
                         response,
                         predictor,
                         rf1,
                         rf2,
                         fire,
                         transect){
  
  modeldat<-data%>%select(hollow = {{response}}, dbh = {{predictor}}, site = {{rf1}}, plot = {{rf2}})
  
  mod<-glmer(as.numeric(hollow)~dbh + (1|site:plot), 
             family = 'binomial',
             data = modeldat)
  
  predmat<-expand.grid(seq(0,300,5), fire, transect)%>%rename(dbh = 1, fire = 2, transect = 3)
  
  predmat$pred<-predict(mod, predmat, type = 'response', ReForm = NA)
  
  ggplot(predmat, aes(x = dbh, y = pred*100, color = fire))+
    geom_line()+
    geom_point()
  
  return(predmat)
  
}

#pre

pre_low_mxd<-holmod_lme_pre(data = pre_nopot%>%filter(transect == 'lowlands'),
                response = 'hollowpres_fix',
                predictor = 'dbh',
                rf = 'plot',
                fire = 'pre',
                transect = 'lowlands')

pre_mid_mxd<-holmod_lme_pre(data = pre_nopot%>%filter(transect == 'mid-hills'),
                response = 'hollowpres_fix',
                predictor = 'dbh',
                rf = 'plot',
                fire = 'pre',
                transect = 'mid-hills')

pre_high_mxd<-holmod_lme_pre(data = pre_nopot%>%filter(transect == 'high-country'),
                 response = 'hollowpres_fix',
                 predictor = 'dbh',
                 rf = 'plot',
                 fire = 'pre',
                 transect = 'high-country')


#post

post_low_mxd<-holmod_lme_post(data = post_model_siteID%>%filter(transect == 'lowlands'),
                 response = 'hollows',
                 predictor = 'dbh',
                 rf1 = 'site_no',
                 rf2 = 'sampling_point',
                 fire = 'post',
                 transect = 'lowlands')

post_mid_mxd<-holmod_lme_post(data = post_model_siteID%>%filter(transect == 'mid-hills'),
                 response = 'hollows',
                 predictor = 'dbh',
                 rf1 = 'site_no',
                 rf2 = 'sampling_point',
                 fire = 'post',
                 transect = 'mid-hills')

post_high_mxd<-holmod_lme_post(data = post_model_siteID%>%filter(transect == 'high-country'),
                  response = 'hollows',
                  predictor = 'dbh',
                  rf1 = 'site_no',
                  rf2 = 'sampling_point',
                  fire = 'post',
                  transect = 'high-country')

a_mix<-ggplot(rbind( pre_low_mxd,
                 pre_mid_mxd,
                 pre_high_mxd,
                 post_low_mxd, 
                 post_mid_mxd, 
                 post_high_mxd,
                 hol_dbh_pred_pre_mixed%>%mutate(transect = 'all'),
                 hol_dbh_pred_post_mixed%>%mutate(transect = 'all')), aes(x = dbh, y = pred*100, color = fire))+
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

diff_all_mixed<-rbind(diff_fun(pre_low_mxd, post_low_mxd),
                diff_fun(pre_mid_mxd, post_mid_mxd),
                diff_fun(pre_high_mxd, post_high_mxd),
                diff_fun(hol_dbh_pred_pre_mixed%>%mutate(transect = 'all')%>%select(1,2,4,3), 
                         hol_dbh_pred_post_mixed%>%mutate(transect = 'all')%>%select(1,2,4,3)))

b_mix<-ggplot(diff_all_mixed, aes(x = dbh, y = diff*100))+
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

a_mix/b_mix + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

ggsave('dbh_holprob_rev.svg',path = 'figures/review/', width = 30, height = 15.6, units = 'cm', dpi = 600)


diff_all_mixed%>%
  group_by(transect)%>%
  filter(diff == max(diff))


# dbh dist HBT and all in one graph ---------------------------------------

### using only burnt sites


post_model_pre_only<-post_model%>%filter(grepl("T",site_ID), firesev > 2)


pre_model_post_only<-pre_nopot%>%filter(plot %in% unique(post_model_pre_only$site_ID))

## post fire - transect (pre burnt only)

post_compare_both<-ggplot()+
  geom_density(data = post_model_pre_only, aes(dbh, color = as.factor(hollows)),
               bins = 15, show.legend = T, lwd = 1.5)+
  #geom_histogram(data = post_model_pre_only%>%filter(hollows == 1), aes(dbh, fill = transect), show.legend = F)+
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = 'HBT diameter (cm)', y = 'Density', fill = 'firesev')+
  scale_fill_brewer(palette = 'Set2', name  ="Transect")+
  scale_color_brewer(palette = 'Dark2', name  ="Hollows", labels = c('Absent', 'Present'))+
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
        panel.spacing = unit(0.8, "cm", data = NULL),
        legend.position = 'bottom')

## pre fire - transect (pre burnt only)

pre_compare_both<-ggplot()+
  geom_density(data = pre_model_post_only, aes(dbh, color = as.factor(hollowpres_fix)), 
               bins = 15, show.legend = F, lwd = 1.5)+
  #geom_histogram(data = post_model_pre_only%>%filter(hollows == 1), aes(dbh, fill = transect), show.legend = F)+
  facet_wrap(~transect, ncol = 3, labeller = as_labeller(c('lowlands' = 'Lowlands', 
                                                           'mid-hills' = 'Mid-elevation',
                                                           'high-country' = 'High elevation')))+
  labs(title=NULL,x = NULL, y = 'Density', fill = 'firesev')+
  #scale_fill_brewer(palette = 'Set2', name  ="transect", labels=c("Lowlands", "Mid-hills", 'High elevation'))+
  scale_color_brewer(palette = 'Dark2', name  ="Hollows", labels = c('Absent', 'Present'))+
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

pre_compare_both / post_compare_both + plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

ggsave('hol_dbh_dist_compare_rev.svg',path = 'figures/review/', width = 25, height = 15.6, units = 'cm', dpi = 600)

# glider and n_spec significance tests ------------------------------------

#n GG

gg_mod_1<-glmer(gg_pres ~ prop_hbt + (1 | site_ID), data = post_holha_av, family = binomial)
gg_mod_1_1<-glmer(gg_pres ~ prop_lhbt + (1 | site_ID), data = post_holha_av, family = binomial)
gg_mod_2<-glmer(gg_pres ~ firesev + (1 | site_ID), data = post_holha_av, family = binomial)
gg_mod_3<-glmer(gg_pres ~ transect + (1 | site_ID), data = post_holha_av, family = binomial)

summary(gg_mod_1)
summary(gg_mod_1_1)
summary(gg_mod_2)
summary(gg_mod_3)

gg_mod_4<-glm(GG ~ prop_hbt + (1 | plot), data = pre_holha_nopot, family = binomial)
gg_mod_4_1<-glmer(GG ~ prop_lhbt + (1 | plot), data = pre_holha_nopot, family = binomial)
gg_mod_5<-glmer(GG ~ transect + (1 | plot), data = pre_holha_nopot, family = binomial)

summary(gg_mod_4)
summary(gg_mod_4_1)
summary(gg_mod_5)

#n spec

spec_mod_1<-glmer(n_spec ~ prop_hbt + (1 | site_ID), data = post_holha_av_spec, family = poisson)
spec_mod_2<-glmer(n_spec ~ prop_lhbt + (1 | site_ID), data = post_holha_av_spec, family = poisson)
spec_mod_3<-glmer(n_spec ~ sev_simp + (1 | site_ID), data = post_holha_av_spec, family = poisson)
spec_mod_4<-glmer(n_spec ~ transect + (1 | site_ID), data = post_holha_av_spec, family = poisson)

summary(spec_mod_1)
summary(spec_mod_2)
summary(spec_mod_3)
summary(spec_mod_4)

spec_mod_5<-glmer(n_spec ~ prop_hbt + (1 | plot), data = pre_holha_nopot_spec, family = poisson)
spec_mod_6<-glmer(n_spec ~ prop_lhbt + (1 | plot), data = pre_holha_nopot_spec, family = poisson)
spec_mod_7<-glmer(n_spec ~ transect + (1 | plot), data = pre_holha_nopot_spec, family = poisson)

summary(spec_mod_5)
summary(spec_mod_6)
summary(spec_mod_7)

#by observer

post_count_obs<-post_spot%>%group_by(Site.Name, observer_no, species)%>%
  summarise(n_obs = sum(n_animals))%>%
  filter(!species %in% c("Dingo & Dog (feral)",
                         "Deer",
                         "Domestic Cat (feral)",
                         "Tawny Frogmouth",
                         "Australian Owlet-nightjar"))%>%
  group_by(Site.Name, observer_no)%>%
  summarise(n_spec = length(species))%>%
  rename(site_ID = 1)

pre_count_obs<-pre_spot%>%group_by(plot, species)%>%
  summarise(n_obs = sum(count))%>%
  filter(!species %in% c("Feral Cat",
                         "Tawny frogmouth ",
                         "Owlet nightjar ",
                         "Feral dog",
                         "Victoria smooth froglet "))%>%
  group_by(plot)%>%
  summarise(n_spec = length(species))

#merge

post_holha_av_spec_obs<-left_join(post_holha_av, post_count_obs, by = 'site_ID')%>%
  mutate(n_spec = case_when(is.na(n_spec) ~ 0,
                            TRUE ~ n_spec),
         observer_no = case_when(is.na(observer_no) ~ 1,
                                 TRUE ~ observer_no))

spec_mod_1o<-glmer(n_spec ~ prop_hbt + (1 | observer_no:site_ID), data = post_holha_av_spec_obs, family = poisson)
spec_mod_2o<-glmer(n_spec ~ prop_lhbt + (1 | observer_no:site_ID), data = post_holha_av_spec_obs, family = poisson)
spec_mod_3o<-glmer(n_spec ~ sev_simp + (1 | observer_no:site_ID), data = post_holha_av_spec_obs, family = poisson)
spec_mod_4o<-glmer(n_spec ~ transect + (1 | observer_no:site_ID), data = post_holha_av_spec_obs, family = poisson)

summary(spec_mod_1o)
summary(spec_mod_2o)
summary(spec_mod_3o)
summary(spec_mod_4o)

#for GG

post_count_obs_GG<-post_spot%>%
  mutate(species = as.factor(species))%>%
  group_by(Site.Name, observer_no, species)%>%
  count(species, .drop = F)%>%
  filter(species %in% c('Southern Greater Glider'))%>%
  mutate(gg_occ = case_when(n>0 ~ 1,
                            TRUE ~ 0))%>%
  ungroup()%>%
  select(1,2,5)%>%
  rename(site_ID = 1)

#merge

post_holha_av_spec_obs_GG<-left_join(post_holha_av, post_count_obs_GG, by = 'site_ID')%>%
  mutate(gg_occ = case_when(is.na(gg_occ) ~ 0,
                            TRUE ~ gg_occ),
         observer_no = case_when(is.na(observer_no) ~ 1,
                                 TRUE ~ observer_no))

GG_mod_1o<-glmer(gg_occ ~ prop_hbt + (1 | observer_no:site_ID), data = post_holha_av_spec_obs_GG, family = binomial)
GG_mod_2o<-glmer(gg_occ ~ prop_lhbt + (1 | observer_no:site_ID), data = post_holha_av_spec_obs_GG, family = binomial)
GG_mod_3o<-glmer(gg_occ ~ sev_simp + (1 | observer_no:site_ID), data = post_holha_av_spec_obs_GG, family = binomial)
GG_mod_4o<-glmer(gg_occ ~ transect + (1 | observer_no:site_ID), data = post_holha_av_spec_obs_GG, family = binomial)

summary(GG_mod_1o)
summary(GG_mod_2o)
summary(GG_mod_3o)
summary(GG_mod_4o)

#save post-table per observer (Table S3)


post_nspec_table<-post_spot%>%group_by(Site.Name, observer_no, species)%>%
  summarise(n_obs = sum(n_animals))%>%
  filter(!species %in% c("Dingo & Dog (feral)",
                         "Deer",
                         "Domestic Cat (feral)",
                         "Tawny Frogmouth",
                         "Australian Owlet-nightjar"))%>%
  rename(site = 1)

#add elev band and + firesev


post_elev_firesev<-post_model%>%select(site = site_ID, transect, elevation, firesev)%>%
  mutate(fireclass = case_when(firesev == 2 ~ 'Unburnt',
                               firesev == 3 ~ 'Low canopy scorch',
                               firesev == 4 ~ 'Medium canopy scorch',
                               firesev == 5 ~ 'High canopy scorch',
                               firesev == 6 ~ 'Canopy burnt'))%>%
  select(!firesev)

#combine

post_nspec_transect_table<-left_join(post_nspec_table, unique(post_elev_firesev), by = 'site')

#merge 

write.csv(post_nspec_transect_table, 'outputs/spot_nspec_post_obs.csv', row.names = F)

# combined multivariate model ---------------------------------------------


post_imp<-post_model%>%filter(crown_area>0)%>%
  select(hollows, ageclass_num, layer_num, Total.Nitrogen,
         pH = pH.Level..CaCl2., crown_area, volume, TWI, burnt = firesev)%>%
  mutate(burnt = case_when(burnt>2 ~ 1, 
                           TRUE ~ 0),
         sampling = 'post')


pre_imp<-pre_nopot_vol%>%select(hollows = hollowpres_fix, 
                                ageclass_num, layer_num, Total.Nitrogen, 
                                pH = soil_pH.CaCl2, crown_area, volume, TWI)%>%
  mutate(burnt = 0,
         sampling = 'pre')

important_all<-rbind(post_imp, pre_imp)%>%mutate(hollows = as.factor(hollows))

table(important_all$burnt)

#model

combined_hollowmod<-glmer(hollows ~ ageclass_num + layer_num + Total.Nitrogen + pH + crown_area + volume + TWI + burnt
                            + (1 | sampling), data = important_all, family = binomial)

summary(combined_hollowmod)

t_all<-tbl_regression(combined_hollowmod, intercept = T)%>%as_gt()

gtsave(t_all, 'tables/hollows_prepost_rev.rtf')


effectsize::standardize_parameters(combined_hollowmod,  exp = TRUE)
