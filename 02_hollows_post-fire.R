# hollow models - post fire -----------------------------------------------

library(effects)
library(doParallel)
library(MuMIn)
library(lme4)
library(lavaan)
library(lavaanPlot)
library(corrplot)
library(pdp)
library(broom)
library(tidyverse)
library(ggrepel)
library(patchwork)

struc_hol_post<-read.csv('outputs/structure_fixed_perha.csv')%>%
  filter(grepl("Euc",species))%>%
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
                                                                            TRUE ~ as.numeric(hollow_num)))

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


# fire model --------------------------------------------------------------

holcols<-as.data.frame(colnames(struc_hol_post))

#past fires and fire severity

hollowmod_fire1<-glm(hollows~fires+firesev, 
                family = 'binomial', 
                data = struc_hol_post, 
                na.action = 'na.fail')

summary(hollowmod_fire1)

plot(allEffects(mod = hollowmod_fire1), type = "response", ylim = c(0, 1))

#interaction between fires and fire severity

hollowmod_fire2<-glm(hollows~fires+firesev+fires*firesev, 
                     family = 'binomial', 
                     data = struc_hol_post, 
                     na.action = 'na.fail')

summary(hollowmod_fire2)

plot(allEffects(mod = hollowmod_fire2), type = "response", ylim = c(0, 1))

#dredge

hollowmod_fire2_dredge<-dredge_fix(hollowmod_fire2)

#fire severity only

hollowmod_fire3<-glm(hollows~firesev, 
                     family = 'binomial', 
                     data = struc_hol_post, 
                     na.action = 'na.fail')

summary(hollowmod_fire3)

plot(allEffects(mod = hollowmod_fire3), type = "response", ylim = c(0, 1))

hollowmod_fire3%>%
  partial( pred.var = c('firesev'),  prob = T, rug = T)%>%
  autoplot(smooth = T, smooth.method = 'glm')


#factorial

hollowmod_fire4<-glm(hollows~firesev_fact, 
                     family = 'binomial', 
                     data = struc_hol_post, 
                     na.action = 'na.fail')

summary(hollowmod_fire4)

plot(allEffects(mod = hollowmod_fire4), type = "response", ylim = c(0, 0.5))


#number of hollows
#NEEDS FIXING of doublicate hollows (see data-clean script)

hollowmod_fire5<-glm(scale(hollow_num)~firesev_fact, 
                     data = struc_hol_post, 
                     na.action = 'na.fail')

summary(hollowmod_fire5)

plot(allEffects(mod = hollowmod_fire5), type = "response", ylim = c(-0.5, 0.5))

hollowmod_fire5%>%
  partial(pred.var = c('firesev_fact'),  prob = T, rug = T)%>%
  autoplot()

#ggplot estimate graph

#function for it

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

est_hbt<-estimate_graph(hollowmod_fire4)
est_holnum<-estimate_graph(hollowmod_fire5)

est_hbt_plot<-est_hbt$table%>%mutate(f = c(3,4,5,6))%>%
  ggplot( 
         aes(x = f, y = Coefficient)) +
  geom_hline(yintercept = 0, 
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = f, 
                 y = Coefficient)) + 
  geom_linerange(aes(x = f, 
                     ymin = conf.low_90,
                     ymax = conf.high_90),
                 lwd = 1) +
  geom_linerange(aes(x = f, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                 lwd = 1/2) + 
  scale_x_reverse()+
  coord_flip()+
  labs(x = 'Fire severity class', y = 'Estimates')+
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

est_holnum_plot<-est_holnum$table%>%mutate(f = c(3,4,5,6))%>%
  ggplot( 
    aes(x = f, y = Coefficient)) +
  geom_hline(yintercept = 0, 
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = f, 
                 y = Coefficient)) + 
  geom_linerange(aes(x = f, 
                     ymin = conf.low_90,
                     ymax = conf.high_90),
                 lwd = 1) +
  geom_linerange(aes(x = f, 
                     ymin = conf.low_95,
                     ymax = conf.high_95),
                 lwd = 1/2) + 
  scale_x_reverse()+
  coord_flip()+
  labs(x = NULL, y = 'Estimates')+
  theme_bw()+
  theme(plot.title = element_text(size =18, face='bold'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"))

est_hbt_plot + est_holnum_plot +
  plot_annotation(tag_levels = 'A')

ggsave('fire_effect_hols.svg',path = 'figures/', width = 20, height = 15.6, units = 'cm', dpi = 600)

# Group model - 1: Structure ----------------------------------------------

struc_hol_post_fire<-struc_hol_post%>%filter(firesev>2)

#check corrleation

corrplot(cor(select(na.omit(struc_hol_post), dbh, ba, ageclass_num, layer_num, layers_num, height, fires, firesev)),
         method = "number", type="lower", order="hclust", tl.cex = 0.75, tl.col="black", tl.srt = 45)

#simple using all rows

hollowmod_struc1<-glm(hollows~dbh+ageclass_num+layer_num+fires+firesev, 
                           family = 'binomial', 
                           data = struc_hol_post)

summary(hollowmod_struc1)
plot(allEffects(mod = hollowmod_struc1), type = "response", ylim = c(0, 1))

hollowmod_struc1_dredge<-dredge_fix(hollowmod_struc1)

#fire affected forests only


hollowmod_struc1_fire<-glm(hollows~dbh+ageclass_num+layer_num+fires+firesev, 
                      family = 'binomial', 
                      data = struc_hol_post_fire)

summary(hollowmod_struc1_fire)
plot(allEffects(mod = hollowmod_struc1_fire), type = "response", ylim = c(0, 1))

hollowmod_struc1_fire_dredge<-dredge_fix(hollowmod_struc1_fire)

#using volume and crown area + no fires

struc_hol_post_crown<-struc_hol_post%>%filter(crown_area>0)

hollowmod_struc2<-glm(hollows~crown_area+ageclass_num+layer_num+volume, 
                      family = 'binomial', 
                      data = struc_hol_post_crown)

summary(hollowmod_struc2)
plot(allEffects(mod = hollowmod_struc2), type = "response", ylim = c(0, 1))

hollowmod_struc2_dredge<-dredge_fix(hollowmod_struc2)

#for number of hollows

hollowmod_struc3<-glm(hollow_num~crown_area+ageclass_num+layer_num+volume+firesev, 
                      data = struc_hol_post_crown)

summary(hollowmod_struc3)
plot(allEffects(mod = hollowmod_struc3), type = "response", ylim = c(0, 5))

hollowmod_struc2_dredge<-dredge_fix(hollowmod_struc2)


# Group model - 2: Site ---------------------------------------------------

corrplot(cor(select(struc_hol_post, elevation,slope.1,
                    aspect,AHMI,climsuit, NDVI, TWI, 
                    cond, soilcode_1, soilcode_2,
                    fires, firesev)),
         method = "number", type="lower", order="hclust", tl.cex = 0.75, tl.col="black", tl.srt = 45)

hollowmod_site<-glm(hollows~slope.1+aspect+soilcode_1+soilcode_2+AHMI+climsuit+TWI+cond+NDVI+fires+firesev, 
                     family = 'binomial', 
                     data = struc_hol_post,
                     na.action = 'na.fail')

summary(hollowmod_site)

hollowmod_site_dredge<-dredge_fix(hollowmod_site)

summary(hollowmod_site_dredge$best_model)

#for number of hollows

hollowmod_site2<-glm(hollow_num~slope.1+aspect+soilcode_1+soilcode_2+AHMI+climsuit+TWI+cond+NDVI+fires+firesev, 
                    data = struc_hol_post,
                    na.action = 'na.fail')

summary(hollowmod_site2)

hollowmod_site_dredge2<-dredge_fix(hollowmod_site2)

plot(allEffects(mod = hollowmod_site_dredge2$best_model), type = "response", ylim = c(0, 3))

summary(hollowmod_site_dredge2$best_model)


# Group model - 3: Nutrients ---------------------------------------------

corrplot(cor(select(struc_hol_post, pH.Level..CaCl2.,Total.Nitrogen,Phosphorus.Colwell,
                    Potassium.Colwell,mean_leaf_n)),
         method = "number", type="lower", order="hclust", tl.cex = 0.75, tl.col="black", tl.srt = 45)

hollowmod_nutrients<-glm(hollows~pH.Level..CaCl2.+Total.Nitrogen+
                           Phosphorus.Colwell+Potassium.Colwell, 
                         family = 'binomial', 
                         data = struc_hol_post,
                         na.action = 'na.fail')

summary(hollowmod_nutrients)
plot(allEffects(mod = hollowmod_nutrients), type = "response", ylim = c(0, 1))

hollowmod_nutrients_dredge<-dredge_fix(hollowmod_nutrients)

summary(hollowmod_nutrients_dredge$best_model)

#for number of hollows

hollowmod_nutrients2<-glm(hollow_num~Total.Nitrogen+
                           Phosphorus.Colwell+Potassium.Colwell,
                          data = struc_hol_post,
                          na.action = 'na.fail')

summary(hollowmod_nutrients2)
plot(allEffects(mod = hollowmod_nutrients2), type = "response", ylim = c(0, 2))

hollowmod_nutrients_dredge2<-dredge_fix(hollowmod_nutrients2)

summary(hollowmod_nutrients_dredge2$best_model)


# Best model --------------------------------------------------------------

#pres:abs

summary(hollowmod_struc1_dredge$best_model)
summary(hollowmod_struc2_dredge$best_model)
summary(hollowmod_site_dredge$best_model)
summary(hollowmod_nutrients_dredge$best_model)

hollowmod_final<-glm(hollows~ageclass_num+layer_num+dbh+
                       AHMI+firesev_fact+TWI+
                       Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                         family = 'binomial', 
                         data = struc_hol_post,
                         na.action = 'na.fail')

summary(hollowmod_final)
plot(allEffects(mod = hollowmod_final), type = "response", ylim = c(0, 1))

hollowmod_final_dredge<-dredge_fix(hollowmod_final)

plot(allEffects(mod = hollowmod_final_dredge$best_model), type = "response", ylim = c(0, 1))

summary(hollowmod_final_dredge$best_model)

hollowmod_final_dredge$best_model%>%
  partial( pred.var = c('dbh'),  prob = T, rug = T)%>%
  autoplot(smooth = T, smooth.method = 'glm')

hollowmod_final_dredge$best_model%>%
  partial(pred.var = c('dbh', 'TWI'), chull = F, prob = T, rug = T, progress = T)%>%
  autoplot(contour = T, legend.title = "Hollow prob.")+
  coord_cartesian(expand = F)+
  labs(title=NULL,x = 'DBH (cm)', y = 'TWI')+
  theme_bw()

#three preds

cl <- makeCluster(4) # use 6 workers
registerDoParallel(cl) # register the parallel backend

partial(hollowmod_final_dredge$best_model, pred.var = c('dbh', 'TWI', 'firesev_fact'), 
        plot = TRUE, 
        chull = TRUE, prob = T, 
        parallel = TRUE, progress = T)

stopCluster(cl) # good practice


#fire affected trees only

hollowmod_final_fire<-glm(hollows~ageclass_num+layer_num+dbh+
                       AHMI+firesev_fact+TWI+
                       Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                     family = 'binomial', 
                     data = struc_hol_post_fire,
                     na.action = 'na.fail')

summary(hollowmod_final_fire)
plot(allEffects(mod = hollowmod_final_fire), type = "response", ylim = c(0, 1))

hollowmod_final_fire_dredge<-dredge_fix(hollowmod_final_fire)

plot(allEffects(mod = hollowmod_final_fire_dredge$best_model), type = "response", ylim = c(0, 1))

#reduced dataset with crowns and volume

struc_hol_post_crown

hollowmod_final2<-glm(hollows~layer_num+volume+crown_area+ageclass_num+
                       AHMI+firesev_fact+TWI+fires+
                       Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                     family = 'binomial', 
                     data = struc_hol_post_crown)

summary(hollowmod_final2)

plot(allEffects(mod = hollowmod_final2), type = "response", ylim = c(0, 1))

hollowmod_final2_dredge<-dredge_fix(hollowmod_final2)

plot(allEffects(mod = hollowmod_final2_dredge$best_model), type = "response", ylim = c(0, 1))

summary(hollowmod_final2_dredge$best_model)# used

#for number of hollows

hollowmod_final2_num<-glm(hollow_num~layer_num+volume+crown_area+ageclass_num+
                        AHMI+firesev_fact+TWI+fires+
                        Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                      data = struc_hol_post_crown)

summary(hollowmod_final2_num)

hollowmod_final2_num_dredge<-dredge_fix(hollowmod_final2_num)

plot(allEffects(mod = hollowmod_final2_num_dredge$best_model), type = "response", ylim = c(0, 4))

summary(hollowmod_final2_num_dredge$best_model)# used

#fire affected trees only

struc_hol_post_crown_fire<-struc_hol_post_fire%>%filter(crown_area>0)

hollowmod_final2_fire<-glm(hollows~layer_num+volume+crown_area+ageclass_num+
                        AHMI+firesev_fact+TWI+fires+
                        Total.Nitrogen+Phosphorus.Colwell+pH.Level..CaCl2., 
                      family = 'binomial', 
                      data = struc_hol_post_crown_fire)

summary(hollowmod_final2_fire)

plot(allEffects(mod = hollowmod_final2_fire), type = "response", ylim = c(0, 1))

hollowmod_final2_fire_dredge<-dredge_fix(hollowmod_final2_fire)

plot(allEffects(mod = hollowmod_final2_fire_dredge$best_model), type = "response", ylim = c(0, 1))

summary(hollowmod_final2_fire_dredge$best_model)


# prediction matrix with varying fire severity ----------------------------

summary(hollowmod_final)

summary(struc_hol_post$pH.Level..CaCl2.)

pred_mat_2<-expand.grid(ageclass_num = 1:5,
                      layer_num = 1:3,
                      dbh = seq(0,300, 10), 
                      AHMI = seq(15,40, 3), 
                      TWI = seq(5, 10, 2),
                      Total.Nitrogen = seq(0,1, 0.2),
                      Phosphorus.Colwell = seq(0,80, 10),
                      pH.Level..CaCl2. = seq(3,6, 2),
                      firesev_fact = 2)%>%
  mutate(firesev_fact = as.factor(firesev_fact))

pred_mat_6<-expand.grid(ageclass_num = 1:5,
                        layer_num = 1:3,
                        dbh = seq(0,300, 10), 
                        AHMI = seq(15,40, 3), 
                        TWI = seq(5, 10, 2),
                        Total.Nitrogen = seq(0,1, 0.2),
                        Phosphorus.Colwell = seq(0,80, 10),
                        pH.Level..CaCl2. = seq(3,6, 2),
                        firesev_fact = 6)%>%
  mutate(firesev_fact = as.factor(firesev_fact))


pred_mat_2$hol_pred<-predict(hollowmod_final, pred_mat_2, type = 'response')

pred_mat_6$hol_pred<-predict(hollowmod_final, pred_mat_6, type = 'response')

#plot

pred_mat_all<-rbind(pred_mat_2, pred_mat_6)%>%sample_n(1000)

pred_mat_all%>%
  ggplot(aes(x = dbh, y = hol_pred))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~firesev_fact)










