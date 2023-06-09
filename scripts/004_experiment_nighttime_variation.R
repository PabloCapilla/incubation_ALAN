###
###
#' 
#' Script for:
#' Experimental light at night explains differences in activity onset between urban and forest Great tits
#' Ciara L. O. McGlade, Pablo Capilla-Lasheras, Robyn J. Womack, Barbara Helm, Davide M. Dominoni
#' 
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of relative onset of activity experimental data
#' 
##
##

##
##
##### Script aim: #####
#' Analysis of night time variation in temperature
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, sjPlot, ggtext,
               lubridate, dplyr, tidyr,
               lme4, performance, rptR,
               ggplot2, extrafont)
loadfonts()

##
##
##### data #####
##
##
data <- readRDS("./data/clean_data_analysis.RDS")
data$first_offbout <- dmy_hm(data$first_offbout)
data$last_onbout <- dmy_hm(data$last_onbout)
data$inc_start_aprildays <- data$inc_start_aprildays - 90 # 90 = 1st April
data
head(data)



##
##
##### Sample sizes and data summary #####
##
##
data <- data %>% 
  filter(!is.na(night_var))

##
## table with sample size
## number of boxes
data %>% 
  group_by(box) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  group_by(type, site) %>% 
  count()

# number of observations
data %>% 
  group_by(box) %>% 
  ungroup() %>% 
  group_by(type, site) %>% 
  count()


##
##
##### Summary stats #####
##
##
head(data)
nrow(data) # total observations
length(unique(data$box)) # number of nest-boxes included

range(data$night_var)


##
##
###### Box plots #####
##
##
boxplot_var <- ggplot(data = data, aes(x = type, y = log(night_var), color = type, fill = type)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 10),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text.x = element_text("Arial", size = 10, angle = -45),
        axis.text.y = element_text("Arial", size = 10)) +
  geom_boxplot(alpha = 0.25, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.1), size = 0.75, shape = 1, fill = NA, color = "black") +
labs(x = " ", 
     y = expression(atop("Variation in log night time", 
                         paste("incubation temperature (", C^2, ")")))) +
  scale_x_discrete(labels =  c("Forest control", "Forest ALAN", "Urban")) +
  scale_fill_manual(name = "", labels = c("Forest control", "Forest ALAN", "Urban"), 
                    values = c("#4daf4a", "#cccc33", "#377eb8")) +
  scale_color_manual(name = "", labels = c("Forest control", "Forest ALAN", "Urban"), 
                     values = c("#4daf4a", "#cccc33", "#80b1d3"))

ggsave(filename = "./plots/Figure S1c.png", 
       plot = boxplot_var, 
       device = "png", 
       units = "mm",
       width = 75, 
       height = 100) 
##
##
##### Models for relative var of activity #####
##
##

# initial (full) model fit
model_var <- lmer(log(night_var) ~ 
                    
                    # interactions
                    type:poly(days_to_hatch,2)[,2] + 
                    type:poly(days_to_hatch,2)[,1] +  
                    
                    type:poly(inc_start_aprildays,2)[,2] + 
                    type:poly(inc_start_aprildays,2)[,1] +  
                    
                    # single effect terms
                    poly(days_to_hatch,2)[,2] +
                    poly(days_to_hatch,2)[,1] +
                    poly(inc_start_aprildays,2)[,2] +
                    poly(inc_start_aprildays,2)[,1] +
                    type + 
                    meantemp +
                    clutch_size + 
                    (1|box), 
                  REML = F,
                  na.action = "na.fail",
                  data= data) #full model
summary(model_var)

# model diagnostics
plot(model_var)
hist(residuals(model_var), freq=F, breaks = 100)
lines(density(x = rnorm(n = 10000,         # expected normal distribution
                        mean = 0, 
                        sd = sd(residuals(model_var)))))

##
## test significance of interactions in full model
drop1(model_var, test = "Chisq")

##
## test of significance for interactions

# 1
model_varb <- update(model_var, 
                     .~. - type:poly(inc_start_aprildays, 2)[, 2])
drop1(model_varb, test = "Chisq")

# 2 
model_varc <- update(model_varb, 
                     .~. - poly(days_to_hatch, 2)[, 2]:type)
drop1(model_varc, test = "Chisq")

# 3
model_vard <- update(model_varc, 
                     .~. - poly(days_to_hatch, 2)[, 1]:type)
drop1(model_vard, test = "Chisq")

# 4
model_vare <- update(model_vard, 
                     .~. - poly(inc_start_aprildays, 2)[, 1]:type)
drop1(model_vare, test = "Chisq")



##
## final full model
model_var_final <- lmer(log(night_var) ~
                          # single effect terms
                          poly(days_to_hatch,2)[,2] +
                          poly(days_to_hatch,2)[,1] +
                          poly(inc_start_aprildays,2)[,2] +
                          poly(inc_start_aprildays,2)[,1] +
                          type + 
                          meantemp +
                          clutch_size +
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data= data)
summary(model_var_final)
drop1(model_var_final, test = "Chisq")


##
## Table model results 
tab_model(model_var_final,
          file="./tables/Table S5 - RAW.doc",
          pred.labels = c("Intercept", 
                          "Days to hatching2",
                          "Days to hatching1",
                          "Date incubation start2", 
                          "Date incubation start1", 
                          "Experimental: Forest lighted",
                          "Experimental: Urban control",
                          "Mean ambient temperature",
                          "Clutch size"),
          string.ci = "CI (95%)",
          string.se = "SE",
          show.se = TRUE, 
          show.stat = FALSE,
          show.p = FALSE,
          show.est = TRUE,
          show.intercept = TRUE,
          rm.terms = NULL,
          show.re.var = T,
          show.ngroups = FALSE,
          show.r2 = FALSE,
          show.obs = FALSE,
          ci.hyphen = ",",
          use.viewer = T)

##
##
##### Plotting Model predictions #####
##
##

##
## model for plot
model_var_plot <- lmer(log(night_var) ~
                         # single effect terms
                         I(days_to_hatch^2) +
                         days_to_hatch +
                         I(inc_start_aprildays^2) +
                         inc_start_aprildays +
                         type + 
                         meantemp +
                         clutch_size +
                         (1|box), 
                       REML = F,
                       na.action = "na.fail",
                       data= data)
summary(model_var_plot)

##
## predictions plot
df_pred <- expand.grid(days_to_hatch = seq(min(data$days_to_hatch), 
                                           max(data$days_to_hatch), 1),
                       type = unique(data$type),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                                 max(data$inc_start_aprildays), 1),
                       meantemp = mean(data$meantemp))

# average temps and clutch size for each habitat
df_pred$meantemp <- ifelse(df_pred$type == "urban",
                           mean(data$meantemp[data$type == "urban"]),
                           mean(data$meantemp[data$type != "urban"]))
df_pred$clutch_size <- ifelse(df_pred$type == "urban",
                              mean(data$clutch_size[data$type == "urban"]),
                              mean(data$clutch_size[data$type != "urban"]))

# predictions
df_pred$prediction <- predict(model_var_plot, df_pred, re.form = NA)

# need to adapt code to new R version
# SE for mean predictions
mm <- model.matrix(~ 
                     # single effect terms
                     I(days_to_hatch^2) +
                     days_to_hatch +
                     I(inc_start_aprildays^2) +
                     inc_start_aprildays +
                     type + 
                     meantemp +
                     clutch_size,
                   data = df_pred)

pvar1 <- diag(mm %*% tcrossprod(vcov(model_var_plot),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)


# plot only data in range
data %>% 
  group_by(type) %>% 
  summarise(min_chr = min(inc_start_aprildays),
            max_chr = max(inc_start_aprildays))

remove_city <- which((df_pred$inc_start_aprildays < 28 | 
                        df_pred$inc_start_aprildays > 44) & 
                       df_pred$type == "urban")
remove_forest <- which((df_pred$inc_start_aprildays < 35 | 
                          df_pred$inc_start_aprildays > 45) & 
                         df_pred$type != "urban")
df_pred <- df_pred[-c(remove_city, remove_forest),]

##
## code type as ordered factor for plots
data$type <- as.character(data$type)
df_pred$type <- as.character(df_pred$type)

## plot hatching
var_hatching <- ggplot(data = data,
                       aes(x = days_to_hatch, 
                           y = log(night_var),
                           fill = type)) +
  geom_point(position = position_dodge(width = 0.85),
             shape = 21,
             color = "black",
             alpha = 0.5,
             size = 0.75) +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text("Arial", size = 10),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10)) +
  geom_errorbar(data = df_pred %>% 
                  group_by(type, days_to_hatch) %>% 
                  summarise(mean_plo = mean(plo),
                            mean_phi = mean(phi),
                            mean_pred = mean(prediction)), 
                aes(ymin = mean_plo, ymax = mean_phi, y = mean_pred, 
                    group = type),
                width = 0,
                size = 1,
                color = "black",
                position = position_dodge(width = 0.85)) +
  geom_point(data = df_pred %>% 
               group_by(type, days_to_hatch) %>% 
               summarise(mean_plo = mean(plo),
                         mean_phi = mean(phi),
                         mean_pred = mean(prediction)), 
             aes(y = mean_pred, group = type), 
             size = 2.5, 
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.85)) +
  labs(x = "Days before hatching", 
       y = expression(atop("Variation in log night time", 
                           paste("incubation temperature (", C^2, ")")))) +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("Forest control", "Forest ALAN", "Urban"), 
                    values = c("#4daf4a", "#cccc33", "#377eb8")) +
  scale_color_manual(name = "", labels = c("Forest control", "Forest ALAN", "Urban"), 
                     values = c("#4daf4a", "#cccc33", "#80b1d3"))

ggsave(filename = "./plots/Figure 2.png", 
       plot = var_hatching, 
       device = "png", 
       units = "mm",
       width = 100, 
       height = 85)  


# date from April 1
group_names <- list(
  'forest-control' = 'Forest control', 
  'forest-experimental light' = 'Forest ALAN', 
  'urban' = 'Urban'
)

group_labeller <- function(variable,value){
  return(group_names[value])
}


var_plot_date <- ggplot(data = data, 
                        aes(x = inc_start_aprildays, 
                            y = log(night_var),
                            fill = type,
                            color = type)) +
  geom_point(alpha = 0.35, size = 1) +
  theme_bw() +
  facet_wrap(~type, labeller=group_labeller) +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 10),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 15),
        axis.title = element_text("Arial", size = 12),
        axis.text = element_text("Arial", size = 12)) +
  geom_ribbon(data = df_pred %>% 
                group_by(type, inc_start_aprildays) %>% 
                summarise(mean_plo = mean(plo),
                          mean_phi = mean(phi),
                          mean_pred = mean(prediction)), 
              aes(ymin = mean_plo, ymax = mean_phi, y = mean_pred),
              color = NA,
              alpha = 0.5) +
  geom_line(data = df_pred %>% 
              group_by(type, inc_start_aprildays) %>% 
              summarise(mean_plo = mean(plo),
                        mean_phi = mean(phi),
                        mean_pred = mean(prediction)), 
            aes(y = mean_pred), 
            size = 1.5) +
  labs(x = "Days after April 1", 
       y = expression(atop("Variation in log night time", 
                           paste("incubation temperature (", C^2, ")")))) +
  scale_fill_manual(name = "", labels = c("Non-urban control", "Non-urban lighted", "Urban"), 
                    values = c("#4daf4a", "#cccc33", "#377eb8")) +
  scale_color_manual(name = "", labels = c("Non-urban control", "Non-urban lighted", "Urban"), 
                     values = c("#4daf4a", "#cccc33", "#80b1d3"))


ggsave(filename = "./plots/Figure S4.png", 
       plot = var_plot_date, 
       device = "png", 
       units = "mm",
       width = 180, 
       height = 100)  





