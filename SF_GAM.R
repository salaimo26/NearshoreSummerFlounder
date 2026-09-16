######################## Load in Packages ######################## 
library(mgcv) 
library(ggplot2)
library(tidyverse)
library(mgcViz)
library(rsample)
library(Metrics)
library(gridExtra)
library(png)
library(marginaleffects)
library(patchwork)
library(gam.hp)
library(scales)
df <- read.csv("njdep_sf_gam.csv")

######################## Data Pre-Processing for the GAM ######################## 
df <- df %>%
  mutate(month = lubridate::month(YRMODA),
         week = lubridate::week(YRMODA), 
         day = lubridate::day(YRMODA))

# sub-setting the data to remove some outliers 
trawls <- subset(df,SALBOT>30)
trawls <- subset(trawls, TEMPBOT>1)
trawls <- subset(trawls, x_dist>0)
trawls <- subset(trawls, y_dist>0)
trawls <- subset(trawls, delta_rho>=0)

#creating the categorical variables 
trawls$cruise <- as.factor(trawls$cruise)

######################## Splitting the data into testing and training sets ######################## 
# Put 70% of the data into the training set 
set.seed(726)
df_split <- rsample::initial_split(trawls, prop = 0.7)
# Create data frames for the two sets:
df_train <- rsample::training(df_split)
df_test  <- rsample::testing(df_split)

######################## Term Correlation ########################
df.vars <- trawls %>% select(x_dist,y_dist,TEMPBOT,TEMPSURF,SALBOT,SALSURF,dtdz,drho,delta_rho,STARTDEPTH,year,week,denssurf,densbot) %>%tidyr::drop_na()

## This is a correlation matrix
jpeg(filename = "term_correlation_plot.jpeg",width = 3600, height = 2400,res = 500)
corrplot::corrplot(cor(df.vars, method = 'spearman'),
                   method = 'number',
                   number.cex = 0.65,
                   type = 'upper',
                   col = rev(colorRampPalette(c("red", "white", "blue"))(200)),
                   tl.col = "black")  # Reversed color palette

dev.off()

######################## Final Model Used in the Paper ######################## 
m1.1 <- gam(para_cpue ~ 
            s(x_dist, by=cruise)
            +s(y_dist, by=cruise)
            +s(TEMPBOT,k=7)
            +s(SALBOT,k=5)
            +s(delta_rho,k=6)
            +te(week,year, k=c(7,10),bs=c('cc','tp'))
            ,method="REML",
            data=df_train,
            family=tw(link = 'log'))

par(mfrow=c(2,2))
gam.check(m1.1)
par(mfrow=c(1,1))

concurvity(m1.1)
summary(m1.1)
getViz(m1.1)%>%plot(allTerms=T)%>%print(pages=1,shade=TRUE)
AIC(m1.1)

######################## Relative Contribution of Each Term ######################## 
gam.hp(m1.1)

######################## GAM Predicted vs. Observed ######################## 
models_predict <- predict(m1.1,newdata = df_test,type = "response")

my_data <- data.frame(
  observed = df_test$para_cpue,
  predicted = models_predict
)

ggplot(my_data, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Observed CPUE",
    y = "Predicted CPUE"
  ) +
  theme_classic()

######################## GAM Index ######################## 
cruise_means <- trawls %>%
  group_by(cruise) %>%
  summarise(
    week = mean(week, na.rm = TRUE),
    TEMPBOT = mean(TEMPBOT, na.rm = TRUE),
    SALBOT = mean(SALBOT, na.rm = TRUE),
    delta_rho = mean(delta_rho, na.rm = TRUE),
    x_dist = mean(x_dist, na.rm = TRUE),
    y_dist = mean(y_dist, na.rm = TRUE)
  )

df_test2 <- expand.grid(
  cruise = unique(trawls$cruise),
  year = 1990:2019
) %>%
  left_join(cruise_means, by = "cruise") %>%
  arrange(year, cruise)
view(df_test2)

models_predict <- predict(m1.1,newdata = df_test2,type = "response")
view(models_predict)

df_test2$predicted <- models_predict

ggplot(df_test2, aes(x = year, y = predicted, group = cruise)) +
  geom_line(aes(linetype = cruise)) +
  geom_point(aes(shape = cruise)) +
  labs(
    x = "Year",
    y = "Predicted CPUE",
    linetype = "Cruise",
    shape = "Cruise"
  ) +
  theme_classic()

######################## Cross Validation--Monte Carlo Method ######################## 
cross_val <- mc_cv(trawls,prop = 0.7,times=100)

AIC_values_mc <- numeric(length(cross_val$id))
MAE_values_mc <- numeric(length(cross_val$id))
RMSE_values_mc <- numeric(length(cross_val$id))
r2_values_mc <- numeric(length(cross_val$id))

for (i in seq_along(cross_val$id)) {
  value <- cross_val$id[i]
  resample <- cross_val$splits[cross_val$id == value][[1]]
  training_set <- analysis(resample)
  testing_set <- assessment(resample)
  
  models <- gam(para_cpue ~ 
                s(x_dist, by=cruise)
                +s(y_dist, by=cruise)
                +s(TEMPBOT,k=7)
                +s(SALBOT,k=5)
                +s(delta_rho,k=6)
                +te(week,year, k=c(7,10),bs=c('cc','tp'))
                ,method="REML",
                data=training_set,
                family=tw(link = 'log'))
  
  AIC_values_mc[i] <- AIC(models)
  
  models_predict <- predict(models, newdata=testing_set, type = "response")
  
  RMSE_values_mc[i] <- rmse(testing_set$para_cpue,models_predict)
  MAE_values_mc[i] <- mae(testing_set$para_cpue, models_predict)
  r2_values_mc[i] <- cor(testing_set$para_cpue,models_predict)^2
  print(i)
}

median(RMSE_values_mc)

jpeg(filename = "crossval_AIC.jpeg",width = 1600, height = 1600,res = 300)
##AIC 
boxplot(
  list("AIC"=AIC_values_mc), 
  main = "Monte Carlo Cross Validation's AIC", 
  xlab = "AIC", 
  ylab = "Values"
)
dev.off()

jpeg(filename = "crossval_RMSE.jpeg",width = 1600, height = 1600,res = 300)
##RMSE 
boxplot(
  list("RMSE"=RMSE_values_mc), 
  main = "Monte Carlo Cross Validation's RMSE", 
  xlab = "RMSE", 
  ylab = "Values"
)
abline(h = 0.5, col = "red", lty = 2)
dev.off()

jpeg(filename = "crossval_R2.jpeg",width = 1600, height = 1600,res = 300)
##r2 
boxplot(
  list("r2"=r2_values_mc), 
  main = "Monte Carlo Cross Validation's r2", 
  xlab = "R2", 
  ylab = "Values"
)
abline(h = 0.5, col = "red", lty = 2)
dev.off()

######################## Plotting Partial Effects ######################## 
m2Viz <- getViz(m1.1)

## April Cross Dist
jpeg(filename = "april_cross_dist.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 2) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Cross Shore Distance (km)",y="s(x)",title = "April Cross Shore Distance")+scale_y_continuous(labels = label_number(accuracy = 0.1)) + theme(
  title = element_text(size=12,face="bold"),
  axis.title.x = element_text(size=13),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## August Cross 
jpeg(filename = "aug_cross_dist.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 4) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Cross Shore Distance (km)",y="s(x)",title = "August Cross Shore Distance")+scale_y_continuous(labels = label_number(accuracy = 0.1))+ theme(
  title = element_text(size=12,face="bold"),
  axis.title.x = element_text(size=13),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## April Along 
jpeg(filename = "april_along_dist.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 7) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Along Shore Distance (km)",y="s(x)",title = "April Along Shore Distance") + theme(
  title = element_text(size=12,face="bold"),
  axis.title.x = element_text(size=13),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## August Along 
jpeg(filename = "aug_along_dist.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 9) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Along Shore Distance (km)",y="s(x)",title = "August Along Shore Distance") + theme(
  title = element_text(size=12,face="bold"),
  axis.title.x = element_text(size=13),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## Bottom Temp 
jpeg(filename = "bot_temp.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 11) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Bottom Temperature (°C)",y="s(x)")+ theme(
  axis.title.x = element_text(size=13, face="bold"),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13,face="bold"),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## Bottom Salinity ~"(\u2030)"
jpeg(filename = "bot_sal.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 12) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs(x="Bottom Salinity (\u2030)",y="s(x)") + theme(
  axis.title.x = element_text(size=13, face="bold"),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13,face="bold"),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 12),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()


library(ggtext)

## delta rho  coord_cartesian(ylim = c(-5, 2)
jpeg(filename = "delta_rho.jpeg",width = 2400, height = 1600,res = 650)
p = plot(m2Viz, select = 13) + l_ciPoly() + l_fitLine() + l_rug()
p = p  + labs( x = "Δρ (kg*m^-3)",
               y = "s(x)") + theme(
                 axis.title.x = element_text(size=13, face="bold"),  # Adjust size of x-axis title
                 axis.title.y = element_text(size = 13,face="bold"),  # Adjust size of y-axis title
                 axis.text.x = element_text(size = 12),   # Adjust size of x-axis text
                 axis.text.y = element_text(size = 12),# Adjust size of y-axis text
                 plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()

## year/week
jpeg(filename = "year_week.jpeg",width = 2400, height = 1600,res = 650)
p <- plot(m2Viz, select = 14)

# Customize the x-axis ticks and labels
month_weeks <- c(4, 16, 25, 34, 43)
month_labels <- c("January", "April",  "June",
                  "August","October")

# Customize the x-axis to show months instead of week numbers
p <- p + scale_x_continuous(
  breaks = month_weeks,
  labels = month_labels
)
# Customize the y-axis ticks and labels
p <- p + scale_y_continuous(
  breaks = seq(1990, 2019, 1), # Tick marks for every year
  labels = ifelse(seq(1990, 2019, 1) %% 2 == 0, seq(1990, 2019, 1), "") # Labels for every other year
) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p = p  + labs(x= NULL,y="Year",title = NULL) + theme(
  title = element_text(size=13,face="bold"),
  axis.title.x = element_text(size=13),  # Adjust size of x-axis title
  axis.title.y = element_text(size = 13),  # Adjust size of y-axis title
  axis.text.x = element_text(size = 8),   # Adjust size of x-axis text
  axis.text.y = element_text(size = 8),# Adjust size of y-axis text
  plot.margin = margin(l=15,t=10,r=15))

print(p)
dev.off()
