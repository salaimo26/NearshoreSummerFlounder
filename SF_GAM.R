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
trawls <- subset(trawls, dtdz>0)
trawls <- subset(trawls, TEMPBOT>1)
trawls <- subset(trawls, x_dist>0)
trawls <- subset(trawls, y_dist>0)
trawls <- subset(trawls, drho<0.5)
trawls <- subset(trawls, drho>=0)
trawls <- subset(trawls, delta_rho>=0)

#creating the categorical variables 
trawls$cruise <- as.factor(trawls$cruise)
trawls$location <- as.factor(trawls$location)
trawls$n_s <- as.factor(trawls$n_s)
trawls$decade <- as.factor(trawls$decade)

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

gam.check(m1.1)
concurvity(m1.1)
summary(m1.1)
getViz(m1.1)%>%plot(allTerms=T)%>%print(pages=1,shade=TRUE)
AIC(m1.1)

######################## GAM Index ######################## 
models_predict <- predict(m1.1, newdata=df_test, type="response")

plot(models_predict,df_test$para_cpue)

df_test

####---- Relative Contribution of Each Term ----####
gam.hp(m1.1)

####---- Cross Validation--Monte Carlo Method----####
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

mean(RMSE_values_mc)

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



