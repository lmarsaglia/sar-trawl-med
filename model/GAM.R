library(dplyr)
library(plyr)
library(vegan)
library(mgcv)
library(voxel)
library(gridExtra)
library(Metrics)
library(ggplot2)

#load data for all med
x = readRDS(here::here("data","df_final.RData")) %>% filter(distance_to_coastline>5556)
x$id = as.character(x$id)

#load matching rate data
matching = readRDS(here::here("data","sar_ais_matching.RData"))
matching$id = as.character(matching$id)

thr = 90

#filter our final df to only these cells with 100% coverage
x = x[which(x$id %in% matching$id[which(matching$percentage <= thr )]),]

#aggregate data by model variables
colnames(x)[which(colnames(x) == "n")] = "nAIS"


xvars = unique(x[, c("id", "Gj", "depth", "distance_to_coastline", "SECT_COD")])
xy = aggregate(data = x[, c("id", "sar_n","nAIS")], 
               . ~ id,  FUN = "sum", na.rm = T)
xy = left_join(xy, xvars)

xy$Gj = as.numeric(decostand(xy$Gj, method = "standardize"))
xy = aggregate(data = xy, . ~ id + depth  + SECT_COD + Gj , FUN = "sum", na.rm = T)

k = 4 #https://rdrr.io/cran/mgcv/man/choose.k.html

# Model selection based on cross-validation: https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1557

# https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf
# https://stats.stackexchange.com/questions/418832/why-the-absence-of-probability-distribution-for-using-quasi-likelihood
# https://www.researchgate.net/post/How_to_select_the_best_model_fit_with_quasipoisson_family

vars = c("sar_n", "Gj", "depth", "distance_to_coastline")
model_tab = expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1))[ -1,]
nfolders = 10 # Number of groups for cross-validation
MSEs = data.frame(Model = numeric(nfolders * nrow(model_tab)),
                  MSE = numeric(nfolders * nrow(model_tab)))

df.fit = data.frame(Observed = numeric(0), Predicted = numeric(0)) # Store the values for graphical purposes

xy_folders = sample(1:10, nrow(xy), rep = T) #Generating the folders

irow = 0
for(i in 1:nrow(model_tab)){
  
  # Generating the formula
  iformula = "nAIS ~ "
  for(j in 1:ncol(model_tab)){
    if(model_tab[i, j] == 1)
      iformula = paste0(iformula, "s(", vars[j], ", k = k) +")
  }
  iformula = as.formula(substr(iformula, 1, nchar(iformula) - 1))
  
  # Cross validation: computing MSE for each combination...
  for(w in 1:nfolders){
    
    # Defining the Test dataset
    test_set = xy[which(xy_folders == w),]
    
    # Defining the Training dataset
    training_set = xy[which(xy_folders != w),]
    
    # Fit 
    i.gam.sar.ais = gam(data = training_set, formula = iformula,  family = quasipoisson(link = "sqrt"))
    
    # Predict on the Test set
    predicted.values = as.numeric(predict(i.gam.sar.ais, newdata = test_set, type = 'response', se.fit = TRUE)$fit)
    
    df.fit = rbind(df.fit,
                   data.frame(Observed = test_set$sar_n, Predicted = predicted.values))
    
    ik.mse = mse(actual = test_set$nAIS, predicted = predicted.values)
    
    irow = irow + 1
    
    MSEs$Model[irow] = i
    MSEs$MSE[irow] = ik.mse
  }
}


MSEs_agg = aggregate(data = MSEs, MSE ~ Model, FUN = "mean")

i.best = MSEs_agg$Model[which.min(MSEs_agg$MSE)]


# Generating and fitting the best model
iformula = "nAIS ~ "
for(j in 1:ncol(model_tab)){
  if(model_tab[i.best, j] == 1)
    iformula = paste0(iformula, "s(", vars[j], ", k = k) +")
}

iformula = as.formula(substr(iformula, 1, nchar(iformula) - 1))

gam.sar.ais <- gam(data = xy,
                   formula = iformula, 
              family = quasipoisson(link = "sqrt"))

summary(gam.sar.ais)


df.res = data.frame(Res = residuals(gam.sar.ais), GSA = xy$SECT_COD)


gres = ggplot(data = df.res, aes(x = Res)) + 
  geom_histogram() +
  theme_test() +
  xlab("Residuals (differences in monthly hours fishing)") + 
  xlim(-250, 250) + 
  ylab("Number of cells") + 
  ggtitle("B - Distribution of the residuals")

xy$predicted = gam.sar.ais$fitted.values

gfit = ggplot(data = xy, aes(x = nAIS, y = predicted)) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "lm") + 
  theme_test() +
  xlab("Observed fishing activity (Hours fishing)") +
  ylab("Predicted fishing activity (Hours fishing)") +
  scale_x_continuous(trans = "log10", limits = c(10, 10^4)) + 
  scale_y_continuous(trans = "log10", limits = c(10, 10^4)) +
  # geom_abline(slope = 1, intercept = 0, col = "red") +
  annotate("text", x = 30, y = 8000, 
           label = expression(paste(R^2, "= 0.58"))) +
  annotate("text", x = 30, y = 5000, 
           label = paste("% Var Exp =", 100*round(summary(gam.sar.ais)$dev.expl, 3))) + 
  ggtitle("A - Scatterplot comparing predicted and observed effort")

gfit
fgam1 = grid.arrange(gfit, gres, 
                     layout_matrix = matrix(c(1,1,2), nrow = 1))

#FIGURE 4
ggsave(plot = fgam1, 
       file = here::here("figures", "F4 - predicted VS observed.png"),
       width = 30, height = 20, units = "cm")

ggsave(plot = fgam1, 
       file = here::here("figures", "F4 - predicted VS observed.tiff"),
       device = "jpeg", width = 30, height = 20, units = "cm",dpi=700)




g1 = plotGAM(gamFit = gam.sar.ais, smooth.cov = "sar_n") +
  geom_segment(data =test_set, aes(x =sar_n, xend = sar_n, yend = 5, y = 0)) +
  ggtitle("Effect of SAR detections")+ xlim(0,80)
  
g2 = plotGAM(gam.sar.ais, smooth.cov = "Gj") +
  geom_segment(data =test_set, aes(x =Gj, xend = Gj, yend = -25, y = -30)) +
  ggtitle("Effect of G* index")

g3 = plotGAM(gam.sar.ais, smooth.cov = "depth") +
  geom_segment(data =test_set, aes(x =depth, xend = depth, yend = 5, y = 0)) +
  ggtitle("Effect of Depth") + xlim(-1500,5)


g4 = plotGAM(gam.sar.ais, smooth.cov = "distance_to_coastline") +
  geom_segment(data =test_set, aes(x =distance_to_coastline, xend = distance_to_coastline, yend = 5, y = 0)) +
  scale_x_continuous(trans = "log10") + 
  ggtitle("Effect of Distance from the coastline")

#FIGURE 5
fgam2 = grid.arrange(g1, g2, g3, g4, nrow = 2)
ggsave(plot = fgam2, 
       file = here::here("figures", "F5 - effect of smoothers.png"),
       device = "jpeg", width = 30, height = 20, units = "cm")
ggsave(plot = fgam2, 
       file = here::here("figures", "F5 - effect of smoothers.tiff"),
       device = "jpeg", width = 30, height = 20, units = "cm",dpi=700)


#PREDICT ALL MEDITERRANEAN
y = readRDS(here::here("data","df_final.RData")) %>% filter(distance_to_coastline>5556)
colnames(y)[which(colnames(y) == "n")] = "nAIS"

y$Gj = as.numeric(decostand(y$Gj, method = "standardize"))

y = y %>% dplyr::group_by(lon,lat,id,Gj,SMU_NAME,SECT_COD,distance_to_coastline,depth) %>% 
  dplyr::summarise(sar_n=sum(sar_n,na.rm = T),nAIS= sum(nAIS,na.rm = TRUE))


predicted_fishing = predict(gam.sar.ais, newdata = y, type = 'response', se.fit = TRUE)


y$predicted_fishing = predicted_fishing$fit
y$se = predicted_fishing$se.fit

saveRDS(y, here::here('data','gridf_pred.RData'))

