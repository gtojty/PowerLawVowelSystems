# @author Tao Gong (gtojty@gmail.com)
# 2021/11/6

# required packages
require(stats)
require(lme4)
require(lmerTest)
require(EMAtools)
require(party)
require(lattice)
require(MuMIn)
require(boot)
require(ggplot2)
require(plotrix)

# for drawing figures
geom.text.size = 5; theme.size = 3*geom.text.size
sizex <- 10; sizey <- 8; hjust <- 0.5; angle <- 45
# number of digits
nsmall <- 4

savefig <- function(fileName, dpi, width, height, units, type){
  #' Function to save figures
  #' @param fileName file name of the figure
  #' @param dpi resolution, e.g., 300
  #' @param width figure width in inches (units)
  #' @param height figure height in inches (units)
  #' @param units unit name for width and height, e.g., 'in'
  #' @param type figure type: 'png', png figure; 'pdf', pdf figure; 'both', both png and pdf figures
  
  if(type=="png"){ file <- paste(fileName, ".png", sep=""); ggsave(file, dpi = dpi, width = width, height = height, units = units) }
  if(type=="pdf"){ file <- paste(fileName, ".pdf", sep=""); ggsave(file, dpi = dpi, width = width, height = height, units = units) }
  if(type=="both"){
    file <- paste(fileName, ".png", sep=""); ggsave(file, dpi = dpi, width = width, height = height, units = units)
    file <- paste(fileName, ".pdf", sep=""); ggsave(file, dpi = dpi, width = width, height = height, units = units)
  }
}


direct <- '.'
# read in data in excel
# require(readxl)
# df <- read_excel('../Supplementary Dataset S1_532 language samples.xlsx')
# df <- df[,c(1:8,10:15)]
# names(df) <- c('IDX',	'ISO693_Code', 'Language', 'Dialect', 'Genus', 'Recording_Method', 'Quantity', 'Vowel_Number', 
#                'Dispersion', 'Focalization', 'Articulatory_Space', 'Effective_DE', 'Geographic_Region', 'Cultural_Society')

# read in data in csv
df <- read.csv(file.path(direct, 'Supplementary Dataset S1_532 language samples.csv'))
df <- df[,c(1:8,10:15)]
names(df) <- c('IDX',	'ISO693_Code', 'Language', 'Dialect', 'Family', 'Recording_Method', 'Quantity', 'Vowel_Number', 
               'Dispersion', 'Focalization', 'Articulatory_Space', 'Effective_DE', 'Geographic_Region_simple', 'Geographic_Region')
write.csv(df, file.path(direct, 'Dataset_S1_Processed.csv'), row.names = FALSE)


#############################
# global level lme and figure
#############################
resdir <- file.path(direct, 'All'); dir.create(resdir, showWarnings = FALSE)

# read in data
df <- read.csv(file.path(direct, 'Dataset_S1_Processed.csv'))

# global level lm
lm_all <- lm(log10(Focalization) ~ log10(Effective_DE), data = df)
summary(lm_all)
# Call:
#   lm(formula = log10(Focalization) ~ log10(Effective_DE), data = df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37677 -0.07650  0.00047  0.07726  0.44183 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -1.08294    0.01811   -59.8   <2e-16 ***
#   log10(Effective_DE) -0.83491    0.02747   -30.4   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1231 on 530 degrees of freedom
# Multiple R-squared:  0.6355,	Adjusted R-squared:  0.6348 
# F-statistic:   924 on 1 and 530 DF,  p-value: < 2.2e-16

# calculate r square
r.squaredGLMM(lm_all)
# R2m        R2c
# 0.6350449 0.6350449


# global level lme
lme_all <- lmer(log10(Focalization) ~ log10(Effective_DE) 
                 + (1|Geographic_Region:ISO693_Code) + (1|ISO693_Code) + (1|Family), data = df, control=lmerControl(optCtrl=list(maxfun=20000)))
summary(lme_all)
# REML criterion at convergence: -873.1
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.2490 -0.5724 -0.0050  0.5713  2.9492 
# 
# Random effects:
#   Groups                        Name        Variance Std.Dev.
# Geographic_Region:ISO693_Code (Intercept) 0.003121 0.05587 
# ISO693_Code                   (Intercept) 0.001613 0.04016 
# Family                        (Intercept) 0.006045 0.07775 
# Residual                                  0.007422 0.08615 
# Number of obs: 532, groups:  Geographic_Region:ISO693_Code, 225; ISO693_Code, 221; Family, 39
# 
# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)          -1.05699    0.02300 115.00246  -45.95   <2e-16 ***
#   log10(Effective_DE)  -0.76178    0.02819 522.11357  -27.02   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# lg10(Ef_DE) 0.690 

# calculate r square
r.squaredGLMM(lme_all)
# R2m       R2c
# 0.5467063 0.8151603

#### So, lmer model is better than simple lm model!

# fitting curve
# raw power law
x <- df$Effective_DE; y <- df$Focalization
z <- nls(y ~ a*x^(-b) + c, start = list(a=1, b=1, c=1))
summary(z)
# Formula: y ~ a * x^(-b) + c
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a  0.11529    0.04467   2.581   0.0101 *  
#   b  0.74332    0.13528   5.495 6.09e-08 ***
#   c -0.04748    0.06315  -0.752   0.4525    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1084 on 529 degrees of freedom
# 
# Number of iterations to convergence: 8 
# Achieved convergence tolerance: 3.147e-06

a <- coef(summary(z))[, 'Estimate'][1]; b <- coef(summary(z))[, 'Estimate'][2]; c <- coef(summary(z))[, 'Estimate'][3]
pwr <- function(x){ a*x^(-b) + c }

lab <- paste("Power-law formula:", 
             "\nFE = ",  format(round(a, nsmall), nsmall), "*", "(Effective DE)^(-", format(round(b, nsmall), nsmall), ")", format(round(c, nsmall), nsmall), 
             "\nAdj. R^2 = ", format(round(r.squaredGLMM(lm_all)[,'R2c'], nsmall), nsmall), sep="")
ggplot(df, aes(x=Effective_DE, y=Focalization)) + 
  geom_point(aes(color = Geographic_Region), size=2) + 
  geom_line(aes(x=Effective_DE, y=pwr(x)), size = 1) + 
  geom_text(aes(x=0.5, y=1, label=lab), hjust=0, size=geom.text.size) + 
  xlab('Effective DE') + ylab('FE') + 
  theme_bw() + 
  theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
fileName <- file.path(resdir, 'lme_all_raw')
savefig(fileName, 300, 10, 4, "in", "both")

# log ~ log linear
fit_intercept <- coef(summary(lme_all))[, 'Estimate'][1]; fit_slope <- coef(summary(lme_all))[,'Estimate'][2]
lab <- paste("Slope K = ", format(round(fit_slope, nsmall), nsmall), 
             "\nIntercept = ", format(round(fit_intercept, nsmall), nsmall), 
             "\np-value = ", round(coef(summary(lme_all))['log10(Effective_DE)','Pr(>|t|)'], abs(floor(log10(coef(summary(lme_all))['log10(Effective_DE)','Pr(>|t|)'])-nsmall))),
             "\nAdj. R^2 = ", format(round(r.squaredGLMM(lme_all)[,'R2c'], nsmall), nsmall), sep="")
ggplot(df, aes(x=log10(Effective_DE), y=log10(Focalization))) +
  geom_point(aes(color=Geographic_Region), size=2) + 
  geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
  geom_text(aes(x=-0.4, y=-0.3, label=lab), hjust=0, size=geom.text.size) + 
  xlab('log10(Effective DE)') + ylab('log10(FE)') + 
  theme_bw() + 
  theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
fileName <- file.path(resdir, 'lme_all_loglog')
savefig(fileName, 300, 10, 4, "in", "both")


##############################
# average level lm and figure
##############################
resdir <- file.path(direct, 'Average'); dir.create(resdir, showWarnings = FALSE)

# read in data
df <- read.csv(file.path(direct, 'Dataset_S1_Processed.csv'))

# create summary data
# by family
df_avg <- data.frame(matrix(ncol = 17, nrow = 0))
colnames(df_avg) <- c('Family', 'NoSample', 
                      'Vowel_Number_mean', 'Vowel_Number_median', 'Vowel_Number_se',  
                      'Dispersion_mean', 'Dispersion_median', 'Dispersion_se', 
                      'Focalization_log_mean', 'Focalization_log_median', 'Focalization_log_se',
                      'Articulatory_Space_mean', 'Articulatory_Space_median', 'Articulatory_Space_se',  
                      'Effective_DE_log_mean', 'Effective_DE_log_median', 'Effective_DE_log_se')
for(family in unique(df$Family)){
  subdf <- df[df$Family==family, ]
  df_avg <- rbind(df_avg, data.frame(Family = family, 
                                     NoSample = nrow(subdf), 
                                     Vowel_Number_mean = mean(subdf$Vowel_Number, na.rm=TRUE), Vowel_Number_median = median(subdf$Vowel_Number, na.rm=TRUE), Vowel_Number_se = std.error(subdf$Vowel_Number, na.rm=TRUE), 
                                     Dispersion_mean = mean(subdf$Dispersion, na.rm=TRUE), Dispersion_median = median(subdf$Dispersion, na.rm=TRUE), Dispersion_se = std.error(subdf$Dispersion, na.rm=TRUE),
                                     Focalization_log_mean = mean(log10(subdf$Focalization), na.rm=TRUE), Focalization_log_median = median(log10(subdf$Focalization), na.rm=TRUE), Focalization_log_se = std.error(log10(subdf$Focalization), na.rm=TRUE),
                                     Articulatory_Space_mean = mean(subdf$Articulatory_Space, na.rm=TRUE), Articulatory_Space_median = median(subdf$Articulatory_Space, na.rm=TRUE), Articulatory_Space_se = std.error(subdf$Articulatory_Space, na.rm=TRUE),
                                     Effective_DE_log_mean = mean(log10(subdf$Effective_DE), na.rm=TRUE), Effective_DE_log_median = median(log10(subdf$Effective_DE), na.rm=TRUE), Effective_DE_log_se = std.error(log10(subdf$Effective_DE), na.rm=TRUE)))
}
# sort by NoSample
df_avg <- df_avg[order(-df_avg$NoSample),]
# write to df
write.csv(df_avg, file.path(resdir, 'data_avg_family.csv'), row.names = FALSE)

# by geographical region
df_avg <- data.frame(matrix(ncol = 17, nrow = 0))
colnames(df_avg) <- c('Region', 'NoSample', 
                      'Vowel_Number_mean', 'Vowel_Number_median', 'Vowel_Number_se',  
                      'Dispersion_mean', 'Dispersion_median', 'Dispersion_se', 
                      'Focalization_log_mean', 'Focalization_log_median', 'Focalization_log_se',
                      'Articulatory_Space_mean', 'Articulatory_Space_median', 'Articulatory_Space_se',  
                      'Effective_DE_log_mean', 'Effective_DE_log_median', 'Effective_DE_log_se')
for(geo in unique(df$Geographic_Region)){
  subdf <- df[df$Geographic_Region==geo, ]
  df_avg <- rbind(df_avg, data.frame(Region = geo, 
                                     NoSample = nrow(subdf), 
                                     Vowel_Number_mean = mean(subdf$Vowel_Number, na.rm=TRUE), Vowel_Number_median = median(subdf$Vowel_Number, na.rm=TRUE), Vowel_Number_se = std.error(subdf$Vowel_Number, na.rm=TRUE), 
                                     Dispersion_mean = mean(subdf$Dispersion, na.rm=TRUE), Dispersion_median = median(subdf$Dispersion, na.rm=TRUE), Dispersion_se = std.error(subdf$Dispersion, na.rm=TRUE),
                                     Focalization_log_mean = mean(log10(subdf$Focalization), na.rm=TRUE), Focalization_log_median = median(log10(subdf$Focalization), na.rm=TRUE), Focalization_log_se = std.error(log10(subdf$Focalization), na.rm=TRUE),
                                     Articulatory_Space_mean = mean(subdf$Articulatory_Space, na.rm=TRUE), Articulatory_Space_median = median(subdf$Articulatory_Space, na.rm=TRUE), Articulatory_Space_se = std.error(subdf$Articulatory_Space, na.rm=TRUE),
                                     Effective_DE_log_mean = mean(log10(subdf$Effective_DE), na.rm=TRUE), Effective_DE_log_median = median(log10(subdf$Effective_DE), na.rm=TRUE), Effective_DE_log_se = std.error(log10(subdf$Effective_DE), na.rm=TRUE)))
}
# write to df
write.csv(df_avg, file.path(resdir, 'data_avg_geo.csv'), row.names = FALSE)


# draw lme figure based on average values
df_avg <- read.csv(file.path(resdir, 'data_avg_family.csv'))
# only use 9 biggest families
df_avg <- df_avg[c(1:5,7:10),]
lm_mean <- lm(Focalization_log_mean ~ Effective_DE_log_mean, data = df_avg)
summary(lm_mean)
# lm(formula = Focalization_log_mean ~ Effective_DE_log_mean, data = df_avg)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.123987  0.002356  0.008981  0.023637  0.141116 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -1.3595     0.1755  -7.744 0.000112 ***
#   Effective_DE_log_mean  -1.2722     0.2561  -4.968 0.001623 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.08711 on 7 degrees of freedom
# Multiple R-squared:  0.7791,	Adjusted R-squared:  0.7475 
# F-statistic: 24.68 on 1 and 7 DF,  p-value: 0.001623

r.squaredGLMM(lm_mean)
# R2m       R2c
# 0.755216 0.755216

# log ~ log linear with two error bars of se
fit_intercept <- coef(summary(lm_mean))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_mean))[,'Estimate'][2]
lab <- paste("Slope K = ", format(round(fit_slope, nsmall), nsmall), 
             "\nIntercept = ", format(round(fit_intercept, nsmall), nsmall), 
             "\np-value = ", round(coef(summary(lm_mean))['Effective_DE_log_mean','Pr(>|t|)'], abs(floor(log10(coef(summary(lm_mean))['Effective_DE_log_mean','Pr(>|t|)'])-nsmall))),
             "\nAdj. R^2 = ", format(round(r.squaredGLMM(lm_mean)[,'R2c'], nsmall), nsmall), sep="")
ggplot(df_avg, aes(x=Effective_DE_log_mean, y=Focalization_log_mean)) +
  geom_point(aes(color=Family), size=5, shape=18) + 
  geom_errorbar(aes(ymin = Focalization_log_mean - Focalization_log_se, ymax = Focalization_log_mean + Focalization_log_se), width=0.01, size=0.5) + 
  geom_errorbarh(aes(xmin = Effective_DE_log_mean - Effective_DE_log_se, xmax = Effective_DE_log_mean + Effective_DE_log_se), height=0.02, size=0.5) + 
  geom_smooth(method='lm') + 
  geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
  geom_text(aes(x=-0.55, y=-0.375, label=lab), hjust=0, size=geom.text.size) + 
  xlab('log10(Effective DE)') + ylab('log10(FE)') + 
  theme_bw() + 
  theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
fileName <- file.path(resdir, 'Family_avg_loglog')
savefig(fileName, 300, 10, 4, "in", "both")


# by geographical region
df_avg <- read.csv(file.path(resdir, 'data_avg_geo.csv'))
lm_mean <- lm(Focalization_log_mean ~ Effective_DE_log_mean, data = df_avg)
summary(lm_mean)
# lm(formula = Focalization_log_mean ~ Effective_DE_log_mean, data = df_avg)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.080345 -0.031656  0.009923  0.026740  0.075122 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -1.1661     0.1463  -7.972 0.000208 ***
#   Effective_DE_log_mean  -0.9591     0.2499  -3.838 0.008574 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05939 on 6 degrees of freedom
# Multiple R-squared:  0.7106,	Adjusted R-squared:  0.6624 
# F-statistic: 14.73 on 1 and 6 DF,  p-value: 0.008574

r.squaredGLMM(lm_mean)
# R2m       R2c
# 0.6779237 0.6779237

# log ~ log linear
fit_intercept <- coef(summary(lm_mean))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_mean))[,'Estimate'][2]
lab <- paste("Slope K = ", format(round(fit_slope, nsmall), nsmall), 
             "\nIntercept = ", format(round(fit_intercept, nsmall), nsmall), 
             "\np-value = ", round(coef(summary(lm_mean))['Effective_DE_log_mean','Pr(>|t|)'], abs(floor(log10(coef(summary(lm_mean))['Effective_DE_log_mean','Pr(>|t|)'])-nsmall))),
             "\nAdj. R^2 = ", format(round(r.squaredGLMM(lm_mean)[,'R2c'], nsmall), nsmall), sep="")
ggplot(df_avg, aes(x=Effective_DE_log_mean, y=Focalization_log_mean)) +
  geom_point(aes(color=Region), size=4, shape=18) + 
  geom_errorbar(aes(ymin = Focalization_log_mean - Focalization_log_se, ymax = Focalization_log_mean + Focalization_log_se), width=0.01, size=0.5) + 
  geom_errorbarh(aes(xmin = Effective_DE_log_mean - Effective_DE_log_se, xmax = Effective_DE_log_mean + Effective_DE_log_se), height=0.02, size=0.5) + 
  geom_smooth(method='lm', se=TRUE) + 
  geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
  geom_text(aes(x=-0.5, y=-0.55, label=lab), hjust=0, size=geom.text.size) + 
  xlab('log10(Effective DE)') + ylab('log10(FE)') + 
  theme_bw() + 
  theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
fileName <- file.path(resdir, 'Geo_avg_loglog')
savefig(fileName, 300, 10, 4, "in", "both")



############################################
# separate group: by Language Family (Genus)
############################################
resdir <- file.path(direct, 'Family'); dir.create(resdir, showWarnings = FALSE)

# read in data
df <- read.csv(file.path(direct, 'Dataset_S1_Processed.csv'))

familyList <- unique(df$Family)
resDF <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(resDF) <- c('Family', 'NoSample', 'SlopeK', 'Intercept', 'p-value', 'adjR2')
for(family in familyList){
  subdf <- df[df$Family==family, ]
  # calculate column values
  nosample <- nrow(subdf)
  if(nosample==1){
    fit_intercept <- NA
    fit_slope <- NA
    pvalue <- NA
    adjr2 <- NA
  }else{
    # train lme model
    result <- try(lme_family <- lmer(log10(Focalization) ~ log10(Effective_DE) 
                                     + (1|Geographic_Region:ISO693_Code) + (1|ISO693_Code), data = subdf, control=lmerControl(optCtrl=list(maxfun=20000))), silent = TRUE)
    if(class(result)[1]=='try-error'){ 
      # all samples are from the same Geographic_Region, use lm model instead
      result <- try(lm_family <- lm(log10(Focalization) ~ log10(Effective_DE), data = subdf), silent = TRUE)
      if(class(result)[1]=='try-error'){
        fit_intercept <- NA; fit_slope <- NA
        pvalue <- NA
        adjr2 <- NA
      }else{
        # calculate column values  
        fit_intercept <- coef(summary(lm_family))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_family))[,'Estimate'][2]
        pvalue <- coef(summary(lm_family))['log10(Effective_DE)','Pr(>|t|)']
        adjr2 <- r.squaredGLMM(lm_family)[,'R2c']
      }
    }else{
      # calculate column values  
      fit_intercept <- coef(summary(lme_family))[, 'Estimate'][1]; fit_slope <- coef(summary(lme_family))[,'Estimate'][2]
      pvalue <- coef(summary(lme_family))['log10(Effective_DE)','Pr(>|t|)']
      adjr2 <- r.squaredGLMM(lme_family)[,'R2c']
    }
  }
  # add to resDF
  resDF <- rbind(resDF, data.frame(Family=family, NoSample=nosample, SlopeK=fit_slope, Intercept=fit_intercept, p_value=pvalue, adjR2=adjr2))
}
# decreasing sorting by NoSample
resDF <- resDF[order(-resDF$NoSample),]
# save the results
write.csv(resDF, file.path(resdir, 'Family_lme_res.csv'), row.names = FALSE)

# get first 9 family data
resDF <- read.csv(file.path(resdir, 'Family_lme_res.csv'))
resDF <- resDF[c(1:5,7:10),]
familyList <- resDF$Family
for(family in familyList){
  subdf <- df[df$Family==family,]
  nosample <- nrow(subdf)
  # train lme model
  result <- try(lme_family <- lmer(log10(Focalization) ~ log10(Effective_DE) 
                                   + (1|Geographic_Region:ISO693_Code) + (1|ISO693_Code), data = subdf, control=lmerControl(optCtrl=list(maxfun=20000))), silent = TRUE)
  # log ~ log linear
  if(class(result)[1]=='try-error'){ 
    # use lm 
    lm_family <- lm(log10(Focalization) ~ log10(Effective_DE), data = subdf)
    
    fit_intercept <- coef(summary(lm_family))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_family))[,'Estimate'][2]
    ggplot(subdf, aes(x=log10(Effective_DE), y=log10(Focalization))) +
      geom_point(size=2, color='red') + 
      #geom_point(aes(color=Geographic_Region), size=2, color='red') + 
      geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
      xlab('log10(Effective DE)') + ylab('log10(FE)') + 
      ggtitle(paste('Family: ', family, "; No. Samples: ", nosample, '; Slope K: ', format(round(fit_slope, nsmall), nsmall), '; Adj R^2: ',  format(round(r.squaredGLMM(lm_family)[,'R2c'], nsmall), nsmall), sep="")) + 
      theme_bw() + 
      theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
            panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
    fileName <- file.path(resdir, paste('Family_lm_', family, '_loglog', sep=""))
    savefig(fileName, 300, 10, 4, "in", "both")
  }else{
    fit_intercept <- coef(summary(lme_family))[, 'Estimate'][1]; fit_slope <- coef(summary(lme_family))[,'Estimate'][2]
    ggplot(subdf, aes(x=log10(Effective_DE), y=log10(Focalization))) +
      geom_point(size=2, color='red') + 
      #geom_point(aes(color=Geographic_Region), size=2, color='red') + 
      geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
      xlab('log10(Effective DE)') + ylab('log10(FE)') + 
      ggtitle(paste('Family: ', family, "; No. Samples: ", nosample, '; Slope K: ', format(round(fit_slope, nsmall), nsmall), '; Adj R^2: ',  format(round(r.squaredGLMM(lme_family)[,'R2c'], nsmall), nsmall), sep="")) + 
      theme_bw() + 
      theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
            panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
    fileName <- file.path(resdir, paste('Family_lme_', family, '_loglog', sep=""))
    savefig(fileName, 300, 10, 4, "in", "both")
  }
}


###################################### 
# separate group: by Geographic_Region
######################################
resdir <- file.path(direct, 'Region'); dir.create(resdir, showWarnings = FALSE)

# read in data
df <- read.csv(file.path(direct, 'Dataset_S1_Processed.csv'))

geoList <- unique(df$Geographic_Region)
resDF <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(resDF) <- c('Geographic_Region', 'NoSample', 'SlopeK', 'Intercept', 'p-value', 'adjR2')
for(geo in geoList){
  subdf <- df[df$Geographic_Region==geo, ]
  # calculate column values
  nosample <- nrow(subdf)
  if(nosample==1){
    fit_intercept <- NA
    fit_slope <- NA
    pvalue <- NA
    adjr2 <- NA
  }else{
    # train lme model
    result <- try(lme_geo <- lmer(log10(Focalization) ~ log10(Effective_DE) 
                                  + (1|ISO693_Code) + (1|Family), data = subdf, control=lmerControl(optCtrl=list(maxfun=20000))), silent = TRUE)
    if(class(result)[1]=='try-error'){ 
      # all samples are from the same Geographic_Region, use lm model instead
      result <- try(lm_geo <- lm(log10(Focalization) ~ log10(Effective_DE), data = subdf), silent = TRUE)
      if(class(result)[1]=='try-error'){
        fit_intercept <- NA; fit_slope <- NA
        pvalue <- NA
        adjr2 <- NA
      }else{
        # calculate column values  
        fit_intercept <- coef(summary(lm_geo))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_geo))[,'Estimate'][2]
        pvalue <- coef(summary(lm_geo))['log10(Effective_DE)','Pr(>|t|)']
        adjr2 <- r.squaredGLMM(lm_geo)[,'R2c']
      }
    }else{
      # calculate column values  
      fit_intercept <- coef(summary(lme_geo))[, 'Estimate'][1]; fit_slope <- coef(summary(lme_geo))[,'Estimate'][2]
      pvalue <- coef(summary(lme_geo))['log10(Effective_DE)','Pr(>|t|)']
      adjr2 <- r.squaredGLMM(lme_geo)[,'R2c']
    }
  }
  # add to resDF
  resDF <- rbind(resDF, data.frame(Geographic_Region=geo, NoSample=nosample, SlopeK=fit_slope, Intercept=fit_intercept, p_value=pvalue, adjR2=adjr2))
}
# save the results
write.csv(resDF, file.path(resdir, 'Geo_lme_res.csv'), row.names = FALSE)

# draw figure
resDF <- read.csv(file.path(resdir, 'Geo_lme_res.csv'))
geoList <- resDF$Geographic_Region
for(geo in geoList){
  subdf <- df[df$Geographic_Region==geo,]
  nosample <- nrow(subdf)
  # train lme model
  result <- try(lme_geo <- lmer(log10(Focalization) ~ log10(Effective_DE) 
                                + (1|ISO693_Code) + (1|Family), data = subdf, control=lmerControl(optCtrl=list(maxfun=20000))), silent = TRUE)
  # log ~ log linear
  if(class(result)[1]=='try-error'){ 
    # use lm 
    lm_geo <- lm(log10(Focalization) ~ log10(Effective_DE), data = subdf)
    
    fit_intercept <- coef(summary(lm_geo))[, 'Estimate'][1]; fit_slope <- coef(summary(lm_geo))[,'Estimate'][2]
    ggplot(subdf, aes(x=log10(Effective_DE), y=log10(Focalization))) +
      geom_point(size=2, color='red') +
      #geom_point(aes(color=Family), size=1, color='red') + 
      geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
      #geom_text(aes(x=-0.4, y=-0.4, label=lab), hjust=0, size=geom.text.size) + 
      xlab('log10(Effective DE)') + ylab('log10(FE)') + 
      ggtitle(paste('Region: ', geo, '; No. Samples: ', nosample, '; Slope K: ', format(round(fit_slope, nsmall), nsmall), '; Adj R^2: ',  format(round(r.squaredGLMM(lm_geo)[,'R2c'], nsmall), nsmall), sep="")) + 
      theme_bw() + 
      theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
            panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
    fileName <- file.path(resdir, paste('Geo_', geo, '_loglog', sep=""))
    savefig(fileName, 300, 10, 4, "in", "both")
  }else{
    fit_intercept <- coef(summary(lme_geo))[, 'Estimate'][1]; fit_slope <- coef(summary(lme_geo))[,'Estimate'][2]
    ggplot(subdf, aes(x=log10(Effective_DE), y=log10(Focalization))) +
      geom_point(size=2, color='red') + 
      #geom_point(aes(color=Family), size=1, color='red') + 
      geom_abline(aes(slope=fit_slope, intercept=fit_intercept), size=1) + 
      #geom_text(aes(x=-0.4, y=-0.4, label=lab), hjust=0, size=geom.text.size) + 
      xlab('log10(Effective DE)') + ylab('log10(FE)') + 
      ggtitle(paste('Region: ', geo, '; No. Samples: ', nosample, '; Slope K: ', format(round(fit_slope, nsmall), nsmall), '; Adj R^2: ',  format(round(r.squaredGLMM(lme_geo)[,'R2c'], nsmall), nsmall), sep="")) + 
      theme_bw() + 
      theme(plot.title=element_text(size=theme.size), axis.text=element_text(size=theme.size), text=element_text(size=theme.size), 
            panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 1))
    fileName <- file.path(resdir, paste('Geo_', geo, '_loglog', sep=""))
    savefig(fileName, 300, 10, 4, "in", "both")
  }
}
