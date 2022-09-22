##################################
#### CRP level analysis
##################################
# laod libraries ----

if(!require("lmerTest")){install.packages("lmerTest")}
require("lmerTest") 

if(!require("groupdata2")){install.packages("groupdata2")}
require("groupdata2") # model crossvalidation

if(!require("performance")){install.packages("performance")}
require("performance") # GLMModel performance comparison

if(!require("bayestestR")){install.packages("bayestestR")}
require("bayestestR") 

if(!require("PerformanceAnalytics")){install.packages("PerformanceAnalytics")}
require("PerformanceAnalytics") 

if(!require("tidyverse")){install.packages("tidyverse")}
require("tidyverse") 

if(!require("gapminder")){install.packages("gapminder")}
require("gapminder") 

if(!require("ggplot2")){install.packages("ggplot2")}
require("ggplot2") 

if(!require("plotly")){install.packages("plotly")}
require("plotly") 

if(!require("lme4")){install.packages("lme4")}
require("lme4") 

if(!require("car")){install.packages("car")}
require("car") 

if(!require("cowplot")){install.packages("cowplot")}
require("cowplot") 

if(!require("plotROC")){install.packages("plotROC")}
require("plotROC") 

if(!require("cowplot")){install.packages("cowplot")}
require("cowplot") 

if(!require("rcompanion")){install.packages("rcompanion")}
require("rcompanion") 

if(!require("ggpubr")){install.packages("ggpubr")}
require("ggpubr") 

if(!require("reportROC")){install.packages("reportROC")}
require("reportROC") 

if(!require("GGally")){install.packages("GGally")}
require("GGally") 



rm(list = ls())

# set directory
setwd("~/Margaret")
# read CRP file  made by Jiayi
CRP<-read.csv(file="CRP_Plates1-7_olink.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, na.strings = c('',' ', 'NA'))
olink_data_bc<-read.csv(file="olink_data_bc.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, na.strings = c('',' ', 'NA'))
parti<-read.csv(file="~/File_for_theAnalysis/participants_2020-05-04.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, na.strings = c('',' ', 'NA'))

###### Check for batch effect ##############################################
CRPnames<-"CRP"
colnames(CRP)

# the good thing is that CRP is not associated with ELISA plate
temp_formula <- as.formula(paste(CRPnames, " ~ as.factor(ELISA.Plate)"))
tempbc <-glm(formula = temp_formula, data = CRP)
summary(tempbc)
Anovtempbc<-Anova(tempbc)

### therefore, we take mean value for each subject among ELISA.plate
# select "consortium_id",

#"redcap_event_name" and "CRP" columns
colnames(CRP)
CRP<-CRP[,c(3,4,8,9,5)]
head(CRP)
MeanCRP <- aggregate (CRP ~ X_consortium_id+X_redcap_event_name , CRP, mean)
head(MeanCRP)
names(MeanCRP)[3] <- "MeanCRP"
CRP <- merge(CRP, MeanCRP, by=c("X_consortium_id", "X_redcap_event_name"))
head(CRP)


# read the anti-TNF file sent by Jiayi (new entry)
setwd("~/Margaret/File_for_theAnalysis")
AntiChange<-read.csv(file="anti_tnf_current.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, na.strings = c('',' ', 'NA'))
NROW(AntiChange)
table(AntiChange$antitnf_current, useNA = "always")

# replace space with underscore
AntiChange$redcap_event_name <- gsub(" ", "_", AntiChange$redcap_event_name)

# replace space with underscore
head(olink_data_bc)
CRP$X_redcap_event_name <- gsub(" ", "_", CRP$X_redcap_event_name)
head(CRP)

# merge by columns of interest from the olink batch corrected file (CRP and olink_data_bc)
CRPMerge<-merge(olink_data_bc, CRP, by.x=c("consortium_id","redcap_event_name"), by.y=c("X_consortium_id","X_redcap_event_name"))
NROW(CRPMerge)#210

head(CRPMerge)
colnames(CRPMerge)

# merge for antiTND
CRPMerge<-merge(CRPMerge, AntiChange, by.x=c("consortium_id","redcap_event_name"), by.y=c("consortium_id","redcap_event_name"))
NROW(CRPMerge)#210
head(CRPMerge)

# merge file with participant file includes sex
CRPMerge<-merge(CRPMerge, parti[,c(1,9)], by.x="consortium_id", by.y="consortium_id")
NROW(CRPMerge)
head(CRPMerge)

# remove duplicated rows 
CRPMerge<-CRPMerge[!duplicated(CRPMerge[ , c("consortium_id","redcap_event_name")]),]
NROW(CRPMerge)
head(CRPMerge)

CRPMerge$antitnf_current[CRPMerge$antitnf_current == "TRUE"]<-1
CRPMerge$antitnf_current[CRPMerge$antitnf_current == "FALSE"]<-0
NROW(CRPMerge)

# remove NAs from rutgeerts score and Inf values introduced by log transform
CRPMerge<-CRPMerge[which(CRPMerge$rutgeerts_score_bin !="NA"),] 
NROW(CRPMerge)# 190


# remove NAs from anti_TNF columns
CRPMerge<-CRPMerge[which(CRPMerge$antitnf_current !="NA"),] 
NROW(CRPMerge)# 186
head(CRPMerge)

### saving the data for CrossValidation
#write.csv(CRPMerge, file="~/Documents/Mount_Sinai_Hospital_documents/Margaret/CRP_analysis_main/CRP_crossValidation/CrossValidation_mainCrossVal2.csv", row.names=FALSE, quote=FALSE)


################################################################################
### Quality control of the data
################################################################################
# log10 tranform the CRP data to make the data normal
CRPMerge$MeanCRP<-log(CRPMerge$MeanCRP+1)

#plot density
ggplot(CRPMerge, aes(x=MeanCRP, fill=as.factor(rutgeerts_score_bin))) +
  geom_histogram(binwidth=1, alpha=1, position="dodge")

# boxplot---
ggbetweenstats(CRPMerge, rutgeerts_score_bin, MeanCRP, outlier.tagging = TRUE) +
  xlab("Rugeerts score binarized") +
  ylab("CRP levels")


###### TEST OF NORMALITY ###########
NROW(CRPMerge)#186

# test of normality Shapiro.test

# CRP
# CRP follow a non-normal distribution pvalue<0.05
shapiro.test(CRPMerge$MeanCRP) 
#Shapiro-Wilk normality test
#data:  CRPMerge$CRP
#W = 0.9365, p-value = 9.592e-08

# cxcl9
# cxcl9 does not follow a normal distribution pvalue<0.05
shapiro.test(CRPMerge$cxcl9) 
#Shapiro-Wilk normality test
#data:  CRPMerge$cxcl9
#W = 0.98273, p-value = 0.0127

# mmp1
# mmp1 does not follow a normal distribution pvalue<0.05
shapiro.test(CRPMerge$mmp1) 
#Shapiro-Wilk normality test
#data:  CRPMerge$mmp1
#W = 0.98168, p-value = 0.009513

# il6
# il6 does not follow a normal distribution pvalue<0.05
shapiro.test(CRPMerge$il6) 
#Shapiro-Wilk normality test
#data:  CRPMerge$il6
#W = 0.91271, p-value = 1.427e-09

#il5
#non normal distribution
shapiro.test(CRPMerge$il5) 
#data:  CRPMerge$il5
#W = 0.74171, p-value < 2.2e-16

#ST1A1
# follow non-normal distribution
shapiro.test(CRPMerge$st1a1) 
#data:  CRPMerge$st1a1
#W = 0.97931, p-value = 0.004328

### QQ plot 
# CRP
q1<-ggqqplot(CRPMerge$MeanCRP, ylab = "CRP", xlab = "QQ plot")
# cxcl9
q2<-ggqqplot(CRPMerge$cxcl9, ylab = "cxcl9", xlab = "QQ plot")

plot_grid(q1, q2, nrow = 2, label_size = 3 , hjust = 0.5)


# plot all correlations ----
# here we would focus on CRP, CXCL9, MMP1 and il6, il5 and st1a1
# convert integer variables
CRPMerge$MeanCRP<-as.numeric(CRPMerge$MeanCRP)
CRPMerge$cxcl9<-as.numeric(CRPMerge$cxcl9)
CRPMerge$il6<-as.numeric(CRPMerge$il6)
CRPMerge$mmp1<-as.numeric(CRPMerge$mmp1)
CRPMerge$il5<-as.numeric(CRPMerge$il5)
CRPMerge$st1a1<-as.numeric(CRPMerge$st1a1)
CRPMerge$rutgeerts_score_bin<-as.factor(CRPMerge$rutgeerts_score_bin)
CRPMerge$sample_age<-as.numeric(CRPMerge$sample_age)
CRPMerge$anti_tnf<-as.factor(CRPMerge$anti_tnf)
CRPMerge$sex<-as.factor(CRPMerge$sex)

DataCRPCor<-subset(CRPMerge, select=c("cxcl9","il6","mmp1","MeanCRP","il5","st1a1"))

chart.Correlation <- function (R, histogram = TRUE, labels = TRUE, method = c("pearson", "kendall", 
                                                               "spearman"), ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = "pearson", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 1/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***","**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.75, 0.75, Signif, cex = cex, col = 2)
  }
  
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "white", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "dark red",lwd = 1.5)
    rug(x)
  }
  
  if (histogram) 
    pairs(x, gap = 0,  lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel, ...)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
}



chart.Correlation(DataCRPCor, 
                  method="spearman",
                  histogram=TRUE,
                  pch = 19)


#####################################################################
##### DEFINE Cutt off value for CRP and its significance ############
#####################################################################

# Plot CRP in two groups (rutgeerts score 0 and 1) ----
D=CRPMerge
D$rutgeerName<-NA
D$rutgeerName[D$rutgeerts_score_bin=="0"]<-"Endsocopic remission\n(i0+i1)"
D$rutgeerName[D$rutgeerts_score_bin=="1"]<-"Endsocopic recurrence\n(i2+i3+i4)"
my_comparisons <- list( c("Endsocopic remission\n(i0+i1)", "Endsocopic recurrence\n(i2+i3+i4)"))
plot1<-ggboxplot(D, x = "rutgeerName", y = "MeanCRP",
                 fill="grey", color = "black", palette = "jco", add = "jitter", width=0.2, notch = FALSE)+stat_compare_means(comparisons =my_comparisons)

plot1+labs(title="", x ="", y = "CRP (mg/l)")

# Stratified sampling
################################################################################
table(CRPMerge$rutgeerts_score_bin)
#  0   1 
# 120  66 

CRPMerge<-stratified(CRPMerge, "rutgeerts_score_bin", 0.8, replace = TRUE)
NROW(CRPMerge)
head(CRPMerge)
table(CRPMerge$rutgeerts_score_bin)
# 0  1 
# 96 53 


################################################
# logistic regression with and without anti_tnf
###############################################

# compare models to find the best fit without random effect (GLM)
modelCRP.1=glm(rutgeerts_score_bin ~ cxcl9, na.action =na.omit, 
            data=CRPMerge, family=binomial(link = "logit"))
modelCRP.2=glm(rutgeerts_score_bin ~ cxcl9+antitnf_current, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.3=glm(rutgeerts_score_bin ~ MeanCRP, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.4=glm(rutgeerts_score_bin ~ MeanCRP+antitnf_current, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.5=glm(rutgeerts_score_bin ~ MeanCRP+cxcl9, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.6=glm(rutgeerts_score_bin ~ cxcl9+MeanCRP+antitnf_current, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.7=glm(rutgeerts_score_bin ~ cxcl9 + MeanCRP+ mmp1, na.action =na.omit,
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.8=glm(rutgeerts_score_bin ~ cxcl9 + MeanCRP + mmp1+antitnf_current, na.action =na.omit, 
               data=CRPMerge, family=binomial(link = "logit"))
modelCRP.9=glm(rutgeerts_score_bin ~  cxcl9 +MeanCRP+mmp1+il5+st1a1, na.action =na.omit,
                data=CRPMerge, family=binomial(link = "logit"))
modelCRP.10=glm(rutgeerts_score_bin ~  cxcl9 +MeanCRP+mmp1+il5+st1a1+antitnf_current, na.action =na.omit,
                data=CRPMerge, family=binomial(link = "logit"))


compareGLM(modelCRP.1, modelCRP.2, modelCRP.3, modelCRP.4,modelCRP.5, modelCRP.6, modelCRP.7,modelCRP.8,modelCRP.9, modelCRP.10)
#compare and rank GLM models
compModel<-compare_performance(modelCRP.1, modelCRP.2, modelCRP.3, modelCRP.4,modelCRP.5, modelCRP.6, modelCRP.7,modelCRP.8, modelCRP.9, modelCRP.10, rank = TRUE)
#knitr::kable(compModel, digits = 2)

# Models 1, 5 and 10 are the best fits 
###############################################################################################
###### Plot ROC for CRP and CXCL9 il5 and st1a1 , mmp1 without anti_TNF #######################
# define names
cytokines<-c("CXCL9", "CRP", "CXCL9/CRP", "CXCL9/CRP/IL5/ST1A1/MMP1")
cytokines<-as.list(cytokines)

# ROC function
ROCResult<-function(x,df){
  # plot size
  par(pty= "s")
  # predict ROC
  df_Predict<-roc(CRPMerge$rutgeerts_score_bin ~ df$fitted.values, legacy.axes = TRUE, 
                     col="black", lwd=2, main="", lty=2, cex.main=1)
  # AUC CI
  AUC<-ci.auc( df_Predict, method="bootstrap")
  print(AUC)
  # plot ROC
  plot(df_Predict, col="black", lwd=2, main=paste("ROC Curve:", cytokines[[x]]), cex.main=0.8, 
       legacy.axes = TRUE, print.auc=TRUE) #cex=0.8)
  #legend("bottom", label=paste(""))
  
  # Youden TEst
  Ydentest<-coords(df_Predict, x="best", input="threshold", best.method="youden", transpose = TRUE)
  #print coordinates
  coords(df_Predict, "best", best.method="youden", ret=c("threshold", "specificity", "sensitivity", "ppv",
                                                            "npv", "recall"), transpose = FALSE)
}

CRPMerge$MeanCRP


# ROC CURVE WITHOUT ANTI_TNF ----
#####################################################################
### CXCL9
ROCResult(1, modelCRP.1)
#95% CI: 0.6-0.7596 (2000 stratified bootstrap replicates)
#threshold specificity sensitivity      ppv   npv    recall
#threshold  0.334502        0.55   0.7878788 0.490566 0.825 0.7878788

### CRP
ROCResult(2, modelCRP.3)
#95% CI: 0.5485-0.7203 (2000 stratified bootstrap replicates)
#threshold specificity sensitivity ppv  npv    recall
#threshold 0.3340158       0.675   0.5909091 0.5 0.75 0.5909091

### CRP/CXCL9
ROCResult(3, modelCRP.5)
#95% CI: 0.616-0.7676 (2000 stratified bootstrap replicates)
#threshold specificity sensitivity       ppv       npv    recall
#threshold 0.3234209   0.5583333   0.7878788 0.4952381 0.8271605 0.7878788

#### CXCL9/CRP/IL5/ST1A1/MMP1
ROCResult(4, modelCRP.10)
#95% CI: 0.6624-0.8102 (2000 stratified bootstrap replicates)
#threshold specificity sensitivity       ppv       npv    recall
#threshold 0.4002163         0.8   0.5909091 0.6190476 0.7804878 0.5909091

########### FInal plot ALL SAMPLES GLM without anti_TNF ####################################################
rocmodel1_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.1$fitted.values, legacy.axes = TRUE, 
                col="black", lwd=2, main="", lty=2, cex.main=1)
rocmodel3_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.3$fitted.values, legacy.axes = TRUE, 
                col="black", lwd=2, main="", lty=2, cex.main=1)
rocmodel10_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.10$fitted.values, legacy.axes = TRUE, 
                col="black", lwd=2, main="", lty=2, cex.main=1)

plot(rocmodel1_glm, col="black", lwd=2, lty=1, main = "All Samples n=186 \n ROC curve GLM model without anti-TNF", legacy.axes = TRUE,cex.main=1)
plot(rocmodel3_glm, col="darkgreen", lwd=2, lty=2, main="", add=TRUE, legacy.axes = TRUE)
plot(rocmodel10_glm, col="black", lwd=2, lty=3, main="", add=TRUE, legacy.axes = TRUE)
legend("bottom", legend=c("CXCL9 (AUC:0.683)", "CRP (AUC:0.644)\n", "CXCL9/CRP/MMP1/IL5/ST1A1\n (AUC:0.741)\n"), col=c("black",  "darkgreen", "black"), lty=c(1,2,3), lwd=2, cex=0.8,box.lty=2, bty="n")

########### FInal plot STRATIFIED SAMPLING GLM without anti_TNF ####################################################
rocmodel1_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.1$fitted.values, legacy.axes = TRUE, 
                   col="black", lwd=2, main="", lty=2, cex.main=1)
rocmodel3_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.3$fitted.values, legacy.axes = TRUE, 
                   col="black", lwd=2, main="", lty=2, cex.main=1)
rocmodel10_glm<-roc(CRPMerge$rutgeerts_score_bin ~ modelCRP.10$fitted.values, legacy.axes = TRUE, 
                    col="black", lwd=2, main="", lty=2, cex.main=1)

plot(rocmodel1_glm, col="black", lwd=2, lty=1, main = "Stratified Sampling n=149 (0=96; 1=53) \n ROC curve GLM model without anti-TNF", legacy.axes = TRUE,cex.main=1)
plot(rocmodel3_glm, col="darkgreen", lwd=2, lty=2, main="", add=TRUE, legacy.axes = TRUE)
plot(rocmodel10_glm, col="black", lwd=2, lty=3, main="", add=TRUE, legacy.axes = TRUE)
legend("bottom", legend=c("CXCL9 (AUC:0.724)", "CRP (AUC:0.679)\n", "CXCL9/CRP/MMP1/IL5/ST1A1\n (AUC:0.801)\n"), col=c("black",  "darkgreen", "black"), lty=c(1,2,3), lwd=2, cex=0.8,box.lty=2, bty="n")

