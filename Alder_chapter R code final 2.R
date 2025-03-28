#Load all packages
library(Hmisc)
library(dplyr)
library(knitr)
library(lme4)
library(rrBLUP)
library(tidyverse)
library(adegenet)
library(ggplot2)
library(vcfR)
library(PerformanceAnalytics)
library(AICcmodavg)
library(emmeans)
library(graph4lg)
library(gplots)
library(pheatmap)

## Check data distribution and descriptive analysis
#The original data comprise of height collected from planted year-2008, 2010, 
#2015, 2017, 2020, DBH collected from 2015, 2017 and 2020. 
#Both height and DBH were converted to relative height growth and relative DBH growth 
#(This was completed using Excel).
corcovety<-read.csv("Corcovety_GS_edited.csv",header = TRUE)

#check structure of data
str(corcovety)
#need to convert Rep from integer to factor;DBH12 to numeric
corcovety$Block<-as.factor(corcovety$Block)
corcovety$Rep<-as.factor(corcovety$Rep)
corcovety$Family<-as.factor(corcovety$Family)


#double check the family and rep
table(corcovety$Family)
table(corcovety$Rep)

summary(corcovety)
library(dplyr)
ht0<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_HT0 = mean(HT0, na.rm=TRUE))
ht2<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_HT2 = mean(HT2, na.rm=TRUE))
ht7<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_HT7 = mean(HT7, na.rm=TRUE))
ht9<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_HT9 = mean(HT9, na.rm=TRUE))
ht12<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_HT.12 = mean(HT.12, na.rm=TRUE))
dbh7<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_DBH7 = mean(DBH7, na.rm=TRUE))
dbh9<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_DBH9 = mean(DBH9, na.rm=TRUE))
dbh12<-corcovety %>%
  group_by(Family) %>%
  summarise(mean_DBH12 = mean(DBH12, na.rm=TRUE))
overall_performance<-cbind(ht0,ht2,ht7,ht9,ht12,dbh7,dbh9,dbh12)
write.csv(overall_performance,"overall_performance.csv")


#plot the relative growth control vs plus-trees
library(tidyr)
library(dplyr)
library(ggplot2)
corcovety_2<-corcovety[,-c(1,2,3,5,6,15,16,17,18,19,20,21)]
h<-corcovety_2[,-c(7,8,9)]
long_corcovety2<- corcovety_2 %>%
  pivot_longer(cols = contains("H"),
    names_to = "trait",
               values_to = "value")

ggplot(long_corcovety2,aes(x = trait, y = value))+
     labs(title = "Overall performance of each family in the progeny trial", x = "Trait", y = "Growth")+
     geom_jitter(data = long_corcovety2,aes(x = trait, y = value, color = ifelse(Family == "control", "control","nocontrol")))+
  geom_boxplot(alpha=0.7)+
     scale_color_manual(values = c("control"="red","nocontrol"="blue"))


long_corcovety2 %>% mutate(trait2 = case_when(trait == "HT0" ~ "Tree Height - 2008",
                                                trait == "HT2" ~ "Tree Height - 2010",
                                                trait == "HT7" ~ "Tree Height - 2015",
                                                trait == "HT9" ~ "Tree Height - 2017",
                                                trait == "HT.12" ~ "Tree Height - 2020",
                                                trait == "DBH7" ~ "DBH - 2015",
                                                trait == "DBH9" ~ "DBH - 2017",
                                                trait == "DBH12" ~ "DBH - 2020")) %>%
  ggplot(aes(x = trait2, y = value))+
  labs(title = "Overall performance of each family in the progeny trial", x = "Trait", y = "Growth")+
  geom_jitter(aes(x = trait2, y = value, color = ifelse(Family == "control", "control","nocontrol")))+
  theme(legend.title=element_blank()) +
  geom_boxplot(alpha=0.7)+
  scale_color_manual(values = c("control"="red","nocontrol"="blue"))


##separate height and dbh
corcovety_2<-corcovety[,-c(1,2,3,5,6,15,16,17,18,19,20,21)]
h<-corcovety_2[,-c(7,8,9)]
h_long<- h %>%
  pivot_longer(cols = contains("HT"),
               names_to = "trait",
               values_to = "value")

#ggplot(h_long,aes(x = trait, y = value))+
  #labs(title = "Overall performance of each family in the progeny trial", x = "Trait", y = "Growth")+
  #geom_jitter(data = h_long,aes(x = trait, y = value, color = ifelse(Family == "control", "Control","Selected Plus-trees")))+
  #geom_boxplot(alpha=0.7)+
  #scale_color_manual(values = c("Control"="red","Selected Plus-trees"="blue"))


hplot<-h_long %>% mutate(trait2 = case_when(trait == "HT0" ~ "2008",
                                              trait == "HT2" ~ "2010",
                                              trait == "HT7" ~ "2015",
                                              trait == "HT9" ~ "2017",
                                              trait == "HT.12" ~ "2020")) %>%
  ggplot(aes(x = trait2, y = value))+
  labs(title = "", x = "Year", y = "Tree height (m)")+
  geom_jitter(aes(x = trait2, y = value, color = ifelse(Family == "control", "Control","Selected Plus-trees")))+
  theme(legend.title=element_blank()) +
  geom_boxplot(alpha=0.7)+
  scale_color_manual(values = c("Control"="red","Selected Plus-trees"="blue"))+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank())


dbh<-corcovety_2[,-c(2,3,4,5,6)]
dbh_long<- dbh %>%
  pivot_longer(cols = contains("DBH"),
               names_to = "trait",
               values_to = "value")

dbhplot<-dbh_long %>% mutate(trait2 = case_when(trait == "DBH7" ~ "2015",
                                     trait == "DBH9" ~ "2017",
                                     trait == "DBH12" ~ "2020")) %>%
  ggplot(aes(x = trait2, y = value))+
  labs(title = "", x = "Year", y = "DBH (cm)")+
  geom_jitter(aes(x = trait2, y = value, color = ifelse(Family == "control", "Control","Selected Plus-trees")))+
  geom_boxplot(alpha=0.7)+
  scale_color_manual(values = c("Control"="red","Selected Plus-trees"="blue"))+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "right",
        legend.title=element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5)))


library(cowplot)
plot_grid(hplot,dbhplot,rel_widths = c(1,1),
          labels = c("A","B"))

ggsave("Figure2.png", dpi=300, h = 6, w = 9)



#a total of 89 families (including control) found in the dataset
#family 33/56 has only 9 trees (appear only in 1 plot)
#family 51/90 has 18 trees(appear only in 2 plots)
#control has 18 trees in every plots
#these families needs to be removed as they would cause unbalance of data

#remove family 33/56, 51/90, control as including them will cause the unbalance of data
corcovety_1 <- corcovety %>%
  filter(Family!= 33 & Family!= 56 & Family!= 51 & Family!= 90 & Family!= "control")
#check the new dataset before running analysis
table(corcovety_1$Family)
table(corcovety_1$Rep)
summary(corcovety_1)
#check the distribution of data
alnus_descript<-describe(corcovety_1)
alnus_descript[c("Family", "Rep","Rep", "X2Yrelativeheight","X7Yrelativeheight",
                 "X9Yrelativeheight","X12Yrelativeheight","X9Yrelativedbh",
                 "X12Yrelativedbh")] %>% html()
#relative growth of height and DBH
#corcovety_1 %>% mutate(case_when(colnames == "X2Yrelativeheight" ~ "RG_H(2008-2010)",
                                          #colnames == "X7Yrelativeheight" ~ "RG_H(2008-2015)",
                                         # colnames == "X9Yrelativeheight" ~ "RG_H(2008-2017)",
                                         # colnames == "X12Yrelativeheight" ~ "RG_H(2008-2020)",
                                         # colnames == "X9Yrelativedbh" ~ "RG_DBH(2015-2017)",
                                # colnames == "X12Yrelativedbh" ~ "RG_DBH(2015-2020)"))

chart.Correlation(corcovety_1[,16:21], histogram = TRUE,pch="+")
#These data seems to be normally distributed and some traits are highly correlated
#a total of 84 families retained 

#calculate the sd of each trait
sd(na.omit(corcovety_1$X2Yrelativeheight))
sd(na.omit(corcovety_1$X7Yrelativeheight))
sd(na.omit(corcovety_1$X9Yrelativeheight))
sd(na.omit(corcovety_1$X12Yrelativeheight))
sd(na.omit(corcovety_1$X9Yrelativedbh))
sd(na.omit(corcovety_1$X12Yrelativedbh))

## Heritability calculation: variance component, adjusted number of trees in each plot

#variance component was extracted using Rep as fixed effect, Family and Rep:Family interaction as random effect
#Variance component result see 'Random effects' section in the summary
mod_h2<-lmer(X2Yrelativeheight~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_h2)

mod_h7<-lmer(X7Yrelativeheight~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_h7)

mod_h9<-lmer(X9Yrelativeheight~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_h9)

mod_h12<-lmer(X12Yrelativeheight~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_h12)

mod_dbh9<-lmer(X9Yrelativedbh~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_dbh9)

mod_dbh12<-lmer(X12Yrelativedbh~Rep+(1|Family)+(1|Rep:Family),REML = TRUE,data = corcovety_1)
summary(mod_dbh12)

var_h2<-as.data.frame(VarCorr(mod_h2),comp=c("Variance","Std.Dev."))
var_h2 <- var_h2[, -which(names(var_h2) %in% c("var1","var2","sdcor"))]
colnames(var_h2)<-c("Effect","x2h")
var_h7<-as.data.frame(VarCorr(mod_h7),comp=c("Variance","Std.Dev."))
var_h7 <- var_h7[, -which(names(var_h7) %in% c("var1","var2","sdcor"))]
colnames(var_h7)<-c("Effect","x7h")
var_h9<-as.data.frame(VarCorr(mod_h9),comp=c("Variance","Std.Dev."))
var_h9 <- var_h9[, -which(names(var_h9) %in% c("var1","var2","sdcor"))]
colnames(var_h9)<-c("Effect","x9h")
var_h12<-as.data.frame(VarCorr(mod_h12),comp=c("Variance","Std.Dev."))
var_h12 <- var_h12[, -which(names(var_h12) %in% c("var1","var2","sdcor"))]
colnames(var_h12)<-c("Effect","x12h")
var_dbh9<-as.data.frame(VarCorr(mod_dbh9),comp=c("Variance","Std.Dev."))
var_dbh9 <- var_dbh9[, -which(names(var_dbh9) %in% c("var1","var2","sdcor"))]
colnames(var_dbh9)<-c("Effect","x9dbh")
var_dbh12<-as.data.frame(VarCorr(mod_dbh12),comp=c("Variance","Std.Dev."))
var_dbh12 <- var_dbh12[, -which(names(var_dbh12) %in% c("var1","var2","sdcor"))]
colnames(var_dbh12)<-c("Effect","x12dbh")
var_alltraits<-data.frame(var_h2,var_h7,var_h9,var_h12,var_dbh9,var_dbh12)
var_alltraits<-var_alltraits[, -which(names(var_alltraits) %in% c("Effect.1","Effect.2","Effect.3","Effect.4","Effect.5"))]
write.csv(var_alltraits,"variance_component_alltraits.csv")

#calculate adjusted numbers of trees in a plot for heritability calculation
#count the number of trees of each family in a plot by trait
treenumber_family_plot<-corcovety_1 %>% 
  group_by(Family, Rep) %>%
  summarise(across(.fns = ~sum(!is.na(.))))
treenumber_family_plot_df<-as.data.frame(treenumber_family_plot)
treenumber_family_plot_df<-treenumber_family_plot_df[,-c(3,4,5,6,7,8,9,10,11,12,13,14,15)]
treenumber_family_plot_df

nh <- treenumber_family_plot_df %>%
  mutate(nh_x2h = 1/X2Yrelativeheight,
         nh_x7h = 1/X7Yrelativeheight,
         nh_x9h = 1/X9Yrelativeheight,
         nh_x12h = 1/X12Yrelativeheight,
         nh_x9dbh = 1/X9Yrelativedbh,
         nh_x12dbh = 1/X12Yrelativedbh)  %>% 
  dplyr::select(nh_x2h,nh_x7h,nh_x9h,nh_x12h,nh_x9dbh,nh_x12dbh) %>% 
  colSums() %>%
  enframe(name = "trait", value = "sum_of_nij") %>%
  mutate(nh = 84*3/sum_of_nij)
nh


## Genomic prediction: first step: estimated means;second step: calculate marker effect

#estimated means of each family by trait extracted using Rep and Family as fixed effect, Family and Rep:Family interaction as random effect
mod_h2_emmean<-lmer(X2Yrelativeheight~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)
mod_h7_emmean<-lmer(X7Yrelativeheight~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)
mod_h9_emmean<-lmer(X9Yrelativeheight~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)
mod_h12_emmean<-lmer(X12Yrelativeheight~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)
mod_dbh9_emmean<-lmer(X9Yrelativedbh~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)
mod_dbh12_emmean<-lmer(X12Yrelativedbh~Rep+Family+(1|Rep:Family),REML = TRUE,data = corcovety_1)

emm_X2h_1<-emmeans(mod_h2_emmean,~Family)
emm_X7h_1<-emmeans(mod_h7_emmean,~Family)
emm_X9h_1<-emmeans(mod_h9_emmean,~Family)
emm_X12h_1<-emmeans(mod_h12_emmean,~Family)
emm_X9dbh_1<-emmeans(mod_dbh9_emmean,~Family)
emm_X12dbh_1<-emmeans(mod_dbh12_emmean,~Family)

emm_all<-data.frame(emm_X2h_1,emm_X7h_1,emm_X9h_1,emm_X12h_1,emm_X9dbh_1,emm_X12dbh_1)
emm_all<-emm_all[,-c(3,4,5,6,7,9,10,11,12,13,15,16,17,18,19,21,22,23,24,25,27,28,29,30,31,33,34,35,36)]
emm_all$Family<-sub("^","A",emm_all$Family)
colnames(emm_all)<-c('Genotype','x2h','x7h','x9h','x12h','x9dbh','x12dbh')
write.csv(emm_all,"emm_all.csv")

#further filter 7 more genotypes those genetic data are not available
emm_all_77 <- emm_all %>%
  filter(Genotype!= 'A2' & Genotype!= 'A38' & Genotype!= 'A39' & Genotype!= 'A73' & Genotype!= 'A79' & Genotype!= 'A9' & Genotype!= 'A91')
write.csv(emm_all_77,"emm_all_77.csv")

#a total of 103 genotypes were sequenced
vcf103<-read.csv("vcf103_2_filter_ch.csv")
vcf103<-as.data.frame(vcf103)
rownames(vcf103)<-vcf103[,1]
vcf103<-vcf103[,-1]
dim(vcf103)

#check how much missing data are in the dataset
na_checking<-vcf103 %>% 
  summarise_all(funs(sum(is.na(.))))
na_checking<-as.data.frame(na_checking)
na_checking<59677*0.1
#the max.missing function in a.mat only remove the marker that does not fullfill the setting
#we want to remove the individual tree if too much missing data appeared rather than marker
#need to check before imputation before transposing
#yes the missing data percentage of each tree are less than 10%

#filter the genotypes those phenotypic data are not available
#a total of 77 families were retained for genomic selection
vcf77 <- vcf103[, -which(names(vcf103) %in% c("A103","A104","A105","A107","A113","A116","A117","A118","A133","A134","A137","A139","A141","A143","A146","A149","A158","A160","A163","A167","A177","A24","A33","A51","A56","A85"))]

# transpose data set for markers in columns and trees in rows
AlderTran <- t(vcf77)
dim(AlderTran) 

# impute missing marker data using EM rather than mean - more appropriate for GBS data
vcf77_Impute <- A.mat(AlderTran, max.missing=0.5, impute.method="EM", return.imputed=T)
#save(vcf103_Impute,file = "vcf103_Impute_77.RData")

#extract markers that have been imputed
#imputed returns the imputed marker matrix
#A returns the additive relationship matrix
markers_impute<-vcf77_Impute$imputed
dim(markers_impute)
#everything seems fine, no need to filter any markers
#setting 70% of data for training and 30% for validation populations
pheno<-as.data.frame(emm_all_77)
dim(pheno)
str(pheno)
print(pheno$Genotype)

#convert the first column to row names
pheno<-pheno %>%
  column_to_rownames(var = "Genotype")
pheno<-as.matrix(pheno)
markers_impute<-as.matrix(markers_impute)
#reorder marker data
markers_impute1<-markers_impute[order(row.names(markers_impute)),]
str(pheno)
str(markers_impute)
dim(pheno)
dim(markers_impute1)

#setting training and test population
train<-as.matrix(sample(1:77,54))
train

test<-setdiff(1:77, train)
test

## prediction of 2 year height growth
#specify the x2h as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy1<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x2h<-pheno_train[,1]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x2h)<-"numeric"
  x2h_answer<-mixed.solve(x2h,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x2h_markereffect<-x2h_answer$u
  e<-as.matrix(x2h_markereffect)
  pred_x2h_valid<-m_valid %*% e
  pred_x2h<-(pred_x2h_valid[,1])+x2h_answer$beta
  pred_x2h
  x2h_valid<-pheno_valid[,1]
  class(pred_x2h_valid)<-"numeric"
  class(x2h_valid)<-"numeric"
  accuracy1[r,1]<-cor(pred_x2h_valid,x2h_valid,use = "complete")
}
mean(accuracy1)

## prediction of 7 year height growth
#specify the x7h as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy2<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x7h<-pheno_train[,2]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x7h)<-"numeric"
  x7h_answer<-mixed.solve(x7h,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x7h_markereffect<-x7h_answer$u
  e<-as.matrix(x7h_markereffect)
  pred_x7h_valid<-m_valid %*% e
  pred_x7h<-(pred_x7h_valid[,1])+x7h_answer$beta
  pred_x7h
  x7h_valid<-pheno_valid[,1]
  class(pred_x7h_valid)<-"numeric"
  class(x7h_valid)<-"numeric"
  accuracy2[r,1]<-cor(pred_x7h_valid,x7h_valid,use = "complete")
}
mean(accuracy2)


## prediction of 9 year height growth
#specify the x9h as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy3<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x9h<-pheno_train[,3]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x9h)<-"numeric"
  x9h_answer<-mixed.solve(x9h,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x9h_markereffect<-x9h_answer$u
  e<-as.matrix(x9h_markereffect)
  pred_x9h_valid<-m_valid %*% e
  pred_x9h<-(pred_x9h_valid[,1])+x9h_answer$beta
  pred_x9h
  x9h_valid<-pheno_valid[,1]
  class(pred_x9h_valid)<-"numeric"
  class(x9h_valid)<-"numeric"
  accuracy3[r,1]<-cor(pred_x9h_valid,x9h_valid,use = "complete")
}
mean(accuracy3)

## prediction of 12 year height growth
#specify the x12h as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy4<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x12h<-pheno_train[,4]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x12h)<-"numeric"
  x12h_answer<-mixed.solve(x12h,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x12h_markereffect<-x12h_answer$u
  e<-as.matrix(x12h_markereffect)
  pred_x12h_valid<-m_valid %*% e
  pred_x12h<-(pred_x12h_valid[,1])+x12h_answer$beta
  pred_x12h
  x12h_valid<-pheno_valid[,1]
  class(pred_x12h_valid)<-"numeric"
  class(x12h_valid)<-"numeric"
  accuracy4[r,1]<-cor(pred_x12h_valid,x12h_valid,use = "complete")
}
mean(accuracy4)

## prediction of 9 year DBH growth
#specify the x9dbh as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy5<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x9dbh<-pheno_train[,5]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x9dbh)<-"numeric"
  x9dbh_answer<-mixed.solve(x9dbh,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x9dbh_markereffect<-x9dbh_answer$u
  e<-as.matrix(x9dbh_markereffect)
  pred_x9dbh_valid<-m_valid %*% e
  pred_x9dbh<-(pred_x9dbh_valid[,1])+x9dbh_answer$beta
  pred_x9dbh
  x9dbh_valid<-pheno_valid[,1]
  class(pred_x9dbh_valid)<-"numeric"
  class(x9dbh_valid)<-"numeric"
  accuracy5[r,1]<-cor(pred_x9dbh_valid,x9dbh_valid,use = "complete")
}
mean(accuracy5)

## prediction of 12 year DBH growth
#specify the x12dbh as response variable
#u is the marker effects
#correlation accuracy with 500 iterations
traits<-1
cycles<-500
accuracy6<-matrix(nrow = cycles,ncol = traits)
for (r in 1:cycles) {
  train<-as.matrix(sample(1:77,54))
  test<-setdiff(1:77, train)
  pheno_train<-pheno[train,]
  m_train<-markers_impute1[train,]
  pheno_valid<-pheno[test,]
  m_valid<-markers_impute1[test,]
  
  x12dbh<-pheno_train[,6]
  class(m_train) <- "numeric"
  class(m_valid) <- "numeric"
  class(x12dbh)<-"numeric"
  x12dbh_answer<-mixed.solve(x12dbh,Z=m_train,K=NULL,SE=FALSE,return.Hinv = FALSE)
  x12dbh_markereffect<-x12dbh_answer$u
  e<-as.matrix(x12dbh_markereffect)
  pred_x12dbh_valid<-m_valid %*% e
  pred_x12dbh<-(pred_x12dbh_valid[,1])+x12dbh_answer$beta
  pred_x12dbh
  x12dbh_valid<-pheno_valid[,1]
  class(pred_x12dbh_valid)<-"numeric"
  class(x12dbh_valid)<-"numeric"
  accuracy6[r,1]<-cor(pred_x12dbh_valid,x12dbh_valid,use = "complete")
}
mean(accuracy6)

#combine result of 500 iterations of each trait
accuracy_all<-data.frame(accuracy1,accuracy2,accuracy3,accuracy4,accuracy5,accuracy6)
write.csv(accuracy_all,"accuracy_all.csv")
#calculate the mean of 500 iterations
colMeans(accuracy_all)

#rename column name
colnames(accuracy_all)<-c("x2h", "x7h","x9h", "x12h", "x9dbh", "x12dbh")
accuracy_all_long<-accuracy_all %>%
  pivot_longer(cols = c('x2h','x7h','x9h','x12h','x9dbh','x12dbh'),
               names_to = 'trait',
               values_to = 'accuracy')
write.csv(accuracy_all_long,"accuracy_all_long.csv")

#plot the GS accuracy distribution
library(dplyr)
library(ggplot2)
order = c("x2h", "x7h","x9h", "x12h", "x9dbh", "x12dbh")
accuracy_all_long %>% mutate(trait2 = case_when(trait == "x2h" ~ "RG_H\n(2008-2010)",
                                                trait == "x7h" ~ "RG_H\n(2008-2015)",
                                                trait == "x9h" ~ "RG_H\n(2008-2017)",
                                                trait == "x12h" ~ "RG_H\n(2008-2020)",
                                                trait == "x9dbh" ~ "RG_DBH\n(2015-2017)",
                                                trait == "x12dbh" ~ "RG_DBH\n(2015-2020)")) %>%
ggplot(aes(factor(trait2, levels = c("RG_H\n(2008-2010)", "RG_H\n(2008-2015)","RG_H\n(2008-2017)", "RG_H\n(2008-2020)", 
                                     "RG_DBH\n(2015-2017)", "RG_DBH\n(2015-2020)")), y=accuracy, fill=trait2))+
  geom_hline(yintercept = 0)+
  geom_violin()+
  #geom_point()+
  geom_boxplot(width=0.25, fill="white")+
  geom_jitter(width = 0.1, alpha=0.1) +
  scale_fill_brewer(palette = "Pastel2") +
  #scale_fill_manual(breaks = c("x2h", "x7h","x9h", "x12h", "x9dbh", "x12dbh"),
   #                 values = c("palegreen","palegreen1","palegreen2","palegreen3","salmon","salmon2"))+
  labs(title = "", x="", y="Predictability")+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank())
#ggsave("Figure4.png", dpi=300, h = 6, w = 9)



#separate height and dbh
library(tidyr)
accuracy_all_height<-accuracy_all[,-c(5,6)]
accuracy_all_height_long<-accuracy_all_height %>%
  pivot_longer(cols = c('x2h','x7h','x9h','x12h'),
               names_to = 'trait',
               values_to = 'accuracy')
order = c("x2h", "x7h","x9h", "x12h")
height<-accuracy_all_height_long %>% mutate(trait2 = case_when(trait == "x2h" ~ "RG_H\n(2008-2010)",
                                                trait == "x7h" ~ "RG_H\n(2008-2015)",
                                                trait == "x9h" ~ "RG_H\n(2008-2017)",
                                                trait == "x12h" ~ "RG_H\n(2008-2020)")) %>%
  ggplot(aes(factor(trait2, levels = c("RG_H\n(2008-2010)", "RG_H\n(2008-2015)","RG_H\n(2008-2017)", "RG_H\n(2008-2020)")), 
             y=accuracy, fill=trait2))+
  geom_hline(yintercept = 0)+
  geom_violin()+
  ylim(-1,1)+
  #geom_point()+
  geom_boxplot(width=0.25, fill="white")+
  geom_jitter(width = 0.1, alpha=0.1) +
  scale_fill_brewer(palette = "Pastel2") +
  #scale_fill_manual(breaks = c("x2h", "x7h","x9h", "x12h"),
                   #values = c("palegreen","palegreen1","palegreen2","palegreen3"))+
  labs(title = "", x="", y="Predictability")+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank())
#ggsave("Figure4a.png", dpi=300, h = 6, w = 9)

accuracy_all_dbh<-accuracy_all[,-c(1,2,3,4)]
accuracy_all_dbh_long<-accuracy_all_dbh %>%
  pivot_longer(cols = c('x9dbh','x12dbh'),
               names_to = 'trait',
               values_to = 'accuracy')
order = c("x9dbh", "x12dbh")
dbh<-accuracy_all_dbh_long %>% mutate(trait2 = case_when(trait == "x9dbh" ~ "RG_DBH\n(2015-2017)",
                                                trait == "x12dbh" ~ "RG_DBH\n(2015-2020)")) %>%
  ggplot(aes(factor(trait2, levels = c("RG_DBH\n(2015-2017)", "RG_DBH\n(2015-2020)")), y=accuracy, fill=trait2))+
  geom_hline(yintercept = 0)+
  geom_violin()+
  ylim(-1,1)+
  #geom_point()+
  geom_boxplot(width=0.25, fill="white")+
  geom_jitter(width = 0.1, alpha=0.1) +
  scale_fill_brewer(palette = "Pastel1") +
  #scale_fill_manual(breaks = c("x9dbh", "x12dbh"),
                   #values = c("salmon","salmon2"))+
  labs(title = "", x="", y="")+
  theme_bw(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank())
#ggsave("Figure4b.png", dpi=300, h = 6, w = 9)

library(cowplot)
plot_grid(height,dbh,rel_widths = c(1,0.5),
          labels = c("A","B"))
ggsave("Figure4.png", dpi=300, h = 6, w = 9)


#GRM of 103 samples
# transpose data set for markers in columns and trees in rows
AlderTran_103 <- t(vcf103)
dim(AlderTran_103)  

# impute missing marker data using EM rather than mean - more appropriate for GBS data
vcf103_Impute <- A.mat(AlderTran_103, max.missing=0.5, impute.method="EM", return.imputed=T)

# extract additive relationship matrix
vcf103_GRM <- vcf103_Impute$A
# add back in the names of clones
vcf103_GRM[1:5, 1:5]

#reordering the samples by population
order103<-read.csv("103alder_location.csv",header = TRUE)
order103<-order103[,-1]
row.names(order103)<-order103[,1]
order103<-as.data.frame(order103)
order103<-order103$Genotype
vcf103_GRM_1<-reorder_mat(mat = as.matrix(vcf103_GRM), order = order103)
vcf103_GRM_1<-as.data.frame(vcf103_GRM_1)
vcf103_GRM_1<-as.matrix(vcf103_GRM_1)

pheatmap(vcf103_GRM_1,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 8,
         fontsize_col = 8)

