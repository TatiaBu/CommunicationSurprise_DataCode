# 1st_data is the one that is used for the final analysis. In the data suprise is categorized with the following rule: The 33.33rd percentile is set as the lowThreshold. Values equal to or below this threshold are considered 'low'.
#The 80.67th percentile  is set as the highThreshold. Values above this threshold are considered 'high'.

## 1) import packages ----
library(parameters)
library(insight)
library(ggeffects)
library(lme4)
library(easystats)
library(splines)
library(patchwork)
library(jtools)
library(lmerTest)
library(effectsize)
library(emmeans)
library(dplyr)

## 2) define colors for ploting
col1 <- '#82c0cc'
col2 <- '#489fb5'
col3 <- '#16697a'

## 3)set working directory and  import pdr data
setwd("C:/Users/user/Desktop/CommunicationSurprise_DataCode") #you need to change your working directory to the location of "CommunicationSurprise_DataCode" folder
FolderD = "Model-based PDR Analysis/pdr_data.csv"
TCGdata <- read.table(FolderD, sep=",", header=TRUE) #read the table

## 4) check the distribution ----
dis<-check_distribution(TCGdata$surprise)
plot(dis)

## 5) aggregate data
TCGdata_aggregated <- TCGdata %>%
  group_by(sub, surprise_levels) %>%
  summarize(Pupil_base = mean(Pupil_base, na.rm = TRUE))


## 6) define factors
TCGdata_aggregated$sub <- as.factor(TCGdata_aggregated$sub)#factor 1
TCGdata_aggregated$surprise_levels <- as.factor(TCGdata_aggregated$surprise_levels)#factor 2


## 7) run the lmer() evarage pdr predicted response  (linear) for different surprise types
m1 <- lmer(Pupil_base ~ surprise_levels-1 + (1  | sub), data = TCGdata_aggregated)

## 8) Estimate marginal means for surprise levels
emm <- emmeans(m1, specs = "surprise_levels")

## 9) Apply the contrasts
con=contrast(emm, method = list("3-2" = c(0, -1, 1), "3-1" = c(-1, 0, 1),"2-1" = c(-1, 1, 0)))
sum=summary(con,adjust = "bonferroni")
print(sum)

## 10) plot and print suprise levels
emm1=ggemmeans(m1,c('surprise_levels'))
plot(emm1)
print(emm1)

## 11) summarize m1
summary(m1)
print(m1)
fixef(m1)
anova(m1)
parameters::p_value(m1)
eta_squaered=eta_squared(anova(m1), partial = FALSE)

## 12) save output
write.csv(x=sum,"lmer_output/sum_model.csv");
write.csv(x=emm1,"lmer_output/pred_model1.csv");


## 13) plot and save save the plots for the paper
saveF = "lmer_output/lmer_1st_data.pdf"
pdf(saveF, width=6, height=5)
p1 <- plot(emm1,ci = TRUE,line.size = 2,show.title = FALSE, show.x.title = TRUE, show.y.title = TRUE,colors = c(col1,col2,col3))
print(plots(p1))
dev.off()
