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
library(ggplot2)     # For plotting
col1 <- '#82c0cc'

## 2) import pdr data
setwd("C:/Users/user/Desktop/CommunicationSurprise_DataCode/") #you need to change your working directory to the location of "CommunicationSurprise_DataCode" folder
FolderD = "Model-based PDR Analysis/pdr_data.csv"
TCGdata <- read.table(FolderD, sep=",", header=TRUE) #read the table

## 3) check the distribution ----
dis<-check_distribution(TCGdata$surprise)
plot(dis)

## 4) Create the step column
TCGdata <- TCGdata %>%
  group_by(sub, trial) %>%
  mutate(step = ceiling(row_number() / 120)) %>%
  ungroup()

## 5) Now you can aggregate the data over the steps
# For example, to get the mean of Pupil_base and surprise for each step:
TCGdata_aggregated <- TCGdata %>%
  group_by(sub, trial, step,surprise_levels) %>%
  summarise(
    Pupil_base = mean(Pupil_base, na.rm = TRUE),
    surprise = mean(surprise, na.rm = TRUE)
  )

## 6) View the aggregated data
print(aggregated_data)


## 7) define factors
TCGdata_aggregated$sub <- as.factor(TCGdata_aggregated$sub)#factor 1
TCGdata$sub <- as.factor(TCGdata$sub)#factor 1

## 8) run the lmer() evarage pdr predicted response  (NONlinear) for surprise as the continous variable
m1 <- lmer(Pupil_base ~ surprise-1 + (1 + surprise | sub), data = TCGdata) # with radnom slope and intercept 

## Check residuals for normality
residuals_m1 <- residuals(m1)
check_distribution(residuals_m1)
plot(check_distribution(residuals_m1))


## Check for homoscedasticity (equal variances)
library(performance)
check_heteroscedasticity(m1)

## 9) estimate the predicted means
emm1 <- ggemmeans(m1, 'surprise[all]')


plot(emm1)
print(emm1)

## 10) summarize m1
sum1=summary(m1)

print(m1)
fixef(m1)
anova(m1)
parameters::p_value(m1)
eta_squaered=eta_squared(anova(m1), partial = FALSE)

## 11) save output
write.csv(x=sum1,"/lmer_output/sum_model.csv");
write.csv(x=emm1,"/lmer_output/pred_model.csv");


## 11) plot and save save the plots for the paper
saveF = "C:/Users/user/Desktop/CommunicationSurprise_DataCode/Model-based PDR Analysis/lmer_output/lmer_1st_data.pdf"
pdf(saveF, width=4, height=5)
p1 <- plot(emm1,ci = TRUE,line.size = 2,show.title = FALSE, show.x.title = TRUE, show.y.title = TRUE,colors = c(col1))
print(plots(p1))
dev.off()
