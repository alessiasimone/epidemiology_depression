library(readxl)
library(skimr)
library(DataExplorer)
library(corrplot)
library(GGally)
library(mediation)
library(MASS)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)
library(survival)
library(epiR)

#----------------------------------------------------------
#DESCRIPTIVE ANALYSIS
df <- read_excel('C:/Users/aless/OneDrive - Università degli Studi di Catania/2 Sem/data analysis for Public Health/report/data.xlsx')

names(df)[9] <- "das28_crp"
names(df)[10] <- "ra_stage"
names(df)[11] <- "haq_di"
df$ra_stage <- as.integer(df$ra_stage)
df$Biologics_use <- as.integer(df$Biologics_use)
df$MTX_use <- as.integer(df$MTX_use)
df$Prednisolone_use <- as.integer(df$Prednisolone_use)
df[, 15:34] <- lapply(df[, 15:34], as.integer)
df <- df[, -1]
df <- df[, -3]
#creating binary depression
df$depression_binary <- ifelse(df$depression_score >= 0 & df$depression_score <= 9, 0,
                        ifelse(df$depression_score >= 10 & df$depression_score <= 21, 1, NA))

#creating healthy variable
df$healthy <- NA
healthy_criteria <- c("Staple_food_breakfast", "Staple_food_lunch", "Staple_food_dinner",
                      "Meat", "Fish", "Tofu", "Vegetable", "Fruits")

unhealthy_criteria <- c("Fried_food", "Cake", "Juice", "Snack_food", "Sweets",
                        "Miso_soup", "Pickles", "Ham", "Frozen_food", "Alcohol")

healthy_count <- apply(df[, healthy_criteria], 1, function(x) length(unique(na.omit(x))))
unhealthy_count <- apply(df[, unhealthy_criteria], 1, function(x) length(unique(na.omit(x))))

df$healthy <- ifelse(healthy_count > unhealthy_count, 0, 1)
df$healthy <- as.integer(df$healthy)

convert_labels <- function(value) {
  if (value == 1 | value == 2) {
    return(0)
  } else if (value == 3 | value == 4) {
    return(1)
  } else {
    return("Unknown")
  }
}

df$Stage_label <- sapply(df$ra_stage, convert_labels)

plot_intro(df)
skim(df)

df <- na.omit(df)

plot_bar(df)
plot_density(df)

correl <- cor(df)
corrplot(correl, method="color", order='hclust', addrect=2, tl.col='black', tl.cex=0.5)

plot(df[, 1:12])

#--------------------------------------------------------------
#PREVALENCE ANALYSIS

(prevalence <- table(df$Stage_label) / nrow(df))
(prevalence <- table(df$depression_binary) / nrow(df))
(prevalence <- table(df$healthy) / nrow(df))
#46% of patients have arthrite on early stage 
#and 54% of patients have arthrite on late stage

#-----------------------------------------------------------------
#RISK FACTOR ANALYSIS

model <- glm(depression_score ~ + Egg + Meat + Fried_food + Juice + Snack_food + Sweets
             + Milk + Cake + Vegetable + Fruits + Tofu + Fish + Alcohol + Miso_soup + Pickles 
             + Staple_food_breakfast + Staple_food_dinner + Staple_food_lunch, data = df)
model.opt <- stepAIC(model, direction='both')
summary(model.opt)
plot(model.opt)

#Measure of association 
epi.2by2(table(-df$depression_binary,-df$healthy))

#----------------------------------------------------------------
#DIETARY ANALYSIS

df$healthy <- as.factor(df$healthy)
# Calculate mean CRP values for healthy and unhealthy diets
mean_crp <- df %>%
  group_by(healthy) %>%
  summarize(mean_crp = mean(CRP))

# Plot mean CRP values
ggplot(mean_crp, aes(x = healthy, y = mean_crp, fill = healthy)) +
  geom_bar(stat = "identity") +
  labs(x = "Diet", y = "Mean CRP") +
  scale_fill_manual(values = c("green", "red"), labels = c("Healthy", "Unhealthy")) +
  theme_minimal()

#-------------------------------------------------------------------
#MEDIATION ANALYSIS

#fit the mediator
fit.mediator <- glm(healthy ~ depression_score, df, family=binomial)

#effect of mediator on dependent variable
fit.dv <- lm(ra_stage ~ depression_score + healthy, df)

results = mediate(fit.mediator, fit.dv, treat='depression_score', mediator='healthy', boot=T)
summary(results)

#--------------------------------------------------------------------
#SURVIVAL ANALYSIS

surv.model <- coxph(Surv(time = df$Disease_duration, event = df$Stage_label) ~ VAS + haq_di + das28_crp + Biologics_use + MTX_use + Prednisolone_use, data = df)
summary(surv.model)
