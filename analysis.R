# ------------------------------------------------------------------
# TB prevalence in wild boar and red deer (Spain) - BIOL 520C project
# Author: Md Tareque Hossain
# Run from the repository root:  source("R/analysis.R")
# Packages: corrplot, xtable, MuMIn, DescTools, ResourceSelection
# ------------------------------------------------------------------
# install.packages(c("corrplot","xtable","MuMIn","DescTools","ResourceSelection"))

library(corrplot); library(xtable); library(MuMIn)
library(DescTools); library(ResourceSelection)

dir.create("figures", showWarnings = FALSE)

# ---- 1. Load and clean ------------------------------------------
# fileEncoding removes the BOM that otherwise renames the first column to "ï..Farm"
deertb <- read.csv("data/Deer_Tb.csv", fileEncoding = "UTF-8-BOM")
deertb <- na.omit(deertb)                       # 23 of 31 farms have complete data
deertb <- subset(deertb, select = -c(DeerPosCervi, DeerSampledCervi, CerviPropPos))
str(deertb)

# ---- 2. Exploratory analysis ------------------------------------
# Spearman correlation matrix (Fig 4, Table 1)
x <- cor(deertb, method = "spearman")
png("figures/fig4_correlation_plot.png", 1344, 960, res = 150)
corrplot(x, type = "upper", order = "hclust")
dev.off()
mcor <- round(x, 2)
print(xtable(mcor), type = "latex")

barfig <- function(file, rows, labels, cols, ylab, main) {
  m <- rbind(deertb$Farm, rows); colnames(m) <- m[1, ]; m <- m[-1, , drop = FALSE]
  rownames(m) <- labels
  png(file, 1344, 960, res = 150)
  barplot(m, col = cols, border = "black", density = 30,
          xlab = "Farm ID", ylab = ylab, main = main)
  legend("topleft", legend = labels, col = cols, lty = 1)
  dev.off()
}
barfig("figures/fig3_proportion_positive_by_farm.png",
       rbind(deertb$BoarPropPos, deertb$DeerPropPos), c("BoarPropPos", "DeerPropPos"),
       c("orange1", "skyblue"), "Proportion of TB pos Boar vs red deer",
       "Proportion of TB positive Boar vs deer in each farm")
barfig("figures/fig1_boar_sampled_vs_positive.png",
       rbind(deertb$BoarSampledTB, deertb$BoarPosTB), c("BoarSampledTB", "BoarPosTB"),
       c("blue", "red"), "Number of Boar sampled vs positive",
       "Number of TB positive boar from sampled in each farm")
barfig("figures/fig2_deer_sampled_vs_positive.png",
       rbind(deertb$DeerSampledTB, deertb$DeerPosTB), c("DeerSampledTB", "DeerPosTB"),
       c("blue", "red"), "Number of deer sampled vs positive",
       "Number of TB positive deer from sampled in each farm")

# ---- 3. Model 1: boar TB proportion -----------------------------
model1 <- glm(BoarPropPos ~ DeerPosTB + DeerSampledTB + DeerPropPos + OpenLand +
                ScrubLand + PinePlantation + QuercusPlants + QuercusTrees +
                WildBoarIndex + ReedDeerIndex + EstateSize + Fenced,
              family = binomial(link = "logit"), data = deertb, na.action = na.fail)
deertb_dredge <- dredge(model1, beta = FALSE, evaluate = TRUE, rank = AICc)
options(na.action = "na.omit")
head(deertb_dredge); nrow(deertb_dredge)        # 4096 candidate models

top_model <- get.models(deertb_dredge, subset = 1)[[1]]
summary(top_model)
bestfit <- glm(BoarPropPos ~ DeerPropPos + 1, family = binomial("logit"),
               data = deertb, na.action = na.fail)
png("figures/fig7_residuals_best_model.png", 1344, 960, res = 150)
plot(fitted(bestfit), resid(bestfit)); abline(0, 0); dev.off()

model_avg <- model.avg(deertb_dredge, subset = delta <= 2)
summary(model_avg)
Average.model <- glm(BoarPropPos ~ DeerPosTB + DeerPropPos + QuercusTrees,
                     family = binomial("logit"), data = deertb, na.action = na.fail)
png("figures/fig9_residuals_average_model.png", 1344, 960, res = 150)
plot(fitted(Average.model), resid(Average.model)); abline(0, 0); dev.off()

PseudoR2(bestfit); PseudoR2(Average.model)
hoslem.test(bestfit$y, fitted(bestfit), g = 10)
hoslem.test(Average.model$y, fitted(Average.model), g = 10)

# ---- 4. Model 2: deer TB proportion -----------------------------
model2 <- glm(DeerPropPos ~ BoarPosTB + BoarSampledTB + BoarPropPos + OpenLand +
                ScrubLand + PinePlantation + QuercusPlants + QuercusTrees +
                WildBoarIndex + ReedDeerIndex + EstateSize + Fenced,
              family = binomial(link = "logit"), data = deertb, na.action = na.fail)
deertb_dred <- dredge(model2, beta = FALSE, evaluate = TRUE, rank = AICc)
options(na.action = "na.omit")
head(deertb_dred); nrow(deertb_dred)
top_model1 <- get.models(deertb_dred, subset = 1)[[1]]
summary(top_model1)                             # intercept-only model is best
bestfit1 <- glm(DeerPropPos ~ 1, family = binomial("logit"),
                data = deertb, na.action = na.fail)
PseudoR2(bestfit1)
