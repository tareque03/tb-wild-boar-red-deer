###############################################################################
# Bovine TB in wild boar and red deer: exploratory analysis and
# model selection (AICc, MuMIn) for farm-level TB prevalence
#
# Author : Md Tareque Hossain
# Course : BIOL 520C
#
# Contents
#   Part 1  Exploratory analysis (correlations, bar plots)
#   Part 2  Boar model : BoarPropPos ~ candidate predictors
#   Part 3  Deer model : DeerPropPos ~ candidate predictors
#
# Input   : Deer_Tb.csv  (place in a "data/" folder next to this script,
#           or change `data_path` below)
###############################################################################

# ---- 0. Setup ---------------------------------------------------------------

# Install any missing packages
pkgs <- c("corrplot", "xtable", "MuMIn", "DescTools", "ResourceSelection")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(corrplot)
library(xtable)
library(MuMIn)
library(DescTools)
library(ResourceSelection)

# Path to the data file (relative path so the script runs on any machine)
data_path <- file.path("data", "Deer_Tb.csv")

# ---- 1. Load and clean data -------------------------------------------------

# fileEncoding = "UTF-8-BOM" avoids the odd "ï..Farm" column name
deertb <- read.csv(data_path, fileEncoding = "UTF-8-BOM")
deertb <- na.omit(deertb)
deertb

# Drop the Cervi (non-TB) columns
deertb <- subset(deertb, select = -c(DeerPosCervi, DeerSampledCervi, CerviPropPos))
str(deertb)

# Farm identifier column (first column of the file)
farm_id <- deertb[[1]]


###############################################################################
# PART 1: EXPLORATORY ANALYSIS
###############################################################################

# ---- 1.1 Correlation plot (Spearman) ----------------------------------------
x <- cor(deertb, method = "spearman")
corrplot(x, type = "upper", order = "hclust")

# ---- 1.2 Correlation table --------------------------------------------------
mcor <- round(cor(deertb, method = "spearman"), 2)
mcor

# LaTeX version of the table
print(xtable(mcor), type = "latex")

# ---- 1.3 Proportion of TB-positive boar and deer in each farm ---------------
df <- rbind(BoarPropPos = deertb$BoarPropPos,
            DeerPropPos = deertb$DeerPropPos)
colnames(df) <- farm_id

barplot(df, col = c("orange1", "skyblue"), border = "black",
        density = 30, xlab = "Farm ID",
        ylab = "Proportion of TB pos Boar vs red deer",
        main = "Proportion of TB positive Boar vs deer in each farm")
legend("topleft", legend = c("BoarPropPos", "DeerPropPos"),
       col = c("orange1", "skyblue"), lty = 1)

# ---- 1.4 Number of TB-positive boar out of those sampled, by farm -----------
df1 <- rbind(BoarSampledTB = deertb$BoarSampledTB,
             BoarPosTB     = deertb$BoarPosTB)
colnames(df1) <- farm_id

barplot(df1, col = c("blue", "red"), border = "black",
        density = 30, xlab = "Farm ID",
        ylab = "Number of Boar sampled vs positive",
        main = "Number of TB positive boar from sampled in each farm")
legend("topleft", legend = c("BoarSampledTB", "BoarPosTB"),
       col = c("blue", "red"), lty = 1)

# ---- 1.5 Number of TB-positive deer out of those sampled, by farm -----------
df2 <- rbind(DeerSampledTB = deertb$DeerSampledTB,
             DeerPosTB     = deertb$DeerPosTB)
colnames(df2) <- farm_id

barplot(df2, col = c("blue", "red"), border = "black",
        density = 30, xlab = "Farm ID",
        ylab = "Number of deer sampled vs positive",
        main = "Number of TB positive deer from sampled in each farm")
legend("topleft", legend = c("DeerSampledTB", "DeerPosTB"),
       col = c("blue", "red"), lty = 1)


###############################################################################
# PART 2: BOAR MODEL (response: BoarPropPos)
###############################################################################

# ---- 2.1 Global model -------------------------------------------------------
# na.action = na.fail is required by dredge()
model1 <- glm(BoarPropPos ~ DeerPosTB + DeerSampledTB + DeerPropPos +
                OpenLand + ScrubLand + PinePlantation + QuercusPlants +
                QuercusTrees + WildBoarIndex + ReedDeerIndex + EstateSize +
                Fenced,
              family = binomial(link = "logit"),
              data = deertb, na.action = na.fail)

# ---- 2.2 Fit all subsets and rank by AICc -----------------------------------
deertb_dredge <- dredge(model1, beta = FALSE, evaluate = TRUE, rank = AICc)

options(na.action = "na.omit")
head(deertb_dredge)

# Number of model combinations
nrow(deertb_dredge)

# ---- 2.3 Best model ---------------------------------------------------------
top_model <- get.models(deertb_dredge, subset = 1)[[1]]
top_model

# Refit of the best model (BoarPropPos ~ DeerPropPos) for diagnostics
bestfit <- glm(formula = BoarPropPos ~ DeerPropPos + 1,
               family = binomial(link = "logit"),
               data = deertb, na.action = na.fail)

# Residual plot of the best-fit model
plot(fitted(bestfit), resid(bestfit))
abline(0, 0)

# Summary of the best model
summary(top_model)

# ---- 2.4 Model averaging (models within delta AICc <= 2) --------------------
model_avg <- model.avg(deertb_dredge, subset = delta <= 2)
summary(model_avg)

# Model containing the averaged terms
Average.model <- glm(formula = BoarPropPos ~ DeerPosTB + DeerPropPos + QuercusTrees,
                     family = binomial(link = "logit"),
                     data = deertb, na.action = na.fail)

# Residual plot of the averaged model
plot(fitted(Average.model), resid(Average.model))
abline(0, 0)

# ---- 2.5 Predictive power (pseudo R-squared) --------------------------------
PseudoR2(bestfit)
PseudoR2(Average.model)

# ---- 2.6 Goodness of fit: Hosmer-Lemeshow test ------------------------------
hoslem.test(bestfit$y, fitted(bestfit), g = 10)
hoslem.test(Average.model$y, fitted(Average.model), g = 10)


###############################################################################
# PART 3: DEER MODEL (response: DeerPropPos)
###############################################################################

# ---- 3.1 Global model -------------------------------------------------------
options(na.action = "na.fail")   # reset, dredge() requires na.fail

model2 <- glm(DeerPropPos ~ BoarPosTB + BoarSampledTB + BoarPropPos +
                OpenLand + ScrubLand + PinePlantation + QuercusPlants +
                QuercusTrees + WildBoarIndex + ReedDeerIndex + EstateSize +
                Fenced,
              family = binomial(link = "logit"),
              data = deertb, na.action = na.fail)

# ---- 3.2 Fit all subsets and rank by AICc -----------------------------------
deertb_dred <- dredge(model2, beta = FALSE, evaluate = TRUE, rank = AICc)

options(na.action = "na.omit")
head(deertb_dred)

# Number of model combinations
nrow(deertb_dred)

# ---- 3.3 Best model ---------------------------------------------------------
top_model1 <- get.models(deertb_dred, subset = 1)[[1]]
top_model1
summary(top_model1)

# The top model is the intercept-only model
bestfit1 <- glm(DeerPropPos ~ 1,
                family = binomial(link = "logit"),
                data = deertb, na.action = na.fail)

# ---- 3.4 Predictive power (pseudo R-squared) --------------------------------
PseudoR2(bestfit1)

# ---- End of script ----------------------------------------------------------
