# tb-wild-boar-red-deer
Prevalences of Tuberculosis in wild boar and red deer

# Tuberculosis in Wild Boar and Red Deer in Spain

Statistical analysis of TB-like lesion prevalence in wild boar (*Sus scrofa*) and red deer (*Cervus elaphus*) across Spanish estates, using logistic regression with AICc model selection and model averaging in R.

*BIOL 520C course project, UBC Okanagan, December 2021. Author: Md Tareque Hossain*

📄 **[Read the full report (PDF)](report/BIO520C_TB_wild_boar_red_deer_report.pdf)**

## Question

Bovine TB (*Mycobacterium bovis*) persists in wildlife and can re-infect livestock. Do habitat, host-density and management variables, or the other host species' prevalence, explain TB prevalence on a farm in (a) wild boar and (b) red deer?

## Data

`data/Deer_Tb.csv`: farm-level TB prevalence and habitat/management covariates (see the [data dictionary](data/README.md)). The course dataset comes from the national survey in [Vicente et al. (2006)](https://doi.org/10.1051/vetres:2005044). After dropping farms with missing values and the *E. cervi* columns (no boar data), **n = 23 farms** remained.

## Methods

1. Spearman rank correlations among all variables.
2. Global logistic (binomial, logit) GLMs: boar TB proportion, and deer TB proportion, each on the remaining covariates.
3. All-subsets selection with `MuMIn::dredge()` (4,096 models each), ranked by AICc.
4. Model averaging over models with ΔAICc ≤ 2.
5. Pseudo-R² (`DescTools`) and Hosmer–Lemeshow goodness of fit (`ResourceSelection`).

## Results

### Exploratory plots

**TB-positive boar out of boar sampled, by farm**
![Boar sampled vs positive](figures/fig1_boar_sampled_vs_positive.png)

**TB-positive deer out of deer sampled, by farm**
![Deer sampled vs positive](figures/fig2_deer_sampled_vs_positive.png)

**Proportion positive in boar vs deer, by farm**
![Proportion positive](figures/fig3_proportion_positive_by_farm.png)

**Spearman correlation matrix** (full values in the [report](report/BIO520C_TB_wild_boar_red_deer_report.pdf), Table 1)
![Correlation plot](figures/fig4_correlation_plot.png)

### Model 1: boar TB proportion

Top models by AICc (best 6 of 4,096):

| Model terms | df | logLik | AICc | ΔAICc | Weight |
|---|---|---|---|---|---|
| DeerPropPos | 2 | -13.25 | 31.1 | 0.00 | 0.348 |
| DeerPosTB | 2 | -14.01 | 32.6 | 1.53 | 0.162 |
| DeerPropPos + QuercusTrees | 3 | -12.73 | 32.7 | 1.63 | 0.154 |
| DeerPropPos + Fenced | 3 | -13.02 | 33.3 | 2.20 | 0.116 |
| DeerPropPos + ScrubLand | 3 | -13.02 | 33.3 | 2.21 | 0.115 |
| DeerPropPos + EstateSize | 3 | -13.12 | 33.5 | 2.41 | 0.105 |

Best model: `BoarPropPos ~ DeerPropPos`

| Term | Estimate | Std. error | z | p |
|---|---|---|---|---|
| Intercept | -0.344 | 0.544 | -0.633 | 0.527 |
| DeerPropPos | 6.382 | 5.068 | 1.259 | 0.208 |

Model-averaged coefficients (ΔAICc ≤ 2, full average): DeerPropPos 4.78 (p = 0.379), DeerPosTB 0.027 (p = 0.739), QuercusTrees -0.030 (p = 0.778). No predictor was statistically significant.

| Model | McFadden pseudo-R² | Hosmer–Lemeshow (g = 10) |
|---|---|---|
| Best model | 0.162 | χ² = 0.861, df = 8, p = 0.999 |
| Averaged model | 0.198 | χ² = 1.472, df = 8, p = 0.993 |

| Best model residuals | Averaged model residuals |
|---|---|
| ![Residuals, best model](figures/fig7_residuals_best_model.png) | ![Residuals, averaged model](figures/fig9_residuals_average_model.png) |

### Model 2: deer TB proportion

The intercept-only model (`DeerPropPos ~ 1`) ranked best (AICc = 5.7); every other covariate gave ΔAICc ≥ 2.41. Intercept = -2.526 (SE 0.796, p = 0.0015), so none of the covariates added explanatory power for deer prevalence.

## Conclusions and limitations (from the report)

- Boar and deer prevalence were positively correlated on farms, and deer prevalence was the most supported predictor of boar prevalence, but the effect was not statistically significant.
- Limited power: only 23 complete farms, with a large candidate set of 12 predictors.
- No animal age data, which is an important driver of infection.
- Samples may not represent the full wild boar and red deer populations in Spain.
- Results are consistent with both species being potential TB reservoirs, but this analysis cannot establish transmission between them.

## Repository structure

```
├── README.md
├── data/        Deer_Tb.csv + data dictionary
├── R/           analysis.R, reproduces all figures and models
├── figures/     PNG plots used above
├── report/      Full project report (PDF)
└── docs/        Original R Markdown appendix outputs (HTML)
```

## Reproduce

```r
install.packages(c("corrplot", "xtable", "MuMIn", "DescTools", "ResourceSelection"))
source("R/analysis.R")   # run from the repository root
```

## Reference

Vicente J, Höfle U, Garrido JM, et al. Wild boar and red deer display high prevalences of tuberculosis-like lesions in Spain. *Vet Res.* 2006;37(1):107-19.
