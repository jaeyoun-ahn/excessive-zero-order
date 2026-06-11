# Excessive-Zero Order

R code and reproducibility notes for the manuscript **Counting models with excessive zeros ensuring stochastic monotonicity**.

This project studies claim-count models for insurance data with many zero outcomes. The main actuarial question is whether a model preserves the **credibility order**: after observing a larger claim history, the posterior predictive distribution for future claims should not move downward.

The manuscript shows that standard excessive-zero models with separate random effects for the zero and positive-count components can violate this monotonicity property. It then proposes comonotonic random-effect hurdle models that preserve stochastic monotonicity while retaining competitive out-of-sample performance.

## Repository structure

Current code structure:

```text
.
├── README.md
├── model/
│   ├── Model2_MVN_Poisson_hurdle_nimble_f.R
│   ├── Model4_Comonotonic_Poisson_hurdle_nimble_f.R
│   └── Model5_Comonotonic_NB_hurdle_nimble_f.R
└── Benchmark_model/
│    ├── BM1_2_Poisson_GLMM.R
│    └── benchmark_models_234.R
└── Data/
     ├── data.RData
     ├── dataout.RData
     └── Data_for_monotonicity_f.R

```


## Code files

| Path | Description |
|---|---|
| `model/Model2_MVN_Poisson_hurdle_nimble_f.R` | Fits Model 2: bivariate normal random-effect Poisson-hurdle model with separate random effects for the hurdle and count components. This model is flexible but does **not** guarantee the credibility order. |
| `model/Model4_Comonotonic_Poisson_hurdle_nimble_f.R` | Fits Model 4: proposed comonotonic random-effect Poisson-hurdle model. This model guarantees the credibility order. |
| `model/Model5_Comonotonic_NB_hurdle_nimble_f.R` | Fits Model 5: proposed comonotonic random-effect Negative Binomial-hurdle model. This model guarantees the credibility order. |
| `Benchmark model/BM1_2_Poisson_GLMM.R` | Fits the Poisson GLMM benchmark model, corresponding to `Benchmark Model 1` / `BM 1`. |
| `Benchmark model/benchmark_models_234.R` | Fits the remaining benchmark models: Poisson GLM, Poisson hurdle, and Poisson zero-inflated models, corresponding to `Benchmark Model 2`-`Benchmark Model 4` / `BM 2`-`BM 4`. |

## Background

Insurance claim-frequency data often contain a high proportion of zeros. Standard Poisson and Negative Binomial models may fit such data poorly because they do not separately model the zero-generating mechanism. Zero-inflated and hurdle models address this issue by adding a zero component, but random-effect versions of these models can create undesirable posterior credibility behavior.

The paper focuses on the following principle:

```text
A larger observed claim history should not imply a lower future claim-risk prediction.
```

This is formalized through stochastic monotonicity / credibility order. Violating this order can be problematic for experience-based pricing because an additional claim could, in some cases, reduce a future premium-related predictive functional.

## Models

| Label in results | Manuscript label | Model | Random effects | Credibility order |
|---|---|---|---|---|
| `Model 2` | `Model 2` | Bivariate RE Poisson-hurdle | Separate hurdle and count random effects | Not guaranteed |
| `Model 4` | `Model 4` | Comonotonic RE Poisson-hurdle | One shared latent effect | Guaranteed |
| `Model 5` | `Model 5` | Comonotonic RE Negative Binomial-hurdle | One shared latent effect | Guaranteed |
| `Benchmark Model 1` | `BM 1` | Poisson GLMM | One random effect | Guaranteed |
| `Benchmark Model 2` | `BM 2` | Poisson GLM | None | Trivial / no credibility update |
| `Benchmark Model 3` | `BM 3` | Poisson hurdle | None | Trivial / no credibility update |
| `Benchmark Model 4` | `BM 4` | Poisson zero-inflated | None | Trivial / no credibility update |


## Empirical data

The empirical study uses claim data from the Wisconsin Local Government Property Insurance Fund (LGPIF). The manuscript focuses on collision coverage for new and old vehicles.

- **Training period:** 2006-2010
- **Validation period:** 2011
- **Response:** claim frequency
- **Main covariates:** entity type and log coverage level
- **Validation target:** one-step-ahead predictive mean for 2011

Raw data are not included in this repository. Place the data in the location expected by each R script, or edit the data-path variables at the top of the scripts before running them.

## Requirements

The code is written in R. The random-effect models are fitted by Bayesian MCMC using `nimble`, while benchmark models are fitted by classical maximum likelihood where possible.

```r
install.packages(c(
  "nimble",
  "coda"
))
```

Depending on the local version of the scripts, additional packages may be required. Check the `library(...)` calls at the top of each `.R` file.

## How to run

From the repository root:

```bash
Rscript "Benchmark model/BM1_2_Poisson_GLMM.R"
Rscript "Benchmark model/benchmark_models_234.R"

Rscript "model/Model2_MVN_Poisson_hurdle_nimble_f.R"
Rscript "model/Model4_Comonotonic_Poisson_hurdle_nimble_f.R"
Rscript "model/Model5_Comonotonic_NB_hurdle_nimble_f.R"
```


## Credibility-order diagnostic

The manuscript evaluates whether posterior predictive means are monotone with respect to the most recent claim history. For each policyholder, the diagnostic compares counterfactual histories with

```text
Y_{i,5} = 0  versus  Y_{i,5} = 1.
```

The following premium functionals are checked:

| Functional | Diagnostic inequality |
|---|---|
| Full coverage | `E[Y_{i,6} | history, Y_{i,5}=0] <= E[Y_{i,6} | history, Y_{i,5}=1]` |
| Deductible coverage | `E[(Y_{i,6}-d)_+ | history, Y_{i,5}=0] <= E[(Y_{i,6}-d)_+ | history, Y_{i,5}=1]`, for `d = 1, 2` |
| Limited coverage | `E[min(Y_{i,6}, d) | history, Y_{i,5}=0] <= E[min(Y_{i,6}, d) | history, Y_{i,5}=1]`, for `d = 1, 2` |

A violation means that adding a claim to the most recent history lowers the corresponding future predictive premium functional.

## Representative results

### Credibility-order violation rates

In the manuscript's LGPIF empirical analysis, the bivariate random-effect Poisson-hurdle model has nonzero credibility-order violations, while the comonotonic random-effect models and the Poisson GLMM benchmark have no observed violations in the reported diagnostic.

| Model | Base `Y` | `(Y-d)+`, `d=1` | `(Y-d)+`, `d=2` | `min(Y,d)`, `d=1` | `min(Y,d)`, `d=2` |
|---|---:|---:|---:|---:|---:|
| Bivariate RE Poisson-hurdle | 9.29% | 10.76% | 12.72% | 0.00% | 0.98% |
| Comonotonic RE Poisson-hurdle | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| Comonotonic RE NB-hurdle | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| Poisson GLMM | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |

### Out-of-sample validation

The following results summarize validation performance on the 2011 test set. 

| Model | Test log likelihood | MSE | MAE | Guaranteed credibility order |
|---|---:|---:|---:|---|
| Benchmark Model 1 | -208.6458 | 1.1171 | 0.6063 | Yes |
| Benchmark Model 2 | -266.5104 | 2.0241 | 0.8575 | Trivial / no random effects |
| Benchmark Model 3 | -301.3224 | 2.4239 | 0.9158 | Trivial / no random effects |
| Benchmark Model 4 | -256.4331 | 1.9762 | 0.8396 | Trivial / no random effects |
| Model 2 | -208.9472 | 1.0489 | 0.5990 | No |
| Model 4 | -208.5055 | 0.9382 | 0.5994 | Yes |
| Model 5 | -214.7405 | 0.9624 | 0.6060 | Yes |

Key points:

- `Model 4`, the comonotonic RE Poisson-hurdle model, has the best test log likelihood and the lowest MSE among the listed models while guaranteeing the credibility order.
- `Model 5`, the comonotonic RE NB-hurdle model, also guarantees the credibility order and has MSE close to `Model 4`.
- The benchmark models without random effects have noticeably worse MSE and MAE in this validation table.

