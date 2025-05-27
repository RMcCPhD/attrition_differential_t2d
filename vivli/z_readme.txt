
## 1. Description of export request ###############################################

This export request contains summary statistics, model outputs and aggregated data
from an analysis evaluating the effect of treatment, age and sex on attrition 
in 92 randomised controlled trials for type 2 diabetes mellitus.

The 92 trials are from seven sponsors:
- AstraZeneca: 5 (5.43%)
- Boehringer Ingelheim: 33 (35.9%)
- Eli lilly: 11 (12%)
- GlaxoSmithKline: 9 (9.78%)
- Janssen: 9 (10.9%)
- Sanofi: 13 (14.1%)
- Takeda: 11 (12%)


## 2. Analysis description ########################################################

Outcome: Attrition, defined as participant non-completion of a trial for any reason
	 (e.g. adverse events, lost to follow-up)

The goal of this analysis is to determine whether attrition is different between 
treatment arms (i.e. differential attrition), and whether there are differences 
by age or sex. The analysis set of 92 trials compared treatment with a newer 
antidiabetic class to a placebo or active comparator (see section 3).

We will subsequently meta-analyse the results of these trial-level models 
to produce overall estimates.

Cox proportional hazard and logistic regression models were fit to each trial:
- Cox: coxph() with strata(id) to stratify by each trial
- Logistic regression: glm() with family = "binomial"

We found the outputs from cox and logistic regression models to be similar, and have fit both
to present comparisons in our results. We are using logistic regression outputs in our 
analysis because it simplifies the analysis.

The following models were specified using cox and logreg:
* cox_unadj | Surv(time, event) ~ trt_class
* lr_unadj  | event ~ trt_class
* cox_adj   | Surv(time, event) ~ sex + age10
* lr_adj    | event ~ sex + age10
* cox_int   | Surv(time, event) ~ trt_class * (sex + age10)
* lr_int    | event ~ trt_class * (sex + age10)

Variables included in models were:
* time      | Duration of randomised participation in weeks
* event     | Attrition status (0 = completed the trial, 1 = non-completion)
* trt_class | Antidiabetic treatment class as per ATC classification (see below)
* age10     | Age, measured in decades (i.e. divided by 10)
* sex       | Participant sex (0 = female, 1 = male)

Treatment class was classified using the WHO anatomical therapeutic classification (ATC)
system:
* Insulin - A10A
* Biguanides - A10BA
* Sulfonylureas - A10BB
* Combinations of antidiabetics - A10BD
* Alpha-glucosidases - A10BF
* Thiazolidinediones - A10BG
* Dipeptidyl peptidase-4 inhibitors - A10BH
* Glucagon-like peptide-1 - A10BJ
* Sodium glucose cotransporter 2 - A10BK
* Other antidiabetics - A10BX


## Contents of export request #####################################################

# Main folder
1. agg_n92.csv - Attrition counts, one row per trial/arm combination.
* nct_id - ClinicalTrials.gov identifier
* min_dur, max_dur, med_dur - Spread of trial duration
* arm_lvl - Treatment arm and dose (coded using WHO ATC classification)
* trttype - Treatment type (experimental, active comparator or placebo)
* n - Number of participants in treatment arm
* attr - Number of attrition events in treatment arm

2. diag_n92.csv - Diagnostics for fitted models, one row per model/trial.
* nct_id
* modeltype - Type of model
* statistic.log - Likelihood ratio test statistic
* p.value.log - P-value associated with likelihood ratio test
* statistic.sc - Log-rank test statistic
* p.value.sc - P-value associated with log-rank test
* statistic.wald - Wald test statistic
* p.value.wald - P-value associated with Wald test
* statistic.robust - Robust likelihood ratio test statistic
* p.value.robust - P-value associated with robust likelihood ratio test
* r.squared - R-squared describing proportion of variation explained by model
* r.squared.max - Maximum R-squared
* concordance - Harrel's concordance (C-index)
* std.error.concordance - Standard error of C-index
* logLik - Log-likelihood
* AIC - Akaike information criterion
* BIC - Bayesian information criterion
* nobs - Number of observations used in fitted model
* null.deviance - Intercept-only model deviance
* df.null - Degrees of freedom for intercept-only model
* deviance - Residual deviance of fitted model
* df.residual - Degrees of freedom of residuals (n observations - n parameters)


3. res_n92.csv - Model coefficients and standard errors, one row per term/model/trial.
* nct_id
* modeltype
* spec - Model specification
* term - Model term (i.e. treatment class, age, sex, interactions)
* ref - Model reference treatment class
* estimate - Model log-estimates
* std.error - Log-standard errors
* statistic - Z-statistic
* p.value - P value
* lci - Lower 95% confidence interval
* uci - Upper 95% confidence interval

4. coef_cor_n92.csv
Correlation matrix for model coefficients (this is NOT the correlation between variables)
3158 rows, one row per pair of terms (i.e. each correlation) per model per trial.
This was obtained by taking upper triangle for the variance-covariance matrix from each
model and reshaping the data into a long format.
* nct_id
* modeltype
* row, col - first and second coefficient terms
* value - correlation

# Scripts
00_config.R - Configuration (i.e. disable scientific notation, clear cache)
01_prepare_data.R - Data preparation for analysis
02_mdl.R - Fit models and extract outputs
03_summarise.R - Prepare summary statistics
04_export.R - Prepare export materials


