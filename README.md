# 🫀 Clinical Survival Analysis: Parametric Accelerated Failure Time (AFT) Modeling

A clinical biostatistical research study modeling survival times and risk factors for heart failure patients. This project evaluates parametric **Accelerated Failure Time (AFT)** models (Weibull, Log-Normal, and Log-Logistic) against semi-parametric alternatives to quantify how clinical covariates accelerate or decelerate patient survival time.

---

## 📌 Project Deliverables
* **Research Poster:** [View PDF](docs/research_poster.pdf)
* **Presentation Deck:** [View Slides](docs/presentation_deck.pdf)
* **Analysis Code:** [R Script](scripts/survival_analysis.R)

---

## 🔬 Methodology & Statistical Framework

### 1. Exploratory Survival Analysis
* **Non-Parametric Baseline:** Estimated baseline survival probabilities using **Kaplan-Meier (KM)** estimators.
* **Hypothesis Testing:** Conducted **Log-Rank tests** across clinical categorical subgroups (e.g., ejection fraction tiers, serum creatinine levels, age brackets) to identify statistically significant divergence in survival distributions.

### 2. Parametric AFT Modeling
Unlike Cox Proportional Hazards models that assume proportional hazard ratios over time, AFT models assess the direct effect of covariates on survival time via an acceleration factor ($\gamma = \exp(\beta)$):
$$\ln(T) = \mu + \mathbf{x}'\boldsymbol{\beta} + \sigma \epsilon$$

* Evaluated distributions: **Weibull**, **Log-Normal**, and **Log-Logistic**.
* **Model Selection:** Compared goodness-of-fit using Akaike Information Criterion (**AIC**), Bayesian Information Criterion (**BIC**), and residual diagnostics (Cox-Snell residuals).

---

## 📊 Key Findings & Clinical Insights
* **Significant Predictors:** Identified key clinical covariates—including low ejection fraction and elevated serum creatinine—as significant drivers that accelerate the time to an event.
* **Parametric Fit:** Parametric AFT formulations provided robust baseline hazard estimation and intuitive time-ratio interpretations for clinical prognostic use.

---

## 🛠️ Tools & Libraries
* **Language:** R
* **Core Packages:** `survival`, `survminer`, `ggplot2`, `dplyr`, `flexsurv`
* **Presentation & Typesetting:** LaTeX / Beamer & Overleaf

---

## 📂 Repository Structure
```text
├── data/                  # Clinical dataset
├── scripts/               # R modeling scripts and data cleaning
├── docs/                  # Academic poster and slide deck PDFs
├── figures/               # Generated KM curves and model diagnostic plots
└── README.md
