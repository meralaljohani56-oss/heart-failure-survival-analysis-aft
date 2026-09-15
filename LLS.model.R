## HEART FAILURE SURVIVAL ANALYSIS 
## Log-Location-Scale Models: Weibull, Log-Normal, Log-Logistic

library(car)
library(survival)
library(readxl)


## 1-READ DATA 
hf.data <- read.table(file="clipboard", sep="\t", header=TRUE)

## 2-RENAME COLUMNS
colnames(hf.data)[5]  <- "ejection"
colnames(hf.data)[6]  <- "blood.prus"
colnames(hf.data)[8]  <- "ser.creatnine"
colnames(hf.data)[9]  <- "ser.sodium"
colnames(hf.data)

## 3-FIT FULL MODEL
all_predictors <- c("age","anaemia","diabetes","ejection","blood.prus","platelets","ser.creatnine","ser.sodium","sex","smoking")

full_formula <- as.formula(paste("Surv(time, DEATH_EVENT) ~",paste(all_predictors, collapse=" + ")))
print(full_formula)

full_model <- survreg(full_formula, data=hf.data, dist="weibull")

## 4-BACKWARD SELECTION 
best.model.selection <- step(full_model, direction="backward")
final_formula        <- formula(best.model.selection)
print(final_formula)

## 5-FIT MODELS
fit.wei   <- survreg(final_formula, data=hf.data, dist="weibull")
fit.lnorm <- survreg(final_formula, data=hf.data, dist="lognormal")
fit.llog  <- survreg(final_formula, data=hf.data, dist="loglogistic")

summary(fit.wei)
summary(fit.lnorm)
summary(fit.llog)

## 6-COEFFICIENT TABLES 
coef.table.wei   <- summary(fit.wei)$table;   print(coef.table.wei)
coef.table.lnorm <- summary(fit.lnorm)$table; print(coef.table.lnorm)
coef.table.llog  <- summary(fit.llog)$table;  print(coef.table.llog)

## 7-AIC / BIC  / CAIC = BIC + k  (stronger penalty)
AIC_vals <- AIC(fit.wei, fit.lnorm, fit.llog)
BIC_vals <- BIC(fit.wei, fit.lnorm, fit.llog)
k_wei   <- length(coef(fit.wei))   + 1
k_lnorm <- length(coef(fit.lnorm)) + 1
k_llog  <- length(coef(fit.llog))  + 1
n       <- nrow(hf.data)

logLik.wei   <- as.numeric(logLik(fit.wei))
logLik.lnorm <- as.numeric(logLik(fit.lnorm))
logLik.llog  <- as.numeric(logLik(fit.llog))

CAIC_wei   <- -2*logLik.wei   + k_wei   * (log(n) + 1)
CAIC_lnorm <- -2*logLik.lnorm + k_lnorm * (log(n) + 1)
CAIC_llog  <- -2*logLik.llog  + k_llog  * (log(n) + 1)

model.comparison <- data.frame(
  Model      = c("Weibull","Log-Normal","Log-Logistic"),
  LogLik     = round(c(logLik.wei, logLik.lnorm, logLik.llog), 3),
  k          = c(k_wei, k_lnorm, k_llog),
  AIC        = round(AIC_vals$AIC, 3),
  BIC        = round(BIC_vals$BIC, 3),
  CAIC       = round(c(CAIC_wei, CAIC_lnorm, CAIC_llog), 3))
print(model.comparison)

## 8-VIF 
vif_result <- vif(lm(update(final_formula, time ~ .), data=hf.data))
print(round(vif_result, 3))


## SECTION A — DESCRIPTIVE STATISTICS 
## Summary table
cat("\n--- Descriptive Statistics ---\n")
desc_vars <- c("age","ejection","ser.creatnine","ser.sodium",
               "platelets","creatinine_phosphokinase","time")

summary(hf.data[, desc_vars])

## Event counts
cat("\nEvent distribution:\n")
table(hf.data$DEATH_EVENT)
prop.table(table(hf.data$DEATH_EVENT))

#################################FIGURE#############################################

## FIGURE 1 — KM by Ejection Fraction 
hf.data$ef_group <- cut(hf.data$ejection,
                        breaks = c(0, 29, 45, 100),
                        labels = c("Low (<30%)","Medium (30-45%)","High (>45%)"),
                        right  = TRUE)

km_ef <- survfit(Surv(time, DEATH_EVENT) ~ ef_group, data=hf.data)

par(mar=c(5,5,4,2))
plot(km_ef,
     xlab  = "Time",
     ylab  = "Survival Probability",
     main  = "Figure 1: Kaplan-Meier by Ejection Fraction ",
     col   = c("#c0392b","#e67e22","#2a7f7f"),
     lwd   = 2.5,
     conf.int = FALSE,
     mark.time = TRUE,
     cex.lab  = 1.2,
     cex.axis = 1.1,
     cex.main = 1.1)
legend("bottomleft",
       legend = c("Low EF (<30%)","Medium EF (30-45%)","High EF (>45%)"),
       col    = c("#c0392b","#e67e22","#2a7f7f"),
       lwd    = 2.5,
       bty    = "n", cex=0.7)

#test
logrank_ef <- survdiff(Surv(time, DEATH_EVENT) ~ ef_group, data=hf.data)
print(logrank_ef)

## FIGURE 2 — KM by Serum Creatinine 
hf.data$sc_group <- cut(hf.data$ser.creatnine,
                        breaks = c(0, 1.0, 1.5, Inf),
                        labels = c("Low (≤1.0)","Medium (1.0-1.5)","High (>1.5)"),
                        right  = TRUE)

km_sc <- survfit(Surv(time, DEATH_EVENT) ~ sc_group, data=hf.data)

par(mar=c(5,5,4,2))
plot(km_sc,
     xlab  = "Time",
     ylab  = "Survival Probability",
     main  = "Figure 2: Kaplan-Meier by Serum Creatinine ",
     col   = c("#2a7f7f","#e67e22","#c0392b"),
     lwd   = 2.5,
     conf.int = FALSE,
     mark.time = TRUE,
     cex.lab  = 1.2,
     cex.axis = 1.1,
     cex.main = 1.1)
legend("bottomleft",
       legend = c("Low SC (≤1.0)","Medium SC (1.0-1.5)","High SC (>1.5)"),
       col    = c("#2a7f7f","#e67e22","#c0392b"),
       lwd    = 2.5,
       bty    = "n", cex=0.7)

#test
logrank_sc <- survdiff(Surv(time, DEATH_EVENT) ~ sc_group, data=hf.data)
print(logrank_sc)

## FIGURE 3— Time Distribution & Age-Death

#Death rate by age group
hf.data$age_group <- cut(hf.data$age,
                         breaks = c(39,49,59,69,79,100),
                         labels = c("40-49","50-59","60-69","70-79","80+"))

km_age <- survfit(Surv(time, DEATH_EVENT) ~ age_group, data=hf.data)

par(mar=c(5,5,4,2))
plot(km_age,
     xlab      = "Time",
     ylab      = "Survival Probability",
     main      = "Figure 3: Kaplan-Meier by Age Group",
     col       = c("#2a7f7f","#5aabab","#e67e22","#d35400","#c0392b"),
     lwd       = 2.5,
     conf.int  = FALSE,
     mark.time = TRUE,
     cex.lab   = 1.2,
     cex.axis  = 1.1,
     cex.main  = 1.1)
legend("bottomleft",
       legend = c("40-49","50-59","60-69","70-79","80+"),
       col    = c("#2a7f7f","#5aabab","#e67e22","#d35400","#c0392b"),
       lwd    = 2.5, bty="n", cex=0.7,
       title  = "Age Group")

#test
logrank_age <- survdiff(Surv(time, DEATH_EVENT) ~ age_group, data=hf.data)
print(logrank_age)

## FIGURE 4 — Empirical - Hazard Function 
#Empirical cumulative hazard 
km_na <- survfit(Surv(time, DEATH_EVENT) ~ 1, data=hf.data, type="kaplan-meier")

#Nelson-Aalen: H(t) = -log(S(t))
emp_time <- km_na$time
emp_haz  <- -log(km_na$surv)   

#Fitted cumulative hazard for each model
mean_x_wei <- colMeans(model.matrix(fit.wei)[,-1])
mean_x_lnorm <- colMeans(model.matrix(fit.lnorm)[,-1])
mean_x_llog <- colMeans(model.matrix(fit.llog)[,-1])

#Fitted linear predictor at mean covariates
lp_wei   <- sum(coef(fit.wei)   * c(1, mean_x_wei))
lp_lnorm <- sum(coef(fit.lnorm) * c(1, mean_x_lnorm))
lp_llog  <- sum(coef(fit.llog)  * c(1, mean_x_llog))

t_grid <- sort(unique(hf.data$time))

#Weibull cumulative hazard: H(t) = exp((log(t) - mu) / sigma)
sigma_wei   <- fit.wei$scale
H_wei <- exp((log(t_grid) - lp_wei) / sigma_wei)

#Log-Normal cumulative hazard: H(t) = -log(1 - Phi(z))
sigma_lnorm <- fit.lnorm$scale
z_lnorm     <- (log(t_grid) - lp_lnorm) / sigma_lnorm
H_lnorm     <- -log(1 - pnorm(z_lnorm))

#Log-Logistic cumulative hazard: H(t) = -log(1 - plogis(z))
sigma_llog <- fit.llog$scale
z_llog     <- (log(t_grid) - lp_llog) / sigma_llog
H_llog     <- -log(1 - plogis(z_llog))

# Plot
ylim_all <- range(c(emp_haz, H_wei, H_lnorm, H_llog), na.rm=TRUE)
ylim_all[2] <- min(ylim_all[2], 4)   # cap y-axis for readability

par(mar=c(5,5,4,2))
plot(emp_time, emp_haz,
     type = "s",
     lwd  = 2.5,
     col  = "black",
     xlab = "Time (Days)",
     ylab = "Cumulative Hazard H(t)",
     main = "Figure 4: Fitted vs. Empirical Cumulative Hazard",
     ylim = ylim_all,
     cex.lab  = 1.2,
     cex.axis = 1.1,
     cex.main = 1.0)

lines(t_grid, H_wei,   col="#e67e22", lwd=2.5, lty=1)
lines(t_grid, H_lnorm, col="#2a7f7f", lwd=2.5, lty=2)
lines(t_grid, H_llog,  col="#c0392b", lwd=2.5, lty=4)

legend("bottomright",
       legend = c("Empirical ",
                  "Weibull ",
                  "Log-Normal ",
                  "Log-Logistic "),
       col    = c("black","#e67e22","#2a7f7f","#c0392b"),
       lwd    = 2.5,
       lty    = c(1,1,2,4),
       bty    = "n", cex=0.7)
#test
logrank_age <- survdiff(Surv(time, DEATH_EVENT) ~ age_group, data=hf.data)
print(logrank_age)

## FIGURE 5 — Cox-Snell Residuals 
cox_snell_plot <- function(fit, model_name, dist) {
  
  if (dist == "weibull") {
    z  <- (log(fit$y[,1]) - fit$linear.predictors) / fit$scale
    cs <- exp(z)                         
  } else if (dist == "lognormal") {
    z  <- (log(fit$y[,1]) - fit$linear.predictors) / fit$scale
    cs <- -log(1 - pnorm(z))
  } else if (dist == "loglogistic") {
    z  <- (log(fit$y[,1]) - fit$linear.predictors) / fit$scale
    cs <- -log(1 - plogis(z))
  }
  
  status <- fit$y[,2]
  
  ## KM estimate of cumulative hazard of residuals
  km_cs  <- survfit(Surv(cs, status) ~ 1)
  cum_haz <- -log(km_cs$surv)
  resid_t <- km_cs$time
  
  plot(resid_t, cum_haz,
       xlab = "Cox-Snell Residual ",
       ylab = "Cumulative Hazard ",
       main = paste0("Figure 5: Cox-Snell — ", model_name),
       pch  = 19, col=adjustcolor("#1a2744",0.6), cex=0.7,
       cex.lab=1.0, cex.main=0.95)
  abline(0, 1, lwd=2, col="#c0392b", lty=2)   
  legend("topleft", legend="Reference",
         col="#c0392b", lwd=2, lty=2, bty="n", cex=0.9)
}

par(mfrow=c(1,3), mar=c(5,5,4,2))
cox_snell_plot(fit.wei,   "Weibull",       "weibull")
cox_snell_plot(fit.lnorm, "Log-Normal",    "lognormal")
cox_snell_plot(fit.llog,  "Log-Logistic",  "loglogistic")
par(mfrow=c(1,1))

## FIGURE 6 — Scatter: Serum Creatinine vs Ejection Fraction
colors_outcome <- ifelse(hf.data$DEATH_EVENT==1,
                         adjustcolor("#c0392b",0.6),
                         adjustcolor("#2a7f7f",0.6))
pch_outcome <- ifelse(hf.data$DEATH_EVENT==1, 17, 16)  

par(mar=c(5,5,4,2))
plot(hf.data$ser.creatnine, hf.data$ejection,
     col  = colors_outcome,
     pch  = pch_outcome,
     cex  = 1.0,
     xlab = "Serum Creatinine (mg/dL)",
     ylab = "Ejection Fraction (%)",
     main = "Figure 6: Serum Creatinine vs. Ejection Fraction ",
     cex.lab=1.2, cex.main=1.0)
legend("bottomright",
       legend = c("Survived","Deceased"),
       col    = c(adjustcolor("#2a7f7f",0.8), adjustcolor("#c0392b",0.8)),
       pch    = c(16, 17), pt.cex=1.2,
       bty    = "n", cex=0.7)

##Figure 7:Overall Kaplan-Meier Survival Estimate
km_overall <- survfit(Surv(time, DEATH_EVENT) ~ 1, data = hf.data, conf.type = "log-log")
plot(km_overall, 
     xlab = "Time", 
     ylab = "Survival Probability",
     main = "Figure 7: Overall Kaplan-Meier Survival Estimate",
     col = "#2a7f7f", 
     lwd = 2.5, 
     conf.int = FALSE,    
     mark.time = FALSE,    
     bty = "l")

lines(km_overall$time, km_overall$upper, col = "#1a2744", lwd = 1.2, lty = 2)
lines(km_overall$time, km_overall$lower, col = "#1a2744", lwd = 1.2, lty = 2)
censored_idx <- which(hf.data$DEATH_EVENT == 0)  
censored_times <- hf.data$time[censored_idx]
for(t in censored_times) {
  idx <- max(which(km_overall$time <= t))
  surv_val <- km_overall$surv[idx]
  points(t, surv_val, pch = 3, col = "#c0392b", cex = 1.2, lwd = 1.5)
}
abline(h = 0.5, lty = 2, col = "grey50", lwd = 1.5)
legend("bottomleft", 
       legend = c("S(t)", "95% CI", "Censored", "Median (50%)"),  
       pch = c(NA, NA, 3, NA), 
       col = c("#2a7f7f", adjustcolor("#2a7f7f", alpha.f = 0.5), "#c0392b", "grey50"),
       lwd = c(2.5, 5, NA, 1.5), 
       bty = "n",
       cex = 0.85)







