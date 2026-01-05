# ==============================================
# Ordinal Mediation Analysis in R (Tailored to ftis_R.csv)
# Author: 
# ----------------------------------------------
# Dataset columns detected:
#   ID, Env_Median, Layout_Median, Staff_Median, Program_Median, Satis_Median, Loyal_Median
# This script will:
#  1) Load ftis_R.csv
#  2) Rename columns to analysis-friendly names
#  3) Coerce to ordered Likert (1-4) if applicable
#  4) Fit ordinal models (mediator & outcome) with MASS::polr
#  5) Check proportional odds (brant)
#  6) Mediation with lavaan (WLSMV, ordered), with bootstrap CIs
#  7) Report fit indices & save plots
# ==============================================

cat("===== Ordinal Mediation Analysis (ftis_R.csv) =====\\n")

# ----------------------
# 0) User options
# ----------------------
DATA_PATH <- "ftis_R.csv"
LIKERT_LEVELS <- c(1,2,3,4)
BOOT_N <- 1000

# ----------------------
# 1) Install & load packages
# ----------------------
need <- c("MASS","psych","car","lavaan","brant","semTools","ggplot2","reshape2")
to_install <- need[!need %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(MASS)
  library(psych)
  library(car)
  library(lavaan)
  library(brant)
  library(semTools)
  library(ggplot2)
  library(reshape2)
})

section <- function(title) {
  cat("\\n\\n=============================\\n")
  cat(title, "\\n")
  cat("=============================\\n")
}

# ----------------------
# 2) Load & rename
# ----------------------
section("2) Load data & rename")
if (!file.exists(DATA_PATH)) stop("Cannot find ftis_R.csv in the working directory.")
data_raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# Required columns in uploaded file:
required_src <- c("Env_Median","Layout_Median","Staff_Median","Program_Median","Satis_Median","Loyal_Median")
missing_src <- setdiff(required_src, names(data_raw))
if (length(missing_src) > 0) {
  stop(paste0("Missing expected columns in ftis_R.csv: ", paste(missing_src, collapse = ", ")))
}

# Rename to analysis names
data <- data_raw[, required_src]
names(data) <- c("environment","layout","staff","program","satisfaction","loyalty")

print(head(data))

# ----------------------
# 3) Coerce to ordered Likert
# ----------------------
section("3) Coerce to ordered Likert (1-4)")
for (v in names(data)) {
  if (all(na.omit(data[[v]]) %in% LIKERT_LEVELS)) {
    data[[v]] <- ordered(data[[v]], levels = LIKERT_LEVELS)
  } else if (is.numeric(data[[v]])) {
    # Try rounding and clipping to 1..4
    z <- pmin(pmax(round(data[[v]]), min(LIKERT_LEVELS)), max(LIKERT_LEVELS))
    unique_vals <- sort(unique(na.omit(z)))
    if (all(unique_vals %in% LIKERT_LEVELS)) {
      data[[v]] <- ordered(z, levels = LIKERT_LEVELS)
      message(paste("Rounded/coerced", v, "to Likert 1..4."))
    } else {
      # Fallback: quantile bin into 4 ordered categories
      q <- quantile(data[[v]], probs = c(.25,.5,.75), na.rm = TRUE)
      data[[v]] <- cut(data[[v]], breaks = c(-Inf, q, Inf), labels = LIKERT_LEVELS, ordered_result = TRUE)
      message(paste("Binned", v, "into 4 ordered categories via quartiles."))
    }
  } else {
    stop(paste("Variable", v, "is not numeric or 1..4 integers; please recode to 1..4 before running."))
  }
}
print(sapply(data, function(x) table(x, useNA = "ifany")))

# ----------------------
# 4) Descriptives
# ----------------------
section("4) Descriptives")
print(psych::describe(as.data.frame(lapply(data, as.numeric))))

# ----------------------
# 5) Mediator model
# ----------------------
section("5) Mediator model (polr): satisfaction ~ environment + layout + staff + program")
model_M <- MASS::polr(satisfaction ~ environment + layout + staff + program,
                      data = data, Hess = TRUE)
summary(model_M)
ct_M <- coef(summary(model_M))
p_M  <- pnorm(abs(ct_M[, "t value"]), lower.tail = FALSE) * 2
ct_M <- cbind(ct_M, "p value" = p_M)
print(ct_M)

# ----------------------
# 6) Outcome model
# ----------------------
section("6) Outcome model (polr): loyalty ~ environment + layout + staff + program + satisfaction")
model_Y <- MASS::polr(loyalty ~ environment + layout + staff + program + satisfaction,
                      data = data, Hess = TRUE)
summary(model_Y)
ct_Y <- coef(summary(model_Y))
p_Y  <- pnorm(abs(ct_Y[, "t value"]), lower.tail = FALSE) * 2
ct_Y <- cbind(ct_Y, "p value" = p_Y)
print(ct_Y)
cat("\\nAIC (Outcome model): ", AIC(model_Y), "\\n")

# ----------------------
# 7) Proportional odds assumption (Brant)
# ----------------------
section("7) Brant test (proportional odds)")
bt <- brant::brant(model_Y)
print(bt)

# ----------------------
# 8) Multicollinearity proxy (VIF on numeric Y ~ X linear model)
# ----------------------
section("8) Multicollinearity (VIF proxy)")
lm_proxy <- lm(as.numeric(loyalty) ~ environment + layout + staff + program + satisfaction, data = data)
print(car::vif(lm_proxy))

# ----------------------
# 9) Mediation via lavaan (WLSMV, ordered), with bootstrap CIs
# ----------------------
section("9) Mediation via lavaan (WLSMV)")

model_med <- '
  # outcome
  loyalty ~ c1*environment + c2*layout + c3*staff + c4*program + b*satisfaction

  # mediator
  satisfaction ~ a1*environment + a2*layout + a3*staff + a4*program

  # indirect effects
  ind_env     := a1*b
  ind_layout  := a2*b
  ind_staff   := a3*b
  ind_program := a4*b

  # total effects
  tot_env     := c1 + ind_env
  tot_layout  := c2 + ind_layout
  tot_staff   := c3 + ind_staff
  tot_program := c4 + ind_program
'

fit <- lavaan::sem(model_med,
                   data = data,
                   ordered = c("satisfaction","loyalty"),
                   estimator = "WLSMV",
                   std.lv = TRUE)

section("9a) lavaan summary (fit, standardized, R2)")
print(summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE))

section("9b) Key fit indices")
print(fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr")))

section("9c) Parameter estimates (standardized + CI)")
print(parameterEstimates(fit, standardized = TRUE, ci = TRUE))

section("9d) Bootstrap CIs for indirect & total effects")
# Use the global BOOT_N if it exists; otherwise default to 1000
BOOT_N <- if (exists("BOOT_N")) BOOT_N else 1000

run_boot <- function(nboot) {
  lavaan::sem(
    model_med,
    data = data,
    ordered = c("satisfaction","loyalty"),
    estimator = "DWLS",            # <-- DWLS required for bootstrap with categorical
    se = "bootstrap",
    bootstrap = nboot,
    std.lv = TRUE,
    test = "mean.var.adjusted",    # similar spirit to WLSMV
    missing = "listwise"
  )
}

# Try bootstrapping with fallback sizes
fit.boot <- try(run_boot(BOOT_N), silent = TRUE)
if (inherits(fit.boot, "try-error")) {
  warning(paste("DWLS bootstrap failed at", BOOT_N, "resamples. Retrying with 300..."))
  fit.boot <- try(run_boot(300), silent = TRUE)
}
if (inherits(fit.boot, "try-error")) {
  warning("Retry with 300 failed. Retrying with 100...")
  fit.boot <- try(run_boot(100), silent = TRUE)
}

med_labels <- c("ind_env","ind_layout","ind_staff","ind_program",
                "tot_env","tot_layout","tot_staff","tot_program")

if (!inherits(fit.boot, "try-error")) {
  # Successful bootstrap: show percentile CI table
  pe.boot <- parameterEstimates(
    fit.boot,
    standardized = TRUE,
    boot.ci.type = "perc"
  )
  med_tbl <- subset(pe.boot, label %in% med_labels)
  if (nrow(med_tbl) == 0) {
    message("Note: mediation labels not found. Check that model_med defines ind_* and tot_* labels.")
  }
  print(med_tbl)
} else {
  # Fallback: non-bootstrap DWLS (still valid; just not resampled)
  section("9d-alt) Bootstrap failed repeatedly – falling back to non-bootstrap SEs (DWLS)")
  fit.noboot <- lavaan::sem(
    model_med,
    data = data,
    ordered = c("satisfaction","loyalty"),
    estimator = "DWLS",
    se = "standard",
    std.lv = TRUE,
    test = "mean.var.adjusted",
    missing = "listwise"
  )
  pe <- parameterEstimates(fit.noboot, standardized = TRUE, ci = TRUE)
  med_tbl <- subset(pe, label %in% med_labels)
  if (nrow(med_tbl) == 0) {
    message("Note: mediation labels not found. Check that model_med defines ind_* and tot_* labels.")
  }
  print(med_tbl)
}


#fit.boot <- lavaan::sem(model_med,
                        #data = data,
                        #ordered = c("satisfaction","loyalty"),
                        #estimator = "WLSMV",
                        #se = "bootstrap",
                        #bootstrap = BOOT_N,
                        #std.lv = TRUE)
#indirect <- parameterEstimates(fit.boot, boot.ci.type = "perc", standardized = TRUE)
#print(indirect[indirect$label %in% c("ind_env","ind_layout","ind_staff","ind_program",
                                     #"tot_env","tot_layout","tot_staff","tot_program"), ])

# ----------------------
# 10) Plots
# ----------------------
section("10) Plots (saved to ./plots)")
if (!dir.exists("plots")) dir.create("plots")

# Satisfaction vs Loyalty (numeric Y)
p1 <- ggplot(data, aes(x = satisfaction, y = as.numeric(loyalty))) +
  geom_jitter(width = 0.15, height = 0.05, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Satisfaction vs Loyalty", y = "Loyalty (numeric)", x = "Satisfaction (ordinal)") +
  theme_minimal()
ggsave(filename = "plots/satisfaction_vs_loyalty.png", plot = p1, width = 7, height = 5, dpi = 150)
cat("Saved: plots/satisfaction_vs_loyalty.png\\n")

# Predicted probabilities for varying satisfaction (others at medians)
med_env    <- median(as.numeric(data$environment), na.rm = TRUE)
med_layout <- median(as.numeric(data$layout), na.rm = TRUE)
med_staff  <- median(as.numeric(data$staff), na.rm = TRUE)
med_program<- median(as.numeric(data$program), na.rm = TRUE)

newdat <- data.frame(
  environment = ordered(rep(med_env, length(LIKERT_LEVELS)), levels = LIKERT_LEVELS),
  layout      = ordered(rep(med_layout, length(LIKERT_LEVELS)), levels = LIKERT_LEVELS),
  staff       = ordered(rep(med_staff, length(LIKERT_LEVELS)), levels = LIKERT_LEVELS),
  program     = ordered(rep(med_program, length(LIKERT_LEVELS)), levels = LIKERT_LEVELS),
  satisfaction= ordered(LIKERT_LEVELS, levels = LIKERT_LEVELS)
)

pred <- predict(model_Y, newdata = newdat, type = "probs")
pred_df <- cbind(newdat["satisfaction"], as.data.frame(pred))
names(pred_df) <- c("satisfaction", paste0("P(Y=", LIKERT_LEVELS, ")"))
print(pred_df)

pred_long <- melt(pred_df, id.vars = "satisfaction", variable.name = "category", value.name = "prob")
p2 <- ggplot(pred_long, aes(x = satisfaction, y = prob, group = category)) +
  geom_line(aes(linetype = category)) + geom_point() +
  labs(title = "Predicted P(Y = k) across Satisfaction levels", y = "Probability", x = "Satisfaction") +
  theme_minimal()
ggsave(filename = "plots/predicted_probs_by_satisfaction.png", plot = p2, width = 7, height = 5, dpi = 150)
cat("Saved: plots/predicted_probs_by_satisfaction.png\\n")

section("DONE")
cat("Outputs:\\n- polr mediator/outcome tables\\n- Brant test\\n- lavaan fit + bootstrap CIs\\n- plots in ./plots\\n")
