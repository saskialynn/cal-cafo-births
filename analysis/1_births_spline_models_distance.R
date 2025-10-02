
#..................................................
# Script:   1_births_spline_models_df.R
# 
# Purpose:  Model PTB as a function of distance to nearest CAFO 
#           using natural splines, with model selection based on AIC. 
#           Models include: 
#           - Overall spline model
#           - Race-stratified spline model
#           - Sex-stratified spline model
#           - Comparison to binned exposure logistic regression
#   Predictions are computed on a fine grid of distances with corresponding 
#   log-odds ratios, odds ratios, and confidence intervals relative to unexposed 
#   reference groups, stratified where appropriate.
#
# Inputs:
#   - receiver_dat_name: cleaned births dataset with nearest CAFO exposure distances
#   - receiver_filter_name: subset of births to define the study population
#
# Outputs:
#   - .rds file with AIC results (best df, knots, model comparison)
#   - .rds file with spline predictions, stratified predictions, and binned results
#
# Helper Functions: 
#   - compute_OR_CIs()
#
# Dependencies: 
#   - R packages: here, dplyr, splines, broom
#   - Custom scripts: 
#       - code/setup_header.R
#       - code/functions/spline_helper_fxns.R
#
#..................................................

##########
# Setup
##########
library(here)
source(file.path(here(), "cal-cafo-births/code/setup_header.R"))
source(file.path(here(), "cal-cafo-births/code/analysis/functions/spline_helper_fxns.R"))
library(dplyr)
library(splines)
library(broom)

#################################
# Specify Parameters
#################################
model_name <- "spline"

# Specify outcome and exposure variables
outcome_var_raw  <- "prem3s"
exposure_var_raw <- "nearest_cafo_km" 
ptb_level        <- 2 # PTB level to model (very/early PTB (1), moderate (2))

if(ptb_level == 1){ptb_name = "early"}
if(ptb_level == 2){ptb_name = "moderate"}

# Specify covariates and stratifying variables
covariates_continuous <- c("mothage")
covariates_cat        <- c("edum")
strat_var             <- c("racem", "sex")

covariates    <- c(covariates_continuous, covariates_cat)
covariates_eq <- paste0(covariates, collapse = "+")

# Specify spline parameters
max_df       <- 10
ref_distance <- 10    # units: kilometers
quantile_num <- 3     # for binned exposure - how many quantiles below ref_distance?
bk_lower     <- "min" # lower boundary knot
bk_upper     <- "max" # upper boundary knot
spline_optional_name <- NULL 

#################################
# Load data
#################################
receiver_dat_name    <- "births_main_analysis_df_cafo20250212_wedgeangle90_maxdist25000_NEARESTCAFO_ALLRECEIVERS_WINDEXP.csv"
receiver_filter_name <- "births_main_analysis_df.csv"
save_results_folder  <- "births_nullparous_dropduplicates_scorege50/10kmref"
data_optional_name   <- NULL  # NULL #"_NC"

receiver        <- fread(file.path(server_cleandata_path, receiver_dat_name))
receiver_filter <- fread(file.path(server_cleandata_path, receiver_filter_name))
receiver        <- receiver %>% dplyr::filter(birth_id %in% receiver_filter$birth_id)

################
# Prepare Data
################
# drop obs that are NA for outcome of interest or any covariates
receiver <- receiver %>% filter(!is.na(!!sym(outcome_var_raw)))
receiver <- receiver %>% drop_na(all_of(covariates))
receiver <- receiver %>% drop_na(all_of(strat_var))

# Factor variables, scale continuous variables
receiver <- receiver %>%
  mutate(case = factor(case, levels = c(0, 1)),
         across(all_of(covariates_cat), as.factor), 
         across(all_of(outcome_var_raw), as.factor), 
         across(all_of(covariates_continuous), ~ as.vector(scale(.x)))) 

# Ensure reference levels for categorical covariates
receiver <- receiver %>%
  mutate(
    racem = relevel(factor(racem), ref = levels(factor(racem))[1]),
    sex   = relevel(factor(sex),   ref = levels(factor(sex))[1]),
    edum  = relevel(factor(edum),  ref = levels(factor(edum))[1])
  )

# Subset: term births and ptb category of interest
receiver_subset <- receiver %>% filter(case == 0 | prem3s == ptb_level)

# Set boundary knots
dist_range <- range(receiver_subset[[exposure_var_raw]], na.rm = TRUE) 
boundary_knots <- c(
  ifelse(bk_lower == "min", dist_range[1], bk_lower),
  ifelse(bk_upper == "max", dist_range[2], bk_upper)
)

#############################################
# Select spline df using AIC
#############################################
aic_list <- list()
model_list <- list()

for (df in 1:max_df) { # Loop over degrees of freedom
  print(paste0("df = ", df))
  
  # Fit spline basis with specified df
  spline_ns_basis <- ns(receiver_subset[[exposure_var_raw]], 
                        df = df,
                        Boundary.knots = boundary_knots)
  colnames(spline_ns_basis) <- paste0("ns_basis.", 1:ncol(spline_ns_basis))
  
  #  Add spline basis to dataset
  receiver_subset_spline <- cbind(receiver_subset, spline_ns_basis)
  
  # Build formula string
  formula_string <- paste0("case_numeric", " ~ ", 
                           paste0("ns_basis.", 1:ncol(spline_ns_basis), collapse = " + "),
                           " + ", covariates_eq)
  formula <- as.formula(formula_string)
  
  # Fit logistic regression 
  model <- glm(formula, family = binomial, data = receiver_subset_spline)
  
  aic_list[[df]] <- list(
    df             = df,
    aic            = AIC(model),
    interior_knots = attr(spline_ns_basis, "knots"),
    boundary_knots = attr(spline_ns_basis, "Boundary.knots")
  )
}

aic_df <- bind_rows(lapply(aic_list, function(x) {
  tibble(
    df             = x$df,
    aic            = x$aic,
    interior_knots = list(x$interior_knots),
    boundary_knots = list(x$boundary_knots)
  )
}))

# Linear model comparison
linear_formula_string <- paste0("case_numeric ~ ", 
                                exposure_var_raw, "+",
                                covariates_eq )
linear_formula <- as.formula(linear_formula_string)
model_linear_glm <- glm(linear_formula, family = binomial, data = receiver_subset)

# Save results
all_aic_results <- list(
  boundary_knots       = boundary_knots,
  aic_df               = aic_df, 
  aic_model_linear_glm = AIC(model_linear_glm)
)

save_name <- paste0("results_MODEL_", model_name, "_AICvsDF",
                    "_OUTCOME_", outcome_var_raw, "_", ptb_name,  
                    "_EXP_", exposure_var_raw, 
                    "_COVAR_", paste0(covariates, collapse = "_"), 
                    spline_optional_name,
                    data_optional_name, ".rds")

saveRDS(all_aic_results, 
        file.path(transfer_modelresults_path, save_results_folder, save_name))

# Re-Load final AIC results
all_aic_results <- readRDS(file.path(transfer_modelresults_path, 
                                     save_results_folder, 
                                     save_name))
aic_df <- all_aic_results$aic_df

###############################################################################
# Fit Spline Models on Observed Data 
###############################################################################
# Select best df and corresponding knots
aic_best <- aic_df[which.min(aic_df$aic), ]
df_best <- aic_best$df
df_best_knots <- aic_best[["interior_knots"]][[1]]

#.................Compute spline basis on the OBSERVED data...............#
spline_ns_basis <- ns(receiver_subset[[exposure_var_raw]], 
                      knots = df_best_knots, 
                      Boundary.knots = boundary_knots)
spline_column_names <-  paste0("ns_basis.", 1:ncol(spline_ns_basis))
colnames(spline_ns_basis) <- spline_column_names

receiver_subset_final_spline <- cbind(receiver_subset, spline_ns_basis)

#.............Fit FULL MODEL logistic regression using spline basis............#
formula_full <- paste0(
  "case_numeric", " ~ ", 
  paste0("ns_basis.", 1:ncol(spline_ns_basis), collapse = " + "),
  " + ", covariates_eq)

model_ns_glm <- glm(as.formula(formula_full), 
                    family = binomial, 
                    data = receiver_subset_final_spline )

#............Fit RACE STRATIFIED logistic regression using spline..............#
# Add race main effects, spline main effects, race*spline interactions
formula_modelmat_race <- paste0("~ racem * (", paste(spline_column_names, collapse = " + "), ")")

spline_race_interactions <- model.matrix(as.formula(formula_modelmat_race), 
                                         data = receiver_subset_final_spline)

receiver_subset_final_spline_race <- cbind(receiver_subset, spline_race_interactions)

# Allow R to drop reference race level to avoid multicollinearity
formula_race <- paste0("case_numeric ~ racem * (", 
                       paste0(colnames(spline_ns_basis), collapse = " + "), ") + ", 
                       covariates_eq)

model_ns_glm_race <- glm(as.formula(formula_race), 
                         family = binomial, 
                         data = receiver_subset_final_spline_race)

#............Fit SEX STRATIFIED logistic regression using spline..............#
# Add sex main effects, spline main effects, sex*spline interactions
formula_modelmat_sex <- paste0("~ sex * (", paste(spline_column_names, collapse = " + "), ")")

spline_sex_interactions <- model.matrix(as.formula(formula_modelmat_sex), 
                                        data = receiver_subset_final_spline)

receiver_subset_final_spline_sex <- cbind(receiver_subset, spline_sex_interactions)

formula_sex <- paste0("case_numeric ~ sex * (", 
                      paste0(colnames(spline_ns_basis), collapse = " + "), ") + ", 
                      covariates_eq)

model_ns_glm_sex <- glm(as.formula(formula_sex), 
                        family = binomial, 
                        data = receiver_subset_final_spline_sex)

###############################################################################
# PREDICT on fine grid of exposure values using fitted model
#   Model was fit on observed distances in previous step.
#   Now want to observe how predicted log odds changes continuously over range
#   of distance values (not just observed points).
###############################################################################
#..................Define reference covariate values.................#
# To ensure proper OR calculation, covariates must be held constant between
# the prediction grid and the unexposed reference group
mothage_ref <- mean(receiver_subset$mothage, na.rm = TRUE)
edum_ref    <- factor("2", levels = levels(receiver_subset$edum))

#.................Build prediction dataset over distance grid.................#
distance_seq <- seq(0, boundary_knots[2], length.out = 400)

# Compute spline basis for the grid using knots from the fitted model
spline_pred_basis <- ns(distance_seq, 
                        knots = df_best_knots, 
                        Boundary.knots = boundary_knots) %>% data.frame()
colnames(spline_pred_basis) <- paste0("ns_basis.", 1:ncol(spline_pred_basis))

# Combine spline basis with reference covariates for prediction
pred_data <- spline_pred_basis %>%
  mutate(distance_km = distance_seq, 
         mothage = mothage_ref, 
         edum = edum_ref)

#..........Build reference matrix using OBSERVED UNEXPOSED population..........#
# Subset to "unexposed" individuals 
unexposed_data <- receiver_subset %>% filter(!!sym(exposure_var_raw) > ref_distance)

# Compute spline basis for observed unexposed distances
unexposed_basis <- ns(unexposed_data[[exposure_var_raw]], 
                      knots = df_best_knots, 
                      Boundary.knots = boundary_knots) %>% data.frame()
colnames(unexposed_basis) <- paste0("ns_basis.", 1:ncol(unexposed_basis))

# Add reference covariates to unexposed basis
unexposed_model_mat <- unexposed_basis %>%
  mutate(racem = unexposed_data$racem,
         sex = unexposed_data$sex, 
         mothage = mothage_ref, 
         edum = edum_ref)

#####################################################
# Prediction: Overall Spline Model
#####################################################
#................. Model predicted log odds over distance grid.................#
preds <- predict(model_ns_glm, newdata = pred_data, type = "link", se.fit = TRUE)
pred_data$log_odds <- preds$fit

#............Reference log odds for OBSERVED UNEXPOSED population..............#
unexposed_preds <- predict(model_ns_glm,  newdata = unexposed_model_mat, type = "link", se.fit = TRUE)
mean_log_odds_baseline <- mean(unexposed_preds$fit) # Point estimate of mean log-odds for reference group

#.........Compute log-odds ratio, odds ratio, and confidence intervals.........#
pred_data <- compute_OR_CIs(fitted_model   = model_ns_glm, 
                            unexposed_data = unexposed_model_mat, 
                            pred_data      = pred_data, 
                            log_odds_ref   = mean_log_odds_baseline,
                            strat_var      = NULL)

#####################################################
# Prediction: Race Stratified Spline Model
#####################################################
#.................Build prediction dataset over distance grid.................#
# Grid: every distance × every race
race_levels    <- levels(receiver_subset$racem)
pred_data_race <- expand.grid(distance_km = distance_seq, racem = race_levels)

# Expand spline basis to match the stratified prediction grid
pred_data_race <- cbind(pred_data_race, 
                        spline_pred_basis[rep(1:nrow(spline_pred_basis), 
                                              times = length(race_levels)), ])
# Add covariates at reference values
pred_data_race <- pred_data_race %>% 
  mutate(
    racem = factor(racem, levels = race_levels), 
    mothage = mothage_ref, 
    edum = edum_ref)

#..............Predict log-odds over the grid of distance x race...............#
preds_race <- predict(model_ns_glm_race, newdata = pred_data_race, 
                      type = "link", se.fit = TRUE)
pred_data_race$log_odds <- preds_race$fit

#.......Compute reference log odds for OBSERVED UNEXPOSED population.........#
# Use the single shared spline basis for all unexposed individuals, regardless of race.
# This is consistent with how the model was originally fit:
# Model was fit using one set of global spline basis variables, 
# with race effects (race*global spline interactions) added on

# Predict log-odds for each unexposed person by race
unexposed_preds_race <- predict(model_ns_glm_race, newdata = unexposed_model_mat,
                                type = "link", se.fit = TRUE)

#.........Compute log-odds ratio, odds ratio, and confidence intervals.........#
pred_data_race <- compute_OR_CIs(fitted_model   = model_ns_glm_race, 
                                 unexposed_data = unexposed_model_mat, 
                                 pred_data      = pred_data_race, 
                                 log_odds_ref   = unexposed_preds_race$fit,
                                 strat_var      = "racem")

#####################################################
# Prediction: Sex-Stratified Spline Model
#####################################################
# Grid: every distance × sex
sex_levels    <- levels(receiver_subset$sex)
pred_data_sex <- expand.grid(distance_km = distance_seq, sex = sex_levels)

# Expand spline basis to match the stratified prediction grid
pred_data_sex <- cbind(pred_data_sex, 
                       spline_pred_basis[rep(1:nrow(spline_pred_basis), 
                                             times = length(sex_levels)), ])
pred_data_sex <- pred_data_sex %>% 
  mutate(
    sex = factor(sex, levels = sex_levels), 
    mothage = mothage_ref, 
    edum = edum_ref)

#..............Predict log-odds over the grid of distance x sex...............#
preds_sex <- predict(model_ns_glm_sex, newdata = pred_data_sex,
                     type = "link", se.fit = TRUE)
pred_data_sex$log_odds <- preds_sex$fit


#.......Compute reference log odds for OBSERVED UNEXPOSED population.........#
unexposed_preds_sex <- predict(model_ns_glm_sex, newdata = unexposed_model_mat,
                               type = "link", se.fit = TRUE)

#.........Compute log-odds ratio, odds ratio, and confidence intervals.........#
pred_data_sex <- compute_OR_CIs(fitted_model = model_ns_glm_sex, 
                                unexposed_data = unexposed_model_mat, 
                                pred_data = pred_data_sex, 
                                log_odds_ref = unexposed_preds_sex$fit,
                                strat_var = "sex")

#######################################################
# Comparison to binned exposure logistic regression
#######################################################
#............................Define Exposure Bins..............................#
# Subset exposed values within reference distance
exposed_values <- receiver[[exposure_var_raw]][receiver[[exposure_var_raw]] <= ref_distance]

# Compute quantiles of nonzero exposures
quantile_probs <- seq(1 / quantile_num, (quantile_num - 1) / quantile_num, by = 1 / quantile_num)
quantile_cutoffs <- quantile(exposed_values, probs = quantile_probs, na.rm = TRUE)

# Define bin breaks: 0 = no exposure, then quantile cutoffs, then max
breaks <- c(0, quantile_cutoffs, max(exposed_values, na.rm = TRUE))

# Define factor levels: 0 = no exposure (reference), then low > medium > high
bin_labels <- paste0("(", breaks[-length(breaks)], ",", breaks[-1], "]")
bin_levels <- c("0", rev(bin_labels))
bin_labels_ordered <- c(paste0(">", ref_distance), rev(bin_labels))
print(bin_labels_ordered)

#......................Assign births to exposure categories....................#
receiver_subset <- receiver_subset %>%
  mutate(
    exposure_cat = case_when(
      .data[[exposure_var_raw]] > ref_distance ~ "0",
      TRUE ~ {
        bins <- as.character(cut(
          .data[[exposure_var_raw]],
          breaks = breaks,
          labels = bin_labels,
          include.lowest = FALSE,
          right = TRUE
        ))
      }
    )
  )
table(receiver_subset$exposure_cat, useNA = "always")

# Convert to factor and order correctly for modeling
receiver_subset <- receiver_subset %>%
  mutate(
    exposure_cat = factor(exposure_cat,
                          levels = bin_levels,         # original numeric/factor levels
                          labels = bin_labels_ordered, # human-readable labels
                          ordered = FALSE              # Avoid polynomial contrasts
    )
  )

# Check distribution
table(receiver_subset$exposure_cat, useNA = "always")
levels(receiver_subset$exposure_cat)
# Reference: unexposed within threshold distance
# low exposure: farther minimum distance
# high exposure: closest minimum distance 

#....................Fit Logistic Regression Model.............................#
# Define exposure variable for modeling
exposure_var <- rlang::sym("exposure_cat")
exposure_var_char <- "exposure_cat"

# Adjusted logistic regression
adjusted_model <- glm(
  reformulate(c(exposure_var_char, covariates), response = "case_numeric"), # "case"
  data = receiver_subset,
  family = binomial(link = "logit")
)
adjusted_results <- tidy(adjusted_model)

#...........................Summarize exposure counts.........................#
# Count number of observations for PTB and term births by exposure quantile
exposure_counts <- receiver_subset %>%
  group_by(.data[[exposure_var]], case) %>%
  summarise(n = n(), .groups = 'drop')

# Count total number of PTB and term births
ptb_term_counts <- receiver_subset %>%
  group_by(case) %>%
  summarise(n = n(), .groups = 'drop')

#...........Combine and save all results......#
all_spline_results <- list(
  outcome_var_raw  = outcome_var_raw,
  ptb_level        = ptb_level,
  exposure_var_raw = exposure_var_raw,
  ref_distance     = ref_distance,
  boundary_knots   = boundary_knots,
  spline_results   = list(
    df_best        = df_best, 
    df_best_knots  = df_best_knots,
    pred_data      = pred_data
  ), 
  spline_results_race = list(
    df_best        = df_best, 
    df_best_knots  = df_best_knots, 
    pred_data_race = pred_data_race
  ),
  spline_results_sex = list(
    df_best        = df_best, 
    df_best_knots  = df_best_knots, 
    pred_data_sex  = pred_data_sex
  ),
  binned_exp_comparison = list(
    quantile_num     = quantile_num,
    breaks           = breaks,
    bin_labels       = bin_labels,
    ptb_term_counts  = ptb_term_counts,
    exposure_counts  = exposure_counts,
    adjusted_results = adjusted_results
  ) 
)

save_name <- paste0("results_MODEL_", model_name,
                    "_OUTCOME_", outcome_var_raw, "_", ptb_name,  
                    "_EXP_", exposure_var_raw, 
                    "_COVAR_", paste0(covariates, collapse = "_"), 
                    "_REFDIST_", ref_distance, "km", 
                    spline_optional_name,
                    data_optional_name, ".rds")

saveRDS(all_spline_results, 
        file.path(transfer_modelresults_path, save_results_folder, save_name))


# STOP

