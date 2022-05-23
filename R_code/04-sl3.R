## ----setup-handbook-utils-noecho, echo = FALSE--------------------------------
library(knitr)
library(kableExtra)
library(data.table)


## ----load-data----------------------------------------------------------------
library(data.table)
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)


## ----show-data-normal-noeval, eval = FALSE------------------------------------
## head(washb_data)

## ----show-data-handbook, echo = FALSE-----------------------------------------
if (knitr::is_latex_output()) {
  head(washb_data) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(washb_data) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----install-sl3, eval = FALSE------------------------------------------------
## library(devtools)
## install_github("tlverse/sl3")


## ----load-sl3-----------------------------------------------------------------
library(sl3)


## ----task---------------------------------------------------------------------
# create the task (i.e., use washb_data to predict outcome using covariates)
task <- make_sl3_Task(
  data = washb_data,
  outcome = "whz",
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)

# let's examine the task
task


## ----list-properties----------------------------------------------------------
sl3_list_properties()


## ----list-learners------------------------------------------------------------
sl3_list_learners(properties = "continuous")


## ----learners-----------------------------------------------------------------
lrn_glm <- Lrnr_glm$new()
lrn_mean <- Lrnr_mean$new()


## ----more-learners------------------------------------------------------------
# penalized regressions:
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)


## ----more-learners-np---------------------------------------------------------
# spline regressions:
lrn_polspline <- Lrnr_polspline$new()
lrn_earth <- Lrnr_earth$new()

# fast highly adaptive lasso (HAL) implementation
lrn_hal <- Lrnr_hal9001$new(max_degree = 2, num_knots = c(3,2), nfolds = 5)

# tree-based methods
lrn_ranger <- Lrnr_ranger$new()
lrn_xgb <- Lrnr_xgboost$new()


## ----more-learners-final------------------------------------------------------
lrn_gam <- Lrnr_gam$new()
lrn_bayesglm <- Lrnr_bayesglm$new()


## ----stack--------------------------------------------------------------------
stack <- Stack$new(
  lrn_glm, lrn_mean, lrn_ridge, lrn_lasso, lrn_polspline, lrn_earth, lrn_hal, 
  lrn_ranger, lrn_xgb, lrn_gam, lrn_bayesglm
)
stack


## ----make-sl------------------------------------------------------------------
sl <- Lrnr_sl$new(learners = stack, metalearner = Lrnr_nnls$new())


## ----train-sl-----------------------------------------------------------------
set.seed(4197)
sl_fit <- sl$train(task)


## ----sl-predictions-----------------------------------------------------------
sl_preds <- sl_fit$predict(task)
head(sl_preds)


## ----lrnr-predictions---------------------------------------------------------
glm_preds <- sl_fit$learner_fits$Lrnr_glm_TRUE$predict(task)
head(glm_preds)


## ----predvobs-----------------------------------------------------------------
# table of observed and predicted outcome values and arrange by observed values
df_plot <- data.table(
  Obs = washb_data[["whz"]], SL_Pred = sl_preds, GLM_Pred = glm_preds,
  Mean_Pred = sl_fit$learner_fits$Lrnr_mean$predict(task)
)
df_plot <- df_plot[order(df_plot$Obs), ] 

## ----predvobs-head, eval = FALSE----------------------------------------------
## head(df_plot)

## ----predvobs-head-handbook, echo = FALSE-------------------------------------
if (knitr::is_latex_output()) {
  head(df_plot) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(df_plot) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}

## ----predobs-plot-------------------------------------------------------------
# melt the table so we can plot observed and predicted values
df_plot$id <- seq(1:nrow(df_plot))
df_plot_melted <- melt(
  df_plot, id.vars = "id",
  measure.vars = c("Obs", "SL_Pred", "GLM_Pred", "Mean_Pred")
)

library(ggplot2)
ggplot(df_plot_melted, aes(id, value, color = variable)) + 
  geom_point(size = 0.1) + 
  labs(x = "Subject id (ordered by increasing whz)", 
       y = "Weight-for-height z-score (whz)") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


## ----cv-predictions-----------------------------------------------------------
cv_preds <- sl_fit$fit_object$cv_fit$predict_fold(task, "validation")

## ----cv-predictions-head, eval = FALSE----------------------------------------
## head(cv_preds)

## ----cv-predictions-head-handbook, echo = FALSE-------------------------------
if (knitr::is_latex_output()) {
  head(cv_preds) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(cv_preds) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----predictions-new-task, eval = FALSE---------------------------------------
## washb_data_new$whz <- rep(NA, nrow(washb_data_new)) # create a fake outcome
## pred_task <- make_sl3_Task(
##   data = washb_data_new,
##   outcome = "whz",
##   outcome_type = "continuous",
##   covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu",
##                  "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec",
##                  "floor", "walls", "roof", "asset_wardrobe", "asset_table",
##                  "asset_chair", "asset_khat", "asset_chouki", "asset_tv",
##                  "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach",
##                  "asset_mobile")
## )
## sl_preds_new_task <- sl_fit$predict(pred_task)


## ----sl-coefs-simple----------------------------------------------------------
round(sl_fit$coefficients, 3)


## ----metalearner-fit----------------------------------------------------------
metalrnr_fit <- sl_fit$fit_object$cv_meta_fit$fit_object
round(metalrnr_fit$coefficients, 3)


## ----sl-summary---------------------------------------------------------------
cv_risk_table <- sl_fit$cv_risk(eval_fun = loss_squared_error)

## ----cv-risk-summary, eval = FALSE--------------------------------------------
## cv_risk_table[,c(1:3)]

## ----cv-risk-summary-handbook, echo = FALSE-----------------------------------
if (knitr::is_latex_output()) {
  cv_risk_table[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  cv_risk_table[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----sl-summary-plot, eval = F------------------------------------------------
## 
## # Column "se" in the CV risk table is the standard error across all losses for
## # a learner, i.e., se = sd(loss)/sqrt(n), where loss is an n length vector of
## # validation set predictions across all folds, and n is the number of
## # validation set observations across all folds. We can use this to
## cv_risk_table[, "lower" := MSE - qnorm(.975)*se]
## cv_risk_table[, "upper" := MSE + qnorm(.975)*se]
## 
## ggplot(cv_risk_table,
##        aes_string(x = "learner", y = "MSE", ymin = "lower", ymax = "upper")) +
##   geom_pointrange() +
##   coord_flip() +
##   ylab("V-fold CV Risk Estimate") +
##   xlab("Learner")


## ----cvsl---------------------------------------------------------------------
set.seed(569)
cv_sl_fit <- cv_sl(lrnr_sl = sl_fit, task = task, eval_fun = loss_squared_error)

## ----cvsl-risk-summary, eval = FALSE------------------------------------------
## cv_sl_fit$cv_risk[,c(1:3)]

## ----cvsl-risk-summary-handbook, echo = FALSE---------------------------------
if (knitr::is_latex_output()) {
  cv_sl_fit$cv_risk[,c(1:3)] %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  cv_sl_fit$cv_risk[,c(1:3)] %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----make-sl-discrete---------------------------------------------------------
cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
stack_with_sl <- Stack$new(stack, sl)
dSL <- Lrnr_sl$new(learners = stack_with_sl, metalearner = cv_selector)


## ----task-with-warning, warning=TRUE------------------------------------------
# create the task (i.e., use washb_data to predict outcome using covariates)
task <- make_sl3_Task(
  data = washb_data,
  outcome = "whz",
  covariates = c("tr", "fracode", "month", "aged", "sex", "momage", "momedu", 
                 "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin", "elec", 
                 "floor", "walls", "roof", "asset_wardrobe", "asset_table", 
                 "asset_chair", "asset_khat", "asset_chouki", "asset_tv", 
                 "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
                 "asset_mobile")
)


## ----which-data-missing-------------------------------------------------------
# which columns have missing values, and how many observations are missing?
colSums(is.na(washb_data))


## ----data-missing-------------------------------------------------------------
some_rows_with_missingness <- which(!complete.cases(washb_data))[31:33]
# note: we chose 31:33 because missingness in momage & momheight is there
washb_data[some_rows_with_missingness, c("momage", "momheight")]


## ----task-data-imputed--------------------------------------------------------
task$data[some_rows_with_missingness, 
          c("momage", "momheight", "delta_momage", "delta_momheight")]
colSums(is.na(task$data))


## ----fruit--------------------------------------------------------------------
cats <- c("calico", "tabby", "cow", "ragdoll", "mancoon", "dwarf", "calico")
cats <- factor(cats)
cats_onehot <- factor_to_indicators(cats)
cats_onehot


## ----show-X, eval = FALSE-----------------------------------------------------
## head(task$X)

## ----show-X-handbook, echo = FALSE--------------------------------------------
if (knitr::is_latex_output()) {
  head(task$X) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(task$X) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----stack-names--------------------------------------------------------------
stack


## ----name-glm, eval = FALSE---------------------------------------------------
## lrn_glm <- Lrnr_glm$new(name = "GLM")


## ----stack-pretty-------------------------------------------------------------
learners_pretty_names <- c(
  "GLM" = lrn_glm, "Mean" = lrn_mean, "Ridge" = lrn_ridge, 
  "Lasso" = lrn_lasso, "Polspline" = lrn_polspline, "Earth" = lrn_earth, 
  "HAL" = lrn_hal, "RF" = lrn_ranger, "XGBoost" = lrn_xgb, "GAM" = lrn_gam, 
  "BayesGLM" = lrn_bayesglm
)
stack_pretty_names <- Stack$new(learners_pretty_names)
stack_pretty_names


## ----lrnr-grid-diy------------------------------------------------------------
grid_params <- list(
  max_depth = c(3, 5, 8),
  eta = c(0.001, 0.1, 0.3),
  nrounds = 100
)
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_xgboost$new, as.list(tuning_params))
})
xgb_learners


## ----lrnr-grid-diy-names------------------------------------------------------
names(xgb_learners) <- c(
  "XGBoost_depth3_eta.001", "XGBoost_depth5_eta.001", "XGBoost_depth8_eta.001", 
  "XGBoost_depth3_eta.1", "XGBoost_depth5_eta.1", "XGBoost_depth8_eta.1", 
  "XGBoost_depth3_eta.3", "XGBoost_depth5_eta.3", "XGBoost_depth8_eta.3"
)


## ----lrnr-grid-caret, eval = FALSE--------------------------------------------
## lrnr_nnet_autotune <- Lrnr_caret$new(method = "nnet", name = "NNET_autotune")


## ----interaction-learner------------------------------------------------------
lrnr_glm_interaction <- Lrnr_glm$new(formula = "~.^2")


## ----screener-properties------------------------------------------------------
sl3_list_learners(properties = "importance")


## ----screener-importance------------------------------------------------------
ranger_with_importance <- Lrnr_ranger$new(importance = "impurity_corrected")
RFscreen_top10 <- Lrnr_screener_importance$new(
  learner = ranger_with_importance, num_screen = 10
)
RFscreen_top10_glm <- Pipeline$new(RFscreen_top10, lrn_glm)


## ----screener-importance-stack------------------------------------------------
RFscreen_top10_stack <- Pipeline$new(RFscreen_top10, stack)


## ----screener-coefs-----------------------------------------------------------
lasso_screen <- Lrnr_screener_coefs$new(learner = lrn_lasso, threshold = 0)
lasso_screen_glm <- Pipeline$new(lasso_screen, lrn_glm)


## ----screener-stack-----------------------------------------------------------
lasso_screen_stack <- Pipeline$new(lasso_screen, stack)


## ----screener-corr------------------------------------------------------------
# select top 10 most correlated covariates
corRank_screen <- Lrnr_screener_correlation$new(
  type = "rank", num_screen = 10
)
corRank_screen_stack <- Pipeline$new(corRank_screen, stack)

# select covariates with correlation p-value below 0.05, and a minimum of 3
corP_screen <- Lrnr_screener_correlation$new(
  type = "threshold", pvalue_threshold = 0.05, min_screen = 3
)
corP_screen_stack <- Pipeline$new(corP_screen, stack)


## ----screener-augment---------------------------------------------------------
keepme <- c("aged", "momage")
# using corRank_screen as an example, but any instantiated screener can be 
# supplied as screener.
corRank_screen_augmented <- Lrnr_screener_augment$new(
  screener = corRank_screen, default_covariates = keepme
)
corRank_screen_augmented_glm <- Pipeline$new(corRank_screen_augmented, lrn_glm)


## ----screeners-stack----------------------------------------------------------
screeners_stack <- Stack$new(stack, corP_screen_stack, corRank_screen_stack, 
                             lasso_screen_stack, RFscreen_top10_stack)


## ----varimp-------------------------------------------------------------------
assets <- c("asset_wardrobe", "asset_table", "asset_chair", "asset_khat",
            "asset_chouki", "asset_tv", "asset_refrig", "asset_bike", 
            "asset_moto", "asset_sewmach", "asset_mobile", "Nlt18", "Ncomp", 
            "watmin", "elec", "floor", "walls", "roof")
set.seed(983)
washb_varimp <- importance(
  fit = sl_fit, eval_fun = loss_squared_error, type = "permute", 
  covariate_groups = list("assets" = assets)
)


## ----varimp-print, eval = FALSE-----------------------------------------------
## washb_varimp

## ----varimp-print-handbook, echo = FALSE--------------------------------------
if (knitr::is_latex_output()) {
  washb_varimp %>%
    kable(digits = 4, format = "latex")
} else if (knitr::is_html_output()) {
  washb_varimp %>%
    kable(digits = 4) %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}


## ----varimp-plot, out.width = "100%"------------------------------------------
# plot variable importance
importance_plot(x = washb_varimp)


## ----cde_using_locscale, eval = FALSE-----------------------------------------
## # TODO: fix code block and flesh out into demo
## 
## # semiparametric density estimator with homoscedastic errors (HOSE)
## hose_hal_lrnr <- Lrnr_density_semiparametric$new(
##   mean_learner = Lrnr_hal9001$new()
## )
## # semiparametric density estimator with heteroscedastic errors (HESE)
## hese_rf_glm_lrnr <- Lrnr_density_semiparametric$new(
##   mean_learner = Lrnr_ranger$new()
##   var_learner = Lrnr_glm$new()
## )
## 
## # SL for the conditional treatment density
## sl_dens_lrnr <- Lrnr_sl$new(
##   learners = list(hose_hal_lrnr, hese_rf_glm_lrnr),
##   metalearner = Lrnr_solnp_density$new()
## )


## ----cde_using_pooledhaz, eval = FALSE----------------------------------------
## # TODO: fix code block and flesh out into demo
## 
## # learners used for conditional densities for (g_n)
## haldensify_lrnr <- Lrnr_haldensify$new(
##   n_bins = c(5, 10)
## )


## ----ex-setup-----------------------------------------------------------------
# load the data set
db_data <- url(
  paste0(
    "https://raw.githubusercontent.com/benkeser/sllecture/master/",
    "chspred.csv"
  )
)
chspred <- read_csv(file = db_data, col_names = TRUE)


## ----ex-head-handbook, echo = FALSE-------------------------------------------
if (knitr::is_latex_output()) {
  head(chspred) %>%
    kable(format = "latex")
} else if (knitr::is_html_output()) {
  head(chspred) %>%
    kable() %>%
    kableExtra:::kable_styling(fixed_thead = TRUE) %>%
    scroll_box(width = "100%", height = "300px")
}

## ----ex-head, eval = FALSE----------------------------------------------------
## head(chspred)


## ----ex-key, eval=FALSE-------------------------------------------------------
## db_data <- url(
##   "https://raw.githubusercontent.com/benkeser/sllecture/master/chspred.csv"
## )
## chspred <- read_csv(file = db_data, col_names = TRUE)
## data.table::setDT(chspred)
## 
## # make task
## chspred_task <- make_sl3_Task(
##   data = chspred,
##   covariates = colnames(chspred)[-1],
##   outcome = "mi"
## )
## 
## # make learners
## glm_learner <- Lrnr_glm$new()
## lasso_learner <- Lrnr_glmnet$new(alpha = 1)
## ridge_learner <- Lrnr_glmnet$new(alpha = 0)
## enet_learner <- Lrnr_glmnet$new(alpha = 0.5)
## # curated_glm_learner uses formula = "mi ~ smoke + beta"
## curated_glm_learner <- Lrnr_glm_fast$new(covariates = c("smoke", "beta"))
## mean_learner <- Lrnr_mean$new() # That is one mean learner!
## glm_fast_learner <- Lrnr_glm_fast$new()
## ranger_learner <- Lrnr_ranger$new()
## svm_learner <- Lrnr_svm$new()
## xgb_learner <- Lrnr_xgboost$new()
## 
## # screening
## screen_cor <- make_learner(Lrnr_screener_correlation)
## glm_pipeline <- make_learner(Pipeline, screen_cor, glm_learner)
## 
## # stack learners together
## stack <- make_learner(
##   Stack,
##   glm_pipeline, glm_learner,
##   lasso_learner, ridge_learner, enet_learner,
##   curated_glm_learner, mean_learner, glm_fast_learner,
##   ranger_learner, svm_learner, xgb_learner
## )
## 
## # make and train SL
## sl <- Lrnr_sl$new(
##   learners = stack
## )
## sl_fit <- sl$train(chspred_task)
## sl_fit$cv_risk(loss_squared_error)

