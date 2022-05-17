## ---- fig.cap="Illustration of a Dynamic Treatment Regime in a Clinical Setting", echo=FALSE, eval=TRUE, out.width='60%'----
knitr::include_graphics(path = "img/png/DynamicA_Illustration.png")


## ----setup-mopttx, message=FALSE, warning=FALSE-------------------------------
library(here)
library(data.table)
library(sl3)
library(tmle3)
library(tmle3mopttx)
library(devtools)
set.seed(116)


## ----load sim_bin_data--------------------------------------------------------
data("data_bin")


## ----data_nodes2-mopttx-------------------------------------------------------
# organize data and nodes for tmle3
data <- data_bin
node_list <- list(
  W = c("W1", "W2", "W3"),
  A = "A",
  Y = "Y"
)


## ----mopttx_sl3_lrnrs2--------------------------------------------------------
# Define sl3 library and metalearners:
lrn_xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm_fast$new()

## Define the Q learner:
Q_learner <- Lrnr_sl$new(
  learners = list(
    lrn_xgboost_50, lrn_xgboost_100,
    lrn_xgboost_500, lrn_mean, lrn_glm
  ),
  metalearner = Lrnr_nnls$new()
)

## Define the g learner:
g_learner <- Lrnr_sl$new(
  learners = list(lrn_xgboost_100, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

## Define the B learner:
b_learner <- Lrnr_sl$new(
  learners = list(
    lrn_xgboost_50, lrn_xgboost_100,
    lrn_xgboost_500, lrn_mean, lrn_glm
  ),
  metalearner = Lrnr_nnls$new()
)


## ----mopttx_make_lrnr_list----------------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----mopttx_spec_init_complex-------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3"), type = "blip1",
  learners = learner_list,
  maximize = TRUE, complex = TRUE,
  realistic = FALSE
)


## ----mopttx_fit_tmle_auto_blip_revere_complex, eval=T-------------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----load sim_cat_data--------------------------------------------------------
data("data_cat_realistic")


## ----data_nodes-mopttx--------------------------------------------------------
# organize data and nodes for tmle3
data <- data_cat_realistic
node_list <- list(
  W = c("W1", "W2", "W3", "W4"),
  A = "A",
  Y = "Y"
)


## ----sl3_lrnrs-mopttx---------------------------------------------------------
# Initialize few of the learners:
lrn_xgboost_50 <- Lrnr_xgboost$new(nrounds = 50)
lrn_xgboost_100 <- Lrnr_xgboost$new(nrounds = 100)
lrn_xgboost_500 <- Lrnr_xgboost$new(nrounds = 500)
lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm_fast$new()

## Define the Q learner, which is just a regular learner:
Q_learner <- Lrnr_sl$new(
  learners = list(lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_500, lrn_mean, lrn_glm),
  metalearner = Lrnr_nnls$new()
)

# Define the g learner, which is a multinomial learner:
# specify the appropriate loss of the multinomial learner:
mn_metalearner <- make_learner(Lrnr_solnp,
  loss_function = loss_loglik_multinomial,
  learner_function = metalearner_linear_multinomial
)
g_learner <- make_learner(Lrnr_sl, list(lrn_xgboost_100, lrn_xgboost_500, lrn_mean), mn_metalearner)

# Define the Blip learner, which is a multivariate learner:
learners <- list(lrn_xgboost_50, lrn_xgboost_100, lrn_xgboost_500, lrn_mean, lrn_glm)
b_learner <- create_mv_learners(learners = learners)


## ----cat_learners-------------------------------------------------------------
# See which learners support multi-class classification:
sl3_list_learners(c("categorical"))


## ----make_lrnr_list-mopttx----------------------------------------------------
# specify outcome and treatment regressions and create learner list
learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----spec_init----------------------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W1", "W2", "W3", "W4"), type = "blip2",
  learners = learner_list, maximize = TRUE, complex = TRUE,
  realistic = FALSE
)


## ----fit_tmle_auto, eval=T----------------------------------------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit

# How many individuals got assigned each treatment?
table(tmle_spec$return_rule)


## ----mopttx_spec_init_noncomplex----------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = FALSE, realistic = FALSE
)


## ----mopttx_fit_tmle_auto_blip_revere_noncomplex, eval=T----------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit


## ----mopttx_spec_init_realistic-----------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_blip_revere(
  V = c("W4", "W3", "W2", "W1"), type = "blip2",
  learners = learner_list,
  maximize = TRUE, complex = TRUE, realistic = TRUE
)


## ----mopttx_fit_tmle_auto_blip_revere_realistic, eval=T-----------------------
# fit the TML estimator
fit <- tmle3(tmle_spec, data, node_list, learner_list)
fit

# How many individuals got assigned each treatment?
table(tmle_spec$return_rule)


## ----data_vim-nodes-mopttx----------------------------------------------------
# bin baseline covariates to 3 categories:
data$W1 <- ifelse(data$W1 < quantile(data$W1)[2], 1, ifelse(data$W1 < quantile(data$W1)[3], 2, 3))

node_list <- list(
  W = c("W3", "W4", "W2"),
  A = c("W1", "A"),
  Y = "Y"
)


## ----mopttx_spec_init_vim-----------------------------------------------------
# initialize a tmle specification
tmle_spec <- tmle3_mopttx_vim(
  V = c("W2"),
  type = "blip2",
  learners = learner_list,
  contrast = "multiplicative",
  maximize = FALSE,
  method = "SL",
  complex = TRUE,
  realistic = FALSE
)


## ----mopttx_fit_tmle_auto_vim, eval=TRUE--------------------------------------
# fit the TML estimator
vim_results <- tmle3_vim(tmle_spec, data, node_list, learner_list,
  adjust_for_other_A = TRUE
)

print(vim_results)


## ----load-washb-data, message=FALSE, warning=FALSE, cache=FALSE, eval=FALSE----
## washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv", stringsAsFactors = TRUE)
## washb_data <- washb_data[!is.na(momage), lapply(.SD, as.numeric)]
## head(washb_data, 3)


## ----washb-data-npsem-mopttx, message=FALSE, warning=FALSE, cache=FALSE, eval=FALSE----
## node_list <- list(
##   W = names(washb_data)[!(names(washb_data) %in% c("whz", "tr"))],
##   A = "tr", Y = "whz"
## )


## ----summary_WASH, eval=FALSE-------------------------------------------------
## # V1, V2 and V3:
## table(washb_data$momedu)
## table(washb_data$floor)
## table(washb_data$asset_refrig)
## 
## # A:
## table(washb_data$tr)
## 
## # Y:
## summary(washb_data$whz)


## ----sl3_lrnrs-WASH-mopttx, eval=FALSE----------------------------------------
## # Initialize few of the learners:
## grid_params <- list(
##   nrounds = c(100, 500),
##   eta = c(0.01, 0.1)
## )
## grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
## xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
##   do.call(Lrnr_xgboost$new, c(as.list(params_tune)))
## })
## lrn_mean <- Lrnr_mean$new()
## 
## ## Define the Q learner, which is just a regular learner:
## Q_learner <- Lrnr_sl$new(
##   learners = list(
##     xgb_learners[[1]], xgb_learners[[2]], xgb_learners[[3]],
##     xgb_learners[[4]], lrn_mean
##   ),
##   metalearner = Lrnr_nnls$new()
## )
## 
## # Define the g learner, which is a multinomial learner:
## # specify the appropriate loss of the multinomial learner:
## mn_metalearner <- make_learner(Lrnr_solnp,
##   loss_function = loss_loglik_multinomial,
##   learner_function = metalearner_linear_multinomial
## )
## g_learner <- make_learner(Lrnr_sl, list(xgb_learners[[4]], lrn_mean), mn_metalearner)
## 
## # Define the Blip learner, which is a multivariate learner:
## learners <- list(
##   xgb_learners[[1]], xgb_learners[[2]], xgb_learners[[3]],
##   xgb_learners[[4]], lrn_mean
## )
## b_learner <- create_mv_learners(learners = learners)
## 
## learner_list <- list(Y = Q_learner, A = g_learner, B = b_learner)


## ----spec_init_WASH, eval=FALSE-----------------------------------------------
## ## Question 2:
## 
## # Initialize a tmle specification
## tmle_spec <- tmle3_mopttx_blip_revere(
##   V = c("momedu", "floor", "asset_refrig"), type = "blip2",
##   learners = learner_list, maximize = TRUE, complex = TRUE,
##   realistic = FALSE
## )
## 
## # Fit the TML estimator.
## fit <- tmle3(tmle_spec, data = washb_data, node_list, learner_list)
## fit
## 
## # Which intervention is the most dominant?
## table(tmle_spec$return_rule)


## ----spec_init_WASH_simple_q3, eval=FALSE-------------------------------------
## ## Question 3:
## 
## # Initialize a tmle specification with "realistic=TRUE":
## tmle_spec <- tmle3_mopttx_blip_revere(
##   V = c("momedu", "floor", "asset_refrig"), type = "blip2",
##   learners = learner_list, maximize = TRUE, complex = TRUE,
##   realistic = TRUE
## )
## 
## # Fit the TML estimator.
## fit <- tmle3(tmle_spec, data = washb_data, node_list, learner_list)
## fit
## 
## table(tmle_spec$return_rule)

