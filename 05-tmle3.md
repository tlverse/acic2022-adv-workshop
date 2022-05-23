# The TMLE Framework (Brief Review) {#tmle3}

_Jeremy Coyle_ and _Nima Hejazi_

Based on the [`tmle3` `R` package](https://github.com/tlverse/tmle3).

## Learning Objectives {#learn-tmle}

By the end of this chapter, you will be able to

1. Use `tmle3` to estimate an Average Treatment Effect (ATE).
2. Understand how to use `tmle3` "Specs" objects.

## Introduction {#tmle-intro}

Mark and Alan introduced the core concepts associated with TMLE in their intro talk. Today, we'll be focused on some more advanced applications of `tmle3`, but we'd like to review the basics of how to use the package. Before we do that, are there any conceptual clarifications on TMLE?

The following sections describe a simple way of
specifying and estimating a TMLE in the `tlverse`. In designing `tmle3`, we
sought to replicate as closely as possible the very general estimation framework
of TMLE, and so each theoretical object relevant to TMLE is encoded in a
corresponding software object/method. More information on this design can be found in the [handbook](http://tlverse.org/tlverse-handbook/tmle3.html#tmle3-components).

## Easy-Bake Example: `tmle3` for ATE

We'll illustrate the most basic use of TMLE using the WASH Benefits data
introduced earlier and estimating an average treatment effect. Similar specifications will be relevant during the later sections on advanced `tmle3` usage.

### Load the Data

We'll use the same WASH Benefits data as the earlier chapters:


```r
library(data.table)
library(dplyr)
library(tmle3)
library(sl3)
washb_data <- fread(
  paste0(
    "https://raw.githubusercontent.com/tlverse/tlverse-data/master/",
    "wash-benefits/washb_data.csv"
  ),
  stringsAsFactors = TRUE
)
```

### Define the variable roles

We'll use the common $W$ (covariates), $A$ (treatment/intervention), $Y$
(outcome) data structure. `tmle3` needs to know what variables in the dataset
correspond to each of these roles. We use a list of character vectors to tell
it. We call this a "Node List" as it corresponds to the nodes in a Directed
Acyclic Graph (DAG), a way of displaying causal relationships between variables.


```r
node_list <- list(
  W = c(
    "month", "aged", "sex", "momage", "momedu",
    "momheight", "hfiacat", "Nlt18", "Ncomp", "watmin",
    "elec", "floor", "walls", "roof", "asset_wardrobe",
    "asset_table", "asset_chair", "asset_khat",
    "asset_chouki", "asset_tv", "asset_refrig",
    "asset_bike", "asset_moto", "asset_sewmach",
    "asset_mobile"
  ),
  A = "tr",
  Y = "whz"
)
```

### Handle Missingness

Currently, missingness in `tmle3` is handled in a fairly simple way:

* Missing covariates are median- (for continuous) or mode- (for discrete)
  imputed, and additional covariates indicating imputation are generated, just
  as described in [the `sl3` chapter](#sl3).
* Missing treatment variables are excluded -- such observations are dropped.
* Missing outcomes are efficiently handled by the automatic calculation (and
  incorporation into estimators) of _inverse probability of censoring weights_
  (IPCW); this is also known as IPCW-TMLE and may be thought of as a joint
  intervention to remove missingness and is analogous to the procedure used with
  classical inverse probability weighted estimators.

These steps are implemented in the `process_missing` function in `tmle3`:


```r
processed <- process_missing(washb_data, node_list)
washb_data <- processed$data
node_list <- processed$node_list
```

### Create a "Spec" Object

`tmle3` is general, and allows most components of the TMLE procedure to be
specified in a modular way. However, most end-users will not be interested in
manually specifying all of these components. Therefore, `tmle3` implements a
`tmle3_Spec` object that bundles a set of components into a _specification_
("Spec") that, with minimal additional detail, can be run by an end-user.

We'll start with using one of the specs, and then work our way down into the
internals of `tmle3`.


```r
ate_spec <- tmle_ATE(
  treatment_level = "Nutrition + WSH",
  control_level = "Control"
)
```

### Define the learners

Currently, the only other thing a user must define are the `sl3` learners used
to estimate the relevant factors of the likelihood: Q and g.

This takes the form of a list of `sl3` learners, one for each likelihood factor
to be estimated with `sl3`:


```r
# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_rf <- make_learner(Lrnr_ranger)

# define metalearners appropriate to data types
ls_metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_rf),
  metalearner = ls_metalearner
)
sl_A <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_rf),
  metalearner = mn_metalearner
)
learner_list <- list(A = sl_A, Y = sl_Y)
```

Here, we use a Super Learner as defined in the previous chapter. In the future,
we plan to include reasonable defaults learners.

### Fit the TMLE

We now have everything we need to fit the tmle using `tmle3`:


```r
tmle_fit <- tmle3(ate_spec, washb_data, node_list, learner_list)
print(tmle_fit)
A tmle3_Fit that took 1 step(s)
   type                                    param   init_est  tmle_est       se
1:  ATE ATE[Y_{A=Nutrition + WSH}-Y_{A=Control}] -0.0031624 0.0077013 0.050351
       lower   upper psi_transformed lower_transformed upper_transformed
1: -0.090985 0.10639       0.0077013         -0.090985           0.10639
```

### Evaluate the Estimates

We can see the summary results by printing the fit object. Alternatively, we
can extra results from the summary by indexing into it:

```r
estimates <- tmle_fit$summary$psi_transformed
print(estimates)
[1] 0.0077013
```

## Summary

`tmle3` is a general purpose framework for generating TML estimates. The easiest
way to use it is to use a predefined spec, allowing you to just fill in the
blanks for the data, variable roles, and `sl3` learners. In the next sections,
we'll see how this framework can be used to estimate advanced parameters such as
optimal treatments and stochastic shift interventions.

There are no exercises for this brief chapter, but you may find the exercises in the corresponding handbook chapter helpful.
