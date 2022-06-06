# Stochastic Treatment Regimes (optional)

_Nima Hejazi_

Based on the [`tmle3shift` `R` package](https://github.com/tlverse/tmle3shift).

Updated: 2022-06-06

## Learning Objectives

1. Differentiate stochastic treatment regimes from static, dynamic, and optimal
   treatment regimes.
2. Describe how estimating causal effects of stochastic interventions informs a
   real-world data analysis.
3. Contrast a population level stochastic intervention policy from a modified
   treatment policy.
4. Estimate causal effects under stochastic treatment regimes with the
   `tmle3shift` `R` package.
5. Specify a grid of counterfactual shift interventions to be used for defining
   a set of stochastic intervention policies.
6. Interpret a set of effect estimates from a grid of counterfactual shift
   interventions.
7. Construct marginal structural models to measure variable importance in terms
   of stochastic interventions, using a grid of shift interventions.
8. Implement a shift intervention at the individual level, to facilitate
   shifting each individual to a value that's supported by the data.
9. Define novel shift intervention functions to extend the `tmle3shift` `R`
   package.

## Introduction

In this section, we examine a simple example of stochastic treatment regimes in
the context of a continuous treatment variable of interest, defining an
intuitive causal effect through which to examine stochastic interventions more
generally. As a first step to using stochastic
treatment regimes in practice, we present the [`tmle3shift` R
package](https://github.com/tlverse/tmle3shift), which features an
implementation of a recently developed algorithm for computing targeted minimum
loss-based estimates of a causal effect based on a stochastic treatment regime
that shifts the natural value of the treatment based on a shifting function
$d(A,W)$. We will also use `tmle3shift` to construct marginal structural models
for variable importance measures, implement shift interventions at the
individual level, and define novel shift intervention functions.

## Stochastic Interventions

* Present a relatively simple yet extremely flexible manner by which _realistic_
  causal effects (and contrasts thereof) may be defined.
* May be applied to nearly any manner of treatment variable -- continuous,
  ordinal, categorical, binary -- allowing for a rich set of causal effects to
  be defined through this formalism.
* Arguably the most general of the classes of interventions through which causal
  effects may be defined, and are conceptually simple.

* We may consider stochastic interventions in two ways:

  1. The equation $f_A$, which produces $A$, is replaced by a probabilistic
     mechanism $g_{A,\delta}(A \mid W)$ that differs from the original $g_A(A
     \mid W)$. The _stochastically modified_ value of the treatment $A_{\delta}$
     is drawn from a user-specified distribution $g_{A,\delta}(A \mid W)$, which
     may depend on the original distribution $g_A(A \mid W)$ and is indexed by a
     user-specified parameter $\delta$. In this case, the stochastically
     modified value of the treatment $A_{\delta} \sim g_{A,\delta}(\cdot \mid
     W)$.

  2. The observed value $A$ is replaced by a new value $A_{d(A,W)}$ based on
     applying a user-defined function $d(A,W)$ to $A$. In this case, the
     stochastic treatment regime may be viewed as an intervention in which $A$
     is set equal to a value based on a hypothetical regime $d(A, W)$, where
     regime $d$ depends on the treatment level $A$ that would be assigned in the
     absence of the regime as well as the covariates $W$. Stochastic
     interventions of this variety may be referred to as depending on the
     _natural value of treatment_ or as _modified treatment policies_
     [@haneuse2013estimation; @young2014identification].

### Identifying the Causal Effect of a Stochastic Intervention

* The stochastic intervention generates a counterfactual random variable
  $Y_{d(A,W)} := f_Y(d(A,W), W, U_Y) \equiv Y_{g_{A,\delta}} := f_Y(A_{\delta},
  W, U_Y)$, where $Y_{d(A,W)} \sim P_0^{\delta}$, for $P_0^{\delta}$ being the
  counterfactual distribution implied by the intervention on $A$.

* The target causal estimand of our analysis is $\psi_{0, \delta} :=
  \mathbb{E}_{P_0^{\delta}}\{Y_{d(A,W)}\}$, the mean of the counterfactual
  outcome variable $Y_{d(A, W)}$. The statistical target parameter may also be
  denoted $\Psi(P_0) = \mathbb{E}_{P_0}{\overline{Q}_{Y,0}(d(A, W), W)}$, where
  $\overline{Q}_{Y,0}(d(A, W), W)$ is the counterfactual outcome value of a
  given individual under the stochastic intervention distribution
  [@diaz2018stochastic].

* In prior work, @diaz2012population showed that the causal quantity of interest
  $\mathbb{E}_{P_0} \{Y_{d(A, W)}\}$ is identified by a functional of the
  distribution of the observed data $O$:
  \begin{align*}\label{eqn:identification2012}
    \psi_{0,\delta} = \int_{\mathcal{W}} \int_{\mathcal{A}} & \mathbb{E}_{P_0}
     \{Y \mid A = d(a, w), W = w\} \\ &g_{A,0}(a \mid W = w) q_{W,0}(w)
     d\mu(a)d\nu(w).
  \end{align*}

* Two standard assumptions are necessary in order to establish identifiability
  of the causal parameter from the observed data via the statistical functional
  above. These were

  ::: {.definition name="Strong Ignorability"}
  $A_i \perp \!\!\! \perp Y^{d(a_i, w_i)}_i \mid W_i$, for $i = 1, \ldots, n$.
  :::

  ::: {.definition name="Treatment Positivity (or Overlap)"}
  $a_i \in \mathcal{A} \implies d(a_i, w_i) \in \mathcal{A}$ for all
  $w \in \mathcal{W}$, where $\mathcal{A}$ denotes the support of
  $A \mid W = w_i \quad \forall i = 1, \ldots n$.
  :::

* With the identification assumptions satisfied, @diaz2018stochastic provide an
  efficient influence function with respect to the nonparametric model
  $\mathcal{M}$ as
  \begin{equation*}\label{eqn:eif}
    D(P_0)(o) = H(a, w)({y - \overline{Q}_{Y,0}(a, w)}) +
    \overline{Q}_{Y,0}(d(a, w), w) - \Psi(P_0),
  \end{equation*}
  where the auxiliary covariate $H(a,w)$ may be expressed
  \begin{equation*}\label{eqn:aux_covar_full}
    H(a,w) = \mathbb{I}(a + \delta < u(w)) \frac{g_{A,0}(a - \delta \mid w)}
      {g_{A,0}(a \mid w)} + \mathbb{I}(a + \delta \geq u(w)),
  \end{equation*}
  which may be reduced to
  \begin{equation*}\label{eqn:aux_covar_simple}
    H(a,w) = \frac{g_{A,0}(a - \delta \mid w)}{g_{A,0}(a \mid w)} + 1
  \end{equation*}
  in the case that the treatment is in the limits that arise from conditioning
  on $W$, i.e., for $A_i \in (u(w) - \delta, u(w))$.



### Interpreting the Causal Effect of a Stochastic Intervention

<div class="figure" style="text-align: center">
<img src="img/gif/shift_animation.gif" alt="Animation of how a counterfactual outcome changes as the natural treatment distribution is subjected to a simple stochastic intervention" width="80%" />
<p class="caption">(\#fig:unnamed-chunk-1)Animation of how a counterfactual outcome changes as the natural treatment distribution is subjected to a simple stochastic intervention</p>
</div>



## Estimating the Causal Effect of a Stochastic Intervention with `tmle3shift`

We use `tmle3shift` to construct a targeted maximum likelihood (TML) estimator of
of a causal effect of a stochastic treatment regime that shifts the natural
value of the treatment based on a shifting function $d(A,W)$. We will follow
the recipe provided by @diaz2018stochastic, tailored to the `tmle3` framework:

1. Construct initial estimators $g_{A,n}$ of $g_{A,0}(A, W)$ and
   $\overline{Q}_{Y,n}$ of $\overline{Q}_{Y,0}(A, W)$, perhaps using
   data-adaptive regression techniques.
2. For each observation $i$, compute an estimate $H_n(a_i, w_i)$ of the
   auxiliary covariate $H(a_i,w_i)$.
3. Estimate the parameter $\epsilon$ in the logistic regression model
   $$ \text{logit}\overline{Q}_{Y, n, \epsilon}(a, w) =
   \text{logit}\overline{Q}_{Y,n}(a, w) + \epsilon H_n(a, w),$$
   or an alternative regression model incorporating weights.
4. Compute TML estimator $\Psi_n$ of the target parameter, defining update
   $\overline{Q}_{Y,n}^{\star}$ of the initial estimate
   $\overline{Q}_{Y, n, \epsilon_n}$:
   \begin{equation*}\label{eqn:tmle}
     \psi_n^{\star} = \Psi(P_n^{\star}) = \frac{1}{n} \sum_{i = 1}^n
     \overline{Q}_{Y,n}^{\star}(d(A_i, W_i), W_i).
   \end{equation*}

To start, let's load the packages we'll use and set a seed for simulation:


```r
library(tidyverse)
library(data.table)
library(sl3)
library(tmle3)
library(tmle3shift)
set.seed(429153)
```

**1. Construct initial estimators $g_{A,n}$ of $g_{A,0}(A, W)$ and
   $\overline{Q}_{Y,n}$ of $\overline{Q}_{Y,0}(A, W)$.**

We need to estimate two components of the likelihood in order to construct a
TML estimator.

1. The outcome regression, $\overline{Q}_{Y,n}$, which is a simple regression
   of the form $\mathbb{E}[Y \mid A,W]$.


```r
# learners used for conditional expectation regression
mean_learner <- Lrnr_mean$new()
fglm_learner <- Lrnr_glm_fast$new()
xgb_learner <- Lrnr_xgboost$new(nrounds = 200)
sl_regression_learner <- Lrnr_sl$new(
  learners = list(mean_learner, fglm_learner, xgb_learner)
)
```

2. The second of these is an estimate of the treatment mechanism, $g_{A,n}$,
   i.e., the _generalized propensity score_. In the case of a continuous
   intervention node $A$, such a quantity takes the form $p(A \mid W)$, which is
   a conditional density of treatment, given covariates. Generally speaking,
   conditional density estimation is a challenging problem that has received
   much attention in the literature. To estimate the treatment mechanism, we
   must make use of learning algorithms specifically suited to conditional
   density estimation; a list of such learners may be extracted from `sl3` by
   using `sl3_list_learners()`:


```r
sl3_list_learners("density")
[1] "Lrnr_density_discretize"     "Lrnr_density_hse"           
[3] "Lrnr_density_semiparametric" "Lrnr_haldensify"            
[5] "Lrnr_solnp_density"         
```

To proceed, we'll select two of the above learners, `Lrnr_haldensify` for using
the highly adaptive lasso for conditional density estimation, based on an
algorithm given by @diaz2011super and implemented in @hejazi2020haldensify, and
semiparametric location-scale conditional density estimators implemented in the
[`sl3` package](https://github.com/tlverse/sl3). A Super Learner may be
constructed by pooling estimates from each of these modified conditional
density estimation techniques.


```r
# learners used for conditional densities (i.e., generalized propensity score)
haldensify_learner <- Lrnr_haldensify$new(
  n_bins = c(3, 5),
  lambda_seq = exp(seq(-1, -10, length = 200))
)
# semiparametric density estimator based on homoscedastic errors (HOSE)
hose_learner_xgb <- make_learner(Lrnr_density_semiparametric,
  mean_learner = xgb_learner
)
# semiparametric density estimator based on heteroscedastic errors (HESE)
hese_learner_xgb_fglm <- make_learner(Lrnr_density_semiparametric,
  mean_learner = xgb_learner,
  var_learner = fglm_learner
)
# SL for the conditional treatment density
sl_density_learner <- Lrnr_sl$new(
  learners = list(
    haldensify_learner, hose_learner_xgb,
    hese_learner_xgb_fglm
  ),
  metalearner = Lrnr_solnp_density$new()
)
```

Finally, we construct a `learner_list` object for use in constructing a TML
estimator of our target parameter of interest:


```r
Q_learner <- sl_regression_learner
g_learner <- sl_density_learner
learner_list <- list(Y = Q_learner, A = g_learner)
```

### Simulate Data


```r
# simulate simple data for tmle-shift sketch
n_obs <- 1000 # number of observations
tx_mult <- 2 # multiplier for the effect of W = 1 on the treatment

## baseline covariates -- simple, binary
W <- replicate(2, rbinom(n_obs, 1, 0.5))

## create treatment based on baseline W
A <- rnorm(n_obs, mean = tx_mult * W, sd = 1)

## create outcome as a linear function of A, W + white noise
Y <- rbinom(n_obs, 1, prob = plogis(A + W))

# organize data and nodes for tmle3
data <- data.table(W, A, Y)
setnames(data, c("W1", "W2", "A", "Y"))
node_list <- list(W = c("W1", "W2"), A = "A", Y = "Y")
head(data)
   W1 W2        A Y
1:  1  1  3.58065 1
2:  1  0  3.20718 1
3:  1  1  1.03584 1
4:  0  0 -0.65785 1
5:  1  1  3.01990 1
6:  1  1  2.78031 1
```

We now have an observed data structure (`data`) and a specification of the role
that each variable in the data set plays as the nodes in a _directed acyclic
graph_ (DAG) via _nonparametric structural equation models_ (NPSEMs).

To start, we will initialize a specification for the TMLE of our parameter of
interest (a `tmle3_Spec` in the `tlverse` nomenclature) simply by calling
`tmle_shift`. We specify the argument `shift_val = 0.5` when initializing the
`tmle3_Spec` object to communicate that we're interested in a shift of $0.5$ on
the scale of the treatment $A$ -- that is, we specify $\delta = 0.5$.


```r
# initialize a tmle specification
tmle_spec <- tmle_shift(
  shift_val = 0.5,
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)
```

As seen above, the `tmle_shift` specification object (like all `tmle3_Spec`
objects) does _not_ store the data for our specific analysis of interest. Later,
we'll see that passing a data object directly to the `tmle3` wrapper function,
alongside the instantiated `tmle_spec`, will serve to construct a `tmle3_Task`
object internally (see the `tmle3` documentation for details).

<!--
Note that in the initialization of the `tmle3_Spec`, we specified a shifting
function `shift_additive_bounded` (and its inverse). This shifting function
corresponds to a stochastic regime slightly more complicated than that
initially considered in @diaz2018stochastic. In particular,
`shift_additive_bounded` is encapsulates a procedure that determines an
acceptable set of shifting values for the shift $\delta$, allowing for the
observed treatment value of a given observation to be shifted if the auxiliary
covariate $H_n$ is bounded by a constant and not shifting the given observation
if this criterion does not hold. We discuss this in greater detail in the
sequel.
-->

### Targeted Estimation of Stochastic Interventions Effects


```r
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)

Iter: 1 fn: 1342.0479	 Pars:  0.56750948 0.00000282 0.43248770
Iter: 2 fn: 1342.0479	 Pars:  0.5675095031 0.0000005883 0.4324899086
solnp--> Completed in 2 iterations
tmle_fit
A tmle3_Fit that took 1 step(s)
   type         param init_est tmle_est       se   lower   upper
1:  TSM E[Y_{A=NULL}]  0.79815  0.79755 0.012914 0.77224 0.82286
   psi_transformed lower_transformed upper_transformed
1:         0.79755           0.77224           0.82286
```

The `print` method of the resultant `tmle_fit` object conveniently displays the
results from computing our TML estimator.

## Stochastic Interventions over a Grid of Counterfactual Shifts

* Consider an arbitrary scalar $\delta$ that defines a counterfactual outcome
  $\psi_n = \mathbb{E}_P \overline{Q}_{Y,n}(d(A, W), W)$, where, for simplicity,
  let $d(A, W) = A + \delta$.  A simplified expression of the auxiliary
  covariate for the TML estimator of $\psi$ is $H_n := \frac{g^{\star}_{A,0}(a
  \mid w)}{g_{A,0}(a \mid w)}$, where $g^{\star}_{A,0}(a \mid w)$ defines the
  treatment mechanism with the stochastic intervention implemented.  In this
  manner, we can specify a _grid_ of shifts $\delta$ to define a set of
  stochastic intervention policies in an _a priori_ manner.

* To ascertain whether a given choice of the shift $\delta$ is acceptable, let
  there be a bound $C(\delta) = \frac{g^{\star}_{A,0}(a \mid w)}{g_{A,0}(a \mid
  w)} \leq M$, where $g^{\star}_{A,0}(a \mid w)$ is a function of $\delta$ in
  part, and $M$ is a user-specified upper bound of $C(\delta)$. Then,
  $C(\delta)$ is a measure of the influence of a given observation (under a
  bound of the ratio of the conditional densities), which provides a way to
  limit the maximum influence of a given observation through a choice of the
  shift $\delta$.

* For the purpose of using such a shift in practice, the present software
  provides the functions `shift_additive_bounded` and
  `shift_additive_bounded_inv`, which define a variation of this shift:
  \begin{equation}
    \delta(a, w) =
      \begin{cases}
        \delta, & C(\delta) \leq M \\
        0, \text{otherwise} \\
      \end{cases},
  \end{equation}
  which corresponds to an intervention in which the natural value of treatment
  of a given observational unit is shifted by a value $\delta$ in the case that
  the ratio of the intervened density $g^{\star}_{A,0}(a \mid w)$ to the natural
  density $g_{A,0}(a \mid w)$ (that is, $C(\delta)$) does not exceed a bound
  $M$. In the case that the ratio $C(\delta)$ exceeds the bound $M$, the
  stochastic intervention policy does not apply to the given unit and they
  remain at their natural value of treatment $a$.

### Initializing `vimshift` through its `tmle3_Spec`

To start, we will initialize a specification for the TMLE of our parameter of
interest (called a `tmle3_Spec` in the `tlverse` nomenclature) simply by calling
`tmle_shift`. We specify the argument `shift_grid = seq(-1, 1, by = 1)`
when initializing the `tmle3_Spec` object to communicate that we're interested
in assessing the mean counterfactual outcome over a grid of shifts -1, 0, 1 on the scale of the treatment $A$.


```r
# what's the grid of shifts we wish to consider?
delta_grid <- seq(from = -1, to = 1, by = 1)

# initialize a tmle specification
tmle_spec <- tmle_vimshift_delta(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)
```

### Targeted Estimation of Stochastic Intervention Effects

One may walk through the step-by-step procedure for fitting the TML estimator
of the mean counterfactual outcome under each shift in the grid, using the
machinery exposed by the [`tmle3` R package](https://tlverse.org/tmle3), or
simply invoke the `tmle3` wrapper function  to fit the series of TML estimators
(one for each parameter defined by the grid delta) in a single function call.
For convenience, we choose the latter:


```r
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)

Iter: 1 fn: 1340.1739	 Pars:  0.59196 0.12956 0.27849
Iter: 2 fn: 1340.1739	 Pars:  0.59196 0.12956 0.27848
solnp--> Completed in 2 iterations
tmle_fit
A tmle3_Fit that took 1 step(s)
         type          param init_est tmle_est        se   lower   upper
1:        TSM  E[Y_{A=NULL}]  0.61155  0.61513 0.0140867 0.58752 0.64274
2:        TSM  E[Y_{A=NULL}]  0.73900  0.73900 0.0138950 0.71177 0.76623
3:        TSM  E[Y_{A=NULL}]  0.84892  0.84841 0.0111864 0.82649 0.87034
4: MSM_linear MSM(intercept)  0.73316  0.73418 0.0125521 0.70958 0.75878
5: MSM_linear     MSM(slope)  0.11869  0.11664 0.0044021 0.10802 0.12527
   psi_transformed lower_transformed upper_transformed
1:         0.61513           0.58752           0.64274
2:         0.73900           0.71177           0.76623
3:         0.84841           0.82649           0.87034
4:         0.73418           0.70958           0.75878
5:         0.11664           0.10802           0.12527
```

_Remark_: The `print` method of the resultant `tmle_fit` object conveniently
displays the results from computing our TML estimator.

### Inference with Marginal Structural Models

Since we consider estimating the mean counterfactual outcome $\psi_n$ under
several values of the intervention $\delta$, taken from the aforementioned
$\delta$-grid, one approach for obtaining inference on a single summary measure
of these estimated quantities involves leveraging working marginal structural
models (MSMs). Summarizing the estimates $\psi_n$ through a working MSM allows
for inference on the _trend_ imposed by a $\delta$-grid to be evaluated via a
simple hypothesis test on a parameter of this working MSM. Letting
$\psi_{\delta}(P_0)$ be the mean outcome under a shift $\delta$ of the
treatment, we have $\vec{\psi}_{\delta} = (\psi_{\delta}: \delta)$ with
corresponding estimators $\vec{\psi}_{n, \delta} = (\psi_{n, \delta}: \delta)$.
Further, let $\beta(\vec{\psi}_{\delta}) = \phi((\psi_{\delta}: \delta))$. By a
straightforward application of the delta method (discussed previously), we may
write the efficient influence function of the MSM parameter $\beta$ in terms of
the EIFs of each of the corresponding point estimates. Based on this, inference
from a working MSM is rather straightforward. To wit, the limiting distribution
for $m_{\beta}(\delta)$ may be expressed $$\sqrt{n}(\beta_n - \beta_0) \to N(0,
\Sigma),$$ where $\Sigma$ is the empirical covariance matrix of
$\text{EIF}_{\beta}(O)$.


```r
tmle_fit$summary[4:5, ]
         type          param init_est tmle_est        se   lower   upper
1: MSM_linear MSM(intercept)  0.73316  0.73418 0.0125521 0.70958 0.75878
2: MSM_linear     MSM(slope)  0.11869  0.11664 0.0044021 0.10802 0.12527
   psi_transformed lower_transformed upper_transformed
1:         0.73418           0.70958           0.75878
2:         0.11664           0.10802           0.12527
```

### Directly Targeting the MSM Parameter $\beta$

Note that in the above, a working MSM is fit to the individual TML estimates of
the mean counterfactual outcome under a given value of the shift $\delta$ in
the supplied grid. The parameter of interest $\beta$ of the MSM is
asymptotically linear (and, in fact, a TML estimator) as a consequence of its
construction from individual TML estimators. In smaller samples, it may be
prudent to perform a TML estimation procedure that targets the parameter
$\beta$ directly, as opposed to constructing it from several independently
targeted TML estimates. An approach for constructing such an estimator is
proposed in the sequel.

Suppose a simple working MSM $\mathbb{E}Y_{g^0_{\delta}} = \beta_0 + \beta_1
\delta$, then a TML estimator targeting $\beta_0$ and $\beta_1$ may be
constructed as
$$\overline{Q}_{n, \epsilon}(A,W) = \overline{Q}_n(A,W) + \epsilon (H_1(g),
H_2(g),$$ for all $\delta$, where $H_1(g)$ is the auxiliary covariate for
$\beta_0$ and $H_2(g)$ is the auxiliary covariate for $\beta_1$.

To construct a targeted maximum likelihood estimator that directly targets the
parameters of the working marginal structural model, we may use the
`tmle_vimshift_msm` Spec (instead of the `tmle_vimshift_delta` Spec that
appears above):


```r
# initialize a tmle specification
tmle_msm_spec <- tmle_vimshift_msm(
  shift_grid = delta_grid,
  max_shifted_ratio = 2
)

# fit the TML estimator and examine the results
tmle_msm_fit <- tmle3(tmle_msm_spec, data, node_list, learner_list)

Iter: 1 fn: 1338.2186	 Pars:  0.62111539 0.00000916 0.37887545
Iter: 2 fn: 1338.2186	 Pars:  0.6211155447 0.0000003734 0.3788840819
solnp--> Completed in 2 iterations
tmle_msm_fit
A tmle3_Fit that took 100 step(s)
         type          param init_est tmle_est        se   lower   upper
1: MSM_linear MSM(intercept)  0.73317  0.73345 0.0125938 0.70877 0.75813
2: MSM_linear     MSM(slope)  0.11842  0.11823 0.0043158 0.10978 0.12669
   psi_transformed lower_transformed upper_transformed
1:         0.73345           0.70877           0.75813
2:         0.11823           0.10978           0.12669
```

### Example with the WASH Benefits Data

To complete our walk through, let's turn to using stochastic interventions to
investigate the data from the WASH Benefits trial. To start, let's load the
data, convert all columns to be of class `numeric`, and take a quick look at it


```r
washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data_subset.csv", stringsAsFactors = TRUE)
washb_data <- washb_data[!is.na(momage) & !is.na(momheight), ]
head(washb_data, 3)
     whz          tr fracode month aged    sex momage         momedu momheight
1: -0.94 Handwashing  N06505     7  237   male     21 Primary (1-5y)    146.00
2: -1.13     Control  N06505     8  310 female     26   No education    148.90
3: -1.61     Control  N06524     3  162   male     25 Primary (1-5y)    153.75
       hfiacat Nlt18 Ncomp watmin elec floor walls roof asset_wardrobe
1: Food Secure     1    25      2    1     0     1    1              0
2: Food Secure     1     7      4    1     0     0    1              0
3: Food Secure     0    15      2    0     0     1    1              0
   asset_table asset_chair asset_khat asset_chouki asset_tv asset_refrig
1:           1           0          0            1        0            0
2:           1           1          0            1        0            0
3:           1           0          1            1        0            0
   asset_bike asset_moto asset_sewmach asset_mobile
1:          0          0             0            0
2:          0          0             0            1
3:          0          0             0            0
```

Next, we specify our NPSEM via the `node_list` object. For our example analysis,
we'll consider the outcome to be the weight-for-height Z-score (as in previous
sections), the intervention of interest to be the mother's age at time of
child's birth, and take all other covariates to be potential confounders.


```r
node_list <- list(
  W = names(washb_data)[!(names(washb_data) %in%
    c("whz", "momage"))],
  A = "momage", Y = "whz"
)
```

Were we to consider the counterfactual weight-for-height Z-score under shifts in
the age of the mother at child's birth, how would we interpret estimates of our
parameter?

To simplify our interpretation, consider a shift (up or down) of two years in
the mother's age (i.e., $\delta = \{-2, 0, 2\}$); in this setting, a stochastic
intervention would correspond to a policy advocating that potential mothers
defer or accelerate plans of having a child for two calendar years, possibly
implemented through the deployment of an encouragement design.

First, let's try a simple upward shift of just two years:

```r
# initialize a tmle specification for just a single delta shift
washb_shift_spec <- tmle_shift(
  shift_val = 2,
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)
```

To examine the effect modification approach we looked at in previous chapters,
we'll estimate the effect of this shift $\delta = 2$ while stratifying on the
mother's education level (`momedu`, a categorical variable with three levels).
For this, we augment our initialized `tmle3_Spec` object like so


```r
# initialize effect modification specification around previous specification
washb_shift_strat_spec <- tmle_stratified(washb_shift_spec)
```

Prior to running our analysis, we'll modify the `learner_list` object we had
created to include just one of the semiparametric location-scale conditional
density estimators, as fitting of these estimators is much faster than the more
computationally intensive approach implemented in the
[`haldensify` package](ihttps://CRAN.R-project.org/package=haldensify)
[@hejazi2020haldensify].


```r
# learners used for conditional density regression (i.e., propensity score),
# but we need to turn on cross-validation for this conditional density learner
hose_learner_xgb_cv <- Lrnr_cv$new(
  learner = hose_learner_xgb,
  full_fit = TRUE
)

# modify learner list, using existing SL for Q fit
learner_list <- list(Y = Q_learner, A = hose_learner_xgb_cv)
```

Now we're ready to construct a TML estimate of the shift parameter at
$\delta = 2$, stratified across levels of our variable of interest:


```r
# fit stratified TMLE
strat_node_list <- copy(node_list)
strat_node_list$W <- setdiff(strat_node_list$W, "momedu")
strat_node_list$V <- "momedu"
washb_shift_strat_fit <- tmle3(
  washb_shift_strat_spec, washb_data, strat_node_list,
  learner_list
)
washb_shift_strat_fit
A tmle3_Fit that took 1 step(s)
             type                             param init_est tmle_est       se
1:            TSM                     E[Y_{A=NULL}] -0.57002 -0.59407 0.058813
2: stratified TSM  E[Y_{A=NULL}] | V=Primary (1-5y) -0.60764 -0.70818 0.075504
3: stratified TSM    E[Y_{A=NULL}] | V=No education -0.65186 -0.85500 0.125652
4: stratified TSM E[Y_{A=NULL}] | V=Secondary (>5y) -0.52289 -0.44745 0.094402
      lower    upper psi_transformed lower_transformed upper_transformed
1: -0.70934 -0.47880        -0.59407          -0.70934          -0.47880
2: -0.85616 -0.56020        -0.70818          -0.85616          -0.56020
3: -1.10128 -0.60873        -0.85500          -1.10128          -0.60873
4: -0.63248 -0.26243        -0.44745          -0.63248          -0.26243
```

For the next example, we'll use the variable importance strategy of considering
a grid of stochastic interventions to evaluate the weight-for-height Z-score
under a shift in the mother's age down by two years ($\delta = -2$) through up
by two years ($\delta = 2$), incrementing by a single year between the two. To
do this, we simply initialize a `Spec` `tmle_vimshift_delta` similar to how we
did in a previous example:


```r
# initialize a tmle specification for the variable importance parameter
washb_vim_spec <- tmle_vimshift_delta(
  shift_grid = seq(from = -2, to = 2, by = 1),
  max_shifted_ratio = 2
)
```

Having made the above preparations, we're now ready to estimate the
counterfactual mean of the weight-for-height Z-score under a small grid of
shifts in the mother's age at child's birth. Just as before, we do this through
a simple call to our `tmle3` wrapper function:


```r
washb_tmle_fit <- tmle3(washb_vim_spec, washb_data, node_list, learner_list)
washb_tmle_fit
A tmle3_Fit that took 1 step(s)
         type          param   init_est   tmle_est        se      lower
1:        TSM  E[Y_{A=NULL}] -0.5614487 -0.5606145 0.0459489 -0.6506728
2:        TSM  E[Y_{A=NULL}] -0.5633714 -0.5603987 0.0464541 -0.6514470
3:        TSM  E[Y_{A=NULL}] -0.5652941 -0.5652941 0.0466314 -0.6566901
4:        TSM  E[Y_{A=NULL}] -0.5672168 -0.5688442 0.0468844 -0.6607360
5:        TSM  E[Y_{A=NULL}] -0.5691396 -0.5653575 0.0479076 -0.6592546
6: MSM_linear MSM(intercept) -0.5652941 -0.5641018 0.0466522 -0.6555384
7: MSM_linear     MSM(slope) -0.0019227 -0.0017931 0.0015528 -0.0048366
        upper psi_transformed lower_transformed upper_transformed
1: -0.4705563      -0.5606145        -0.6506728        -0.4705563
2: -0.4693503      -0.5603987        -0.6514470        -0.4693503
3: -0.4738982      -0.5652941        -0.6566901        -0.4738982
4: -0.4769525      -0.5688442        -0.6607360        -0.4769525
5: -0.4714603      -0.5653575        -0.6592546        -0.4714603
6: -0.4726652      -0.5641018        -0.6555384        -0.4726652
7:  0.0012503      -0.0017931        -0.0048366         0.0012503
```

---

## Exercises

::: {.exercise}
Set the `sl3` library of algorithms for the Super Learner to a simple,
interpretable library and use this new library to estimate the counterfactual
mean of mother's age at child's birth (`momage`) under a shift $\delta = 0$.
What does this counterfactual mean equate to in terms of the observed data?
:::

::: {.exercise}
Describe two (equivalent) ways in which the causal effects of stochastic
interventions may be interpreted.
:::

::: {.exercise}
Using a grid of values of the shift parameter $\delta$ (e.g., $\{-1, 0, +1\}$),
repeat the analysis on the variable of interest (`momage`), summarizing the
trend for this sequence of shifts using a marginal structural model.
:::

::: {.exercise}
For either the grid of shifts in the example preceding the exercises or that
estimated in (3) above, plot the resultant estimates against their respective
counterfactual shifts. Graphically add to the scatterplot a line with slope and
intercept equivalent to the MSM fit through the individual TML estimates.
:::

::: {.exercise}
How does the marginal structural model we used to summarize the trend along the
sequence of shifts previously help to contextualize the estimated effect for a
single shift? That is, how does access to estimates across several shifts and
the marginal structural model parameters allow us to more richly interpret our
findings?
:::
