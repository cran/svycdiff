## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(svycdiff)

## ----example------------------------------------------------------------------

#-- Set Seed for Random Number Generation

set.seed(1)

#-- Define Population Parameter Values

#- Population Size

N <- 10000

#- Propensity Model Parameter

tau <- 1

#- Selection Model Parameters

beta0 <- -3

beta1 <- 1 

beta2 <- 1 

#- Outcome Model Parameters

gamma0 <- 1
  
gamma1 <- 1

gamma2 <- 1

gamma3 <- 0.1

#-- Simulate Data

X <- rnorm(N, 1)

p_A <- plogis(tau * X)

A <- rbinom(N, 1, p_A)

p_S <- plogis(beta0 + beta1 * A + beta2 * X + rnorm(N, 0, 0.1))

s_wt <- 1/p_S

aa <- 1; Y1 <- gamma0 + gamma1 * X + gamma2 * aa + gamma3 * X * aa + rnorm(N)
aa <- 0; Y0 <- gamma0 + gamma1 * X + gamma2 * aa + gamma3 * X * aa + rnorm(N)

Y <- A * Y1 + (1 - A) * Y0

dat <- data.frame(Y, A, X, p_A, p_S, s_wt)

true_cdiff <- mean(Y1 - Y0)

S <- rbinom(N, 1, p_S)

samp <- dat[S == 1, ]


## ----fit----------------------------------------------------------------------

#-- Fit Model

y_mod <- Y ~ A * X

a_mod <- A ~ X

s_mod <- p_S ~ A + X

fit <- svycdiff(samp, "OM", a_mod, s_mod, y_mod, "gaussian")

fit


