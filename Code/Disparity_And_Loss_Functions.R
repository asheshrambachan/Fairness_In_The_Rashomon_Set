# This script contains functions that implement different choices of the disparity function and loss functions
#
# For disparities, we implement the following measures for a binary group attribute
#   (1) Statistical Parity
#   (2) Balance for Positive Class
#   (3) Balance for Negative Class
#   (4) Fairness promoting intervention, affirmative action
#   (5) Fairness promoting intervention, qualified affirmative action
# These functions simply characterize the extremes of the fairness frontier.
#
# For losses, we implement 
#   (1) logistic regression loss
#   (2) least squares los
# Each function takes in two numeric values:
#   - y = the outcome 
#   - u = the predicted value.
# 
# These helper functions are then used to define the helper cost function, which implements
#   c(y, z) = N * (l(y, z + alpha/2) - l(y, z - alpha/2))
# where N is the grid size, alpha = 1/N.
#
# The c_function helper function is used to implement the function cost, which takes as input
# an augmented dataset of rows (Y_i, A_i, X_i, z) and adds as a column c(Y, Z).

library(tidyverse)

# DISPARITY FUNCTIONS -----
computeDisparity_riskScore <- function(augmented_data, disparity_measure, protected_class) {
  # Provides a simple wrapper function to compute the disparity measure associated with 
  # user specified choice. 
  #
  #  Input
  #   - augmented_data; n*N dataset with rows (Y_i, X_i, A_i, Z, C, Clambda, riskScore), where riskScore = f(X_i) is the prediction for that observation
  #   -                 or with rows (Y_i, g_i, X_i, A_i, Z, C, Clambda, riskScore) depending on whether disparity_measure has INCLUDE_UNOBS
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFP" for balance for positive class
  #       "BNP" for balance for negative class
  #       "FPI_AA" for FPI affirmative action
  #       "FPI_QAAA" for FPI qualified affirmative action
  #       "BFPC_INCLUDE_UNOBS" for balance for positive class on the full (selectively observed and unobserved) data
  #       "BFNC_INCLUDE_UNOBS" for balance for negative class on the full (selectively observed and unobserved) data
  #       "FPI_QAAA_INCLUDE_UNOBS" for FPI qualified affirmative action on the full (selectively observed and unobserved) data
  #   - protected_class: 0, 1 value that specifies protected class.
  #
  # Output
  #   - empirical estimate of disparity associated with risk score.
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  if (!(disparity_measure %in% c("SP", "BFPC", "BFNC", "FPI_AA", "FPI_QAA",
                                 "BFPC_INCLUDE_UNOBS", "BFNC_INCLUDE_UNOBS", "FPI_QAA_INCLUDE_UNOBS"))) { 
    stop("Invalid input to 'disparity_measure'. Must equal 'SP', 'BFPC', 'BFNC', 'FPI_AA', 'FPI_QAA',
         'BFPC_INCLUDE_UNOBS', 'BFNC_INCLUDE_UNOBS', or 'FPI_QAA_INCLUDE_UNOB'")
  }
  
  augmented_data <- augmented_data %>% mutate(
    threshold_classifier = ifelse(riskScore >= Z, 1, 0)
  )
  
  if (disparity_measure == "SP") {
    result = mean( augmented_data %>% filter(A == protected_class) %>% pull(threshold_classifier) ) - mean( augmented_data %>% filter(A != protected_class) %>% pull(threshold_classifier) )
  } else if ( disparity_measure == "BFPC") {
    result = mean( augmented_data %>% filter(A == protected_class & Y == 1) %>% pull(threshold_classifier) ) - mean( augmented_data %>% filter(A != protected_class & Y == 1) %>% pull(threshold_classifier) )
  } else if ( disparity_measure == "BFNC") {
    result = mean( augmented_data %>% filter(A == protected_class & Y == 0) %>% pull(threshold_classifier) ) -  mean( augmented_data %>% filter(A != protected_class & Y == 0) %>% pull(threshold_classifier) )
  } else if ( disparity_measure == "FPI_AA") {
    result = -mean( augmented_data %>% filter(A == protected_class) %>% pull(threshold_classifier) ) # we flip the sign on this since the researcher would like to maximize
  } else if ( disparity_measure == "FPI_QAA") {
    result = -mean( augmented_data %>% filter(A == protected_class, Y == 1) %>% pull(threshold_classifier) ) # we flip the sign on this since the researcher would like to maximize
  } else if ( disparity_measure == "FPI_QAA_INCLUDE_UNOBS") {
    result = -mean( augmented_data %>% filter(A == protected_class) %>%
                      mutate(temp = threshold_classifier*g) %>% 
                      pull (temp)) / mean( augmented_data %>% 
                                           filter(A == protected_class) %>% 
                                           pull (g)) # we flip the sign on this since the researcher would like to maximize
  } else { # same for both "BFPC_INCLUDE_UNOBS" and  "BFNC_INCLUDE_UNOBS" since the difference is in g (= mu for BFPC and = (1-mu) for BFNC)
    result = mean( augmented_data %>% filter(A == protected_class) %>%
                     mutate(temp = threshold_classifier*g) %>% 
                     pull (temp)) / mean( augmented_data %>% 
                                            filter(A == protected_class) %>% 
                                            pull (g)) - mean( augmented_data %>% filter(A != protected_class) %>%
                                                                mutate(temp = threshold_classifier*g) %>% 
                                                                pull (temp)) / mean( augmented_data %>% 
                                                                                       filter(A != protected_class) %>% 
                                                                                       pull (g))
  }
  return(result)
}

# LOSS FUNCTION HELPERS ----
.loss_logisticRegression_helper <- function(y, u, C_log = 5) {
  # implements the logistic regression loss given in Example 2 of 
  # Agarwal, Dudik and Wu (2019). The constant C_logistic is a pre-specificed
  # parameter that we always set to be equal to 5. 
  
  return( 
    log( ( 1 + exp(-C_log*(2*y-1)*(2*u-1)) ) ) / ( 2*log(1 + exp(C_log)) ) 
    )
}

.loss_leastSquares_helper <- function(y, u) {
  # implements the least squares loss.
  
  return( 
    (y-u)^2/2
  )
}

# COST FUNCTION HELPERS ----
.cost_logisticRegression_cHelper <- function(y, u, N) {
  # implements the helper cost function for logistic regression loss
  
  return(
    N * (.loss_logisticRegression_helper(y, u + 1/(2*N)) - .loss_logisticRegression_helper(y, u - 1/(2*N)))
  )
}

.cost_leastSquares_cHelper <- function(y, u, N) {
  # implements the helper cost function for the least squares loss
  
  return(
    N * (.loss_leastSquares_helper(y, u + 1/(2*N)) - .loss_leastSquares_helper(y, u - 1/(2*N)))
    )
}

# COST FUNCTIONS -----
cost_logisticRegression_cFunction <- function(augmented_data, N) {
  # Constructs the vector C(\underline{Y_i}, Z) for each observation (Y_i, A_i, X_i, Z) in the 
  # augmented dataset associated with the logistic regression loss.
  # Inputs
  #   - augmented_data: n*N dataset with rows (Y_i, A_i, X_i, Z)
  #   - N: integer that contains the grid size
  # Output
  #   n*N dataset with rows (Y_i, A_i, X_i, C) where C = C(\underline{Y}, Z).
  
  return(
    augmented_data %>% 
      mutate(C = .cost_logisticRegression_cHelper(discretize_Y(Y = Y, N = N), Z, N)
      )
  )
}

cost_leastSquares_cFunction <- function(augmented_data, N) {
  # Constructs the vector C(\underline{Y_i}, Z) for each observation (Y_i, A_i, X_i, Z) in the 
  # augmented dataset associated with the least squares loss.
  # Inputs
  #   - augmented_data: n*N dataset with rows (Y_i, A_i, X_i, Z)
  #   - N: integer that contains the grid size
  # Output
  #   n*N dataset with rows (Y_i, A_i, X_i, Z, C) where C = C(\underline{Y}, Z).
  
  return(
    augmented_data %>% 
      mutate(C = .cost_leastSquares_cHelper(discretize_Y(Y = Y, N = N), Z, N))
      )
}

cost_RiskScore <- function(augmented_data) {
  # Constructs an estimate of cost(h_f) in the augmented dataset associated with the any loss.
  # Inputs
  #   - augmented_data: n*N dataset with rows (Y_i, A_i, X_i, Z, C, riskScore), where riskScore = f(X_i) is the prediction for that observation
  # Outputs
  #   - numeric estimate of cost(h_f) = E[ C * 1{F >= Z} ].
  
  augmented_data <- augmented_data %>% mutate(
    threshold_classifier = ifelse(riskScore >= Z, 1, 0)
  )
  
  return( 
    mean( augmented_data$C * augmented_data$threshold_classifier)
  )
}
