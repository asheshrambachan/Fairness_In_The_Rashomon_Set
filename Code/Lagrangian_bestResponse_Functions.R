# This script implements helper functions related to the Lagrangian and computation of 
# the best response functions for the Q, lambda players respectively. 
#
# Part 1 implements the data augmentation function.
#
# Part 2 implements moment functions, which construct estimates of p_0, p_1 
# from the augmented dataset. 
#
# Part 3 implements functions to construct C_lambda(Y, A, Z) and the function 
# g_lambda(Y, A, U), which are used to construct the best response for the Q player.
#
# Part 4 implements the best response functions for the Q, lambda players respectively.

library(tidyverse)

# DATA AUGMENTATION FUNCTION -----
data_augment_grid <- function( data, N ) {
  # Inputs 
  #   - data: dataframe with n rows, each row containing (Y_i, A_i, X_i)
  #   - N: integer that contains the grid size
  # Outputs
  #   - dataframe with n*N rows, each row containing Y_i, A_i, X_i, z)
  
  obs = NROW(data)
  
  return(
    data[rep(seq_len(nrow(data)), each = N), ] %>% mutate(Z = rep(1:N, obs)/N)
  )
}


# MOMENT FUNCTIONS -------
.computeMoments_statisticalParity <- function(augmented_data, protected_class) {
  # Constructs the empirical moments p0, p1 associated with the statistical parity measure of disparity.
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class)]
  #       p1 = E_n[1(A = !protected_class)]

  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class)/nrow(augmented_data),
      p1 = sum(augmented_data$A != protected_class)/nrow(augmented_data)
    )
  )
}

.computeMoments_balanceForPositiveClass <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with the balance for positive class measure of disparity.
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class, Y = 1)]
  #       p1 = E_n[1(A = !protected_class, Y = 1)]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class & augmented_data$Y == 1)/nrow(augmented_data),
      p1 = sum(augmented_data$A != protected_class & augmented_data$Y == 1)/nrow(augmented_data)
    )
  )
}

.computeMoments_balanceForPositiveClass_includeUnobs <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with the balance for positive class measure of disparity when including 
  # unobserved labels
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, g_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[g] = p1
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = mean(augmented_data$g),
      p1 = mean(augmented_data$g)
    )
  )
}

.computeMoments_balanceForNegativeClass <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with the balance for negative class measure of disparity.
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class, Y = 0)]
  #       p1 = E_n[1(A = !protected_class, Y = 0)]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class & augmented_data$Y == 0)/nrow(augmented_data),
      p1 = sum(augmented_data$A != protected_class & augmented_data$Y == 0)/nrow(augmented_data)
    )
  )
}

.computeMoments_balanceForNegativeClass_includeUnobs <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with the balance for negative class measure of disparity when including 
  # unobserved labels
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, g_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[g] = p1
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = mean(augmented_data$g),
      p1 = mean(augmented_data$g)
    )
  )
}

.computeMoments_FPI_affirmativeAction <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with FPI Affimative Action
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class)]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class)/nrow(augmented_data),
    )
  )
}

.computeMoments_FPI_affirmativeAction <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with FPI Affimative Action
  #
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class)]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class)/nrow(augmented_data)
    )
  )
}

.computeMoments_FPI_qualifiedAffirmativeAction <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with FPI Qualified Affimative Action
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[1(A = protected_class, Y = 1)]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = sum(augmented_data$A == protected_class & augmented_data$Y == 1)/nrow(augmented_data)
    )
  )
}

.computeMoments_FPI_qualifiedAffirmativeAction_includeUnobs <- function(augmented_data, protected_class) {
  # Constructs empirical moments p0, p1 associated with FPI Qualified Affimative Action when including unobserved labels
  #
  # Input
  #   - augmented_data: n*N dataset with rows (Y_i, g_i, X_i, A_i, Z)
  #   - protected_class: 0, 1 value that specifies protected class.
  # Output
  #   - list with entries
  #       p0 = E_n[g]
  
  if (protected_class != 0 & protected_class != 1) { 
    stop("Invalid input to 'protected_class'. Must equal 0 or 1")
  }
  
  return(
    list(
      p0 = mean(augmented_data$g)
    )
  )
}


computeDisparityMoments <- function(augmented_data, disparity_measure, protected_class) {
  # Provides a simple wrapper function to compute the disparity measure associated with 
  # user specified choice. 
  #
  #  Input
  #   - augmented_data; n*N dataset with rows (Y_i, X_i, A_i, Z, C)
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFP" for balance for positive class
  #       "BFNC" for balance for negative class
  #       "FPI_AA" for FPI affirmative action
  #       "FPI_QAAA" for FPI qualified affirmative action
  #       "BFPC_INCLUDE_UNOBS" for balance for positive class on the full (selectively observed and unobserved) data
  #       "BFNC_INCLUDE_UNOBS" for balance for negative class on the full (selectively observed and unobserved) data
  #       "FPI_QAAA_INCLUDE_UNOBS" for FPI qualified affirmative action on the full (selectively observed and unobserved) data
  #   - protected_class: 0, 1 value that specifies protected class.
  #
  # Output
  #   - list with with entries
  #       p0, p1 depending on the choice of disparity measure
  
  
  if (!(disparity_measure %in% c("SP", "BFPC", "BFNC", "FPI_AA", "FPI_QAA",
                                 "BFPC_INCLUDE_UNOBS", "BFNC_INCLUDE_UNOBS", "FPI_QAA_INCLUDE_UNOBS"))) { 
    stop("Invalid input to 'disparity_measure'. Must equal 'SP', 'BFPC', 'BFNC', 'FPI_AA', 'FPI_QAA', 'BFPC_INCLUDE_UNOBS', 'BFNC_INCLUDE_UNOBS', or 'FPI_QAA_INCLUDE_UNOB'")
  }
  
  if (disparity_measure == "SP") {
    result = .computeMoments_statisticalParity(augmented_data, protected_class)
  } else if ( disparity_measure == "BFPC") {
    result = .computeMoments_balanceForPositiveClass(augmented_data, protected_class)
  } else if ( disparity_measure == "BFNC") {
    result = .computeMoments_balanceForNegativeClass(augmented_data, protected_class)
  } else if ( disparity_measure == "FPI_AA") {
    result = .computeMoments_FPI_affirmativeAction(augmented_data, protected_class) 
  } else if ( disparity_measure == "FPI_QAA") {
    result = .computeMoments_FPI_qualifiedAffirmativeAction(augmented_data, protected_class)
  } else if ( disparity_measure == "FPI_QAA_INCLUDE_UNOBS") {
    result = .computeMoments_FPI_qualifiedAffirmativeAction_includeUnobs(augmented_data, protected_class)
  } else if ( disparity_measure == "BFNC_INCLUDE_UNOBS") {
    result = .computeMoments_balanceForNegativeClass_includeUnobs(augmented_data, protected_class)
  } else if ( disparity_measure == "BFPC_INCLUDE_UNOBS") {
    result = .computeMoments_balanceForPositiveClass_includeUnobs(augmented_data, protected_class)
  } 
  return(result)
}

# C_lambda, G_lambda FUNCTIONS ------
compute_CLambda_extremes <- function(augmented_data, lagrange_mult, protected_class, 
                            disparity_measure, disparity_moments) {
  # Computes the function c_lambda(Y, A, Z) for an augmented dataset for the extremes version of the fairness frontier
  #
  # Input
  #   - augmented_data; n*N dataset with rows (Y_i, X_i, A_i, Z, C)
  #     or rows rows (Y_i, g_i, X_i, A_i, Z, C) (depending on whether disparity_measure has INCLUDE_UNOBS)
  #   - lagrange_mult: numeric specifying value of lagrange multiplier
  #   - protected_class: 0, 1 value that specifies protected class.
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFPC" for balance for positive class
  #       "BFNC" for balance for negative class
  #       "FPI_AA" for FPI affirmative action
  #       "FPI_QAAA" for FPI qualified affirmative action
  #       "BFPC_INCLUDE_UNOBS" for balance for positive class on the full (selectively observed and unobserved) data
  #       "BFNC_INCLUDE_UNOBS" for balance for negative class on the full (selectively observed and unobserved) data
  #       "FPI_QAAA_INCLUDE_UNOBS" for FPI qualified affirmative action on the full (selectively observed and unobserved) data
  #   - disparity_moments: list containing values p0, p1 based on the user's choice of disparity measure
  #
  # Output
  #   - n*N dataset with rows (Y_i, X_i, A_i, Z, C, Clambda)
  
  if (!(disparity_measure %in% c("SP", "BFPC", "BFNC", "FPI_AA", "FPI_QAA",
                                 "BFPC_INCLUDE_UNOBS", "BFNC_INCLUDE_UNOBS", "FPI_QAA_INCLUDE_UNOBS"))) { 
    stop("Invalid input to 'disparity_measure'. Must equal 'SP', 'BFPC', 'BFNC', 'FPI_AA', 'FPI_QAA',
         'BFPC_INCLUDE_UNOBS', 'BFNC_INCLUDE_UNOBS', or 'FPI_QAA_INCLUDE_UNOB'")
  }

  # Compute g function for each possible choice of the disparity_measure
  if (disparity_measure == "SP") {
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          A == protected_class ~ (1/disparity_moments$p0) + lagrange_mult*C, 
          A != protected_class ~ -(1/disparity_moments$p1) + lagrange_mult*C
        ))
  } else if (disparity_measure == "BFPC") { 
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 1) ~ (1/disparity_moments$p0) + lagrange_mult*C, 
          (A != protected_class & Y == 1) ~ -(1/disparity_moments$p1) + lagrange_mult*C,
          TRUE ~ lagrange_mult*C
        )
      )
  } else if (disparity_measure == "BFNC") {
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 0) ~ (1/disparity_moments$p0) + lagrange_mult*C, 
          (A != protected_class & Y == 0) ~ -(1/disparity_moments$p1) + lagrange_mult*C,
          TRUE ~ lagrange_mult*C
        )
      )
  } else if (disparity_measure == "FPI_AA") {
    augmented_data <- augmented_data %>% 
      mutate(
        Clambda = case_when(
          (A == protected_class) ~ -(1/disparity_moments$p0) + lagrange_mult*C, 
          TRUE ~ lagrange_mult*C
        )
      )
  } else if (disparity_measure == "FPI_QAA") {
    augmented_data <- augmented_data %>% 
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 1) ~ -(1/disparity_moments$p0) + lagrange_mult*C, 
          TRUE ~ lagrange_mult*C
        )
      )
  } else if (disparity_measure == "FPI_QAA_INCLUDE_UNOBS") {
    augmented_data <- augmented_data %>% 
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 1) ~ -(g/disparity_moments$p0) + lagrange_mult*C, 
          TRUE ~ lagrange_mult*C
        )
      )
  }  else if (disparity_measure == "BFNC_INCLUDE_UNOBS") {
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 0) ~ (g/disparity_moments$p0) + lagrange_mult*C, 
          (A != protected_class & Y == 0) ~ -(g/disparity_moments$p1) + lagrange_mult*C,
          TRUE ~ lagrange_mult*C
        )
      )
  } else if (disparity_measure == "BFPC_INCLUDE_UNOBS") { 
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 1) ~ (g/disparity_moments$p0) + lagrange_mult*C, 
          (A != protected_class & Y == 1) ~ -(g/disparity_moments$p1) + lagrange_mult*C,
          TRUE ~ lagrange_mult*C
        )
      )
  }
  
  return(
    augmented_data
  )
}

compute_CLambda_minDisp <- function(augmented_data, lagrange_mult, protected_class, 
                                    disparity_measure, disparity_moments) {
  # Computes the function c_lambda(Y, A, Z) for an augmented dataset for the min disparity version of the fairness frontier
  #
  # Input
  #   - augmented_data; n*N dataset with rows (Y_i, X_i, A_i, Z, C)
  #   - lagrange_mult: numeric specifying value of lagrange multiplier
  #   - protected_class: 0, 1 value that specifies protected class.
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFP" for balance for positive class
  #       "BFNC" for balance for negative class
  #   - disparity_moments: list containing values p0, p1 based on the user's choice of disparity measure
  #
  # Output
  #   - n*N dataset with rows (Y_i, X_i, A_i, Z, C, Clambda)
  
  if (disparity_measure != "SP" & disparity_measure != "BFPC" & disparity_measure != "BFNC") { 
    stop("Invalid input to 'disparity_measure'. Must equal 'SP', 'BFP', 'BFNC', 'FPI_AA', or 'FPI_QAA'")
  }
  
  # Compute g function for each possible choice of the disparity_measure
  if (disparity_measure == "SP") {
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          A == protected_class ~ (lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p0 + lagrange_mult[3]*C, 
          A != protected_class ~ -(lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p1 + lagrange_mult[3]*C
        ))
  } else if (disparity_measure == "BFPC") { 
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 1) ~ (lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p0 + lagrange_mult[3]*C, 
          (A != protected_class & Y == 1) ~ -(lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p1 + lagrange_mult[3]*C,
          TRUE ~ lagrange_mult[3]*C
        )
      )
  } else {
    augmented_data <- augmented_data %>%
      mutate(
        Clambda = case_when(
          (A == protected_class & Y == 0) ~ (lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p0 + lagrange_mult[3]*C, 
          (A != protected_class & Y == 0) ~ -(lagrange_mult[1] - lagrange_mult[2])/disparity_moments$p1 + lagrange_mult[3]*C,
          TRUE ~ lagrange_mult[3]*C
        )
      )
  }
  
  return(
    augmented_data
  )
}

.compute_Glambda_helper <- function( obs_number, augmented_data, N ) {
  # This is a helper function that takes in the number of a unique observation in the data, 
  # an n*N dataset with rows (X, A, Y, Z, C, Clambda) and:
  #   (1) selects the rows (X, A, Y, Z, C, Clambda) associated with the unique obs (i.e., X, A, Y are replicated but Z varies over the grid)
  #   (2) computes glambda at different values of U in the grid. 
  #   (3) returns the row (X, A, Y, U) with the value of U that minimizes glambda.
  #
  # Inputs
  #   - obs_number: integer from 1 to n,
  #   - augmented_data: n*N dataset with rows (X, A, Y, Z, C, Clambda)
  #   - N: the grid size
  # 
  #  Output:
  #  - 1 row dataset with entries (X, Y, Z) where Z = U is the minimizer of glambda as a function of U
  
  # Select rows associated with unique obs
  obs_data = augmented_data[(N*obs_number - N + 1):(N*obs_number), ]
  glambda = cumsum(obs_data$Clambda)/N
  min_index = which.min(glambda)

  # Return minimizing row
  return(
    subset(obs_data[min_index, ], select = -c(A, C, Clambda))
  )
}

compute_Umin_Glambda <- function( augmented_data, N) {
  # This function takes in the n*N augmented dataset with rows (X, A, Y, Z, C, Clambda) and 
  # iterates over the each unique observation to construct the value of U that minimizes 
  # function Glambda. It return a dataset with n rows (X, A, Y, U) that can be used 
  # to construct the best response for the h player.
  #
  # Inputs:
  #  augmented_data: n*N dataset with rows (X, A, Y, Z, C, Clambda)
  #  N: grid size
  #
  # Output:
  #  n row dataset with columns (X, U) that can be used in computing the best response for Q player.
  
  obs = nrow(augmented_data)/N
  
  # Loop over each unique observation and compute (X, A, Y, U = Z) associated with that row.
    reduced_data = purrr::map_dfr(.x = 1:obs, 
                                  .f = ~.compute_Glambda_helper(obs_number = .x, augmented_data = augmented_data, N = N))
  return(
    reduced_data %>% rename(U = Z)
  )
}

