library(tidyverse)
library(ranger)

# This script implements the exponentiated gradient algorithm for the extremes of the fairness frontier. 
run_expgrad_extremes <- function(data, 
                         disparity_measure, protected_class, learner, loss_function,
                         N, eps, n.iters, 
                         B = NULL, nu = NULL, eta = NULL,
                         Bscale = NULL, theta_init = NULL,debug = FALSE) {
  # This script implements the exponentiated gradient algorithm and 
  # returns a risk score that lies on the extremes of the fairness frontier, i.e. the disparity-measure
  # minimizing risk score.
  # 
  # Inputs
  #   - data: dataset with n rows and columns (Y, X, A) if disparity_measure is SP, BFPC, BFNC, FPIA_AA, FPIA_QAA 
  #           or columns (Y, g, X, A) if disparity_measure is BFPC_INCLUDE_UNOBS, BFNC_INCLUDE_UNOBS, or FPIA_QAA_INCLUDE_UNOBS
  #           where g is g(X,Y) per section 5 (i.e. the estimate of P(Y=1 |X))
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFPC" for balance for positive class
  #       "BFNC" for balance for negative class
  #       "BFPC_INCLUDE_UNOBS" for balance for positive class on the full (selectively observed and unobserved) data
  #       "BFNC_INCLUDE_UNOBS" for balance for negative class on the full (selectively observed and unobserved) data
  #       "FPI_AA" for FPI affirmative action
  #       "FPI_QAAA" for FPI qualified affirmative action
  #       "FPI_QAAA_INCLUDE_UNOBS" for FPI qualified affirmative action on the full (selectively observed and unobserved) data
  #   - protected_class: numeric 0, 1 value that specifies the protected class.
  #   - learner: string indicating the learner used to construct risk score. Must take values
  #       "least squares" for least squares regression
  #       "logistic" for logistic regression
  #       "random forest" for random forest
  #   - loss_function: string inicating the loss function to be used. Must take values
  #       "logistic" for logistic regression loss
  #       "least squares" for least squares loss
  #   - N: grid size
  #   - eps: allowed cost constraint violation
  #   - n.iters: maximum number of iterations
  #   - B: max setting of lagrange multiplier, default equals NULL and default setting described below
  #   - nu: convergence threshold, default equals NULL and default setting described below
  #   - eta: learning rate, default equals NULL and default setting described below
  #   - debug: logical value specifying whether the user wants function to return debugging output
  #
  # Outputs:
  #  - List containing following objects:
  #       'risk_scores'  - list of risk score models produced by each iteration of exponentiated gradient algorithm
  #       `costs'        - vector of cost(h_t) associated with each model produced over run of the algorithm
  #       `losses`       - vector of loss(f_t) associated with each model produced over run of the algorithm
  #       `disps`        - vector of disp(h_t) associated with each model produced over run of the algorithm
  #       `lagrange_mult` - vector of lagrange multiplier produced over run of the algorithm
  #       `gaps`          - vector of lagrangian dual gaps produced over run of the algorithm
  #       `param`         - list of parameters used by algorithm
  # 
  # Note:
  #   The parameters nu, eta and B are set automatically by the function following the tuning parameter choices
  #   in the implementation of the fair regression methods by Agarwal et al (2019) (https://github.com/steven7woo/fair_regression_reduction).
  #   
  #   By default, we set 
  #     - B = n^(1/2)
  #     - nu = n^(-1/2)
  #     - eta = 2 / B and then we adaptively shrink the learning rate.
  #     - eps_hat = eps - hat{c_0}, where hat{c_0} is E_n(l(Y, 1/2n)).
  
  # Type checks
  if (learner != "least squares" & learner != "logistic" & learner != "random forest") {
    stop("Invalid input to 'learner'. Must equal 'least squares', 'logistic' or 'random forest'")
  }
  if (loss_function != "logistic" & loss_function != "least squares") {
    stop("Invalid input to 'loss_function'. Must equal 'logistic' or 'least squares'")
  }
  
  # Parameters
  obs = nrow(data)
  
  if (loss_function == "logistic") {
    hatC0 = mean( data %>% mutate( C0 = .loss_logisticRegression_helper(discretize_Y(Y = Y, N = N), 
                                                                        1/(2*N)) ) %>% 
                    pull(C0) )
  } else {
    hatC0 = mean( data %>% mutate( C0 = .loss_leastSquares_helper(discretize_Y(Y = Y, N = N), 
                                                                  1/(2*N)) ) %>% 
                    pull(C0) )
  }
  epshat = eps - hatC0
  
  if (is.null(B)) { B = abs(1/epshat) }
  if (!is.null(Bscale)) { B = B/Bscale }
  if (is.null(nu)) {nu = 1/sqrt(obs) }
  if (is.null(eta)) { eta = 2 }
  
  cat(sprintf("Parameters: epshat = %.2f, B = %.2f, nu = %.2f, eta = %.2f \n", epshat, B, nu, eta))
  
  # Parameters for adaptive shrinking of eta
  SHRINK_REGRET = 0.8
  SHRINK_ETA = 0.8
  ETA_CHECK_INCR_T = 1.6
  last_eta_checked = 5
  last_gap = Inf
  
  # Preallocate vectors
  theta_vec  = rep(0, n.iters) # store vector of all theta value
  if (!is.null(theta_init)) { theta_vec[1] = theta_init }
  lambda_vec = rep(0, n.iters) # store vector of all lagrange multipliers
  nu_vec     = rep(0, n.iters) # store vector of nu_t gaps
  
  ht_cost_vec = rep(0, n.iters) # store vector of all costs of each risk score
  ht_disp_vec = rep(0, n.iters) # store vector of all disparities of each risk score
  cost_Qval   = 0               # store running sum of cost associated with Q that places equal weight on each h.
  disp_Qval   = 0               # store running sum of disparity associated with Q that places equal weight on each h.
  
  riskScore_list = vector(mode = "list", length = n.iters) # list to contain all risk score models
  
  # Create augmented dataframe with n*N rows containing (X, A, Y, Z)
  augmented_data <- data_augment_grid(data, N)
  
  # Construct C(\underline{Y_i}, Z) in the augmented dataset based on user choice of loss. Adds column C to augmented dataframe.
  if (loss_function == "logistic") {
    augmented_data <- cost_logisticRegression_cFunction(augmented_data, N)
  } else {
    augmented_data <- cost_leastSquares_cFunction(augmented_data, N)
  }
  
  # Compute moments p0, p1 associated with user choice of disparity measure
  disparityMoments <- computeDisparityMoments(augmented_data = augmented_data, 
                                              disparity_measure = disparity_measure, 
                                              protected_class = protected_class) # need to validate update for selective condition on y1

  # Begin exponentiated algorithm
  if (debug) {
    cat(sprintf("...Beginning EG Algorithm... \n"))
    cat(sprintf("Total iterations: %i \n", n.iters))
  }
  for (t in 1:n.iters) {
    #### Set lambda value ####
    lambda_vec[t] <- replace_na(B*(exp(theta_vec[t]) / (1 + exp(theta_vec[t])) ), sign(theta_vec[t]))
    
    #### Construct h best response ####
    # Add column Clambda to augmented_data
    augmented_data_alg <- compute_CLambda_extremes(augmented_data = augmented_data, lagrange_mult = lambda_vec[t],
                                      protected_class = protected_class, disparity_measure = disparity_measure,
                                      disparity_moments = disparityMoments) # need to validate update for selective condition on y1
    
    # Construct reduced dataframe with n rows and columns X, U
    reduced_data <- compute_Umin_Glambda(augmented_data = augmented_data_alg, N = N)
   
    # Construct risk score ht by performing least squares of U on X and use ht to predict over augmented data
    if (learner == "least squares") {
      ht <- lm(U ~., data = select(reduced_data, -any_of(c("Y", "g")))) 
      augmented_data_alg$riskScore <- predict(ht, newdata = augmented_data_alg)
    } else if (learner == "logistic") {
      ht <- glm(U ~., data = select(reduced_data, -any_of(c("Y", "g"))), family = "quasibinomial", control = list(maxit = 500))
      augmented_data_alg$riskScore <- predict(ht, newdata = augmented_data_alg, type = "response")
    } else {
      ht <- ranger::ranger(U~., data = select(reduced_data, -any_of(c("Y", "g"))))
      augmented_data_alg$riskScore <- predict(ht, data = augmented_data_alg, type = "response")$predictions
    }
    augmented_data_alg$riskScore[augmented_data_alg$riskScore < 0] = 0
    augmented_data_alg$riskScore[augmented_data_alg$riskScore > 1] = 1
    
    # Compute cost associated with ht and disparity associated with ht and ht to list of models
    cost_ht <- cost_RiskScore(augmented_data_alg)
    disp_ht <- computeDisparity_riskScore(augmented_data_alg, disparity_measure, protected_class) # need to validate update for selective condition on y1
    riskScore_list[[t]] <- ht
    ht_cost_vec[t] <- cost_ht
    ht_disp_vec[t] <- disp_ht
    
    # Update cost_Qval, disp_Qval
    cost_Qval <- cost_Qval + cost_ht
    disp_Qval <- disp_Qval + disp_ht
    
    #### Compute L(Qhat_t, lambdahat_t) ####
    L_t = disp_Qval/t + mean(lambda_vec[1:t])*(cost_Qval/t - epshat)
    
    #### Compute barL and barnu #####
    barL <- ifelse( cost_Qval/t <= epshat, disp_Qval/t, disp_Qval/t + B*(cost_Qval/t - epshat) )
    barnu <- barL - L_t
    
    #### Compute ubarL and ubarnu ####
    augmented_data_alg <- compute_CLambda_extremes(augmented_data = augmented_data, lagrange_mult = mean(lambda_vec[1:t]),
                                          protected_class = protected_class, disparity_measure = disparity_measure,
                                          disparity_moments = disparityMoments) 
    reduced_data <- compute_Umin_Glambda(augmented_data = augmented_data_alg, N = N)
    if (learner == "least squares") {
      ubarL_ht <- lm(U ~., data = select(reduced_data, -any_of(c("Y", "g")))) 
      augmented_data_alg$riskScore <- predict(ubarL_ht, newdata = augmented_data_alg)
    } else if (learner == "logistic") {
      ubarL_ht <- glm(U ~., data = select(reduced_data, -any_of(c("Y", "g"))), family = "quasibinomial", control = list(maxit = 500))
      augmented_data_alg$riskScore <- predict(ubarL_ht, newdata = augmented_data_alg, type = "response")
    } else {
      ubarL_ht <- ranger::ranger(U~., data = select(reduced_data, -any_of(c("Y", "g"))))
      augmented_data_alg$riskScore <- predict(ubarL_ht, data = augmented_data_alg, type = "response")$predictions
    }
    augmented_data_alg$riskScore[augmented_data_alg$riskScore < 0] = 0
    augmented_data_alg$riskScore[augmented_data_alg$riskScore > 1] = 1
    ubarL_cost <- cost_RiskScore(augmented_data_alg)
    ubarL_disp <- computeDisparity_riskScore(augmented_data_alg, disparity_measure, protected_class)
    ubarL <- ubarL_disp + mean(lambda_vec[1:t])*(ubarL_cost - epshat)
    ubarnu <- L_t - ubarL
    
    #### Set nu_t ####
    nu_vec[t] = max(barnu, ubarnu)
    
    if (debug) {
      cat(sprintf("----- Completing iteration: %i ------ \n", t))
      cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
      cat(sprintf("barL: %.2f, Lagrangian: %.2f, ubarL: %.2f \n", barL, L_t, ubarL))
      cat(sprintf("theta_t: %.2f, LM: %.2f, update: %.2f  \n", theta_vec[t], mean(lambda_vec[1:t]), eta*(cost_ht - epshat)))
    }
    
    if (nu_vec[t] <= nu) {
      break
    }
    
    #### Update eta if necessary ####
    if (t >= last_eta_checked*ETA_CHECK_INCR_T) {
      best_gap = min(nu_vec[1:t])
      if (best_gap > last_gap*SHRINK_REGRET) { 
        eta = eta*SHRINK_ETA 
      } else {
        eta = eta/SHRINK_ETA
      }
      last_eta_checked = t
      last_gap = best_gap
    }
    
    #### Update theta ####
    if (t < n.iters) {
      theta_vec[t+1] <- theta_vec[t] + eta*(cost_ht - epshat)
    }
  }
  
  # Warn if no convergence
  if (min(nu_vec) > nu) {
    cat(sprintf("WARNING: Algorithm did not converge! \n"))
    cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
  }
  
  # We flip the sign of the disparity measure if the researcher selected an FPI
  if (disparity_measure == "FPI_AA" | disparity_measure == "FPI_QAA") {
    ht_disp_vec = -ht_disp_vec
  }
  
  # Check if feasible
  if (disparity_measure == "SP" | disparity_measure == "BFPC" | disparity_measure == "BFNC") {
    feasible = (mean(ht_cost_vec[1:t]) <= epshat + (2 + 2*nu)/B)
  } else {
    feasible = (mean(ht_cost_vec[1:t]) <= epshat + (1 + 2*nu)/B)
  }
  
  return(
    list(
      risk_scores = riskScore_list[1:t],
      costs = ht_cost_vec[1:t],
      losses = ht_cost_vec[1:t] + hatC0,
      disps = ht_disp_vec[1:t],
      lagrange_mult = lambda_vec[1:t],
      gaps = nu_vec[1:t],
      feasible = feasible,
      param = list(epshat = epshat, eps = eps, B = B, nu = nu, iters = t, 
                   learner = learner, loss = loss_function, disparity = disparity_measure)
    )
  )
}

run_expgrad_minDisp <- function( data, 
                                 disparity_measure, protected_class, learner, loss_function,
                                 N, eps, n.iters, 
                                 B = NULL, nu = NULL, eta = NULL, 
                                 debug = FALSE, Bscale = NULL, thetacost_init = NULL ) {
  # This script implements the exponentiated gradient algorithm and 
  # returns a risk score that minimizes the absolute deviation from the disparity_measure. 
  #
  # This function should only be used with parity-based disparity measures such as 
  # statistical parity, balance for the positive class, and balance for the negative class.
  # 
  # Inputs
  #   - data: dataset with n rows and columns (Y, X, A)
  #   - disparity_measure: string indicating the disparity measure of interest. Must take values
  #       "SP" for statistical parity
  #       "BFPC" for balance for positive class
  #       "BFNC" for balance for negative class
  #   - protected_class: numeric 0, 1 value that specifies the protected class.
  #   - learner: string indicating the learner used to construct risk score. Must take values
  #       "least squares" for least squares regression
  #       "logistic" for logistic regression
  #       "random forest" for random forest
  #   - loss_function: string inicating the loss function to be used. Must take values
  #       "logistic" for logistic regression loss
  #       "least squares" for least squares loss
  #   - N: grid size
  #   - eps: allowed cost constraint violation
  #   - n.iters: maximum number of iterations
  #   - B: max setting of lagrange multiplier, default equals NULL and default setting described below
  #   - nu: convergence threshold, default equals NULL and default setting described below
  #   - eta: learning rate, default equals NULL and default setting described below
  #   - debug: logical value specifying whether the user wants function to return debugging output
  #
  # Outputs:
  #  - List containing following objects:
  #       'risk_scores'  - list of risk score models produced by each iteration of exponentiated gradient algorithm
  #       `costs'        - vector of cost(h_t) associated with each model produced over run of the algorithm
  #       `losses`       - vector of loss(f_t) associated with each model produced over run of the algorithm
  #       `disps`        - vector of disp(h_t) associated with each model produced over run of the algorithm
  #       `lagrange_mult` - vector of lagrange multiplier produced over run of the algorithm
  #       `gaps`          - vector of lagrangian dual gaps produced over run of the algorithm
  #       `param`         - list of parameters used by algorithm
  # 
  # Note:
  #   The parameters nu, eta and B are set automatically by the function following the tuning parameter choices
  #   in the implementation of the fair regression methods by Agarwal et al (2019) (https://github.com/steven7woo/fair_regression_reduction).
  #   
  #   By default, we set 
  #     - B = n^(1/2)
  #     - nu = n^(-1/2)
  #     - eta = 2 / B and then we adaptively shrink the learning rate.
  #     - eps_hat = eps - hat{c_0}, where hat{c_0} is E_n(l(Y, 1/2n)).
  
  # Type checks
  if (learner != "least squares" & learner != "logistic" & learner != "random forest") {
    stop("Invalid input to 'learner'. Must equal 'least squares', 'logistic' or 'random forest'")
  }
  if (loss_function != "logistic" & loss_function != "least squares") {
    stop("Invalid input to 'loss_function'. Must equal 'logistic' or 'least squares'")
  }
  if (disparity_measure != "SP" & disparity_measure != "BFPC" & disparity_measure != "BFNC") {
    stop("Invalid input to 'disparity_measure'. User should only use 'disparity_measure' to be 'SP', 'BFPC', 'BFNC' with run_expgrad_minDisp().")
  }
  
  # Parameters
  obs = nrow(data)
  if (loss_function == "logistic") {
    hatC0 = mean( data %>% mutate( C0 = .loss_logisticRegression_helper(discretize_Y(Y = Y, N = N),
                                                                        1/(2*N)) ) 
                  %>% pull(C0) )
  } else {
    hatC0 = mean( data %>% mutate( C0 = .loss_leastSquares_helper(discretize_Y(Y = Y, N = N),
                                                                  1/(2*N)) ) %>%
                    pull(C0) )
  }
  epshat = eps - hatC0
  
  if (is.null(B)) { B = abs(1/epshat) }
  if (!is.null(Bscale)) { B = Bscale }
  if (is.null(nu)) {nu = 1/sqrt(obs) }
  if (is.null(eta)) { eta = 2 }
  cat(sprintf("Parameters: epshat = %.2f, B = %.2f, nu = %.2f, eta = %.2f \n", epshat, B, nu, eta))
  
  # Parameters for adaptive shrinking of eta
  SHRINK_REGRET = 0.8
  SHRINK_ETA = 0.8
  ETA_CHECK_INCR_T = 1.6
  last_eta_checked = 5
  last_gap = Inf
  
  # Preallocate vectors
  thetaplus_vec  = rep(0, n.iters) # store vector of all theta_{+} over run of algorithm
  thetaminus_vec = rep(0, n.iters) # store vector of all theta_{-} over run of algorithm
  thetacost_vec  = rep(0, n.iters) # store vector of all theta_{cost} over run of algorithm
  if (!is.null(thetacost_init)) { thetacost_vec[1] = thetacost_init }
  
  lambdaplus_vec  = rep(0, n.iters)  # store vector of all lambda_{+} over run of algorithm
  lambdaminus_vec = rep(0, n.iters)  # store vector of all lambda_{-} over run of algorithm 
  lambdacost_vec  = rep(0, n.iters)  # store vector of all lambda_{cost} over run of algorithm
  
  xi_vec = rep(0, n.iters)      # store list of all slack variables
  nu_vec  = rep(0, n.iters)      # store vector of nu_t gaps
  
  ht_cost_vec    = rep(0, n.iters)    # store vector of all costs of each risk score
  ht_absdisp_vec = rep(0, n.iters)    # store vector of all absolute disparities of each risk score
  cost_Qval      = 0                  # store running sum of cost associated with Q that places equal weight on each h.
  disp_Qval      = 0                  # store running sum of disparity associated with Q that places equal weight on each h.
  
  riskScore_list = vector(mode = "list", length = n.iters) # list to contain all risk score models
  
  # Create augmented dataframe with n*N rows containing (X, A, Y, Z)
  augmented_data <- data_augment_grid(data, N)
  
  # Construct C(Y, Z) in the augmented dataset based on user choice of loss. Adds row C to augmented dataframe.
  if (loss_function == "logistic") {
    augmented_data <- cost_logisticRegression_cFunction(augmented_data, N)
  } else {
    augmented_data <- cost_leastSquares_cFunction(augmented_data, N)
  }
  
  # Compute moments p0, p1 associated with user choice of disparity measure
  disparityMoments <- computeDisparityMoments(augmented_data = augmented_data, 
                                              disparity_measure = disparity_measure, 
                                              protected_class = protected_class)
  

  # Begin exponentiated algorithm
  if (debug) {
    cat(sprintf("...Beginning EG Algorithm... \n"))
    cat(sprintf("Total iterations: %i \n", n.iters))
  }
  for (t in 1:n.iters) {
    #### Set lambda value ####
    lambdaplus_vec[t] = B*( exp(thetaplus_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    lambdaminus_vec[t] = B*( exp(thetaminus_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    lambdacost_vec[t] = B*( exp(thetacost_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    
    #### Construct eta best respose ####
    xi_vec[t] = case_when( (1 - lambdaplus_vec[t] - lambdaminus_vec[t] < 0) ~ 1,
                            (1 - lambdaplus_vec[t] - lambdaminus_vec[t] >= 0) ~ 0 )
    
    #### Construct h best response ####
    # Add column Clambda to augmented_data
    augmented_data_alg <- compute_CLambda_minDisp(augmented_data = augmented_data, lagrange_mult = c(lambdaplus_vec[t], lambdaminus_vec[t], lambdacost_vec[t]),
                                          protected_class = protected_class, disparity_measure = disparity_measure,
                                          disparity_moments = disparityMoments)
    
    # Construct reduced dataframe with n rows and columns X, U
    reduced_data <- compute_Umin_Glambda(augmented_data = augmented_data_alg, N = N)
    
    # Construct risk score ht by performing least squares of U on X and use ht to predict over augmented data
    if (learner == "least squares") {
      ht <- lm(U ~., data = select(reduced_data, -Y))
      augmented_data_alg$riskScore <- predict(ht, newdata = augmented_data_alg)
    } else if (learner == "logistic") {
      ht <- glm(U ~., data = select(reduced_data, -Y), family = "quasibinomial", control = list(maxit = 500))
      augmented_data_alg$riskScore <- predict(ht, newdata = augmented_data_alg, type = "response")
    } else {
      ht <- ranger::ranger(U~., data = select(reduced_data, -Y))
      augmented_data_alg$riskScore <- predict(ht, data = augmented_data_alg, type = "response")$predictions
    }
    augmented_data_alg$riskScore[augmented_data_alg$riskScore < 0] = 0
    augmented_data_alg$riskScore[augmented_data_alg$riskScore > 1] = 1
    
    # Compute cost associated with ht and disparity associated with ht and ht to list of models
    cost_ht <- cost_RiskScore(augmented_data_alg)
    disp_ht <- computeDisparity_riskScore(augmented_data_alg, disparity_measure, protected_class)
    riskScore_list[[t]] <- ht
    ht_cost_vec[t] <- cost_ht
    ht_absdisp_vec[t] <- abs(disp_ht)
    
    # Update cost_Qval, disp_Qval
    cost_Qval <- cost_Qval + cost_ht
    disp_Qval <- disp_Qval + disp_ht
    
    #### Compute L(Qhat_t, lambdahat_t) ####
    L_t = (1 - (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t])))*mean(xi_vec[1:t]) + 
      (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]))*(disp_Qval/t) + mean(lambdacost_vec[1:t])*(cost_Qval/t - epshat)
    
    #### Compute barL and barnu #####
    barL <- case_when( 
      max(disp_Qval/t - mean(xi_vec[1:t]), 
          -disp_Qval/t - mean(xi_vec[1:t]), 
          cost_Qval/t - epshat) <= 0 ~ mean(xi_vec[1:t]), 
      max(disp_Qval/t - mean(xi_vec[1:t]), 
          -disp_Qval/t - mean(xi_vec[1:t]), 
          cost_Qval/t - epshat) > 0 ~ mean(xi_vec[1:t]) + B*max(disp_Qval/t - mean(xi_vec[1:t]), -disp_Qval/t - mean(xi_vec[1:t]), cost_Qval/t - epshat)
      )
    
    #### Compute ubarL and ubarnu ####
    ubarL_xi = case_when( (1 - mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]) < 0) ~ 1,
                           (1 - mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]) >= 0) ~ 0 )
    
    augmented_data_alg <- compute_CLambda_minDisp(augmented_data = augmented_data, 
                                          lagrange_mult = c(mean(lambdaplus_vec[1:t]), mean(lambdaminus_vec[1:t]), mean(lambdacost_vec[1:t])),
                                          protected_class = protected_class, disparity_measure = disparity_measure,
                                          disparity_moments = disparityMoments)
    reduced_data <- compute_Umin_Glambda(augmented_data = augmented_data_alg, N = N)
    if (learner == "least squares") {
      ubarL_ht <- lm(U ~.-Y, data = reduced_data)
      augmented_data_alg$riskScore <- predict(ubarL_ht, newdata = augmented_data_alg)
    } else if (learner == "logistic") {
      ubarL_ht <- glm(U ~.-Y, data = reduced_data, family = "quasibinomial", control = list(maxit = 500))
      augmented_data_alg$riskScore <- predict(ubarL_ht, newdata = augmented_data_alg, type = "response")
    } else {
      ubarL_ht <- ranger::ranger(U~.-Y, data = reduced_data)
      augmented_data_alg$riskScore <- predict(ubarL_ht, data = augmented_data_alg, type = "response")$predictions
    }
    augmented_data_alg$riskScore[augmented_data_alg$riskScore < 0] = 0
    augmented_data_alg$riskScore[augmented_data_alg$riskScore > 1] = 1
    ubarL_cost <- cost_RiskScore(augmented_data_alg)
    ubarL_disp <- computeDisparity_riskScore(augmented_data_alg, disparity_measure, protected_class)
    ubarL <- (1 - (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t])))*mean(ubarL_xi) + 
      (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]))*(ubarL_disp) + mean(lambdacost_vec[1:t])*(ubarL_cost - epshat)
  
    #### Set nu_t ####
    nu_vec[t] = max(barL - L_t, L_t - ubarL)
    
    if (debug) {
      cat(sprintf("----- Completing iteration: %i ------ \n", t))
      cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
      cat(sprintf("barL: %.2f, Lagrangian: %.2f, ubarL: %.2f \n", barL, L_t, ubarL))
      cat(sprintf("thetacost_t: %.2f, LM: %.2f \n", thetacost_vec[t], mean(lambdacost_vec[1:t])))
    }
    
    if (nu_vec[t] <= nu) {
      break
    }
    
    #### Update eta if necessary ####
    if (t >= last_eta_checked*ETA_CHECK_INCR_T) {
      best_gap = min(nu_vec[1:t])
      if (best_gap > last_gap*SHRINK_REGRET) { 
        eta = eta*SHRINK_ETA 
      } else {
        eta = eta/SHRINK_ETA
      }
      last_eta_checked = t
      last_gap = best_gap
    }
    
    #### Update theta ####
    if (t < n.iters) {
      thetaplus_vec[t+1] <- thetaplus_vec[t] + eta*(disp_ht - xi_vec[t])
      thetaminus_vec[t+1] <- thetaminus_vec[t] + eta*(-disp_ht - xi_vec[t])
      thetacost_vec[t+1] <- thetacost_vec[t] + eta*(cost_ht - epshat)
    }
  }
  
  # Warn if no convergence
  if (min(nu_vec) > nu) {
    cat(sprintf("WARNING: Algorithm did not converge! \n"))
    cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
  }
  
  # Check if feasible
  feasible = (mean(ht_cost_vec[1:t]) <= epshat + (1 + 2*nu)/B)
  
  return(
    list(
      risk_scores = riskScore_list[1:t],
      costs = ht_cost_vec[1:t],
      losses = ht_cost_vec[1:t] + hatC0,
      absdisps = ht_absdisp_vec[1:t],
      lagrange_mult = cbind(lambdaplus_vec[1:t], lambdaminus_vec[1:t], lambdacost_vec[1:t]),
      gaps = nu_vec[1:t],
      feasible = feasible,
      param = list(epshat_cost = epshat, eps_cost = eps, B = B, nu = nu, iters = t, 
                   learner = learner, loss = loss_function, disparity = disparity_measure)
    )
  )
}

# BOUNDED GROUP LOSS, EXPONENTIATED GRADIENT ALGORITHM -------------------------
run_expgrad_extremes_BGL <- function(data, protected_class, learner, loss_function,
                                     N, eps, n.iters, B = NULL, nu = NULL, eta = NULL, 
                                     Bscale = NULL, theta_init = NULL,
                                     debug = FALSE) {
  # This script implements the exponentiated gradient algorithm for bounded group loss and 
  # returns a risk score that lies on the extremes of the fairness frontier.
  # minimizing risk score.
  # 
  # Inputs
  #   - data: dataset with n rows and columns (Y, X, A)
  #   - protected_class: numeric 0, 1 value that specifies the protected class.
  #   - learner: string indicating the learner used to construct risk score. Must take values
  #       "least squares" for least squares regression
  #       "logistic" for logistic regression
  #       "random forest" for random forest
  #   - loss_function: string inicating the loss function to be used. Must take values
  #       "logistic" for logistic regression loss
  #       "least squares" for least squares loss
  #   - N: grid size
  #   - eps: allowed cost constraint violation
  #   - n.iters: maximum number of iterations
  #   - B: max setting of lagrange multiplier, default equals NULL and default setting described below
  #   - nu: convergence threshold, default equals NULL and default setting described below
  #   - eta: learning rate, default equals NULL and default setting described below
  #   - debug: logical value specifying whether the user wants function to return debugging output
  #
  # Outputs:
  #  - List containing following objects:
  #       'risk_scores'  - list of risk score models produced by each iteration of exponentiated gradient algorithm
  #       `costs'        - vector of cost(h_t) associated with each model produced over run of the algorithm
  #       `losses`       - vector of loss(f_t) associated with each model produced over run of the algorithm
  #       `disps`        - vector of disp(h_t) associated with each model produced over run of the algorithm
  #       `lagrange_mult` - vector of lagrange multiplier produced over run of the algorithm
  #       `gaps`          - vector of lagrangian dual gaps produced over run of the algorithm
  #       `param`         - list of parameters used by algorithm
  # 
  # Note:
  #   The parameters nu, eta and B are set automatically by the function following the tuning parameter choices
  #   in the implementation of the fair regression methods by Agarwal et al (2019) (https://github.com/steven7woo/fair_regression_reduction).
  #   
  #   By default, we set 
  #     - B = 1/abs(epshat)
  #     - nu = n^(-1/2)
  #     - eta = 2 and then we adaptively shrink the learning rate.
  #     - eps_hat = eps - hat{c_0}, where hat{c_0} is E_n(l(Y, 1/2n)).
  
  # Type checks
  if (learner != "least squares" & learner != "logistic" & learner != "random forest") {
    stop("Invalid input to 'learner'. Must equal 'least squares', 'logistic' or 'random forest'")
  }
  if (loss_function != "logistic" & loss_function != "least squares") {
    stop("Invalid input to 'loss_function'. Must equal 'logistic' or 'least squares'")
  }
  
  # Parameters
  obs = nrow(data)
  p0 = sum(data$A == protected_class)/obs
  p1 = sum(data$A != protected_class)/obs
  epshat = eps
  if (is.null(B)) { B = abs(1/epshat) }
  if (!is.null(Bscale)) { B = B/Bscale }
  if (is.null(nu)) {nu = 1/sqrt(obs) }
  if (is.null(eta)) { eta = 2}
  cat(sprintf("Parameters: epshat = %.2f, B = %.2f, nu = %.2f, eta = %.2f \n", epshat, B, nu, eta))
  
  # Parameters for adaptive shrinking of eta
  SHRINK_REGRET = 0.8
  SHRINK_ETA = 0.8
  ETA_CHECK_INCR_T = 1.6
  last_eta_checked = 5
  last_gap = Inf
  
  # Preallocate vectors
  theta_vec  = rep(0, n.iters) # store vector of all theta value
  if (!is.null(theta_init)) { theta_vec[1] = theta_init }
  lambda_vec = rep(0, n.iters) # store vector of all lagrange multipliers
  nu_vec     = rep(0, n.iters) # store vector of nu_t gaps
  
  ht_cost_vec = rep(0, n.iters) # store vector of all costs of each risk score
  ht_disp_vec = rep(0, n.iters) # store vector of all disparities of each risk score
  cost_Qval   = 0               # store running sum of cost associated with Q that places equal weight on each h.
  disp_Qval   = 0               # store running sum of disparity associated with Q that places equal weight on each h.
  
  riskScore_list = vector(mode = "list", length = n.iters) # list to contain all risk score models
  
  # Begin exponentiated algorithm
  if (debug) {
    cat(sprintf("...Beginning EG Algorithm... \n"))
    cat(sprintf("Total iterations: %i \n", n.iters))
  }
  for (t in 1:n.iters) {
    
    #### Set lambda value ####
    lambda_vec[t] = B*(exp(theta_vec[t]) / (1 + exp(theta_vec[t])) )
    
    #### Construct h best response ####
    # Add column of weights to data
    Wlambda = ifelse(data$A == protected_class, 1/p0 + lambda_vec[t], -1/p1 + lambda_vec[t])
    Wlambda[Wlambda < 0] = 0
    
    # Construct risk score ht by performing least squares of U on X and use ht to predict over augmented data
    if (learner == "least squares") {
      ht <- lm(Y ~.-A, data = data, weights = Wlambda)
      riskScore <- predict(ht, newdata = data)
    } else if (learner == "logistic") {
      ht <- glm(Y ~.-A, data = data, family = "quasibinomial", weights = Wlambda, control = list(maxit = 500))
      riskScore <- predict(ht, newdata = data, type = "response")
    } else {
      ht <- ranger::ranger(Y ~.-A, data = data, case.weights = Wlambda)
      riskScore <- predict(ht, data = data, type = "response")$predictions
    }
    riskScore[riskScore < 0] = 0
    riskScore[riskScore > 1] = 1
    
    # Compute cost associated with ht and disparity associated with ht and ht to list of models
    if (loss_function == "least squares") {
      cost_ht <- mean(.loss_leastSquares_helper(y = data$Y, u = riskScore))
      disp_ht <- mean(.loss_leastSquares_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_leastSquares_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    } else {
      cost_ht <- mean(.loss_logisticRegression_helper(y = data$Y, u = riskScore))
      disp_ht <- mean(.loss_logisticRegression_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_logisticRegression_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    }
    riskScore_list[[t]] <- ht
    ht_cost_vec[t] <- cost_ht
    ht_disp_vec[t] <- disp_ht
    
    # Update cost_Qval, disp_Qval
    cost_Qval <- cost_Qval + cost_ht
    disp_Qval <- disp_Qval + disp_ht
    
    #### Compute L(Qhat_t, lambdahat_t) ####
    L_t = disp_Qval/t + mean(lambda_vec[1:t])*(cost_Qval/t - epshat)
    
    #### Compute barL and barnu #####
    barL <- ifelse( cost_Qval/t <= epshat, disp_Qval/t, disp_Qval/t + B*(cost_Qval/t - epshat) )
    barnu <- barL - L_t
    
    #### Compute ubarL and ubarnu ####
    Wlambda = ifelse(data$A == protected_class, 1/p0 + mean(lambda_vec[1:t]), -1/p1 + mean(lambda_vec[1:t]))
    Wlambda[Wlambda < 0] = 0
    if (learner == "least squares") {
      ht <- lm(Y ~.-A, data = data, weights = Wlambda)
      riskScore <- predict(ht, newdata = data)
    } else if (learner == "logistic") {
      ht <- glm(Y ~.-A, data = data, family = "quasibinomial", weights = Wlambda, control = list(maxit = 500))
      riskScore <- predict(ht, newdata = data, type = "response")
    } else {
      ht <- ranger::ranger(Y ~.-A, data = data, case.weights = Wlambda)
      riskScore <- predict(ht, data = data, type = "response")$predictions
    }
    riskScore[riskScore < 0] = 0
    riskScore[riskScore > 1] = 1
    
    if (loss_function == "least squares") {
      ubarL_cost <- mean(.loss_leastSquares_helper(y = data$Y, u = riskScore))
      ubarL_disp <- mean(.loss_leastSquares_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_leastSquares_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    } else {
      ubarL_cost <- mean(.loss_logisticRegression_helper(y = data$Y, u = riskScore))
      ubarL_disp <- mean(.loss_logisticRegression_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_logisticRegression_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    }
    ubarL <- ubarL_disp + mean(lambda_vec[1:t])*(ubarL_cost - epshat)
    ubarnu <- L_t - ubarL
    
    #### Set nu_t ####
    nu_vec[t] = max(barnu, ubarnu)
    
    if (debug) {
      cat(sprintf("----- Completing iteration: %i ------ \n", t))
      cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
      cat(sprintf("barL: %.2f, Lagrangian: %.2f, ubarL: %.2f \n", barL, L_t, ubarL))
      cat(sprintf("theta_t: %.2f, LM: %.2f, update: %.2f  \n", theta_vec[t], mean(lambda_vec[1:t]), eta*(cost_ht - epshat)))
    }
    
    if (nu_vec[t] <= nu) {
      break
    }
    
    #### Update eta if necessary ####
    if (t >= last_eta_checked*ETA_CHECK_INCR_T) {
      best_gap = min(nu_vec[1:t])
      if (best_gap > last_gap*SHRINK_REGRET) { 
        eta = eta*SHRINK_ETA 
      } else {
        eta = eta/SHRINK_ETA
      }
      last_eta_checked = t
      last_gap = best_gap
    }
    
    #### Update theta ####
    if (t < n.iters) {
      theta_vec[t+1] <- theta_vec[t] + eta*(cost_ht - epshat)
    }
  }
  
  # Warn if no convergence
  if (min(nu_vec) > nu) {
    cat(sprintf("WARNING: Algorithm did not converge! \n"))
    cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
  }
  
  # Check if feasible
  feasible = (mean(ht_cost_vec[1:t]) <= epshat + (2+2*nu)/B)
  
  return(
    list(
      risk_scores = riskScore_list[1:t],
      losses = ht_cost_vec[1:t],
      disps = ht_disp_vec[1:t],
      feasible = feasible,
      lagrange_mult = lambda_vec[1:t],
      gaps = nu_vec[1:t],
      param = list(epshat = epshat, eps = eps, B = B, nu = nu, iters = t, 
                   learner = learner, loss = loss_function)
    )
  )
}

run_expgrad_minDisp_BGL <- function( data, 
                                 protected_class, learner, loss_function,
                                 N, eps, n.iters, 
                                 B = NULL, nu = NULL, eta = NULL, 
                                 Bscale = NULL, thetacost_init = NULL,
                                 debug = FALSE ) {
  # This script implements the exponentiated gradient algorithm and 
  # returns a risk score that minimizes the absolute deviation of bounded group loss
  # 
  # Inputs
  #   - data: dataset with n rows and columns (Y, X, A)
  #   - protected_class: numeric 0, 1 value that specifies the protected class.
  #   - learner: string indicating the learner used to construct risk score. Must take values
  #       "least squares" for least squares regression
  #       "logistic" for logistic regression
  #       "random forest" for random forest
  #   - loss_function: string inicating the loss function to be used. Must take values
  #       "logistic" for logistic regression loss
  #       "least squares" for least squares loss
  #   - N: grid size
  #   - eps: allowed cost constraint violation
  #   - n.iters: maximum number of iterations
  #   - B: max setting of lagrange multiplier, default equals NULL and default setting described below
  #   - nu: convergence threshold, default equals NULL and default setting described below
  #   - eta: learning rate, default equals NULL and default setting described below
  #   - debug: logical value specifying whether the user wants function to return debugging output
  #
  # Outputs:
  #  - List containing following objects:
  #       'risk_scores'  - list of risk score models produced by each iteration of exponentiated gradient algorithm
  #       `costs'        - vector of cost(h_t) associated with each model produced over run of the algorithm
  #       `losses`       - vector of loss(f_t) associated with each model produced over run of the algorithm
  #       `disps`        - vector of disp(h_t) associated with each model produced over run of the algorithm
  #       `lagrange_mult` - vector of lagrange multiplier produced over run of the algorithm
  #       `gaps`          - vector of lagrangian dual gaps produced over run of the algorithm
  #       `param`         - list of parameters used by algorithm
  # 
  # Note:
  #   The parameters nu, eta and B are set automatically by the function following the tuning parameter choices
  #   in the implementation of the fair regression methods by Agarwal et al (2019) (https://github.com/steven7woo/fair_regression_reduction).
  #   
  #   By default, we set 
  #     - B = n^(1/2)
  #     - nu = n^(-1/2)
  #     - eta = 2 / B and then we adaptively shrink the learning rate.
  #     - eps_hat = eps - hat{c_0}, where hat{c_0} is E_n(l(Y, 1/2n)).
  
  # Type checks
  if (learner != "least squares" & learner != "logistic" & learner != "random forest") {
    stop("Invalid input to 'learner'. Must equal 'least squares', 'logistic' or 'random forest'")
  }
  if (loss_function != "logistic" & loss_function != "least squares") {
    stop("Invalid input to 'loss_function'. Must equal 'logistic' or 'least squares'")
  }
  
  # Parameters
  epshat = eps
  obs = nrow(data)
  p0 = sum(data$A == protected_class)/obs
  p1 = sum(data$A != protected_class)/obs
  if (is.null(B)) { B = abs(1/epshat) }
  if (!is.null(Bscale)) { B = B/Bscale }
  if (is.null(nu)) {nu = 1/sqrt(obs) }
  if (is.null(eta)) { eta = 2 }
  cat(sprintf("Parameters: epshat = %.2f, B = %.2f, nu = %.2f, eta = %.2f \n", epshat, B, nu, eta))
  
  # Parameters for adaptive shrinking of eta
  SHRINK_REGRET = 0.8
  SHRINK_ETA = 0.8
  ETA_CHECK_INCR_T = 1.6
  last_eta_checked = 5
  last_gap = Inf
  
  # Preallocate vectors
  thetaplus_vec  = rep(0, n.iters) # store vector of all theta_{+} over run of algorithm
  thetaminus_vec = rep(0, n.iters) # store vector of all theta_{-} over run of algorithm
  thetacost_vec  = rep(0, n.iters) # store vector of all theta_{cost} over run of algorithm
  if (!is.null(thetacost_init)) { thetacost_vec[1] = thetacost_init }
  
  lambdaplus_vec  = rep(0, n.iters)  # store vector of all lambda_{+} over run of algorithm
  lambdaminus_vec = rep(0, n.iters)  # store vector of all lambda_{-} over run of algorithm 
  lambdacost_vec  = rep(0, n.iters)  # store vector of all lambda_{cost} over run of algorithm
  
  xi_vec = rep(0, n.iters)      # store list of all slack variables
  nu_vec  = rep(0, n.iters)      # store vector of nu_t gaps
  
  ht_cost_vec    = rep(0, n.iters)    # store vector of all costs of each risk score
  ht_absdisp_vec = rep(0, n.iters)    # store vector of all absolute disparities of each risk score
  cost_Qval      = 0                  # store running sum of cost associated with Q that places equal weight on each h.
  disp_Qval      = 0                  # store running sum of disparity associated with Q that places equal weight on each h.
  
  riskScore_list = vector(mode = "list", length = n.iters) # list to contain all risk score models
  
  
  # Begin exponentiated algorithm
  if (debug) {
    cat(sprintf("...Beginning EG Algorithm... \n"))
    cat(sprintf("Total iterations: %i \n", n.iters))
  }
  for (t in 1:n.iters) {
    
    #### Set lambda value ####
    lambdaplus_vec[t] = B*( exp(thetaplus_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    lambdaminus_vec[t] = B*( exp(thetaminus_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    lambdacost_vec[t] = B*( exp(thetacost_vec[t])/(1 + exp(thetaplus_vec[t]) + exp(thetaminus_vec[t]) + exp(thetacost_vec[t])) )
    
    #### Construct eta best respose ####
    xi_vec[t] = case_when( (1 - lambdaplus_vec[t] - lambdaminus_vec[t] < 0) ~ 1,
                           (1 - lambdaplus_vec[t] - lambdaminus_vec[t] >= 0) ~ 0 )
    
    #### Construct h best response ####
    # Add column Clambda to augmented_data
    Wlambda = ifelse(data$A == protected_class, 
                     (lambdaplus_vec[t] - lambdamins_vec[t])/p0 + lambdacost_vec[t], 
                     -(lambdaplus_vec[t] - lambdaminus_vec[t])/p1 + lambdacost_vec[t]
                     )
    Wlambda[Wlambda < 0] = 0
    
    if (learner == "least squares") {
      ht <- lm(Y ~.-A, data = data, weights = Wlambda)
      riskScore <- predict(ht, newdata = data)
    } else if (learner == "logistic") {
      ht <- glm(Y ~.-A, data = data, family = "quasibinomial", weights = Wlambda, control = list(maxit = 500))
      riskScore <- predict(ht, newdata = data, type = "response")
    } else {
      ht <- ranger::ranger(Y~.-A, data = data, case.weights = Wlambda)
      riskScore <- predict(ht, data = data, type = "response")$predictions
    }
    riskScore[riskScore < 0] = 0
    riskScore[riskScore > 1] = 1
    
    # Compute cost associated with ht and disparity associated with ht and ht to list of models
    if (loss_function == "least squares") {
      cost_ht <- mean(.loss_leastSquares_helper(y = data$Y, u = riskScore))
      disp_ht <- mean(.loss_leastSquares_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_leastSquares_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    } else {
      cost_ht <- mean(.loss_logisticRegression_helper(y = data$Y, u = riskScore))
      disp_ht <- mean(.loss_logisticRegression_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_logisticRegression_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    }
    riskScore_list[[t]] <- ht
    ht_cost_vec[t] <- cost_ht
    ht_absdisp_vec[t] <- abs(disp_ht)
    
    # Update cost_Qval, disp_Qval
    cost_Qval <- cost_Qval + cost_ht
    disp_Qval <- disp_Qval + disp_ht
    
    #### Compute L(Qhat_t, lambdahat_t) ####
    L_t = (1 - (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t])))*mean(xi_vec[1:t]) + 
      (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]))*(disp_Qval/t) + mean(lambdacost_vec[1:t])*(cost_Qval/t - epshat)
    
    #### Compute barL and barnu #####
    barL <- case_when( 
      max(disp_Qval/t - mean(xi_vec[1:t]), 
          -disp_Qval/t - mean(xi_vec[1:t]), 
          cost_Qval/t - epshat) <= 0 ~ mean(xi_vec[1:t]), 
      max(disp_Qval/t - mean(xi_vec[1:t]), 
          -disp_Qval/t - mean(xi_vec[1:t]), 
          cost_Qval/t - epshat) > 0 ~ mean(xi_vec[1:t]) + B*max(disp_Qval/t - mean(xi_vec[1:t]), -disp_Qval/t - mean(xi_vec[1:t]), cost_Qval/t - epshat)
    )
    
    #### Compute ubarL and ubarnu ####
    ubarL_xi = case_when( (1 - mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]) < 0) ~ 1,
                          (1 - mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]) >= 0) ~ 0 )
    
    Wlambda = ifelse(data$A == protected_class, 
                     (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]))/p0 + mean(lambdacost_vec[1:t]), 
                     -(mean(lambdaplus_vec[1:t])/p1 - mean(lambdaminus_vec[1:t]))/p1 + mean(lambdacost_vec[1:t])
                     )
    Wlambda[Wlambda < 0] = 0
    
    # Construct risk score ht by performing least squares of U on X and use ht to predict over augmented data
    if (learner == "least squares") {
      ht <- lm(Y ~.-A, data = data, weights = Wlambda)
      riskScore <- predict(ht, newdata = data)
    } else if (learner == "logistic") {
      ht <- glm(Y ~.-A, data = data, family = "quasibinomial", weights = Wlambda, control = list(maxit = 500))
      riskScore <- predict(ht, newdata = data, type = "response")
    } else {
      ht <- ranger::ranger(Y~.-A, data = data, case.weights = Wlambda)
      riskScore <- predict(ht, data = data, type = "response")$predictions
    }
    riskScore[riskScore < 0] = 0
    riskScore[riskScore > 1] = 1
    if (loss_function == "least squares") {
      ubarL_cost <- mean(.loss_leastSquares_helper(y = data$Y, u = riskScore))
      ubarL_disp <- mean(.loss_leastSquares_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_leastSquares_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    } else {
      ubarL_cost <- mean(.loss_logisticRegression_helper(y = data$Y, u = riskScore))
      ubarL_disp <- mean(.loss_logisticRegression_helper(y = data$Y[data$A == protected_class], u = riskScore[data$A == protected_class])) -
        mean(.loss_logisticRegression_helper(y = data$Y[data$A != protected_class], u = riskScore[data$A != protected_class]))
    }
    ubarL <- (1 - (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t])))*mean(ubarL_xi) + 
      (mean(lambdaplus_vec[1:t]) - mean(lambdaminus_vec[1:t]))*(ubarL_disp) + mean(lambdacost_vec[1:t])*(ubarL_cost - epshat)
    
    #### Set nu_t ####
    nu_vec[t] = max(barL - L_t, L_t - ubarL)
    
    if (debug) {
      cat(sprintf("----- Completing iteration: %i ------ \n", t))
      cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
      cat(sprintf("barL: %.2f, Lagrangian: %.2f, ubarL: %.2f \n", barL, L_t, ubarL))
      cat(sprintf("theta_t: %.2f, LM: %.2f \n", theta_vec[t], mean(lambda_vec[1:t])))
    }
    
    if (nu_vec[t] <= nu) {
      break
    }
    
    #### Update eta if necessary ####
    if (t >= last_eta_checked*ETA_CHECK_INCR_T) {
      best_gap = min(nu_vec[1:t])
      if (best_gap > last_gap*SHRINK_REGRET) { 
        eta = eta*SHRINK_ETA 
      } else {
        eta = eta/SHRINK_ETA
      }
      last_eta_checked = t
      last_gap = best_gap
    }
    
    #### Update theta ####
    if (t < n.iters) {
      thetaplus_vec[t+1] <- thetaplus_vec[t] + eta*(disp_ht - xi_vec[t])
      thetaminus_vec[t+1] <- thetaminus_vec[t] + eta*(-disp_ht - xi_vec[t])
      thetacost_vec[t+1] <- thetacost_vec[t] + eta*(cost_ht - epshat)
    }
  }
  
  # Warn if no convergence
  if (min(nu_vec) > nu) {
    cat(sprintf("WARNING: Algorithm did not converge! \n"))
    cat(sprintf("nu_t: %.2f, nu: %.2f \n", nu_vec[t], nu))
  }
  
  # Check if feasible
  feasible = (mean(ht_cost_vec[1:t]) <= epshat + (1 + 2*nu)/B)
  
  return(
    list(
      risk_scores = riskScore_list[1:t],
      losses = ht_cost_vec[1:t],
      absdisps = ht_absdisp_vec[1:t],
      feasible = feasible,
      lagrange_mult = cbind(lambdaplus_vec[1:t], lambdaminus_vec[1:t], lambdacost_vec[1:t]),
      gaps = nu_vec[1:t],
      param = list(epshat = epshat, eps = eps, B = B, nu = nu, iters = t, 
                   learner = learner, loss = loss_function)
    )
  )
}

# REDUCE SOLUTION SIZE --------------------------------------------------
reduceSolutions_extremes <- function(costs, disps, epshat) {
  # Returns indices of which models are kept in the reduction based on Cotter et al (2019).
  # We choose the bound on the linear program to be the minimal nu > 0 such that the linear program 
  # has a feasible solution. 
  #
  # We set the nu in the linear program equal to:
  #   if min(costs - epshat) < 0, set rhs = epshat
  #   if min(costs - epshat) > 0, set nu = epshat + min(costs - epshat) 

  rhs = ifelse(min(costs - epshat) < 0, epshat, epshat + min(costs - epshat))
  lpSolution = lpSolve::lp(direction = "min", objective.in = disps, 
                           const.mat = rbind(costs, rep(1, length(disps))), const.dir = c("<=", "=="),
                           const.rhs = c(rhs, 1))
  return(
    lpSolution$solution
  )
}

reduceSolutions_minDisp <- function(costs, disps, epshat) {
  # Returns indices of which models are kept in the reduction based on Cotter et al (2019).
  # We choose the bound on the linear program to be the minimal nu > 0 such that the linear program 
  # has a feasible solution. 
  #
  # We set the nu in the linear program equal to:
  #   if min(costs - epshat) < 0, set rhs = epshat
  #   if min(costs - epshat) > 0, set nu = epshat + min(costs - epshat) 
  
  rhs = ifelse(min(costs - epshat) < 0, epshat, epshat + min(costs - epshat))
  lpSolution = lpSolve::lp(direction = "min", objective.in = disps, 
                           const.mat = rbind(costs, rep(1, length(disps))), const.dir = c("<=", "=="),
                           const.rhs = c(rhs, 1))
  return(
    lpSolution$solution
  )
}

