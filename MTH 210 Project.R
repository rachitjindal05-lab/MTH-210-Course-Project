# ====================================================================
# Project-MTH210: Statistical Computing
# MLE and Bootstrap for Exponentiated Exponential Distribution
# ====================================================================

# 1. Load the Dataset (fort.38)
x_data <- c(0.811, 2.750, 0.952, 0.760, 0.585, 1.477, 0.660, 1.479, 
            1.582, 2.693, 1.568, 0.990, 0.484, 1.819, 2.669, 2.349, 
            0.672, 0.655, 1.241, 0.576, 0.747, 0.852, 1.558, 1.007, 3.216)
n <- length(x_data)

hist(x_data, 
     breaks = 8, 
     freq = FALSE, 
     main = "Histogram of fort.38 Dataset",
     xlab = "Observations (x)", 
     ylab = "Density",
     col = "skyblue", 
     border = "white",
     ylim = c(0, 1.0))

# 2. Define the Newton-Raphson Function
# Returns the MLEs [alpha, lambda] and the number of iterations
newton_raphson_mle <- function(x, tol=1e-6, max_iter=100) {
  n_obs <- length(x)
  
  # Initial Guesses
  alpha <- 1.0
  lambda <- 1.0 / mean(x) 
  theta <- c(alpha, lambda)
  
  for (i in 1:max_iter) {
    a <- theta[1]
    l <- theta[2]
    
    # Pre-compute common terms to keep code clean
    term1 <- exp(-l * x)
    term2 <- 1 - term1
    
    # Prevent division by zero or log(0) in extreme edge cases
    term2[term2 < 1e-10] <- 1e-10 
    
    # Gradient Vector (U)
    dl_da <- n_obs/a + sum(log(term2))
    dl_dl <- n_obs/l - sum(x) + (a - 1) * sum((x * term1) / term2)
    U <- c(dl_da, dl_dl)
    
    # Hessian Matrix (H)
    d2l_da2 <- -n_obs / (a^2)
    d2l_dadl <- sum((x * term1) / term2)
    d2l_dl2 <- -n_obs / (l^2) - (a - 1) * sum((x^2 * term1) / (term2^2))
    
    H <- matrix(c(d2l_da2, d2l_dadl, 
                  d2l_dadl, d2l_dl2), nrow=2, byrow=TRUE)
    
    # Newton-Raphson Update Step
    # theta_new = theta_old - H^(-1) * U
    theta_new <- theta - solve(H) %*% U
    
    # Check Stopping Criterion (Euclidean distance < tolerance)
    if (sqrt(sum((theta_new - theta)^2)) < tol) {
      return(list(mle=as.vector(theta_new), iter=i, converged=TRUE))
    }
    
    theta <- theta_new
  }
  return(list(mle=as.vector(theta), iter=max_iter, converged=FALSE))
}

# 3. Calculate MLEs for Original Data
results <- newton_raphson_mle(x_data)
alpha_hat <- results$mle[1]
lambda_hat <- results$mle[2]

cat(sprintf("Iterations to converge: %d\n", results$iter))
cat(sprintf("MLE of Alpha:   %.4f\n", alpha_hat))
cat(sprintf("MLE of Lambda: %.4f\n\n", lambda_hat))


# 4. Define function to generate random variables from Exponentiated Exponential
# Using the Inverse Transform Method: F(x) = (1 - exp(-lambda*x))^alpha = u
generate_ee_data <- function(n, a, l) {
  u <- runif(n)
  x_gen <- -(1/l) * log(1 - u^(1/a))
  return(x_gen)
}

# 5. Bootstrapping Setup
B <- 1000 # Number of replications

# Arrays to store bootstrap estimates
param_boot_alpha <- numeric(B)
param_boot_lambda <- numeric(B)
nonparam_boot_alpha <- numeric(B)
nonparam_boot_lambda <- numeric(B)

# Set seed for reproducibility (so your verifiable code matches your report)
set.seed(2026) 

# 6. Run Bootstraps
for (i in 1:B) {
  # --- Parametric Bootstrap ---
  # Generate new data from the estimated distribution
  x_param <- generate_ee_data(n, alpha_hat, lambda_hat)
  res_param <- newton_raphson_mle(x_param)
  param_boot_alpha[i] <- res_param$mle[1]
  param_boot_lambda[i] <- res_param$mle[2]
  
  # --- Non-Parametric Bootstrap ---
  # Resample original data with replacement
  x_nonparam <- sample(x_data, n, replace=TRUE)
  res_nonparam <- newton_raphson_mle(x_nonparam)
  nonparam_boot_alpha[i] <- res_nonparam$mle[1]
  nonparam_boot_lambda[i] <- res_nonparam$mle[2]
}

# 7. Calculate 95% Confidence Intervals (2.5th and 97.5th percentiles)
ci_param_alpha <- quantile(param_boot_alpha, probs=c(0.025, 0.975))
ci_param_lambda <- quantile(param_boot_lambda, probs=c(0.025, 0.975))

ci_nonparam_alpha <- quantile(nonparam_boot_alpha, probs=c(0.025, 0.975))
ci_nonparam_lambda <- quantile(nonparam_boot_lambda, probs=c(0.025, 0.975))

# 8. Print Results for the Report
cat("--- 2. PARAMETRIC BOOTSTRAP 95% CI ---\n")
cat(sprintf("Alpha CI:  [%.4f,  %.4f]\n", ci_param_alpha[1], ci_param_alpha[2]))
cat(sprintf("Lambda CI: [%.4f,  %.4f]\n\n", ci_param_lambda[1], ci_param_lambda[2]))

cat("--- 3. NON-PARAMETRIC BOOTSTRAP 95% CI ---\n")
cat(sprintf("Alpha CI:  [%.4f,  %.4f]\n", ci_nonparam_alpha[1], ci_nonparam_alpha[2]))
cat(sprintf("Lambda CI: [%.4f,  %.4f]\n", ci_nonparam_lambda[1], ci_nonparam_lambda[2]))
