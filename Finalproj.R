# Load necessary libraries
library(readxl)
library(ggplot2)
library(mvtnorm)
data <- read_excel("Final dataset.xlsx")

# Select numeric features and scale
numeric_data <- scale(data[, c("Median Aging(in Days)", "Product Price(USD)", "Average Distinct Customers Ordered(per month)", "Average Monthly Profit(USD)")])

# Convert to matrix for computation
data_matrix <- as.matrix(numeric_data)
# Safeguard against small log probabilities
safe_log <- function(x, epsilon = 1e-9) {
  ws<-pmax(x,epsilon)

  return(list(logs=log(pmax(x, epsilon)),ws=ws))  # Replace values close to 0 with epsilon
}

gmm_soft <- function(X, k, max_iter = 10000, tol = 1e-4, epsilon = 1e-6) {
  n <- nrow(X)  # Number of data points
  d <- ncol(X)  # Dimensionality of data

  set.seed(123)
  mu <- X[sample(1:n, k), ]  # Randomly initialize cluster means
  sigma <- array(0, dim = c(d, d, k))
  for (j in 1:k) sigma[,,j] <- diag(d)  # Initialize covariances as identity matrices
  pi_k <- rep(1/k, k)  # Equal cluster weights
  gamma <- matrix(0, n, k)  # Responsibility matrix
  log_likelihoods <- c()

  for (iter in 1:max_iter) {
    # E-Step: Compute responsibilities
    for (j in 1:k) {
      regularized_sigma <- sigma[,,j] + diag(epsilon, d)
      gamma[, j] <- pi_k[j] * dmvnorm(X, mean = mu[j, ], sigma = regularized_sigma)
    }

    # Safeguard against zero probabilities
    gamma[gamma < epsilon] <- epsilon
    gamma <- gamma / rowSums(gamma)  # Normalize responsibilities

    # M-Step: Update parameters
    n_k <- colSums(gamma)  # Effective cluster sizes
    for (j in 1:k) {
      mu[j, ] <- colSums(gamma[, j] * X) / n_k[j]
      centered_data <- sweep(X, 2, mu[j, ], "-")
      sigma_j <- t(centered_data) %*% (gamma[, j] * centered_data) / n_k[j]
      sigma[,,j] <- sigma_j + diag(epsilon, d)  # Regularization for stability
      pi_k[j] <- n_k[j] / n
    }
    cal<-(rowSums(sapply(1:k, function(j) {
      pi_k[j] * dmvnorm(X, mean = mu[j, ], sigma = sigma[,,j])
    })))
    res<-safe_log(cal)
    # Log-Likelihood: Safeguard against NaN values
    log_likelihood <- sum(res$logs)

    # Store the log-likelihood for convergence check
    log_likelihoods <- c(log_likelihoods, log_likelihood)

    
    # Convergence check: Compare log-likelihood from current and previous iteration
    if (iter > 1 && abs(log_likelihoods[iter] - log_likelihoods[iter - 1]) < tol) {
      cat("Converged at iteration:", iter, "\n")
      break
    }
  }
  return(list(means = mu, covariances = sigma, weights = pi_k, responsibilities = gamma,
              log_likelihood = log_likelihood, log_likelihoods = log_likelihoods))
}

# Function to calculate BIC
compute_bic <- function(log_likelihood, n, k, d) {
  num_params <- k * (d + d * (d + 1) / 2 + 1)  # Number of parameters in GMM
  return(-2 * log_likelihood + num_params * log(n))
}
compute_aic <- function(log_likelihood, n, k, d) {
  num_params <- k * (d + d * (d + 1) / 2 + 1)  # Number of parameters in GMM
  return(-2 * log_likelihood + 2 * num_params)
}
k_values <- c(4,5,6, 7, 8,9,10,11)
bic_values <- c()
aic_values <- c()
results <- list()

for (k in k_values) {
  cat("Running GMM for k =", k, "\n")
  gmm_result <- gmm_soft(data_matrix, k)
  bic <- compute_bic(gmm_result$log_likelihood, nrow(data_matrix), k, ncol(data_matrix))
  aic <- compute_aic(gmm_result$log_likelihood, nrow(data_matrix), k, ncol(data_matrix))
  bic_values <- c(bic_values, bic)
  aic_values<- c(aic_values,aic)
  results[[paste0("k_", k)]] <- gmm_result
}

# Determine the best k based on BIC
best_k <- k_values[which.min(bic_values)]

# Create a data frame for BIC and AIC values
bic_aic_df <- data.frame(
  k = k_values,
  BIC = bic_values,
  AIC = aic_values
)
print(bic_aic_df)
cat("Optimal number of clusters based on BIC:", best_k, "\n")
# Convert to long format for ggplot
library(reshape2)
bic_aic_long <- melt(bic_aic_df, id.vars = "k", variable.name = "Metric", value.name = "Value")

# Generate the BIC and AIC plot
ggplot(bic_aic_long, aes(x = k, y = Value, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "BIC and AIC vs. Number of Clusters (k)",
       x = "Number of Clusters (k)",
       y = "Value",
       color = "Metric") +
  theme_minimal()

# PCA scatter plot
library(ggrepel)
pca <- prcomp(data_matrix)
# Prepare PCA data for plotting
pca_data <- data.frame(
  PC1 = pca$x[, 1],  
  PC2 = pca$x[, 2], 
  Cluster = data$Cluster,
  Product = data$Product   
)

pca_plot_with_labels <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +  
  geom_text_repel(aes(label = Product), size = 3, max.overlaps = 10) 
  labs(
    title = paste("PCA Visualization of Clusters with Product Types for k = ",best_k),
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Cluster"
  ) +
  theme_minimal()

ggsave("pca_plot_with_labels.png", pca_plot_with_labels, width = 10, height = 7)

# Display the plot
print(pca_plot_with_labels)
