orthogonal_correction_old <- function(data, variable){
  data_mean <- colMeans(data)
  variable_mean <- mean(variable)
  data_centered <- sweep(data, 2, data_mean, "-")
  variable_centered <- variable - variable_mean
  scores <- as.numeric(t(data_centered) %*% variable_centered / sum(variable_centered^2))
  projection <- outer(variable_centered, scores)
  corrected_data <- data_centered - projection
  corrected_data <- sweep(corrected_data, 2, data_mean, "+")
  return(list(corrected = corrected_data, projection = projection))
}

orthogonal_correction <- function(data, variables){
  # Convert to numeric matrices
  data <- as.matrix(data)
  variables <- as.matrix(variables)

  # Center data
  data_mean <- colMeans(data)
  data_centered <- sweep(data, 2, data_mean, "-")

  # Center external variables
  V_centered <- scale(variables, center = TRUE, scale = FALSE)

  # Compute projection matrix P
  XtX <- crossprod(V_centered)              # k x k
  XtX_inv <- solve(XtX)                     # k x k
  P <- V_centered %*% XtX_inv %*% t(V_centered)  # n x n

  # Projected component (aligned with external variables)
  projection <- P %*% data_centered         # n x p

  # Corrected data (orthogonal part)
  corrected_data <- data_centered - projection
  corrected_data <- sweep(corrected_data, 2, data_mean, "+")

  return(list(corrected = corrected_data, projection = projection))
}

