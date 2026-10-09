#install.packages("./abcWasserstein", repos = NULL, type="source")

library(R.matlab)
library(T4transport)
library(transport)
library(ptw)
library(abc)
library(future.apply)
plan(multisession)  # or multiprocess on older R
library(rhdf5)

# Load the data
data <- readMat('narray_Area=LCT_HumanRange_mu=1.2e-07_M=5000.mat')

# Subset only the summary statistics.
sust <- lapply(data[["n"]], function(x) x[[1]])
length(sust)


# Subset only the parameters.
parameters <- as.data.frame(rbind(as.matrix(data$N),as.matrix(data$s),as.matrix(data$D),as.matrix(data$t)))
parameters<- data.frame(t(parameters))
colnames(parameters)<-c("N","s","D","t")

  # Create a storage data frame for results
  results <- data.frame(
    iteration = 1:500,
    orig_param1 = NA_real_,
    median_param1 = NA_real_,
    ci_low_param1 = NA_real_,
    ci_high_param1 = NA_real_,
    orig_param2 = NA_real_,
    median_param2 = NA_real_,
    ci_low_param2 = NA_real_,
    ci_high_param2 = NA_real_,
    orig_param3 = NA_real_,
    median_param3 = NA_real_,
    ci_low_param3 = NA_real_,
    ci_high_param3 = NA_real_,
    orig_param4 = NA_real_,
    median_param4 = NA_real_,
    ci_low_param4 = NA_real_,
    ci_high_param4 = NA_real_
  )

  # Define a function to pad smaller matrices on the reight side, so that the Wasserstein distance can be computed.
  pad_right <- function(x, target_cols) {
  x <- as.matrix(x)
  nc <- NCOL(x)
  if (nc >= target_cols) return(x)
  cbind(x, matrix(0, nrow = NROW(x), ncol = target_cols - nc))
  }

  # Define a function to calculate the Wasserstein distance.
  wasserstein_pad <- function(x, y, p = 2) {
  tgt <- max(NCOL(x), NCOL(y))
  x2 <- pad_right(x, tgt)
  y2 <- pad_right(y, tgt)
  transport::wasserstein(transport::pp(x2), transport::pp(y2), p = 2)
  }

# Loop through the datasets calculating the distances and using them to perform ABC.
i <- 1
for (p in 1:500) {

  dist <- future_sapply(sust, function(su, tgt) 
    wasserstein_pad(su, tgt, p = 2), tgt = sust[[p]])
  
    # Run abc for the i-th line of sust and parameters
    post.parest <- abc::abc(
      target = 0, 
      param = parameters[-p, ], 
      sumstat = dist[-p], 
      tol = 0.005,
      transf=c("none","none","none","none"), 
      method = "rejection"
    )

    # Extract original values from parameters
    orig_vals <- parameters[p, ]
    
    # Compute confidence intervals (e.g., 95% CI)
    ci_vals <- apply(post.parest$unadj.values, 2, quantile, probs = c(0.025, 0.975))
    
    # Store results
    results$orig_param1[i] <- orig_vals[[1]]
    results$median_param1[i] <- median(post.parest$unadj.values[,1])
    results$ci_low_param1[i] <- ci_vals[1, 1]
    results$ci_high_param1[i] <- ci_vals[2, 1]
    
    results$orig_param2[i] <- orig_vals[[2]]
    results$median_param2[i] <- median(post.parest$unadj.values[,2])
    results$ci_low_param2[i] <- ci_vals[1, 2]
    results$ci_high_param2[i] <- ci_vals[2, 2]
    
    results$orig_param3[i] <- orig_vals[[3]]
    results$median_param3[i] <- median(post.parest$unadj.values[,3])
    results$ci_low_param3[i] <- ci_vals[1, 3]
    results$ci_high_param3[i] <- ci_vals[2, 3]
    
    results$orig_param4[i] <- orig_vals[[4]]
    results$median_param4[i] <- median(post.parest$unadj.values[,4])
    results$ci_low_param4[i] <- ci_vals[1, 4]
    results$ci_high_param4[i] <- ci_vals[2, 4]
    
    i <- i+1
    
  }
# Save to output file.
write.csv(results,paste("Wasserstein_CI_results_mu=1.2e-07.csv",sep=""))

# Load the vector of empirical summary statistics.
emp <- h5read('LCT_SummaryStats_AfricaRange.mat', "ndata")

# Calculate the distances
dist <- future_sapply(sust, function(su, tgt) 
  wasserstein_pad(su, emp, p = 2), tgt = emp)
  
  # Run abc with a simple rejection approach on the Wasserstein distances.
  post.parest <- abc::abc(
    target = 0, 
    param = parameters, 
    sumstat = dist, 
    tol = 0.005,
    transf=c("none","none","none","none","none"), 
    method = "rejection"
  )
  
# Print parameter estimates.
summary(post.parest)

# Save the parameter estimatess.
write.csv(post.parest$unadj.values,paste("Wasserstein_posterior_results_mu=1.2e-07.csv",sep=""))
  

