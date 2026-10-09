library(abc)

# Simulated parameters and summary statistics.
LCT_AfricaRange_mu12e07 <- read.csv("./SampleSummaryStats_Area=LCT_HumanRange_mu=1.2e-07_M=5000.csv")

datasets <- list(
  LCT_AfricaRange_mu12e07 = LCT_AfricaRange_mu12e07
)

# Loop through the different datasets.
counter=0
for (d in datasets) {
  
  # Subset only the parameters.
  parameters<-d[c(1:4)]
  counter=counter+1

  # Subset only the summary statistics.
  sust<-d[c(6,9:12)]

  # Identify rows with NaN values
  parameters <- parameters[is.finite(rowSums(sust)),]
  sust <- sust[is.finite(rowSums(sust)),]
  
  # Perform 100 iterations of cross-validation with different thresholds to select the value with the highest accuracy.
  cv.parest <- NULL
  while( is.null(cv.parest) ) {
    try(  
  invisible(capture.output(cv.parest <- cv4abc(parameters, sust, nval=100, tol =c(.05,.01,.005), transf=c("log","log","log","log"), prior.range = rbind(c(min(parameters$N),max(parameters$N)),c(min(parameters$s),max(parameters$s)),c(min(parameters$D),max(parameters$D)),c(min(parameters$t),max(parameters$t))), method = "neuralnet")))
    )
  }
  print(names(datasets[counter]))
  print(summary(cv.parest))
  

  # Create a storage data frame for results
  results <- data.frame(
    iteration = 1:nrow(sust),
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
  
  for (i in 1:nrow(sust)) {
    # Run abc for the i-th line of sust and parameters
    post.parest <- NULL
    while( is.null(post.parest) ) {
      try( 
        invisible(capture.output(post.parest <- abc(
      target = sust[i, ], 
      param = parameters[-i, ], 
      sumstat = sust[-i, ], 
      tol = .05,
      transf=c("log","log","log","log"), prior.range = rbind(c(min(parameters$N),max(parameters$N)),c(min(parameters$s),max(parameters$s)),c(min(parameters$D),max(parameters$D)),c(min(parameters$t),max(parameters$t))), method = "neuralnet"
    )))
      )
    }
    # Extract original values from parameters
    orig_vals <- parameters[i, ]
    
    # Compute confidence intervals (e.g., 95% CI)
    ci_vals <- apply(post.parest$adj.values, 2, quantile, probs = c(0.025, 0.975))
    
    # Store results
    results$orig_param1[i] <- orig_vals[[1]]
    results$median_param1[i] <- median(post.parest$adj.values[,1])
    results$ci_low_param1[i] <- ci_vals[1, 1]
    results$ci_high_param1[i] <- ci_vals[2, 1]
    
    results$orig_param2[i] <- orig_vals[[2]]
    results$median_param2[i] <- median(post.parest$adj.values[,2])
    results$ci_low_param2[i] <- ci_vals[1, 2]
    results$ci_high_param2[i] <- ci_vals[2, 2]
    
    results$orig_param3[i] <- orig_vals[[3]]
    results$median_param3[i] <- median(post.parest$adj.values[,3])
    results$ci_low_param3[i] <- ci_vals[1, 3]
    results$ci_high_param3[i] <- ci_vals[2, 3]
    
    results$orig_param4[i] <- orig_vals[[4]]
    results$median_param4[i] <- median(post.parest$adj.values[,4])
    results$ci_low_param4[i] <- ci_vals[1, 4]
    results$ci_high_param4[i] <- ci_vals[2, 4]
    
  }
  # Save to output file.
  write.csv(results,paste("NN_CI_results_eta_",names(datasets[counter]),".csv",sep=""))
}  

# Load the vector of empirical summary statistics.
emp<-c(3,1084,0.1479,0.0551,0.0727)

# Run abc with a neural network approach.
post.parest <- abc(
  target = emp, 
  param = parameters, 
  sumstat = sust, 
  tol = .05,
  transf=c("log","log","log","log"), 
  prior.range = rbind(c(min(parameters$N),max(parameters$N)),c(min(parameters$s),max(parameters$s)),c(min(parameters$D),max(parameters$D)),c(min(parameters$t),max(parameters$t))), 
  method = "neuralnet"
)

# Print parameter estimates.
summary(post.parest)

# Save the parameter estimatess.
write.csv(post.parest$adj.values,paste("NN_posterior_",names(datasets[counter]),".csv",sep=""))

