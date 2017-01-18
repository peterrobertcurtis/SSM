
setMethod("show", "SSM",
          function(object){
            cat("Smooth supersaturated model in d=",
                object@dimension,
                " variables with N=",
                object@basis_size,
                " terms over design of n=", object@design_size, " points.\n", sep="")
            if (object@fail) {
              cat("No SSM has been fit.\n")
            } else {
              cat("Smoothness measure Psi_2=", object@theta %*% object@K %*% object@theta, "\n", sep="")
            }
            if (length(object@variances)>0) {
              cat("Main effect Sobol' indices:\n")
              print(round(object@main_sobol, digits = 3))
              cat("Total indices:\n")
              print(round(object@total_sobol, digits = 3))
            }
            if (length(object@LOO_RMSE)>0)
              cat("Standardised LOO RMSE =", object@LOO_RMSE,"\n")
          }
)
