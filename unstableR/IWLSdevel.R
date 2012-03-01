setMethod("IWLS", "PointProcessModel",
          function(model, control, ...) {

            value <- computeMinusLogLikelihood(model)
            reltol <- sqrt(.Machine$double.eps)
            maxit <- 100
            trace <- 1
            
            i <- 1
            while(i < maxit) {
              w <- ppstat:::computeWeights(model)
              z <- ppstat:::computeWorkingResponse(model)
              ## TODO: check when lm.fit.sparse becomes exported.
              ## This is relying on an algorithm in the MatrixModels package still
              ## under development.
              coefficients(model) <- MatrixModels:::lm.fit.sparse(ppstat:::getModelMatrix(model), z, w)
              val <- computeMinusLogLikelihood(model)
              i <- i + 1
              if(trace > 0)
                cat("Value: ", val, "NOTHELLO\n")
              if(val < value && value < reltol * (abs(value) + reltol) + val)
                break
              value <- val
            }
            
            return(model)
          }
          )
