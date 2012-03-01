setMethod("formula", "MultivariatePointProcess",
          function(x, ...) {
            return(lapply(getModels(x), formula))
          }
          )

setMethod("getModels", "MultivariatePointProcess",
          function(model, ...) {
            return(model@models)
          }
          )

setMethod("getAdjacencyMatrix", "MultivariatePointProcess",
          function(model, ...) {
            return(model@adjMat)
          }
          )

setMethod("getLocalIndependenceGraph", "MultivariatePointProcess",
          function(model, ...) {
            return(graph.adjacency(getAdjacencyMatrix(model), add.rownames = "label"))
          }
          )

setMethod("termPlot", "MultivariatePointProcess",
          function(model, alpha = 0.05, layer = geom_line(), trans = NULL, ...) {
            
            if(all(sapply(model@models, function(m) length(getFilterTerms(m)) == 0))){
              print("No filter function terms to plot")
              return(invisible())
            }

            plotData <- lapply(model@models, function(model) {
              pd <- getTermPlotData(model = model, alpha = alpha, trans = trans, ...)
              pd$response <- do.call(function(...) paste(..., sep = "+"), as.list(response(model)))
              return(pd)
            })
            plotData <- do.call(rbind, plotData)

            linearFilterPlot <- ggplot(data = plotData, aes(x = x, y = value)) +
              facet_grid(variable ~ response, scales = "free_y") +
                scale_x_continuous("position") +
                  scale_y_continuous("") + layer
            
            if(!isTRUE(all.equal(alpha, 1)))
              linearFilterPlot <- linearFilterPlot + geom_ribbon(aes(min = cf.lower, max = cf.upper), fill = alpha("blue", 0.2))

            return(linearFilterPlot)
          }
          )

setMethod("ppmFit", "MultivariatePointProcess",
          function(model, control = list(), ...) {
            models <- getModels(model)
            
            setModels(model) <- foreach(i = seq_along(models)) %dopar% {
              ppmFit(models[[i]], control = control, ...)          
            }

            ## TODO:
            ## This can run if multicore architechture is
            ## not available. Is this faster without multicore?
            ## for(i in seq_along(models)) {
            ##   models[[i]] <- ppmFit(models[[i]], control = control, ...)          
            ## }
            ## setModels(model) <- models
            return(model)
          } 
          )

setReplaceMethod("setModels", "MultivariatePointProcess",
                 function(model, value) {
                   model@models <- value

                   ## Computing the adjacency matrix
                   nodes <- unique(unlist(lapply(value, function(m) all.vars(formula(m)))))
                   adjMat <- getAdjacencyMatrix(model) 
                   if(all(nodes %in% colnames(adjMat))) {
                     adjMat[] <- 0L
                   } else {
                     n = length(nodes)
                     adjMat <- matrix(0L, ncol = n, nrow = n)
                     rownames(adjMat) <- colnames(adjMat) <- nodes
                   } 
                                        
                   for(i in seq_along(value)) {
                     to <- response(value[[i]])
                     from <- setdiff(all.vars(formula(value[[i]])), to)
                     adjMat[from, to] <- 1L
                   }
                   
                   model@adjMat <- adjMat
                   
                   return(model)
                 }
                 )

setMethod("stepInformation", "MultivariatePointProcess",
          function(model, direction = "both", trace = 1, steps = 1000, warmStart = TRUE, k = 2, ...) {
            models <- getModels(model)

             if(trace > 0 && getDoParName() != "doSEQ") {
               cat("No tracing information implemented when running in parallel.\n")
               trace <- 0
             }

            setModels(model) <- foreach(i = seq_along(models)) %dopar% {
              
              stepInformation(models[[i]],  direction = direction,
                              trace = trace, steps = steps,
                              warmStart = warmStart, k = k, ...)          
            }
            
            ## for(i in seq_along(models)) {
            ##   if(trace > 0)
            ##     cat("Response ", response(models[[i]]), "\n")
              
            ##   models[[i]] <- stepInformation(models[[i]],  direction = direction,
            ##                                  trace = trace, steps = steps,
            ##                                  warmStart = warmStart, k = k, ...)          
            ## }
            ## setModels(model) <- models
            return(model)
          }
          )

setMethod("summary", "MultivariatePointProcess",
          function(object,...) {
            return(lapply(getModels(object), summary))
          }
          )



