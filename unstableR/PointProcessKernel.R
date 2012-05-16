pointProcessKernel <- function(
                        formula,
                        data,
                        family,
                        support = 1,
                        N = 200,
                        Delta,
                        Omega,
                        coefficients,
                        fixedCoefficients = list(),
                        modelMatrix = TRUE,
                        fit = modelMatrix,
                        varMethod = 'Fisher',
                        basisEnv,
                        kernel = sobolevKernel
                        ...) {
  
  call <- match.call()
  argList <- as.list(call)[-1]
  argList$fit <- FALSE
  argList$modelMatrix <- FALSE
  model <- do.call("pointProcessModel", argList)
  if (class(model) == "MultivariatePointProcess")
    stop("Multivariate models not currently supported with 'pointProcessKernel'.")
  
  if (modelMatrix) {
      model <- computeModelMatrix(as(model, "PointProcessKernel"))
    }
  else {
    model <- updateKernelMatrix(as(model, "PointProcessKernel"),
                                Matrix(),
                                assign = numeric(),
                                form = formula(~0)
                                )
  }

  grid <- model@basisPoints
  support <- model@support
    
  G1 <- outer(grid, grid, kernel, t = support[2], sub = 1)
  G <- outer(grid, grid, kernel, t = support[2])
  G1svd <- svd(G1, nv = 0)
  ### Select computationally non-zero singular values.
  model@d <- G1svd$d/G1svd$d[1] > 5*.Machine$double.eps   
  G1svd$d <- ifelse(model@d, G1svd$d, 1)
  model@U <- G %*% t(t(G1svd$u)/sqrt(G1svd$d))

  if(fit) {
    model <- ppmFit(model, selfStart = selfStart, ...)
  } else {
    ## Initializing the variance matrix without computing it.
    model <- computeVar(model, method = "none")  
    ##  TODO: correct to work with multivariate models
    ##    model@optimResult <- list(value = computeMinusLogLikelihood(model),
    ##                              counts = c(0, 0),
    ##                              convergence = NA)
  }

  return(model)
}


setMethod("computeModelMatrix", "PointProcessKernel",
          function(model, evaluationPositions = NULL, ...){

            ## The 'model' of class PointProcessKernel contains the data
            ## as an object of class MarkedPointProcess and the formula for the
            ## model specification. The 'evalPositions' below corresponding to
            ## the model matrix rows are either given by the 'evaluationPositions'
            ## argument or extracted from the the MarkedPointProcess object (default).

            if(is.null(evaluationPositions)) {
              evalPositions <- tapply(getPosition(processData(model)),
                                      getId(processData(model)), list)
            } else {
              evalPositions <- evaluationPositions
            }

            ## Checks if the model is allowed to be anticipating and sets
            ## the 'zero' accordingly.

            if(ppstat:::anticipating(model)) {
              zero <- which(model@basisPoints == 0) - 1
            } else {
              zero <- 0
            }
              
            ## The observed points ('positions') for the marked point process,
            ## the corresponding 'id' labels and 'marks' are extracted.
            
            processData <- processData(model)
            positions <- getPointPosition(processData)
            r <- length(model@basisPoints)
            
            id <- factor(getPointId(processData))
            idLevels <- levels(id)
            
            marks <- getMarkType(processData, drop = FALSE)
            markLevels <- levels(marks)

            ## The special terms in the formula that encodes the
            ## linear filters that are modeled non-parametrically are
            ## identified and the formula for the remaining model
            ## specification is constructed.

            formula <- formula(model)
            terms <- terms(formula, "k")
            kernels <- attr(delete.response(terms), "specials")$k
            formulaNoKernels <- formula(terms[-kernels])
            kernelVar <- 1 + attr(terms, "response") + kernels
            kernelVar <- sapply(as.list(attr(terms, "variables"))[kernelVar],
                                 all.vars)
            terms <- delete.response(terms)
            termLabels <- attr(terms, "term.labels")

            if(!all(kernelVar %in% markLevels))
              stop("The use of kernel filters is only implemented for point process variables.")

            ## The points where the basis functions are evaluated are extracted
            ## and the list of model matrices ('design') is set up, which holds model
            ## matrices for the different terms. 'assign' will be an attribute to  
            ## the model matrix of length equal to the number of columns, and for
            ## each column pointing to the term number. 
                        
            ## Model matrix computations for the terms involving filters:

            design <- ppstat:::lapplyParallel(kernels,
                                     function(i, ...) {
                                       term <- termLabels[i]
                                       variable <- all.vars(terms[i])
                                       
                                       ## The occurrence matrix is computed by a loop over 
                                       ## each value of 'id' whose result is stored in
                                       ## 'designList'.
                                       
                                       ## TODO: C level computation?
                                       
                                       designList <- list()
                                       
                                       ## Central loop over 'idLevels' and computations of
                                       ## the occurrence matrix as a sparse matrix, bound together
                                       ## in one matrix below
                                       ## and stored in the variable 'localDesign'.
                                       
                                       for(i in idLevels) {
                                         posi <- positions[marks == variable & id == i]
                                         ## posi is sorted for a valid data object. This is
                                         ## assumed in the following computation.                                        

                                         xt <- evalPositions[[i]]
                                         nt <- length(xt)
                                         xs <- posi
                                         ns <- 1
                                         xZ <- numeric(nt*r)
                                         d <- model@Delta
                                         antip <- zero
                                         w <- (r-1)*d
                                                      
                                         for(ii in 1:nt) {
                                           target = xt[ii] + antip;
                                           while(ns < length(xs) && target > xs[ns+1])
                                             ns <- ns + 1;
                                             nss = ns;
                                             diff = target - xs[ns];
                                             if(diff > 0) {
                                               while(diff <= w) 
                                                 { 
                                                   lookupIndex = floor(diff/d + 0.5);
                                                   entry = ii + nt*lookupIndex;
                                                   xZ[entry] <- xZ[entry] + 1;
                                                   ns <- ns - 1;
                                                   if(ns < 1) break;
                                                   diff = target - xs[ns];
                                                 }
                                             }
                                             ns = nss;
                                           }

                                         designList[[i]] <- Matrix(xZ, nrow = nt, sparse = TRUE)
                                       }
                                       localDesign <- do.call("rBind", designList)
                                       colnames(localDesign) <- paste(term, 1:r, sep = "")
                                       localDesign ## The return value
                                     },     
                                              mc.preschedule = FALSE 
                                     ) ## End lapplyParallel
            
            assign <- unlist(lapply(kernels,
                                    function(i) {
                                      rep(i, r)
                                    }
                                    )
                             )
            names(design) <- termLabels[kernels]
            kernelMatrix <- do.call("cBind", design)
            form <- formula(model)
            attr(form, "kernelTerms") <- kernels
            model <- updateKernelMatrix(model, kernelMatrix, assign, form)
            lockEnvironment(model@kernelMatrixEnv, binding = TRUE)

            formula(model) <- formulaNoSpecials
            model <- computeModelMatrix(as(model, "PointProcessModel"),
                                        evaluationPositions = evaluationPositions,
                                        ...)
            formula(model) <- formula
            
            return(model)
          }
          )

setMethod("coefficients", "PointProcessKernel",
          function(object,...){
            list(object@coefficients, object@g)
          }
          )
            
setMethod("computeLinearPredictor", "PointProcessKernel",
          function(model, coefficients = NULL, ...) {
            if(is.null(coefficients)) 
              coefficients <- coefficients(model)
            as.numeric(getModelMatrix(model) %*% coefficients[[1]]) +
              as.numeric(getKernelMatrix(model) %*% coefficients[[2]])             
          }
          )

setMethod("computeDMinusLogLikelihood", "PointProcessKernel",
          function(model, coefficients = NULL, eta = NULL, ...) {
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            if(is.null(eta))
              eta <- computeLinearPredictor(model, coefficients, ...)

            kM <- getKernelMatrix(model)
              
            if(model@family@link == "log") {

              dmll <- (as.vector(t(exp(eta)*model@delta) %*% kM) -
                colSums(kM[getPointPointer(processData(model), response(model)), , drop = FALSE])) %*% model@gram

            } else {
              
              etaP <- eta[getPointPointer(processData(model), response(model))]
              mmP <- kM[getPointPointer(processData(model), response(model)), , drop = FALSE]

              dmll <-  (as.vector(t(model@family@Dphi(eta)*model@delta) %*% kM) -
                as.vector(t(model@family@Dphi(etaP)/model@family@phi(etaP)) %*% mmP)) %*% model@gram
              
            }

            c(dmll, computeDMinusLogLikelihood(as(model, "PointProcessModel"),
                                               eta = eta, ...))
          }
          )

setMethod("getKernelMatrix", c(model = "PointProcessKernel", col = "ANY"),
          function(model, col,...){
            if(missing(col))
              col <- model@kernelMatrixCol
            if(length(col) == 0) {
              kernelMatrix <- model@kernelMatrixEnv$kernelMatrix
            } else {
              kernelMatrix <- model@kernelMatrixEnv$kernelMatrix[, col, drop = FALSE]
            }
            return(kernelMatrix)
          }
          )

setMethod("updateKernelMatrix", "PointProcessKernel",
          function(model, kernelMatrix = getKernelMatrix(model), assign = getAssign(model), form){
            force(kernelMatrix)
            model@kernelMatrixEnv <- new.env(parent = emptyenv())
            model@kernelMatrixEnv$kernelMatrix <- kernelMatrix
            model@kernelMatrixEnv$assign <- assign
            if(!missing(form))
              model@kernelMatrixEnv$formula <- form
            model@kernelMatrixCol <- numeric()
            return(model)
          }
          )
