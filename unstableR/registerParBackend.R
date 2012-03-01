.registerParBackend <- function(backends = "doMC", cores = NULL) {
  if(!(class(cores) %in% c("NULL", "numeric")))
    stop("Variable 'cores' must be numeric or NULL")
  
  ## This checks if the doMC package is installed and if the GUI is
  ## appropriate for using the multicore backend for parallel
  ## computations.
  if("doMC" %in% backends && .Platform$GUI %in% c("X11", "unknown") && require("doMC")){
    registerDoMC(cores)
    cat("Registering the doMC parallel backend. \n")
  } else {
    ## Registering the sequential backend to avoid 1st warning from
    ## foreach when no parallel backend is registered.
    registerDoSEQ()
    cat("Backend set for sequential execution. \n")
  }

  ## The following does in principle work, that is, if doSMP is available
  ## it loads and subsequent foreach calls seem to run in parallel, but
  ## computations seem to require specific packages loaded in the foreach
  ## calls etc. It does not work "out of the box" as doMC does.
  ## if("doSMP" %in% backends && require(doSMP)) {
  ##   if(is.null(cores)) {
  ##     workers <- startWorkers(1) ## Should figure out the number of available cores on a windows machine, for now only 1 is assumed ...
  ##   } else {
  ##     workers <- startWorkers(cores)
  ##   }
  ##   registerDoSMP(workers)
  ##   cat("Registering the doSMP parallel backend. \n")
  ## } 

  ## The following does _not_ work correctly! Clusters seem to have to
  ## be set up and explicitly stopped between every foreach call?
  ## if("doSNOW" %in% backends && require(doSNOW)) {
  ##   if(is.null(cores)) {
  ##     cl <<- makeCluster(1, type='SOCK')  ## Should figure out the number of available cores on the machine, for now only 1 is assumed ...
  ##   } else {
  ##     cl <<- makeCluster(cores, type='SOCK')
  ##   }
  ##   registerDoSNOW(cl)
  ##   cat("Registering the doSNOW parallel backend. \n")
  ## } 
}
