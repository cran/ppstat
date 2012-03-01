.onLoad <- function(libname, pkgname) {

  ## TODO: This is not safe! Default must be not to use multiple cores
  ## perhaps unless specified for the package at install time. 
  
  ## Registration of parallel backend
  if (getDoParRegistered()) {
    warning("The registered parallel backend", getDoParName(), "is about to be changed.\n")
  } else {
    .registerParBackend(backend = "doMC")
  }  
}

.onUnload <- function(libpath)
{
    library.dynam.unload("ppstat", libpath)
}
