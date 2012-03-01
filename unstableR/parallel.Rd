\name{parallel}
\Rdversion{1.1}
\alias{parallelExecution}
\alias{.registerParBackend}

\title{
  Parallel execution of point process model fitting and model selection. 
}

\description{

  Parallel execution of several of the fundamental computations in
  \pkg{ppstat} is supported via the \pkg{foreach} looping construct.
  The function is a convenience wrapper that attempts to register some
  backend supporting multiple cores or CPUs. 

}

\usage{
ppstat:::.registerParBackend(backends = "doMC", cores = NULL)
}

\arguments{
  
  \item{backends}{ a \code{character} vector. The default value
    \code{"doMC"} is the only backend currently
    supported by the function.}
  
  \item{cores}{ a \code{numeric}. The number of cores to use. The
    default value \code{NULL} gives the number of cores available for
    the \code{doMC} backend and one core otherwise.}
  
}

\details{

  Several different fundamental loops in the model building, fitting
  and selection algorithms in \pkg{ppstat} are implemented using a
  \code{foreach} loop with the \code{\%dopar\%} operator, which results
  in parallel evaluation of the expressions in the loop, if a parallel
  backend is registered.

  The default behavior of the \pkg{ppstat} package upon loading is to
  look for the suggested package \pkg{doMC}, and if is is
  available and R does not run in a GUI, it is loaded and the doMC
  backend is registred. Otherwise the sequential backend is registred
  and all computations are done sequentially. At load time \pkg{ppstat}
  checks for an already registered parallel backend using the function
  \code{getDoParRegistered}. If a backend is registered, the
  user is warned that the backend is going to be changed.
  
  The \pkg{doMC} is based on the \pkg{multicore} package, which does not
  run on Windows machines, and should not be used when running a GUI.
  Alternative solutions using e.g. \pkg{doSMP} from Revolution Analytics
  (available on R-Forge) or \pkg{doSNOW} may be developed, but do
  currently not work 'out of the box' as the doMC backend does. This
  effectively limits the use of multiple cores for this package to unix-like
  machines with R not running in a GUI. 

}

\seealso{

  \code{\link{pointProcessModel}}, \code{\link{stepInformation}},
  \code{\link{foreach}}, \code{\link{getDoParRegistered}}.

  }

\author{
  Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
}


\keyword{utilities}
