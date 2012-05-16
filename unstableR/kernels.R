### The reproducing kernel for the 2-Sobolev Hilbert space on [0,t].

sobolevKernel <- function(s, r, t = 1, sub = NULL) {
  smin <- pmin(s, r)/t
  r <- pmax(s, r)/t
  s <- smin
  if (is.null(sub)) {
    1 + s * r / 48 +
      (1 - r) * s * (2 * r - r * r - s * s) / 6
  } else if (sub == 0) {
    1 + s * r / 48
  } else {
    (1 - r) * s * (2 * r - r * r - s * s) / 6
  }
}

### A Gaussian kernel

gaussianKernel <- function(s, r, t = 1, c = 5e3) {
  s <- s/t
  r <- r/t
##  1 + s * r +
    exp( - c * (s - r)^2)
}

grid <- seq(0, 10, 0.05)
support <- c(0, 10)

grid <- model@basisPoints
G1 <- outer(grid, grid, kernel, t = support[2], sub = 1)
G <- outer(grid, grid, kernel, t = support[2])
G1svd <- svd(G1, nv = 0)
d <- G1svd$d/G1svd$d[1] > 5*.Machine$double.eps   ### Computationally non-zero singular values.
G1svd$d <- ifelse(d, G1svd$d, 1)
U <- G %*% t(t(G1svd$u)/sqrt(G1svd$d))





image(Matrix(U))
diag(t(U) %*% G1 %*% U)
t(U) %*% G1 %*% U


G <- outer(seq(0, 10, 0.05), seq(0, 10, 0.05), sobolevKernel, t = 10, sub = 1)
G <- outer(seq(0, 10, 0.05), seq(0, 10, 0.05), gaussianKernel, t = 10)
tmp <- svd(G)
plot(log10(tmp$d))
plot(tmp$u[, 600], type = "l")


OmegaSvd <- svd(Omega[2:208, 2:208])
C <- diag(1, 208)
C[2:208, 2:208] <- t(t(OmegaSvd$u)/sqrt(OmegaSvd$d))

mll <- function(par,...) ppstat:::computeMinusLogLikelihood(model, C %*% par, ...) + as.numeric(par[-1] %*% par[-1])
dmll <- function(par,...) t(C) %*% ppstat:::computeDMinusLogLikelihood(model, C %*% par, ...) + 2*as.numeric(c(0, par[-1]))

initPar <- model@coefficients
initPar[] <- 0
initPar[1] <- 1

tmp <- optim(initPar, fn = mll, gr = dmll, method = "BFGS", control = list(trace = 5))

summary(Tpps)
model <- Tpps
termPlot(model) + coord_cartesian(ylim = c(-0.3, 0.7)) + truth
model@coefficients[] <- as.numeric(C %*% tmp$par)

model@coefficients %*% Omega %*% model@coefficients
Tpps@coefficients %*% Omega %*% Tpps@coefficients


################################################################################
###
################################################################################



### One example of a history based intensity; stochastic integral of exponentials

h0 <- function(t) exp(-(t-0.5)^2) ## True filter function

lambda0 <- function(t, T = numeric(), alpha = 1, beta = 0.5, ...) 
  beta * sum(h0((t - T[T < t]))) + alpha

### For later plotting and integration.

lambda1 <- Vectorize(lambda0, "t")

### Implementation of Ogata's thinning algorithm

OgataThin <- function(n = 1, ...) {
  T <- numeric(n)
  t <- 0
  i <- 0
  K <- lambda0(0)
  while(i < n) {
    S <- rexp(1, K)
    U <- runif(1) 
    t <- t + S
    lambda <- lambda0(t, T[0:i], ...) ## T[0:0] is numeric(0), OK
    if(U < lambda/K) {
      i <- i + 1
      T[i] <- t
      K <- lambda + lambda0(0)
    } else {
      K <- lambda
    }
  }
  T
} 

### Example simulation

N <- 999
T <- OgataThin(N)

library(ppstat)  ## loads the 'processdata' package as well

## Creates the marked point process object with suitable time grid.

Tmpp <- markedPointProcess(T[T < 200], seq(0, 200, 0.05))
summary(Tmpp)

Tppm <- pointProcessModel(point ~ h0(point),
                          data = Tmpp,
                          family = Hawkes(link = "identity"),
                          support = 4
                          )
termPlot(Tppm)

Tppk <- as(Tppm, "PointProcessKernel")
Tppk@formula <- point ~ k(point)
setGeneric("computeModelMatrix", function(model, evaluationPositions=NULL, ...) standardGeneric("computeModelMatrix"))
tmp <- computeModelMatrix(Tppk)
