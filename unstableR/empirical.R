library(inline)
library(Rcpp)
library(RUnit)
library(ppstat)

setGeneric("computeDistanceTable", function(x, H, h, ...) standardGeneric("computeDistanceTable"))
setGeneric("getRanges", function(object, ...) standardGeneric("getRanges"))

setMethod("getRanges", "MarkedPointProcess",
          function(object, ...) {
            splitPositions <- split(getPosition(object), getId(object))
            ranges <- lapply(splitPositions, function(p) c(p[1], p[length(p)]))
            return(ranges)
          }
          )
            
setCMethod("computeDistanceTable",
           signature(x = "numeric", H = "numeric", h = "numeric"),
            body = "
NumericVector xx(x);
int nx = xx.size(), i, ip, ii = 1, j = 0, m, id = 1;
IntegerVector bi(nx), cId(nx), count(nx+1), size(nx+1);
double HH = as<double>(H), hh = as<double>(h), diff;


cId[0] = id;
size[0] = 0;
size[1] = 0;
count[0] = 0;
count[1] = 1;

while(ii < nx & xx[ii] - xx[0] < HH) 
    ii++;

if(ii > 1) {
  bi[0] = ii - 1;
  m = ii - 1;
} else {
  ii = 2;
  m = 0;
}

for(i = 1; i < nx; i++) {
  diff =  xx[i] - xx[i-1];
  if(diff > hh) {
    id++;
    count[id] = 1;
    size[id] = 0;
  } else {
    size[id] = size[id] + diff;
    count[id]++;
  }

  cId[i] = id;

  while(ii < nx & xx[ii] - xx[i] < HH) 
    ii++;
  
  if(ii > i + 1) {
    bi[i] = ii - 1;
    m += ii - i - 1;
  } else {
    ii++;
  }
}


NumericVector dist(m);
IntegerVector left(m);
IntegerVector right(m);
IntegerVector clusterId(m);
IntegerVector clusterCount(m);
NumericVector clusterSize(m);

for(i = 0; i < nx; i++) {
  ip = i + 1;
  for(ii = ip; ii <= bi[i]; ii++) {
    left[j] = ip;
    right[j] = ii + 1;
    dist[j] = xx[ii] - xx[i];
    if(cId[ii] == cId[i])
      clusterId[j] = cId[i];
    clusterSize[j] = size[clusterId[j]];
    clusterCount[j] = count[clusterId[j]];
    j++;
  }
}

List result = List::create(
       _[\"dist\"] = dist,
       _[\"left\"] = left,
       _[\"right\"] = right,
       _[\"clusterId\"] = clusterId,
       _[\"clusterCount\"] = clusterCount,
       _[\"clusterSize\"] = clusterSize);

return(result);
",
            includes = "using namespace Rcpp;",
            Rcpp = TRUE)


## The following function computes the cluster information assuming
## that the pointers are ordered w.r.t. to the left pointer.  !!!!
## REMARK !!! This is a ridiculous way to compute this information,
## the function is obsolete, though it can serve some testing
## purposes. The computation is embedded in the function for distance
## table computations now.

setCMethod("computeCluster",
           signature(leftR = "numeric",
                     rightR = "numeric",
                     distR = "numeric",
                     cluster = "numeric"),
           body =  "
NumericVector left(leftR), right(rightR), dist(distR);
int i, n = left.size();
NumericVector clusterId(n), span(n+1), sspan(n), clusterSize(n+1);
double l = left[0], h = as<double>(cluster);

clusterId[0] = 1;
clusterSize[1] = 2;

for(i = 1; i < n; i++) {
  if(l != left[i]) { 
    l = left[i];
    if(dist[i] > h) {
      clusterId[i] = clusterId[i-1] + 1;
      clusterSize[clusterId[i]] = 1;
      sspan[i] = 0;
      span[clusterId[i]] = 0;
    } else {                    
      if(right[i-1] < l || clusterSize[clusterId[i-1]] == 1) {
        clusterId[i] = clusterId[i-1] + 1;
        clusterSize[clusterId[i]] = 2;
        sspan[i] = 0;
        span[clusterId[i]] = dist[i];
      } else { 
        clusterId[i] = clusterId[i-1];
        clusterSize[clusterId[i]] = clusterSize[clusterId[i]] + 1;
        sspan[i] = span[clusterId[i]];
        span[clusterId[i]] = span[clusterId[i]] + dist[i];
      }
    }
  } else { 
    clusterId[i] = clusterId[i-1];
    sspan[i] = sspan[i-1];
  }
}

NumericVector size(n);

for(i = 0; i < n; i++) {
  if(dist[i] > h & dist[i] > span[clusterId[i]]-sspan[i]) {
    size[i] = 0;
    clusterId[i] = 0;
  } else {
    size[i] = clusterSize[clusterId[i]];
}
}

List result = List::create(
       _[\"clusterId\"] = clusterId,
       _[\"clusterSize\"] = size,
_[\"span\"] = sspan);

return(result);
", includes = "using namespace Rcpp;",
            Rcpp = TRUE)

setMethod("computeDistanceTable", c(x = "MarkedPointProcess", H = "numeric", h = "missing"),
          function(x, H, h, shift = NULL, ...) {
            callGeneric(x = x, H = H, h = H, shift = shift, ...)
          }
          )

setMethod("computeDistanceTable", c(x = "MarkedPointProcess", H = "numeric", h = "numeric"),
          function(x, H, h, shift = NULL, size = NULL, ...) {
            id <- getPointId(x)
            ip <- split(seq_along(id), id)
            z <- getPointPosition(x)
            marks <- getMarkType(x)

            ## Shifting implemented to allow for easy, random shifting of some marks without having to recreate the entire MarkedPointProcess object
            if(is.matrix(shift)) {
              if(!identical(rownames(shift), names(ip)) || !identical(colnames(shift), levels(marks))) 
                stop("Either some row names or column names of 'shift' are invalid.", call = FALSE)
              if(is.null(size)) {
                ranges <- getRanges(x)
                for(i in names(ip)) {
                  z[ip[[i]]] <- (z[ip[[i]]] +
                                 shift[i, marks[ip[[i]]]]) %% ranges[[i]][2] +
                                   ranges[[i]][1]
                }
              } else {
                for(i in names(ip)) {
                  z[ip[[i]]] <- (z[ip[[i]]] %/% size) * size +
                    (z[ip[[i]]] + shift[i, marks[ip[[i]]]]) %% size
                }
              }
                
              ord <- order(id, z)
              z <-  z[ord]
              marks <- marks[ord]
            } else {
              if(!is.null(shift))
                warning("Argument 'shift' is not a 'matrix' as required and is ignored.", call = FALSE)
            }
            z <- sapply(names(ip), function(i) z[ip[[i]]], simplify = FALSE)
            dt <- lapply(z, computeDistanceTable, H = H, h = h)
            dt <- sapply(names(ip), function(i) {
              dt[[i]][["left"]] <- ip[[i]][dt[[i]][["left"]]]
              dt[[i]][["right"]] <- ip[[i]][dt[[i]][["right"]]]
              return(dt[[i]])}, simplify = FALSE)
            dt <- sapply(dt, as.data.frame, simplify = FALSE)
            id <- factor(rep(names(dt),
                             times = sapply(dt, function(k) dim(k)[1])),
                         levels = levels(id))
            dt <- do.call(rbind, dt)
            ## TODO: Should the function return pointers to the
            ## orignal positions of the points, or as now only the
            ## marks?
            ## if(h > 0) {
            ##    dt <- cbind(dt, as.data.frame(computeCluster(dt$left,
            ##                                                 dt$right,
            ##                                                 dt$dist,
            ##                                                 h)))
                ## Deprecated, buggy R code ##
                ##   left <- dt$left
                ##   right <- dt$right
                ##   dist <- dt$dist
                ##   n <- length(left)
                ##   span <- numeric(n)
                ##   clusterId <- numeric(n)
                ##   clusterSize <- numeric(n)
                ##   clusterId[1] <- 1
                ##   clusterSize[1] <- 2
                ##   l <- left[1]
                ##   for(i in 2:n) {
                ##     if(l != left[i]) { ## New left point.
                ##       l <- left[i]
                ##       if(dist[i] > h) { ## New "pseudo"-cluster
                ##         clusterId[i] <- clusterId[i-1] + 1
                ##         clusterSize[i] <- 1
                ##         span[clusterId[i]] <- 0
                ##       } else {                    
                ##         if(right[i-1] < l) { ## New "real"-cluster.
                ##           clusterId[i] <- clusterId[i-1] + 1
                ##           clusterSize[i] <- 2
                ##           span[clusterId[i]] <- dist[i]
                ##         } else { ## right[i-1] == l, same cluster.
                ##           clusterId[i] <- clusterId[i-1]
                ##           clusterSize[i] <- clusterSize[i-1] + 1 
                ##           span[clusterId[i]] <- span[clusterId[i]] + dist[i]
                ##         }
                ##       }
                ##     } else { ## Same left point
                ##       clusterId[i] <- clusterId[i-1]
                ##       clusterSize[i] <- clusterSize[i-1]
                ##     }
                ##   }
                ##   outsideCluster <- which(dist > span[clusterId])
                ##   clusterSize[outsideCluster] <- 0
                ##   clusterId[outsideCluster] <- 0
                ##   clusterSize <- as.vector(tapply(clusterSize, clusterId, max)[as.character(clusterId)])
                ##   dt <- cbind(dt, data.frame(clusterId = clusterId, clusterSize = clusterSize))
                ## }
             # }
            dt$left <- marks[dt$left]
            dt$right <- marks[dt$right]
            dt$id <- id
            dt <- dt[order(dt$dist),]
            return(dt)
          }
          )



shiftTable <- function(object, shiftMarks = levels(getMarkType(object))[-1], integer = FALSE, min = 0, size = NULL) {
  shift <- matrix(0, nrow = length(levels(getId(object))),
                  ncol = length(levels(getMarkType(object))),
                  dimnames = list(levels(getId(object)), levels(getMarkType(object))))
  if(is.null(size)) {
    ranges <- sapply(getRanges(object), function(r) r[2])
  } else {
    ranges <- rep(size, dim(shift)[1])
  }
  
  if(integer)
    runif <- function(n, max) sample.int(max, n)
  
  for(m in shiftMarks) {
    shift[ , m] <- sapply(ranges, function(r) runif(n = 1, max = r-2*min) + min)
  }
  return(shift)
}
                                     
computeMeanCum <- function(x, method = 'mean') {
  ## Assuming 'x' is a data frame with columns 'id', 'cum', and that
  ## the rows can be are ordered according to a 'dist' variable.
  ord <- order(x$dist)
  y <- x[ord, ]
  id <- as.factor(y$id)
  n <- dim(y)[1]
  cumMean <- numeric(n)
  iCum <- integer(length(levels(id)))
  names(iCum) <- levels(id)
  for(i in seq_len(n)) {
    iCum[id[i]] <- y$cum[i]
    cumMean[i] <- do.call(method, list(iCum))
  }
  x$rightRatesMean <- rep(mean(unique(x$rightRates)), dim(x)[1])
  x$cumMean <- numeric(n)
  x$cumMean[ord] <- cumMean
  return(x)
}

summarizeBoot <- function(x, method = function(x) quantile(x, c(0.05, 0.5, 0.95))) {
  ## Assuming 'x' is a data frame with columns 'L1', 'cumMean' and
  ## 'rightRatesMean', and that the rows can be are ordered according
  ## to a 'dist' variable.
  ord <- order(x$dist)
  y <- x[ord, ]
  id <- as.factor(y$L1)
  n <- dim(y)[1]
  z <- vector("list", max(x$dist))
  iz <- integer(length(levels(id)))
  names(iz) <- levels(id)
  for(i in seq_len(n)) {
    iz[id[i]] <- sqrt(y$cumMean[i]/y$rightRatesMean[i]) - sqrt(y$dist[i])
    z[[y$dist[i]+1]] <- c(dist = y$dist[i], do.call(method, list(iz)))
  }
  z <- as.data.frame(do.call(rbind, z))
  z$left <- rep(x$left[1], dim(z)[1])
  z$right <- rep(x$right[1], dim(z)[1])
  return(z)
}


cum <- function(x) {
  x$cum <- seq_len(dim(x)[1])/(counts[x$id[1], as.character(x$left[1])])
  x$rightRates <- rep(rates[as.character(x$id[1]), as.character(x$right[1])], dim(x)[1])
  return(x)
}
