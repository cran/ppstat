source("~/programs/ppstat/pkg/unstableR/empirical.R")

x <- c(1, 4, 7, 35, 44)
checkIdentical(computeDistanceTable(x, 5),
               list(dist = c(3, 3),
                    left = c(1L, 2L),
                    right = c(2L, 3L)))

checkIdentical(computeDistanceTable(x, 10),
               list(dist = c(3, 6, 3, 9),
                    left = c(1L, 1L, 2L, 4L),
                    right = c(2L, 3L, 3L, 5L)))

checkIdentical(computeDistanceTable(x, 50),
               list(dist = c(3, 6, 34, 43, 3, 31, 40, 28, 37, 9),
                    left = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 4L),
                    right = c(2L, 3L, 4L, 5L, 3L, 4L, 5L, 4L, 5L, 5L)))

checkIdentical(computeDistanceTable(x, 2),
               list(dist = numeric(),
                    left = integer(),
                    right = integer()))

x <- cumsum(rexp(1000))
d <- outer(x, x, "-")
d0 <- c(d[d > 0])
checkEquals(computeDistanceTable(x, 10)$dist,
            d0[d0 < 10])
checkEquals(computeDistanceTable(x, 20)$dist,
            d0[d0 < 20])
checkEquals(computeDistanceTable(x, 1000)$dist,
            d0[d0 < 1000])


N <- 1e5
T <- N
x1 <- cumsum(rexp(N))
## x2 <- cumsum(rexp(N)) ##Null case
x2 <- rnorm(length(x1), x1+0.05, sd = 0.01)
lab <- factor(c(rep("A", N), rep("B", N)))
x <- c(x1, x2)
ord <- order(x)
lab <- lab[ord]
x <- x[ord]
lab <- lab[x < T]
x <- x[x < T]
system.time(dt <- as.data.frame(computeDistanceTable(x, 0.1)))
ed <- sort(dt[lab[dt$left] == "A" & lab[dt$right] == "B", "dist"])
plot(ed, seq_along(ed)/T, type = "l")

ed <- sort(dt[lab[dt$left] == "B" & lab[dt$right] == "A", "dist"])
plot(ed, seq_along(ed)/T, type = "l")

x <- split(getPointPosition(encodeData), getPointId(encodeData))
dt <- lapply(x, computeDistanceTable, H = 1000)
ip <- split(seq_along(getPointId(encodeData)), getPointId(encodeData))
dt <- sapply(names(ip), function(i) {
  dt[[i]][["left"]] <- ip[[i]][dt[[i]][["left"]]]
  dt[[i]][["right"]] <- ip[[i]][dt[[i]][["right"]]]
  return(dt[[i]])}, simplify = FALSE)
dt <- sapply(dt, as.data.frame, simplify = FALSE)
id <- rep(names(dt), times = sapply(dt, function(k) dim(k)[1]))
dt <- do.call(rbind, dt)
marks <- getMarkType(encodeData)
dt$left <- marks[dt$left]
dt$right <- marks[dt$right]
dt$id <- id
dt <- dt[order(dt$dist),]
library(ggplot2)
ggplot(dt, aes(x = dist)) + geom_density() + facet_grid(right~left) +
  coord_cartesian(ylim = c(0,0.002))

cum <- function(x) seq_along(x)
dt.cum <- ddply(dt, .(left, right, id), transform,
                 cum = cum(dist))

ggplot(dt.cum, aes(x = dist, y = cum, colour = id, group = id)) + geom_line() + facet_grid(right~left) + coord_cartesian(ylim = c(0, 50)) + opts(legend.position = "none")


ggplot(dt.cum, aes(x = dist, y = cum)) + geom_smooth() + geom_line(aes(colour = id)) + facet_grid(right~left) + coord_cartesian(ylim = c(0, 50)) + opts(legend.position = "none")

 
