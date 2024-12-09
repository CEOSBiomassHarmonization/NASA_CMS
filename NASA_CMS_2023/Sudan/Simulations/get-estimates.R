rm(list=ls())
###############
## libraries ##

############################
## load synthetic dataset ##
load("synthetic-data.RData")

y.true <- mean(y.all)

###############################
## draw simple random sample ##

n <- 250 ## number of plots we want to imagine we collected
sub <- sample(1:nrow(coords.all), n)
y <- y.all[sub]
x <- x.all[sub]
coords <- coords.all[sub,]

plot(coords, pch = 19, col = hcl.colors(50)[cut(x,50)],
     main = "x",xlab = "", ylab = "")
plot(coords, pch = 19, col = hcl.colors(50)[cut(y,50)],
     main = "y",xlab = "", ylab = "")

############################
## design-based estimator ##
get.db.est <- function(y){
  
  n <- length(y)
  
  y.bar.hat <- sum(y)/n
  
  s2 <- sum((y - y.bar.hat)^2)/n
  
  y.se.hat <- sqrt(s2/n)
  
  y.ci <- y.bar.hat + c(-1,1)*1.96*y.se.hat
  
  return(list(y.bar.hat,y.se.hat,y.ci))
}


##############################
## model-assisted estimator ##
get.ma.est <- function(y,x,x.all){
  
  n <- length(y)
  N <- length(x.all)
  
  assisting.model <- lm(sqrt(y) ~ x)
  
  b0.hat <- coef(assisting.model)[1]
  b1.hat <- coef(assisting.model)[2]

  y.hat <- (b0.hat + b1.hat*x)^2
  y.pred <- (b0.hat + b1.hat*x.all)^2

  y.bar.hat <- sum(y.pred)/N + sum(y - y.hat)/n

  s2 <- sum((y - y.hat)^2)/n

  y.se.hat <- sqrt(s2/n)

  y.ci <- y.bar.hat + c(-1,1)*1.96*y.se.hat

  return(list(y.bar.hat,y.se.hat,y.ci))
}

n.iter <- 10000

y.est.db <- numeric(n.iter)
y.se.db <- numeric(n.iter)
y.ci.db <- matrix(ncol = 2, nrow = n.iter)

y.est.ma <- numeric(n.iter)
y.se.ma <- numeric(n.iter)
y.ci.ma <- matrix(ncol = 2, nrow = n.iter)

for(i in 1:n.iter){
    
    sub <- sample(1:nrow(coords.all), n)
    y <- y.all[sub]
    x <- x.all[sub]
    coords <- coords.all[sub,]
    
    db.est <- get.db.est(y)
    
    y.est.db[i] <- db.est[[1]]
    y.se.db[i] <- db.est[[2]]
    y.ci.db[i,] <- db.est[[3]]
    
    
    ma.est <- get.ma.est(y,x,x.all)
    
    y.est.ma[i] <- ma.est[[1]]
    y.se.ma[i] <- ma.est[[2]]
    y.ci.ma[i,] <- ma.est[[3]]

    cat(i, "\r")
}


hist(y.est.db)
abline(v = y.true, col = "blue", lwd = 2)

hist(y.se.db)
abline(v = sd(y.est.db), col = "blue", lwd = 2)

mean(y.ci.db[,1] < y.true & y.true < y.ci.db[,2])


hist(y.est.ma)
abline(v = y.true, col = "blue", lwd = 2)

hist(y.se.ma)
abline(v = sd(y.est.ma), col = "blue", lwd = 2)

mean(y.ci.ma[,1] < y.true & y.true < y.ci.ma[,2])
