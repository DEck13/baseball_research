

## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
numCores <- detectCores() - 1

## Ordinary Ftilde
Ftilde <- function(y, t){
  y <- sort(y)
  n <- length(y)
  ytilde <- rep(0, n + 1)
  ytilde[1] <- y[1] - 1/(y[2] - y[1]+1)
  ytilde[n+1] <- y[n] + 1/(y[n] - y[n-1])
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (y[j]+y[j-1])/2 
  }))
  
  j <- length(which(ytilde < t))
  (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  
}

## Adjusted Ftilde
Ftilde_adj <- function(y, t){
  y <- sort(y)
  n <- length(y)
  rank_d <- rank(diff(y))
  ytilde <- rep(0, n + 1)
  ytilde[1] <- y[1] - 1/(y[2] - y[1]+1)
  ytilde[n+1] <- y[n] + 1/(y[n] - y[n-1])
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    y[j] * (1 - rank_d[j-1]/ (n-1)) + y[j-1] * rank_d[j-1] / (n-1)
  }))
  
  j <- length(which(ytilde < t))
  (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  
}


order_pbino <- function(p = 0, k = 1, n = 1e4){
  pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
}

order_bino_vec <- function(p){
  p <- sort(p) # just in case
  n <- length(p)
  unlist(lapply(1:n, function(j){
    order_pbino(p[j], k = j, n = n)
  }))
}

## Simulation 
set.seed(13)
n <- 1000
yy <- sort(rnorm(n))
yy[n-1] <- 4.99
yy[n] <- 5
yy <- yy *100
us <- unlist(lapply(yy, function(xx) Ftilde(y = yy, t = xx)))
a <- seq(0.01,0.99,0.01) *100 # 99 numbers 
b <- c(seq(0.01,0.883, 0.009), 0.99) *100 # 99 numbers 
## Ftilde is the function that build \widetilde{F}_{Y}(t) 
ra <- unlist(lapply(a, function(xx) Ftilde(y = a, t = xx)))
rb <- unlist(lapply(b, function(xx) Ftilde(y = b, t = xx)))
ra[99]
rb[99]
a <- c(3,5,9,10)
a <- seq(0.01,0.99,0.01) *100 # 99 numbers 
b <- c(seq(0.01,0.874, 0.009), 0.989,0.99) *100 # 99 numbers 
## Ftilde is the function that build \widetilde{F}_{Y}(t) 
ra <- unlist(lapply(b, function(xx) Ftilde(y = b, t = xx)))
rb <- unlist(lapply(b, function(xx) Ftilde_adj(y = b, t = xx)))
ra[98:99]
rb[98:99]
us[999:1000]

a <- batters %>% filter(AB > AB_thresh) %>% filter(yearID == 1959)
c <- sort(a$AVG)
ra <- unlist(lapply(c, function(xx) Ftilde(y = c, t = xx)))
rb <- unlist(lapply(c, function(xx) Ftilde_adj(y = c, t = xx)))
ra[114:115]
rb[114:115]

order_bino_vec(ra)
## order_Pareto_vec function 
# this function transform percentiles from order stats (in increasing order)
# to Pareto values corresponding to the general population 
# of a greater than or equal to size
# default alpha is that of the Pareto principle 80-20
order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
  n <- length(u)
  if(length(npop) == 1) npop <- rep(npop, n)
  unlist(lapply(1:n, function(j){
    qPareto(qbeta(u[j], j + npop[j], n + 1 - j), t = t, alpha = alpha)
  }))
}

## map_Pareto_vals_vec function 
# this function transforms ordered Pareto values corresponding to 
# the general population to percentiles from order stats 
# (in increasing order)
map_Pareto_vals_vec <- function(x, t = 1, alpha = 1.16, npop = 1e4){
  n <- length(x)
  if(length(npop) == 1) npop <- rep(npop, n)  
  unlist(lapply(1:n, function(j){
    pbeta(pPareto(x[j], t = t, alpha = alpha), j + npop[j], n + 1 - j)
  }))
}

## Simulation 
order_Pareto_vec(u = us, npop = 1e7)[999:1000]
max(abs(us - map_Pareto_vals_vec(
  order_Pareto_vec(u = us, npop = 1e7), npop = 1e7)))
us


yy <- sort(yy)
n <- length(yy)
ytilde <- rep(0, n + 1)
ytilde[1] <- yy[1] - 1/(yy[2] - yy[1])
ytilde[n+1] <- yy[n] + 1/(yy[n] - yy[n-1])
ytilde[2:n] <- unlist(lapply(2:n, function(j){
  (yy[j]+yy[j-1])/2 
}))

## Same function as Yempirical

map_Y <- function(u, ytilde){
  n <- length(ytilde)-1
  seqence <- seq(0, 1, 1/n)
  pos <- findInterval(u, seqence)
  out <- (n*u -pos + 1) * (ytilde[(pos+1)] - ytilde[pos]) + ytilde[pos]
  return(out)
}

## Simulation
y_emprical <- sapply(1:length(us), function(x) map_Y(us[x], ytilde = ytilde))

yy - y_emprical

Yempirical <- function(u, ytilde){
  if(length(u) + 1 != length(ytilde)){
    stop("y must be the same dimension")
  } 
  n <- length(u)
  #y <- sort(y)
  #ytilde <- rep(0, n + 1)
  #ytilde[1] <- y[1] - 1/(y[2] - y[1])
  #ytilde[n+1] <- y[n] + 1/(y[n] - y[n-1])
  #ytilde[2:n] <- unlist(lapply(2:n, function(j){
  #  (y[j] + y[j-1])/2 
  #}))
  #qs <- qbeta(u, shape1 = 1:n, shape2 = n:1)
  unlist(lapply(1:n, function(j){
    out <- (n*u[j] - j + 1)*(ytilde[j+1] - ytilde[j]) + ytilde[j]
    #if(j %in% c(1,n)){
    #  print(j)
    #  out <- (n*u[j] - j + 1)*(ytilde[j+1] - ytilde[j])/sqrt(n) + ytilde[j]
    #}
    out
  }))
}

yy - Yempirical(u = us, ytilde = ytilde)


## map the quantile to the predicated sample values
order_qempirical <- function(u){
  n <- length(u)
  a <- qbeta(u, shape1 = 1:n, shape2 = n:1)
  out <- sapply(1:n, function(x) map_Y(a[x], ytilde = ytilde))
  out
}

y_emprical <- order_qempirical(a)
yy - y_emprical


## Atbat threshold
AB_thresh <- 320 # rougly 2 per game (same as Gould)
## years to consider
years <- 1880:2019
ncores <- detectCores() - 2

## load in local csv file
batters <- read_csv("batters_adj.csv") %>%
  select(-pops)

## decenial population counts from Eck (2020)
pops <- c(4.40, 5.010, 5.580, 8.550, 8.930, 9.920, 11.130, 11.590 + 1, # 1950
          19.580, 26.100, 36.370, 40.680, 64.060, 76.190, 98) * 1e6

## interpolated population counts
pops <- approx(x = 1880 + 0:14*10, y = pops, xout = years)$y


## decenial weights (approximate). These weights account for: ## weights to be one 
# 1) AL/NL differences in integration (see Armour), 
# 2) Gallup polling of baseball interest with lag factor 
# 3) decreases in talent pool size due to WWI/WWII
# We make extrapolations for pre-1937 Gallup data and 
# post integration differences in leagues
#weightsNL <- c(0.12, 0.13, 0.14, 0.15, 0.20, 0.23, 0.25, 0.20,
 #             0.35, 0.35, 0.33,
#              0.31, 0.25, 0.22, 0.20, 0.18)
# weightsAL <- c(0.12, 0.13, 0.14, 0.15, 0.20, 0.23, 0.25, 0.20,
#               0.25, 0.25, 0.29,
#               0.30, 0.25, 0.22, 0.20, 0.18)
weightsNL <- rep(1, 16)
weightsAL <- rep(1, 16)

## interpolated population counts
weightsNL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsNL, xout = years)$y
weightsAL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsAL, xout = years)$y
w_popsNL <- pops * weightsNL
w_popsAL <- pops * weightsAL

## add weighted populations to data
pops_data <- rbind(
  data.frame(pops = w_popsNL, lgID = "NL", yearID = 1880 + 0:(14*10 - 1)),
  data.frame(pops = w_popsAL, lgID = "AL", yearID = 1880 + 0:(14*10 - 1)))
batters <- merge(batters, pops_data, by = c("lgID", "yearID"))

batters <- batters %>% filter(AB >= AB_thresh)
xx = 1918
k <- 100

AVG_leader_board <- do.call(rbind, mclapply(years, mc.cores = numCores, function(xx){
batters %>% filter(yearID == xx) %>% filter(AB >= AB_thresh) %>% 
arrange(AVG) %>% mutate(scale_AVG = AVG*k) %>% mutate(
AVG_talent = order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_AVG, function(xx) 
  Ftilde(y = scale_AVG, t = xx)))), npop = pops))  %>% 
    select(playerID, yearID, AVG, scale_AVG, AVG_talent, pops) #%>% 
#filter(OPS_plus_talent > 1e5)
})) %>% arrange(-AVG_talent) 

## adjusted Ftilde
AVG_leader_board <- do.call(rbind, mclapply(years, mc.cores = numCores, function(xx){
  batters %>% filter(yearID == xx) %>% filter(AB > AB_thresh) %>% 
    arrange(AVG) %>% mutate(scale_AVG = AVG*k) %>% mutate(
      AVG_talent = order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_AVG, function(xx) 
        Ftilde_adj(y = scale_AVG, t = xx)))), npop = pops))  %>% 
    select(playerID, yearID, AVG, scale_AVG, AVG_talent, pops) #%>% 
  #filter(OPS_plus_talent > 1e5)
})) %>% arrange(-AVG_talent) 


foo <- AVG_leader_board %>% filter(yearID < 1997)
foo$playerID <- droplevels(as.factor(foo$playerID))
bar <- split(foo, f = foo$playerID)
baz <- (do.call(rbind, mclapply(bar, mc.cores = ncores, function(xx){
  xx[which.max(xx$AVG_talent), ]
})) %>% arrange(-AVG_talent))[1:25, ]

ref_1997 <- AVG_leader_board %>% filter(yearID == 1997)
w_pops_1997 <- unique(ref_1997$pops)

## Ordinary Ftilde
yy <- sort(ref_1997$scale_AVG)
n <- length(yy)
ytilde <- rep(0, n + 1)
ytilde[1] <- yy[1] - 1/(yy[2] - yy[1])
ytilde[n+1] <- yy[n] + 1/(yy[n] - yy[n-1])
ytilde[2:n] <- unlist(lapply(2:n, function(j){
  (yy[j]+yy[j-1])/2 
}))

## Adjusted Ordinary Ftilde
y <- sort(ref_1997$scale_AVG)
n <- length(y)
rank_d <- rank(diff(y))
ytilde <- rep(0, n + 1)
ytilde[1] <- y[1] - 1/(y[2] - y[1])
ytilde[n+1] <- y[n] + 1/(y[n] - y[n-1])
ytilde[2:n] <- unlist(lapply(2:n, function(j){
  y[j] * (1 - rank_d[j-1]/ (n-1)) + y[j-1] * rank_d[j-1] / (n-1)
}))

top25seasons_bf1997_AVG <- as.data.frame(do.call(rbind, mclapply(1:25, mc.cores = ncores, function(j){
  xx <- baz[j, ]
  xx$pops <- w_pops_1997
  xx$playerID <- paste(xx$playerID, "_proj", sep = "")
  int <- rbind(ref_1997, xx) %>% arrange(AVG_talent) %>% 
    mutate(adj_AVG = order_qempirical(map_Pareto_vals_vec(AVG_talent, npop = pops))/100)
  int[grepl("_proj", int$playerID), ]
})))

b <- merge(top25seasons_bf1997_AVG, a, by = c("playerID"))

top25seasons_bf1997_AVG %>% filter(yearID < 1950)
write_csv(top25seasons_bf1997_AVG, "bf1996_interpolation_adjusted.csv")

plot(c(7,9.5), c(1/2,3/4), xlim = c(4,10), ylim = c(0.4,0.85), type = "b", 
     xlab = "tilde{Y}_{(j)}", ylab = "widetilde{F}_{Y}(Y_{(j)})")

abline(h = 0.5, v = 7, lty = 2)
abline(v = 9.5, lty = 2)
abline(v = 9, lty = 4, col="red")
lines(c(5,9.5), c(1/2,3/4), col = "blue")
lines(c(7,9), c(1/2,3/4), col = "blue")
abline(h = 3/4, lty = 2)

abline(v = 5, lty = 2)

x <- c(2.5, 4, 7, 9.5, 11)
y <- c(0, 1/4, 1/2, 3/4, 1)
plot(x, y, type = "b")
abline(h = 0, lty = 2)
abline(h = 1/4, v = 4, lty = 2)
abline(h = 2/4, v = 7, lty = 2)
abline(h = 3/4, v = 9.5, lty = 2)
abline(h = 1, v = 11, lty = 2)
abline( v = 4, lty = 2)
abline(v = 7, lty = 2)
abline( v = 9.5, lty = 2)
abline( v = 11, lty = 2)
abline(v = 3, lty = 4, col="red")
abline(v = 5, lty = 4, col="red")
abline(v = 9, lty = 4, col="red")
abline(v = 10, lty = 4, col="red")


######################################
## Try it on baseball data using different statistics
######################################
k  = 1
xx = 1980
years = seq(1880,2019,1)
WAR_leader_board <- do.call(rbind, mclapply(years, mc.cores = numCores, function(xx){
  a <- batters %>% filter(yearID == xx) %>% filter(AB > AB_thresh) %>% 
    arrange(WAR) %>% mutate(scale_WAR = WAR*k) %>% 
   mutate(WAR_talent = 
            order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_WAR, function(xx) 
     Ftilde(y = scale_WAR, t = xx)))), npop = pops))  %>% 
    select(playerID, yearID, WAR, scale_WAR,WAR_talent, pops)
})) %>% arrange(-WAR_talent)  

ruth1926 <- WAR_leader_board %>% filter(yearID == 1926)
ruth1923 <- WAR_leader_board %>% filter(yearID == 1923)
ruth1921 <- WAR_leader_board %>% filter(yearID == 1921)
a <- batters %>% filter(yearID == 1991)
k = 6
## adjusted Ftilde
WAR_leader_board <- do.call(rbind, mclapply(years, mc.cores = numCores, function(xx){
  batters %>% filter(yearID == xx) %>% filter(AB > AB_thresh) %>% 
    arrange(WAR) %>% mutate(scale_WAR = WAR*k) %>% mutate(
      WAR_talent = order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_WAR, function(xx) 
        Ftilde_adj(y = scale_WAR, t = xx)))), npop = pops))  %>% 
    select(playerID, yearID, WAR, scale_WAR, WAR_talent, pops) #%>% 
  #filter(OPS_plus_talent > 1e5)
})) %>% arrange(-WAR_talent) 

## top 100
sum(WAR_leader_board[1:100,]$yearID <= 1950)
sum(WAR_leader_board[1:50,]$yearID <= 1950)


a <- WAR_leader_board[1:200,]
write.csv(a, "top200_WAR_nonpara.csv")

a <- batters %>% filter(yearID == 1975)

foo <- WAR_leader_board %>% filter(yearID < 1997)
foo$playerID <- droplevels(as.factor(foo$playerID))
bar <- split(foo, f = foo$playerID)
baz <- (do.call(rbind, mclapply(bar, mc.cores = ncores, function(xx){
  xx[which.max(xx$WAR_talent), ]
})) %>% arrange(-WAR_talent))[1:25, ]

ref_1997 <- AVG_leader_board %>% filter(yearID == 1997)
w_pops_1997 <- unique(ref_1997$pops)

## Ordinary Ftilde
yy <- sort(ref_1997$scale_AVG)
n <- length(yy)
ytilde <- rep(0, n + 1)
ytilde[1] <- yy[1] - 1/(yy[2] - yy[1])
ytilde[n+1] <- yy[n] + 1/(yy[n] - yy[n-1])
ytilde[2:n] <- unlist(lapply(2:n, function(j){
  (yy[j]+yy[j-1])/2 
}))

## Adjusted Ordinary Ftilde
y <- sort(ref_1997$scale_AVG)
n <- length(y)
rank_d <- rank(diff(y))
ytilde <- rep(0, n + 1)
ytilde[1] <- y[1] - 1/(y[2] - y[1])
ytilde[n+1] <- y[n] + 1/(y[n] - y[n-1])
ytilde[2:n] <- unlist(lapply(2:n, function(j){
  y[j] * (1 - rank_d[j-1]/ (n-1)) + y[j-1] * rank_d[j-1] / (n-1)
}))

top25seasons_bf1997_AVG <- as.data.frame(do.call(rbind, mclapply(1:25, mc.cores = ncores, function(j){
  xx <- baz[j, ]
  xx$pops <- w_pops_1997
  xx$playerID <- paste(xx$playerID, "_proj", sep = "")
  int <- rbind(ref_1997, xx) %>% arrange(AVG_talent) %>% 
    mutate(adj_AVG = order_qempirical(map_Pareto_vals_vec(AVG_talent, npop = pops))/100)
  int[grepl("_proj", int$playerID), ]
})))


plot(batters$WAR)




