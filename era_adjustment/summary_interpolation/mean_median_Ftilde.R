## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readr)
library(EnvStats)
library(extRemes)
library(gPdtest)
ncores <- detectCores() - 1

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/Rcode")

## input yearID and playID
## years to consider
years <- 1880:2019
# AB_thresh <- 320
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

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/pitching")
pitchers <- read_csv('pitchers_adj.csv')
pitchers_weighted <- read_csv('pitchers_adj_weighted.csv')

## Ordinary Ftilde
Ftilde <- function(y, t, thres){
  y <- sort(y)
  n <- length(y)
  ytilde <- rep(0, n + 1)
  if (thres == "median") {
    ytilde[1] <- y[1] - 1/(median(y) - y[1])
    ytilde[n+1] <- y[n] + 1/(y[n] - median(y))
  }
  if (thres == "mean") {
    ytilde[1] <- y[1] - 1/(mean(y) - y[1])
    ytilde[n+1] <- y[n] + 1/(y[n] - mean(y))
  }
  if (thres == "top"){
    ytilde[1] <- y[1] - 1/(y[3] - y[1])
    ytilde[n+1] <- y[n] + 1/(y[n] - y[n-2])
  }
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (y[j]+y[j-1])/2 
  }))
  
  j <- length(which(ytilde < t))
  (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  
}

## order_pbino function
# converts order stats to their percentiles
order_pbino <- function(p = 0, k = 1, n = 1e4){
  pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
}

## order_pbino_vec function
# this function converts a vector of order stats 
# to their percentiles. This vector should be the entire 
# sample sorted in increasing order
order_bino_vec <- function(p){
  p <- sort(p) # just in case
  n <- length(p)
  unlist(lapply(1:n, function(j){
    order_pbino(p[j], k = j, n = n)
  }))
}

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


map_Y <- function(u, ytilde){
  n <- length(ytilde)-1
  seqence <- seq(0, 1, 1/n)
  pos <- findInterval(u, seqence)
  out <- (n*u -pos + 1) * (ytilde[(pos+1)] - ytilde[pos]) + ytilde[pos]
  return(out)
}

## map the quantile to the predicated sample values
order_qempirical <- function(u, ytilde){
  n <- length(u)
  a <- qbeta(u, shape1 = 1:n, shape2 = n:1)
  out <- sapply(1:n, function(x) map_Y(a[x], ytilde = ytilde))
  out
}

## WAR of batters ##

WAR_batter_median <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  batters %>% filter(yearID == xx)  %>% 
    arrange(WAR) %>%  mutate(WAR_talent = order_Pareto_vec(u = 
    order_bino_vec(unlist(lapply(WAR, function(xx) 
      Ftilde(y = WAR, t = xx, thres = "median")))), npop = pops))
}))

WAR_batter_mean <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  batters %>% filter(yearID == xx) %>% 
    arrange(WAR) %>%  mutate(WAR_talent = order_Pareto_vec(u = 
    order_bino_vec(unlist(lapply(WAR, function(xx) 
    Ftilde(y = WAR, t = xx, thres = "mean")))), npop = pops))
}))
WAR_batter_top <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  batters %>% filter(yearID == xx) %>% 
    arrange(WAR) %>%  mutate(WAR_talent = order_Pareto_vec(u = 
    order_bino_vec(unlist(lapply(WAR, function(xx) 
    Ftilde(y = WAR, t = xx, thres = "top")))), npop = pops))
}))
setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/summary_interpolation")
write.csv(WAR_batter_median, "WAR_batter_median.csv")
write.csv(WAR_batter_mean, "WAR_batter_mean.csv")
write.csv(WAR_batter_top, "WAR_batter_test.csv")

## new Ftilde: yearID == 2016
## If the difference between the first two players is not large enough, 
## but they havea large discrepancy between the third or the forth player. 

Ftilde <- function(y,t){
  y <- sort(y)
  n <- length(y)
  ytilde <- rep(0, n + 1)
  ytilde[1] <- y[1] - 1/(y[3] - y[1])
  ytilde[n+1] <- y[n] + 1/(y[n] - y[n-2])
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (y[j]+y[j-1])/2 
  }))
  
  j <- length(which(ytilde < t))
  (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  
}

WAR_test <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  batters %>% filter(yearID == xx)  %>% 
    arrange(WAR) %>%  mutate(WAR_talent = order_Pareto_vec(u =
    order_bino_vec(unlist(lapply(WAR, function(xx) 
    Ftilde(y = WAR, t = xx)))), npop = pops))
})) %>% arrange(-WAR_talent)

write.csv(WAR_test, "WAR_batter_test.csv")

y <- as.matrix(batters %>% filter(yearID == 1991)  %>% select(WAR) %>% arrange(WAR) )
sapply(y, function(x) Ftilde(y,x, thres = "top") )
sapply(y, function(x) Ftilde(y,x, thres = "mean") )

a <- sapply(years, function(x) 
  as.matrix(batters %>% filter(yearID == x)  %>% 
              select(WAR) %>% arrange(-WAR))[1] - 
    (as.matrix(batters %>% filter(yearID == x)  %>% 
                       select(WAR) %>% arrange(-WAR))[3]))

d <- cbind(years,a)
## WAR of pitchers with unweighted population ##

WAR_pitchers_median_uw <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "median")))), npop = pops))
})) %>% arrange(-WAR_talent)

WAR_pitchers_mean_uw <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "mean")))), npop = pops))
})) %>% arrange(-WAR_talent)


WAR_pitchers_top_uw <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "top")))), npop = pops))
})) %>% arrange(-WAR_talent)

y <- as.matrix(pitchers %>% filter(yearID == 1913)  %>% select(WAR) %>% arrange(WAR) )

write.csv(WAR_pitchers_median_uw, "WAR_pitchers_median_uw.csv")
write.csv(WAR_pitchers_mean_uw, "WAR_pitchers_mean_uw.csv")
write.csv(WAR_pitchers_top_uw, "WAR_pitchers_test_uw.csv")

## WAR of pitchers with weighted population ##

WAR_pitchers_median_w <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers_weighted %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "median")))), npop = pops))
})) %>% arrange(-WAR_talent)

WAR_pitchers_mean_w <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers_weighted %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "mean")))), npop = pops))
})) %>% arrange(-WAR_talent)


WAR_pitchers_top_w <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers_weighted %>% filter(yearID == xx)  %>% arrange(WAR) %>% mutate(WAR_talent = 
  order_Pareto_vec(u = order_bino_vec(unlist(lapply(WAR, function(xx) 
  Ftilde(y = WAR, t = xx, thres = "top")))), npop = pops))
})) %>% arrange(-WAR_talent)

y <- as.matrix(pitchers_weighted %>% filter(yearID == 1913)  %>% select(WAR) %>% arrange(WAR) )

write.csv(WAR_pitchers_median_w, "WAR_pitchers_median_w.csv")
write.csv(WAR_pitchers_mean_w, "WAR_pitchers_mean_w.csv")
write.csv(WAR_pitchers_top_w, "WAR_pitchers_test_w.csv")



