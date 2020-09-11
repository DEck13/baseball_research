setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/pitching")

## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readxl)
library(readr)
years <- 1880:2019
ncores <- detectCores() - 1

## load in local csv file
pitchers <- read_csv('pitchers_adj.csv')
pitchers_weighted <- read_csv('pitchers_adj_weighted.csv')

## Ipouts: Outs Pitched
# Ipout_thresh <- 375

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
  m <- round(n * 2 /3)
  ytilde[2:m] <- unlist(lapply(2:m, function(j){
    (y[j]+y[j-1])/2 
  }))
  ytilde[(m+1):n] <- unlist(lapply((m+1):n, function(j){
    y[j] * (1 - rank_d[j-1]/ (n-1)) + y[j-1] * rank_d[j-1] / (n-1)
  }))
  
  j <- length(which(ytilde <= t))
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
# xx = 1881
k <- 7

WAR_leader_board <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers_weighted %>% filter(yearID == xx)  %>% 
    arrange(WAR) %>% mutate(scale_WAR = WAR*k) %>% mutate(
      WAR_talent = order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_WAR, function(xx) 
        Ftilde(y = scale_WAR, t = xx)))), npop = pops))  
})) %>% arrange(-WAR_talent) 

foo <- WAR_leader_board
foo$playerID <- droplevels(as.factor(foo$playerID))
bar <- split(foo, f = foo$playerID)
baz <- (do.call(rbind, mclapply(bar, mc.cores = ncores, function(xx){
  xx[which.max(xx$WAR_talent), ]
})) %>% arrange(-WAR_talent))

nrow(baz[1:100,] %>% filter(yearID < 1950))

WAR_leader_board_adj <- do.call(rbind, mclapply(years, mc.cores = ncores, function(xx){
  pitchers_weighted %>% filter(yearID == xx)  %>% 
    arrange(WAR) %>% mutate(scale_WAR = WAR*k) %>% mutate(
      WAR_talent = order_Pareto_vec(u = order_bino_vec(unlist(lapply(scale_WAR, function(xx) 
        Ftilde_adj(y = scale_WAR, t = xx)))), npop = pops))  
})) %>% arrange(-WAR_talent) 


foo_adj <- WAR_leader_board_adj
foo_adj$playerID <- droplevels(as.factor(foo_adj$playerID))
bar_adj <- split(foo_adj, f = foo_adj$playerID)
baz_adj <- (do.call(rbind, mclapply(bar_adj, mc.cores = ncores, function(xx){
  xx[which.max(xx$WAR_talent), ]
})) %>% arrange(-WAR_talent))

## clean the code 
## try different methods
setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/pitching/interpolation/weighted")
write.csv(baz, 'pitchers_WAR_nonpara_weighted_single.csv')
write.csv(baz_adj, 'pitchers_WAR_nonpara_adj_single.csv')

write.csv(WAR_leader_board, 'pitchers_WAR_nonpara_weighted.csv')
write.csv(WAR_leader_board_adj, 'pitchers_WAR_nonpara_adj_weighted.csv')




