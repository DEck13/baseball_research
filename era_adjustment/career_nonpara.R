## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readr)
`%notin%` <- Negate(`%in%`)
numCores <- detectCores() - 1
ncores <- detectCores() - 1
#############################################
## Get adjusted career batting averages
#############################################

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/batters_WAR")

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


# careers starting in 1996
foo <- read_csv("WAR_nonpara.csv")
foo$playerID <- droplevels(as.factor(foo$playerID))
bar <- split(foo, f = foo$playerID)
year_start = 1996; year_finish = 2019

# career averages for AB > 3000
index_kAB <- which(unlist(lapply(bar, function(xx){
  ifelse(sum(xx$AB) >= 3e3,1,0)
})) == 1)

# split by players
talent_kAB <- do.call(rbind, lapply(index_kAB, function(j){
  arrange(bar[[j]], yearID)
}))
talent_kAB$playerID <- droplevels(talent_kAB$playerID)

k <- 4
# compute career averages at reference span
career_talent <- function(snippet, year_start = 1996, year_finish = 2019){
  span <- year_start:year_finish
  if(nrow(snippet) < length(span)) span <- span[1:nrow(snippet)]
  snippet <- snippet[1:length(span), ]
  snippet <- snippet %>% mutate(playerID = paste(playerID, "_proj", sep = ""))
  do.call(rbind, lapply(span, function(xx){
    index <- which(span == xx)
    batters_int <- foo %>% filter(yearID == span[index])
    
    ## Ordinary Ftilde
    yy <- sort(batters_int$scale_WAR)
    n <- length(yy)
    ytilde <- rep(0, n + 1)
    ytilde[1] <- yy[1] - 1/(yy[2] - yy[1])
    ytilde[n+1] <- yy[n] + 1/(yy[n] - yy[n-1])
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (yy[j]+yy[j-1])/2 
    }))
    
    batters_int <- rbind(batters_int, snippet[index, ])
    batters_int$pops[nrow(batters_int)] <- batters_int$pops[1]
    batters_int <- batters_int %>% arrange(WAR_talent) %>% 
      mutate(adj_WAR = order_qempirical(map_Pareto_vals_vec(WAR_talent, npop = pops), 
                                        ytilde = ytilde)/k) %>% 
      filter(playerID == unique(snippet$playerID))
    
    batters_int
  })) %>% mutate(span = span)
}

# get career adjustment on season by season basis
career_kAB <- do.call(rbind, mclapply(unique(talent_kAB$playerID), 
function(xx){
int <- career_talent(talent_kAB %>% filter(playerID == xx), 
year_start = year_start, year_finish = year_finish) 
int
}, mc.cores = ncores))


# compute totals
career_totals <- do.call(rbind, mclapply(
  split(career_kAB, f = droplevels(as.factor(career_kAB$playerID))), 
  mc.cores = ncores, FUN = function(xx){
    data.frame(unique(xx$playerID), min(xx$yearID), 
               sum(xx$adj_WAR))
  }))

colnames(career_totals) <- c("name", "rookie_year", "WAR_total")
write.csv((career_totals %>% arrange(-WAR_total)), "career_nonpara.csv")

