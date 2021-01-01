## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readr)
numCores <- detectCores() - 1
setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/summary_interpolation")

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

predicted_WAR <- function(player1, year1, year2, thres, pos){
  if (pos == "batter") {
    if (thres == "median") {
      WAR_nonpara = read_csv("WAR_batter_median.csv")
    }
    if (thres == "mean") {
      WAR_nonpara = read_csv("WAR_batter_mean.csv")
    }
    if (thres == "top") {
      WAR_nonpara = read_csv("WAR_batter_test.csv")
    }
  }
  if (pos == "pitcher_uw") {
    if (thres == "median") {
      WAR_nonpara = read_csv("WAR_pitchers_median_uw.csv")
    }
    if (thres == "mean") {
      WAR_nonpara = read_csv("WAR_pitchers_mean_uw.csv")
    }
    if (thres == "top") {
      WAR_nonpara = read_csv("WAR_pitchers_test_uw.csv")
    }
  }
  if (pos == "pitcher_w") {
    if (thres == "median") {
      WAR_nonpara = read_csv("WAR_pitchers_median_w.csv")
    }
    if (thres == "mean") {
      WAR_nonpara = read_csv("WAR_pitchers_mean_w.csv")
    }
    if (thres == "top") {
      WAR_nonpara = read_csv("WAR_pitchers_test_w.csv")
    }
  }
  
  
  m <- WAR_nonpara %>% filter(yearID == year1) %>% filter(playerID == player1)

  if(nrow(m) == 0 ){
    stop("This player has never played in this season")
  } 
  
  print(paste("WAR of", player1, "in",year1, "season","is" ,m$WAR, sep = " "))
  
  nonpara_tal <- m$WAR_talent

  ref <- WAR_nonpara %>%  filter(yearID == year2)
  w_pops <- mean(ref$pops)
  
  yy <- sort(ref$WAR)
  n <- length(yy)
  ytilde <- rep(0, n + 1)
  
  if (thres == "median") {
    ytilde[1] <- yy[1] - 1/(median(yy) - yy[1])
    ytilde[n+1] <- yy[n] + 1/(yy[n] - median(yy))
  }
  if (thres == "mean") {
    ytilde[1] <- yy[1] - 1/(mean(yy) - yy[1])
    ytilde[n+1] <- yy[n] + 1/(yy[n] - mean(yy))
  }
  if (thres == "top"){
    ytilde[1] <- yy[1] - 1/(yy[3] - yy[1])
    ytilde[n+1] <- yy[n] + 1/(yy[n] - yy[n-2])
  }
  
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (yy[j]+yy[j-1])/2 
  }))
  
  m$playerID <- paste(m$playerID, "_proj", sep = "")
  m$pops <- w_pops
  int <- rbind(ref, m) %>% arrange(WAR_talent) %>% 
    mutate(adj_WAR = order_qempirical(map_Pareto_vals_vec(WAR_talent, npop = pops), 
                                      ytilde = ytilde))
  a <- int[grepl("_proj", int$playerID), ]
  
  return(a$adj_WAR)
}

predicted_WAR("ripkeca01", 1991, 1996, thres = "top", pos = "batter")
