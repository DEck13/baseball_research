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

setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/Reuslt/res_nonpara")

# careers starting in 1996
foo <- read_csv("WAR_nonpara_adj.csv")
foo$playerID <- droplevels(as.factor(foo$playerID))
bar <- split(foo, f = foo$playerID)
year_start = 1996; year_finish = 2019

# career averages for AB > 3000
index_kAB <- which(unlist(lapply(bar, function(xx){
  ifelse(sum(xx$AB) >= 5e3,1,0)
})) == 1)

# split by players
talent_kAB <- do.call(rbind, lapply(index_kAB, function(j){
  arrange(bar[[j]], yearID)
}))
talent_kAB$playerID <- droplevels(talent_kAB$playerID)

AB_thresh <- 320
k <- 6
# compute career averages at reference span
career_talent <- function(snippet, year_start = 1996, year_finish = 2019){
  span <- year_start:year_finish
  if(nrow(snippet) < length(span)) span <- span[1:nrow(snippet)]
  snippet <- snippet[1:length(span), ]
  snippet <- snippet %>% mutate(playerID = paste(playerID, "_proj", sep = ""))
  do.call(rbind, lapply(span, function(xx){
    index <- which(span == xx)
    batters_int <- foo %>% filter(yearID == span[index], AB >= AB_thresh)
    min_int <- min(batters_int$WAR)
    
    ## Adjusted Ordinary Ftilde
    y <- sort(batters_int$scale_WAR)
    n <- length(y)
    rank_d <- rank(diff(y))
    ytilde_adj <- rep(0, n + 1)
    ytilde_adj[1] <- y[1] - 1/(y[2] - y[1])
    ytilde_adj[n+1] <- y[n] + 1/(y[n] - y[n-1])
    ytilde_adj[2:n] <- unlist(lapply(2:n, function(j){
      y[j] * (1 - rank_d[j-1]/ (n-1)) + y[j-1] * rank_d[j-1] / (n-1)
    }))
    
    batters_int <- rbind(batters_int, snippet[index, ])
    batters_int$pops[nrow(batters_int)] <- batters_int$pops[1]
    batters_int <- batters_int %>% arrange(WAR_talent) %>% 
      mutate(adj_WAR = order_qempirical(map_Pareto_vals_vec(WAR_talent, npop = pops), 
                                        ytilde = ytilde_adj)/k) %>% 
      filter(playerID == unique(snippet$playerID))
    batters_int$adj_WAR <- ifelse(batters_int$adj_WAR < min_int, 
                                  min_int, batters_int$adj_WAR)
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

## check mike trout 

m <- (career_totals %>% arrange(-WAR_total))[1:100,]

write.csv((career_totals %>% arrange(-WAR_total))[1:25, ], "career_nonpara_adj.csv")

