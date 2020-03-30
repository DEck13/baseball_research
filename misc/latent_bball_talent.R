
## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)


set.seed(13)
n <- 500
Z <- sort(rnorm(n))
Z[n] <- 1.5*Z[n]

## order_pnorm function
# converts normal order stats to their percentiles
order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
  p <- pnorm(q, mean = mean, sd = sd)
  pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
}

## order_pnorm_vec function
# this function converts a vector of normal order stats 
# to their percentiles. This vector should be the entire 
# sample sorted in increasing order
order_pnorm_vec <- function(q, mean = 0, sd = 1){
  q <- sort(q) # just in case
  n <- length(q)
  unlist(lapply(1:n, function(j){
    order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
  }))
}

## order_pnorm_vec function
# new implementation taken from Order Stat Notes
#order_pnorm_vec <- function(q, mean = 0, sd = 1){
#  q <- sort(q) # just in case
#  pnorm(q, mean = mean, sd = sd)
#}

# try it out
u <- order_pnorm_vec(Z)


## order_Pareto_vec function 
# this function transform percentiles from order stats (in increasing order)
# to Pareto values corresponding to the general population 
# of a greater than or equal to size
order_Pareto_vec <- function(u, t = 1, alpha = 1.5, npop = 1e4){
  n <- length(u)
  if(length(npop) == 1) npop <- rep(npop, n)
  unlist(lapply(1:n, function(j){
    qPareto(qbeta(u[j], j + npop[j], n + 1 - j), t = t, alpha = alpha)
  }))
}

# try it out on the percentiles obtained from order_norm_vec
X <- order_Pareto_vec(u, npop = 1e6)



## map_Pareto_vals_vec function 
# this function transforms ordered Pareto values corresponding to 
# the general population to percentiles from order stats 
# (in increasing order)
map_Pareto_vals_vec <- function(x, t = 1, alpha = 1.5, npop = 1e4){
  n <- length(x)
  if(length(npop) == 1) npop <- rep(npop, n)  
  unlist(lapply(1:n, function(j){
    pbeta(pPareto(x[j], t = t, alpha = alpha), j + npop[j], n + 1 - j)
  }))
}

# try it out on the percentiles obtained from order_Pareto_vec
u2 <- map_Pareto_vals_vec(X, npop = 1e6)
sqrt(crossprod(u2 - u)) # check


## oder_qnorm function
# take the ordered percentiles and convert them to order statistics 
# from a normal distribution
order_qnorm <- function(u, mean = 0, sd = 1){
  n <- length(u)
  qnorm(qbeta(u, shape1 = 1:n, shape2 = n:1), mean = mean, sd = sd)
}


# try it out on the percentiles obtained from map_Pareto_vals_vec
Z2 <- order_qnorm(u2)
sqrt(crossprod(Z - Z2)) # check
#cbind(Z,Z2)[which.max(abs(Z - Z2)) + c(-7:7), ]
#tail(cbind(Z,Z2), 10)




## add fielder information to create positions.
fielders <- read.csv("~/research/baseball/baseballdatabank-2019.2/core/Fielding.csv", 
                     header = TRUE) #%>% select(playerID, yearID, POS)
fielders$InnOuts[is.na(fielders$InnOuts)] <- 1
fielders_by_player <- split(fielders, as.factor(fielders$playerID))
ncores <- detectCores() - 1
unique_POS <- do.call(rbind, mclapply(fielders_by_player, mc.cores = ncores, 
  function(xx){
    #print(xx$playerID)
    do.call(rbind, lapply(split(xx, f = droplevels(as.factor(xx$yearID))), 
         function(yy){
           if(dim(yy)[1] >= 2) yy$POS <- rep(yy$POS[which.max(yy$InnOuts)], dim(yy)[1])
           data.frame(playerID = unique(yy$playerID), 
                      yearID = unique(yy$yearID), 
                      POS = unique(yy$POS))
         }))
}))


## load in batters information, incorporate positions and park factors
batters <- read.csv("~/research/baseball/baseballdatabank-2019.2/core/Batting.csv", 
                    header = TRUE)
batters <- merge(batters, unique_POS)

## years to consider
years <- 1880:2018

## load in teams info and obtain park factors.
# from 1880-1903: we will use runs park factors to adjust for batting averages
#                 that are obtained in different ball parks
Teams <- read.csv("~/research/baseball/baseballdatabank-2019.2/core/Teams.csv", 
                    header = TRUE)
BPF <- Teams %>% select(teamID, yearID, BPF)
batters <- merge(batters, BPF)

## park effect adjusted batting average and OPS
# 2/3 weight to league park factor; 1/2 weight to neutral park
#batters <- batters %>% mutate(AVG = H/AB * (2/3 * 100/BPF + 1/3) ) %>% 
batters <- batters %>% filter(AB > 1) %>% 
  mutate(AVG = H/AB * 100/BPF) %>% 
  mutate(OPS = round( (((BB + H) / (AB + BB)) +
    (((H - X2B - X3B- HR) + 2*X2B + 3*X3B + 4*HR) / AB)) * 100/BPF, 3)) %>% 
  select(playerID, POS, yearID, teamID, lgID, AVG, OPS, HR, AB) %>% 
  filter(yearID %in% years)
batters$playerID <- droplevels(batters$playerID)
  
## remove multiple stints/combine intraseason stats after trade
foo <- split(batters, f = as.factor(batters$playerID))
batters <- do.call(rbind, mclapply(foo, mc.cores = ncores,  function(bar){
  stint_int <- sort(unique(bar$yearID))[which(table(bar$yearID) > 1)]
  print(stint_int); print(any(stint_int))
  if(any(stint_int)){
    for(j in stint_int){
      xx <- bar[bar$yearID %in% j, ]
      yy <- xx[1, ]
      yy$teamID <- xx$teamID[which.max(xx$AB)]
      yy$lgID <- xx$lgID[which.max(xx$AB)]
      yy$HR <- sum(xx$HR)
      yy$AVG <- sum(xx$AVG * xx$AB / sum(xx$AB))
      yy$OPS <- sum(xx$OPS * xx$AB / sum(xx$AB))
      yy$AB <- sum(xx$AB)
      bar <- filter(bar, yearID != j)
      bar <- rbind(bar, yy) %>% arrange(yearID)
    }
  }
  bar
}))


## Atbat threshold
AB_thresh <- 320 # rougly 2 per game (same as Gould)

## decenial population counts from Eck (2020)
pops <- c(4.40, 5.010, 5.580, 8.550, 8.930, 9.920, 11.130, 11.590 + 1, # 1950
          19.580, 26.100, 36.370, 40.680, 64.060, 76.190, 98) * 1e6

## interpolated population counts
pops <- approx(x = 1880 + 0:14*10, y = pops, xout = years)$y


## decenial weights (correction for AL/NL differences in integration)
weightsNL <- c(0.15, 0.15, 0.18, 0.20, 0.22, 0.26, 0.24, 0.23, #1947 
               0.28, 0.28, 0.28, #1956, 1960, 1970
               0.26, 0.20, 0.17, 0.14, 0.11) #1980-
weightsAL <- c(0.15, 0.15, 0.18, 0.20, 0.22, 0.26, 0.24, 0.23, #1947 
               0.24, 0.24, 0.25, #1956, 1960, 1970
               0.25, 0.20, 0.17, 0.14, 0.11) #1980-

## interpolated population counts
weightsNL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsNL, xout = years)$y
weightsAL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsAL, xout = years)$y
w_popsNL <- pops * weightsNL
w_popsAL <- pops * weightsAL

## add weighted populations to data
pops_data <- rbind(
data.frame(pops = w_popsNL, lgID = "NL", yearID = 1880 + 0:(14*10 - 2)),
data.frame(pops = w_popsAL, lgID = "AL", yearID = 1880 + 0:(14*10 - 2)))
batters <- merge(batters, pops_data, by = c("lgID", "yearID"))


## map normal batting averages to Pareto space 
# cap talent below max talent for AB < AB_thresh
AVG_talent <- do.call(rbind, mclapply(years, mc.cores = ncores, function(yy){
  foo <- batters %>% filter(yearID == yy) %>% arrange(AVG)
  bar <- foo %>% filter(AB >= AB_thresh)
  nsys <- nrow(bar)
  mu_AVG <- mean(bar$AVG)
  sd_AVG <- sd(bar$AVG)
  AVG_sorted <- bar$AVG
  pops <- bar$pops
  u_int <- order_pnorm_vec(AVG_sorted, mean = mu_AVG, sd = sd_AVG)
  bar <- bar %>% mutate(talent = order_Pareto_vec(u_int, npop = pops))
  max_talent <- max(bar$talent) - 1
  range <- which(!(foo$playerID %in% bar$playerID))
  bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
    #pops_full <- baz$pops
    #u_int_full <- order_pnorm_vec(baz$AVG, mean = mu_AVG, sd = sd_AVG)
    #baz %>% mutate(talent = order_Pareto_vec(u_int_full, npop = pops_full)) %>% 
    #  filter(AB < AB_thresh) %>% 
    #  mutate(talent = ifelse(talent > max_talent+1, max_talent, talent))
    rbind(bar %>% select(-talent), foo[j, ]) %>% arrange(AVG) %>% 
      mutate(talent = order_Pareto_vec(order_pnorm_vec(AVG, mean = mu_AVG, sd = sd_AVG), 
                                     npop = pops)) %>% 
      filter(AB < AB_thresh) %>% 
      mutate(talent = ifelse(talent > max_talent+1, max_talent, talent))
  }))) %>% arrange(talent)
}))

## map normal batting averages to Pareto space 
#talent_eval <- function(yy){
#  batters_int <- batters %>% filter(yearID == yy & AB >= AB_thresh) %>% 
#    arrange(AVG)
#  nsys <- nrow(batters_int)
#  mu_AVG <- mean(batters_int$AVG)
#  sd_AVG <- sd(batters_int$AVG)
#  AVG_sorted <- batters_int$AVG
#  u_int <- order_pnorm_vec(AVG_sorted, mean = mu_AVG, sd = sd_AVG)
#  #print(cbind(yy, u_int))
#  
  ## years input and pops need to be the same length
  #years_int <- which(years == yy)
  #npop <- w_pops[years_int]
#  batters_int <- batters_int %>% mutate(talent = 
#    order_Pareto_vec(u_int, npop = pops)
#  ) %>% select(yearID, playerID, pops, AB, AVG, talent)
#  
#  batters_int #%>% filter(talent > 1e7) %>% mutate(talent = talent / 1e7)
#}

## compute talent scores for all qualifying seasons
#AVG_talent <- do.call(rbind, lapply(years, talent_eval))

## top 25 talent-scores
top25seasons <- (AVG_talent %>% filter(AB >= AB_thresh) %>% arrange(-talent))[1:25, ]
nrow((AVG_talent %>% filter(AB >= AB_thresh) %>% 
        arrange(-talent))[1:25, ] %>% filter(yearID < 1947))


## get peak performance with respect to 1996
ref_1996 <- AVG_talent %>% filter(yearID == 1996) %>% filter(AB >= AB_thresh)
mean_AVG_1996 <- mean(ref_1996$AVG)
sd_AVG_1996 <- sd(ref_1996$AVG)
w_pops_1996 <- w_popsNL[which(years == 1996)]

## restrict attention to players before 1997 to compare
## with the era-bridging paper and Gould
foo <- AVG_talent %>% filter(AB >= AB_thresh, yearID < 1997)
foo$playerID <- droplevels(foo$playerID)
bar <- split(foo, f = foo$playerID)
baz <- (do.call(rbind, mclapply(bar, mc.cores = ncores, function(xx){
  xx[which.max(xx$talent), ]
})) %>% arrange(-talent))[1:25, ]

# top individual seasons in 1996 for players who came before 1997
top25seasons_bf1996 <- do.call(rbind, mclapply(1:25, mc.cores = ncores, function(j){
  xx <- baz[j, ]
  xx$pops <- w_pops_1996
  xx$playerID <- paste(xx$playerID, "_proj", sep = "")
  int <- rbind(ref_1996, xx) %>% arrange(talent) %>% 
          mutate(adj_AVG = order_qnorm(map_Pareto_vals_vec(
            talent, npop = pops), 
            mean = mean_AVG_1996, sd = sd_AVG_1996))
  int[grepl("_proj", int$playerID), ]
})) %>% select(-AVG,-pops)
top25seasons_bf1996

# top individual seasons in 1996 period
foo <- AVG_talent %>% filter(AB >= AB_thresh)
foo$playerID <- droplevels(foo$playerID)
bar <- split(foo, f = foo$playerID)
baz <- (do.call(rbind, mclapply(bar, mc.cores = ncores, function(xx){
  xx[which.max(xx$talent), ]
})) %>% arrange(-talent))[1:25, ]
top25seasons_1996 <- do.call(rbind, mclapply(1:25, mc.cores = ncores, function(j){
  xx <- baz[j, ]
  xx$pops <- w_pops_1996
  xx$playerID <- paste(xx$playerID, "_proj", sep = "")
  int <- rbind(ref_1996, xx) %>% arrange(talent) %>% 
    mutate(adj_AVG = order_qnorm(map_Pareto_vals_vec(
      talent, npop = pops), 
      mean = mean_AVG_1996, sd = sd_AVG_1996))
  int[grepl("_proj", int$playerID), ]
})) %>% select(-AVG,-pops)
top25seasons_1996


## career averages for AB > 3000 (from a desired baseline)
## AB > 3000 is chosen to be nearly consistent with Brefs list
#foo <- batters %>% filter(POS != "P") 
foo <- AVG_talent #%>% filter(AB >= AB_thresh)
foo$playerID <- droplevels(foo$playerID)
bar <- split(foo, f = foo$playerID)
year_start = 1996; year_finish = 2018

# may check 5000 AB
index_kAB <- which(unlist(lapply(bar, function(xx){
  ifelse(sum(xx$AB) >= 3e3,1,0)
})) == 1)
AVG_talent_kAB <- do.call(rbind, lapply(index_kAB, function(j){
  arrange(bar[[j]], yearID)
}))
AVG_talent_kAB$playerID <- droplevels(AVG_talent_kAB$playerID)


# compute career averages at reference span
career_AVG <- function(snippet, year_start = 1996, year_finish = 2018){
  span <- year_start:year_finish
  if(nrow(snippet) < length(span)) span <- span[1:nrow(snippet)]
  snippet <- snippet[1:length(span), ]
  snippet <- snippet %>% mutate(playerID = paste(playerID, "_proj", sep = ""))
  do.call(rbind, lapply(span, function(yy){
    index <- which(span == yy)
    batters_int <- AVG_talent %>% filter(yearID == span[index], AB >= AB_thresh)
    mean_int <- mean(batters_int$AVG)
    sd_int <- sd(batters_int$AVG)
    #pops_int <- batters_int %>% select(pops)
    batters_int <- rbind(batters_int, snippet[index, ])
    batters_int$pops[nrow(batters_int)] <- batters_int$pops[1]
    batters_int <- batters_int %>% arrange(talent) %>% 
      mutate(adj_AVG = order_qnorm(map_Pareto_vals_vec(
        talent, npop = pops), 
        mean = mean_int, sd = sd_int)) %>% 
      filter(playerID == unique(snippet$playerID))
  })) %>% mutate(span = span)
  
}

career_AVG_kAB <- do.call(rbind, mclapply(unique(AVG_talent_kAB$playerID), 
  function(xx){
    int <- career_AVG(AVG_talent_kAB %>% filter(playerID == xx), 
      year_start = year_start, year_finish = year_finish)%>% 
      #mutate(adj_H = adj_AVG * AB) %>% 
      select(-pops, -talent, -AVG)
  #data.frame(xx, sum(int$adj_H) / sum(int$AB), sum(int$adj_H))
    int
}, mc.cores = ncores))

career_total_AVG <- do.call(rbind, mclapply(
  split(career_AVG_kAB, f = droplevels(as.factor(career_AVG_kAB$playerID))), 
  mc.cores = ncores, FUN = function(xx){
    xx <- xx %>% mutate(adj_H = adj_AVG * AB)
    data.frame(unique(xx$playerID), min(xx$yearID), 
               sum(xx$adj_H) / sum(xx$AB), sum(xx$adj_H))
  }))
colnames(career_total_AVG) <- c("name", "rookie year", "average", "hits")
(career_total_AVG %>% arrange(-hits))[1:25, ]
(career_total_AVG %>% arrange(-average))[1:25, ]



## select old time players
career_AVG_kAB %>% filter(name == "willite01")
career_AVG_kAB %>% filter(name == "hornsro01")
career_AVG_kAB %>% filter(name == "collied01")
career_AVG_kAB %>% filter(name == "wagneho01")
career_AVG_kAB %>% filter(name == "gehrilo01")
career_AVG_kAB %>% filter(name == "ruthba01")
career_AVG_kAB %>% filter(name == "delahed01")
career_AVG_kAB %>% filter(name == "lajoina01")
career_AVG_kAB %>% filter(name == "ansonca01")
career_AVG_kAB %>% filter(name == "oneilti01")
career_AVG_kAB %>% filter(name == "duffyhu01")

## Keeler's career if starting in 1996
career_AVG_kAB %>% filter(name == "keelewi01")



## C
career_AVG_kAB %>% filter(name == "poseybu01")
career_AVG_kAB %>% filter(name == "piazzmi01")
career_AVG_kAB %>% filter(name == "mauerjo01")
career_AVG_kAB %>% filter(name == "rodriiv01")
career_AVG_kAB %>% filter(name == "molinya01")
career_AVG_kAB %>% filter(name == "benchjo01")
career_AVG_kAB %>% filter(name == "cartega01")

## 1B
career_AVG_kAB %>% filter(name == "cabremi01") # or 3B
career_AVG_kAB %>% filter(name == "pujolal01")
career_AVG_kAB %>% filter(name == "martied01") # DH
career_AVG_kAB %>% filter(name == "vottojo01")
career_AVG_kAB %>% filter(name == "thomafr04")
career_AVG_kAB %>% filter(name == "musiast01") # or LF
career_AVG_kAB %>% filter(name == "bagweje01")
career_AVG_kAB %>% filter(name == "gehrilo01")
career_AVG_kAB %>% filter(name == "thomeji01")

## 2B
career_AVG_kAB %>% filter(name == "carewro01")
career_AVG_kAB %>% filter(name == "hornsro01")
career_AVG_kAB %>% filter(name == "molitpa01")
career_AVG_kAB %>% filter(name == "rosepe01")
career_AVG_kAB %>% filter(name == "alomaro01")
career_AVG_kAB %>% filter(name == "biggicr01")
career_AVG_kAB %>% filter(name == "morgajo02")
career_AVG_kAB %>% filter(name == "sandbry01")


## 3B
career_AVG_kAB %>% filter(name == "cabremi01") # or 1B
career_AVG_kAB %>% filter(name == "molitpa01")
career_AVG_kAB %>% filter(name == "brettge01")
career_AVG_kAB %>% filter(name == "rosepe01")
career_AVG_kAB %>% filter(name == "jonesch06")
career_AVG_kAB %>% filter(name == "beltrad01")
career_AVG_kAB %>% filter(name == "schmimi01")

## SS
career_AVG_kAB %>% filter(name == "jeterde01")
career_AVG_kAB %>% filter(name == "garcino01")
career_AVG_kAB %>% filter(name == "trammal01")
career_AVG_kAB %>% filter(name == "ripkeca01")
career_AVG_kAB %>% filter(name == "smithoz01")
career_AVG_kAB %>% filter(name == "wagneho01")
#career_AVG_kAB %>% filter(name == "yountro01")


## LF
career_AVG_kAB %>% filter(name == "musiast01") # or 1B
career_AVG_kAB %>% filter(name == "willite01")
career_AVG_kAB %>% filter(name == "bondsba01")
career_AVG_kAB %>% filter(name == "henderi01")

## CF
career_AVG_kAB %>% filter(name == "troutmi01")
career_AVG_kAB %>% filter(name == "cobbty01")
career_AVG_kAB %>% filter(name == "bettsmo01") # not too many PAs
career_AVG_kAB %>% filter(name == "mantlmi01")
career_AVG_kAB %>% filter(name == "mayswi01")
career_AVG_kAB %>% filter(name == "griffke01")
career_AVG_kAB %>% filter(name == "beltrca01")
career_AVG_kAB %>% filter(name == "jonesan01")

## RF
career_AVG_kAB %>% filter(name == "suzukic01")
career_AVG_kAB %>% filter(name == "clemero01")
career_AVG_kAB %>% filter(name == "guerrvl01")
career_AVG_kAB %>% filter(name == "aaronha01")
career_AVG_kAB %>% filter(name == "robinfr02")
career_AVG_kAB %>% filter(name == "ruthba01")
career_AVG_kAB %>% filter(name == "walkela01")








test <- batters %>% filter(yearID == 1996) %>% 
  filter(AB > AB_thresh) %>% select(playerID, AB, HR) %>% 
  mutate(HRpAB = HR/AB) %>% arrange(-HRpAB)
head(test)
hist(test$HRpAB)
ecdf_HR <- ecdf(test$HRpAB)
ecdf_HR(test$HRpAB)

plot(ecdf_HR)





