setwd("~/Desktop/PhD_UIUC/ProfEck/baseball/pitching")

## average the population ## 

## load in software
rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(readxl)
library(readr)
library(dplyr)
`%notin%` <- Negate(`%in%`)
ncores <- detectCores() - 1
## load in local csv file
pitchers_war <- read_excel("bref_pitchers_war.xlsx")
load('Pitching.rdata')
colnames(pitchers_war)[4] <- "playerID"
colnames(pitchers_war)[5] <- "yearID"
colnames(pitchers_war)[8] <- "lgID"

# m <- Pitching %>% filter(playerID == "santajo01")
# n <- pitchers_war %>% filter(playerID == "santajo02")

table(pitchers_war$lgID)
table(Pitching$lgID)

## check the error in pitcher_war ##

ID_1 <- unique(Pitching$playerID)
ID_2 <- unique(pitchers_war$playerID)
# which(b == 'santajo02')

a <- ID_1[which(ID_1 %notin% ID_2)]
b <- ID_2[which(ID_2 %notin% ID_1)]

s1 = sapply(1:length(a), function(x) unlist(strsplit(a[x], split='0', fixed=TRUE))[1])
s2 = sapply(1:length(b), function(x) unlist(strsplit(b[x], split='0', fixed=TRUE))[1])
# a[which(s1 %in% s2)]
# b[which(s2 %in% s1)]

## rodrijo02 is an error # 11
m <- Pitching %>% filter(playerID == a[which(s1 %in% s2)][k])
n <- pitchers_war %>% filter(playerID == b[which(s2 == s1[which(s1 %in% s2)][k])])
m$yearID == n$yearID

for (k in c(seq(1:20)[-11])) {
  n <- nrow(pitchers_war %>% filter(playerID == 
  b[which(s2 == s1[which(s1 %in% s2)][k])]))
  pitchers_war[pitchers_war$playerID == b[which(s2 == 
  s1[which(s1 %in% s2)][k])],]$playerID = rep(a[which(s1 %in% s2)][k], n)
}

# m <- Pitching %>% filter(playerID == "santajo01")
# n <- pitchers_war %>% filter(playerID == "santajo01")


# ID_1 <- unique(Pitching$playerID)
# ID_2 <- unique(pitchers_war$playerID)

## select players who are from NL and AL 

Pitching <- Pitching %>% 
  select('playerID', 'yearID', 'IPouts', 'lgID') %>% 
  filter(yearID >= 1880) %>% 
  filter(lgID %in% c('NL', 'AL'))
pitchers_war <- pitchers_war %>% 
  select('playerID', 'yearID', 'IPouts', 'lgID', 'WAR') %>% 
  filter(yearID >= 1880) %>% 
  filter(lgID %in% c('NL', 'AL'))

## unweighted population ##

years <- 1880:2019
ncores <- detectCores() - 1

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

## weighted population ##

years <- 1880:2019
ncores <- detectCores() - 1

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
weightsNL <- c(0.12, 0.13, 0.14, 0.15, 0.20, 0.23, 0.25, 0.20,
            0.35, 0.35, 0.33,
             0.31, 0.25, 0.22, 0.20, 0.18)
weightsAL <- c(0.12, 0.13, 0.14, 0.15, 0.20, 0.23, 0.25, 0.20,
              0.25, 0.25, 0.29,
              0.30, 0.25, 0.22, 0.20, 0.18)
# weightsNL <- rep(1, 16)
# weightsAL <- rep(1, 16)

## interpolated population counts
weightsNL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsNL, xout = years)$y
weightsAL <- approx(x = c(1880 + 0:6*10, 1947, 1956, 1960 + 0:6 * 10), 
                    y = weightsAL, xout = years)$y
w_popsNL <- pops * weightsNL
w_popsAL <- pops * weightsAL

## add weighted populations to data
pops_data_weighted <- rbind(
  data.frame(pops = w_popsNL, lgID = "NL", yearID = 1880 + 0:(14*10 - 1)),
  data.frame(pops = w_popsAL, lgID = "AL", yearID = 1880 + 0:(14*10 - 1)))



## combine stats for players were traded ##

Pitching_adj <- do.call(rbind, mclapply(split(Pitching, f = droplevels(as.factor(Pitching$playerID))), 
mc.cores = ncores, FUN = function(xx){
n <- unique(xx$yearID)
IPouts <- sapply(n, function(yy) sum(xx %>% filter(yearID == yy) %>% 
                                       select(IPouts)))
lgID <- unlist(sapply(n, function(yy) (xx %>% filter(yearID == yy) %>% 
                                     select(lgID))[1,]   ))
data.frame(rep(unique(xx$playerID), length(n)),n,IPouts, lgID)
}
))
colnames(Pitching_adj) <- c("playerID", "yearID", "IPouts", 'lgID')
# xx <- split(Pitching, f = droplevels(as.factor(Pitching$playerID)))[[3]]
# a <- as.numeric(pitchers_war$WAR)
# b <- pitchers_war[which(is.na(a) == TRUE),]

pitchers_war$WAR <- as.numeric(pitchers_war$WAR)
pitchers_war <- pitchers_war[complete.cases(pitchers_war),]


Pitchers_adj <- do.call(rbind, mclapply(split(pitchers_war, f = droplevels(as.factor(pitchers_war$playerID))), 
mc.cores = ncores, FUN = function(xx){
n <- unique(xx$yearID)
IPouts <- sapply(n, function(yy) sum(xx %>% filter(yearID == yy) %>% 
                                     select(IPouts)))
WAR <- sapply(n, function(yy) sum(xx %>% filter(yearID == yy) %>% 
                                       select(WAR)))
lgID <- unlist(sapply(n, function(yy) (xx %>% filter(yearID == yy) %>% 
                                         select(lgID))[1,]   ))
data.frame(rep(unique(xx$playerID), length(n)),n,IPouts, lgID, WAR)
}
))
colnames(Pitchers_adj) <- c("playerID", "yearID", "IPouts", "lgID" ,"WAR")

pitchers_adj <- merge(Pitching_adj, Pitchers_adj, by = c('playerID','yearID'))
# a <- merge(Pitching_adj, Pitchers_adj, by = c('playerID','yearID', 'IPouts'))
b <- pitchers_adj[pitchers_adj$lgID.x != pitchers_adj$lgID.y,]

## select the maximam IPouts 
l <- nrow(pitchers_adj)
pitchers_adj <- pitchers_adj %>% mutate(IPouts = 
sapply(1:l, function(zz) max(pitchers_adj$IPouts.x[zz], 
pitchers_adj$IPouts.y[zz]))) %>% 
  select("playerID", "yearID", "IPouts", "WAR", "lgID.x") %>% filter(yearID >= 1880)

colnames(pitchers_adj)[5] <- 'lgID'
pitchers <- merge(pitchers_adj, pops_data, by = c("lgID", "yearID"))
write_csv(pitchers, 'pitchers_adj.csv')

## weighted case ##


## take out the players that play in the different league in a single season ##

p_s1 <- do.call(rbind, mclapply(split(Pitching, f = droplevels(as.factor(Pitching$playerID))), 
mc.cores = ncores, FUN = function(xx){
n <- unique(xx$yearID)
if (sum(duplicated(xx$yearID)) != 0) xx[xx$yearID %in% xx$yearID[duplicated(xx$yearID)],]
}
))

s_s1 <- anti_join(Pitching, p_s1, by = c("playerID", "yearID", "IPouts", "lgID")) 

p_sp1 <- merge(p_s1, pops_data_weighted, by = c("lgID", "yearID")) %>% select(-lgID)
s_sp1 <- merge(s_s1, pops_data_weighted, by = c("lgID", "yearID")) %>% select(-lgID)

p_sf1 <- do.call(rbind, mclapply(split(p_sp1, f = droplevels(as.factor(p_sp1$playerID))), 
mc.cores = ncores, FUN = function(xx){
  a <- sapply(unique(xx$yearID), function(yy) sum(xx[xx$yearID == yy,]$IPouts))
  b <- sapply(unique(xx$yearID), function(yy) 
    weighted.mean(x = xx[xx$yearID == yy,]$pops, 
   w = xx[xx$yearID == yy,]$IPouts/sum(xx[xx$yearID == yy,]$IPouts)))
  data.frame(playerID = unique(xx$playerID), yearID = unique(xx$yearID), 
IPouts = a, pops = b)
  }
))

pitchers1 <- rbind(s_sp1, p_sf1)


p_s2 <- do.call(rbind, mclapply(split(pitchers_war, f = droplevels(as.factor(pitchers_war$playerID))), 
mc.cores = ncores, FUN = function(xx){
n <- unique(xx$yearID)
if (sum(duplicated(xx$yearID)) != 0) xx[xx$yearID %in% xx$yearID[duplicated(xx$yearID)],]
}
))

s_s2 <- anti_join(pitchers_war, p_s2, 
by = c("playerID", "yearID", "IPouts", "lgID", "WAR")) 

p_sp2 <- merge(p_s2, pops_data_weighted, by = c("lgID", "yearID")) %>% select(-lgID)
s_sp2 <- merge(s_s2, pops_data_weighted, by = c("lgID", "yearID")) %>% select(-lgID)

p_sf2 <- do.call(rbind, mclapply(split(p_sp2, f = droplevels(as.factor(p_sp2$playerID))), 
mc.cores = ncores, FUN = function(xx){
a <- sapply(unique(xx$yearID), function(yy) sum(xx[xx$yearID == yy,]$IPouts))
b <- sapply(unique(xx$yearID), function(yy) 
weighted.mean(x = xx[xx$yearID == yy,]$pops, 
w = xx[xx$yearID == yy,]$IPouts/sum(xx[xx$yearID == yy,]$IPouts)))
c <- sapply(unique(xx$yearID), function(yy) sum(xx[xx$yearID == yy,]$WAR))
data.frame(playerID = unique(xx$playerID), yearID = unique(xx$yearID), 
IPouts = a, WAR = c, pops = b)
}
))

pitchers2 <- rbind(s_sp2, p_sf2)

pitchers_adj_weighted <- merge(pitchers1, pitchers2, by = c('playerID','yearID'))
# a <- merge(Pitching_adj, Pitchers_adj, by = c('playerID','yearID', 'IPouts'))

# population is consistent. 
b <- pitchers_adj_weighted[pitchers_adj_weighted$pop.x != pitchers_adj_weighted$pop.y,]

## select the maximam IPouts 
l <- nrow(pitchers_adj_weighted)
pitchers_adj_weighted <- pitchers_adj_weighted %>% mutate(IPouts = 
sapply(1:l, function(zz) max(pitchers_adj_weighted$IPouts.x[zz], 
pitchers_adj_weighted$IPouts.y[zz]))) %>% 
  select("playerID", "yearID", "IPouts", "WAR", "pops.x") %>% filter(yearID >= 1880)

colnames(pitchers_adj_weighted)[5] <- 'pops'

write_csv(pitchers_adj_weighted, 'pitchers_adj_weighted.csv')

Pitching %>% filter(playerID == "abbotpa01")
pitchers_war %>% filter(playerID == "abbotpa01")

