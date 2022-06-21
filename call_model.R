library(plyr)
library(lubridate)

source('model.R')
table <- read.csv(file = "NAAMES_model_data.csv", header = T,stringsAsFactors = F)

naames_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
naames_npp <- lapply(naames_npp, npp_func_var_c_naames)

naames_npp <- data.frame(matrix(unlist(naames_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(naames_npp) <- table[,c("yd")]
naames_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(naames_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cvar <- subset(table_long, depth <= Ez1)
npp_cvar <- aggregate(.~yd,  npp_cvar[,c("yd","value")], FUN = sum)
colnames(npp_cvar) <- c("yd","Ez_NPP_cvar")

## Cmod 

naames_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
naames_npp <- lapply(naames_npp, npp_func_mod_c)

naames_npp <- data.frame(matrix(unlist(naames_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(naames_npp) <- table[,c("yd")]
naames_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(naames_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cmod <- subset(table_long, depth <= Ez1)
npp_cmod <- aggregate(.~yd,  npp_cmod[,c("yd","value")], FUN = sum)
colnames(npp_cmod) <- c("yd","Ez_NPP_cmod")

## Cbbp

naames_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
naames_npp <- lapply(naames_npp, npp_func_bbp_c)

naames_npp <- data.frame(matrix(unlist(naames_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(naames_npp) <- table[,c("yd")]
naames_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(naames_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cbbp <- subset(table_long, depth <= Ez1)
npp_cbbp <- aggregate(.~yd,  npp_cbbp[,c("yd","value")], FUN = sum)
colnames(npp_cbbp) <- c("yd","Ez_NPP_cbbp")

naames_npp <- merge(npp_cbbp,npp_cmod, id = "yd")
naames_npp <- merge(naames_npp,npp_cvar, id = "yd")

npp14c <- read.csv(file = "NAAMES_14C.csv", header = T,stringsAsFactors = F)
npp14c$yd <- npp14c$doy
naames_npp <- join(naames_npp,npp14c[,c("yd","int_pp")], by = "yd")

table$year <-  year(as.POSIXct(paste(table$date), format="%Y%m%d",tz="UTC"))
naames_npp <- join(naames_npp,table[,c("yd","year")], by = "yd")

naames_npp$Cruise <- ifelse(naames_npp$year == 2015,"NAAMES1", NA)
naames_npp$Cruise <- ifelse(naames_npp$year == 2016,"NAAMES2", naames_npp$Cruise)
naames_npp$Cruise <- ifelse(naames_npp$year == 2017,"NAAMES3", naames_npp$Cruise)
naames_npp$Cruise <- ifelse(naames_npp$year == 2018,"NAAMES4", naames_npp$Cruise)
naames_npp$year <- NULL
colnames(naames_npp) <- c("yd","Ez_NPP_cbbp","Ez_NPP_cmod","Ez_NPP_cvar","14C","Cruise")

write.csv(naames_npp, "c:/users/foxja/Documents/OSU/Active projects/VARIABLE_SLOPE/Data/Plotting_files/NAAMES_14C.csv",row.names = F)


####  AMT  ####

rm(list = ls())
source('model.R')
table <- read.csv(file = "AMT_model_data.csv", header = T,stringsAsFactors = F)

amt_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
amt_npp <- lapply(amt_npp, npp_func_var_c)

amt_npp <- data.frame(matrix(unlist(amt_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(amt_npp) <- table[,c("yd")]
amt_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(amt_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cvar <- subset(table_long, depth <= Ez1)
npp_cvar <- aggregate(.~yd,  npp_cvar[,c("yd","value")], FUN = sum)
colnames(npp_cvar) <- c("yd","Ez_NPP_cvar")

## Cmod 

amt_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
amt_npp <- lapply(amt_npp, npp_func_mod_c)

amt_npp <- data.frame(matrix(unlist(amt_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(amt_npp) <- table[,c("yd")]
amt_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(amt_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cmod <- subset(table_long, depth <= Ez1)
npp_cmod <- aggregate(.~yd,  npp_cmod[,c("yd","value")], FUN = sum)
colnames(npp_cmod) <- c("yd","Ez_NPP_cmod")


## Cbbp

amt_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
amt_npp <- lapply(amt_npp, npp_func_bbp_c)

amt_npp <- data.frame(matrix(unlist(amt_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(amt_npp) <- table[,c("yd")]
amt_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(amt_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cbbp <- subset(table_long, depth <= Ez1)
npp_cbbp <- aggregate(.~yd,  npp_cbbp[,c("yd","value")], FUN = sum)
colnames(npp_cbbp) <- c("yd","Ez_NPP_cbbp")

amt_npp <- merge(npp_cbbp,npp_cmod, id = "yd")
amt_npp <- merge(amt_npp,npp_cvar, id = "yd")


npp <- read.csv(file = "AMT_14C.csv", header = T,stringsAsFactors = F)
npp$yd <- npp$doy
amt_npp <- join(amt_npp,npp[,c("yd","int_pp")], by = "yd")
amt_npp$Cruise <- "AMT22"
colnames(amt_npp) <- c("yd","Ez_NPP_cbbp","Ez_NPP_cmod","Ez_NPP_cvar","14C","Cruise")

write.csv(amt_npp, "amt_npp.csv",row.names = F)


####  EXPORTS NPP and mu analysis ####

rm(list = ls())
table <- read.csv(file = "EXPORTS_model_data.csv", header = T,stringsAsFactors = F)
source('model.R')

exp_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
exp_npp <- lapply(exp_npp, npp_func_var_c)

exp_npp <- data.frame(matrix(unlist(exp_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(exp_npp) <- table[,c("yd")]
exp_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(exp_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cvar <- subset(table_long, depth <= Ez1)
npp_cvar <- aggregate(.~yd,  npp_cvar[,c("yd","value")], FUN = sum)
colnames(npp_cvar) <- c("yd","Ez_NPP_cvar")

## Cmod 

exp_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
exp_npp <- lapply(exp_npp, npp_func_mod_c)

exp_npp <- data.frame(matrix(unlist(exp_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(exp_npp) <- table[,c("yd")]
exp_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(exp_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd<- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cmod <- subset(table_long, depth <= Ez1)
npp_cmod <- aggregate(.~yd,  npp_cmod[,c("yd","value")], FUN = sum)
colnames(npp_cmod) <- c("yd","Ez_NPP_cmod")


## Cbbp 

exp_npp <- setNames(split(table, seq(nrow(table))), rownames(table))
exp_npp <- lapply(exp_npp, npp_func_bbp_c)

exp_npp <- data.frame(matrix(unlist(exp_npp), nrow=201, byrow=F),stringsAsFactors=FALSE)

colnames(exp_npp) <- table[,c("yd")]
exp_npp$depth <- seq(0,200,1)

table_long <- reshape2::melt(exp_npp, id = "depth")
table_long$yd <- as.character(table_long$variable)
table_long$yd <- as.numeric(table_long$yd)

table_long <- join(table_long,table[,c("yd","mld","Ez1")],
                by = c("yd"),match = "all")

npp_cbbp <- subset(table_long, depth <= Ez1)
npp_cbbp <- aggregate(.~yd,  npp_cbbp[,c("yd","value")], FUN = sum)
colnames(npp_cbbp) <- c("yd","Ez_NPP_cbbp")

exp_npp <- merge(npp_cbbp,npp_cmod, id = "yd")
exp_npp <- merge(exp_npp,npp_cvar, id = "yd")


npp <- read.csv(file = "EXPORTS_14C.csv", header = T,stringsAsFactors = F)
npp$yd <- npp$doy
exp_npp <- join(exp_npp,npp[,c("yd","int_pp")], by = "yd", match = "first")
exp_npp$Cruise <- "EXP1"
colnames(exp_npp) <- c("yd","Ez_NPP_cbbp","Ez_NPP_cmod","Ez_NPP_cvar","14C","Cruise")

write.csv(exp_npp, "exp_npp.csv",row.names = F)

