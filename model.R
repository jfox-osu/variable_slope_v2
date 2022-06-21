npp_func_mod_c <- function(x){
  
  bbp   = x$bbp # change bbp for modelled carbon
  chl   = x$chl
  irr   = x$par
  kd  = x$kd_M07
  mld   = x$mld
  yd = x$yd
  sf =  x$sf*12287
  fe = x$fe
  a490 = x$ap490 
  lat = x$lat
  
  mld[mld > 200] = 200
  
  kd_chl = 0.0166+0.0773*(chl^0.6715)
  PAR_kd = irr^0.45/kd_chl
  Ig = irr * exp(-0.5*kd_chl*mld)
  bb490 = bbp * (490 / 470)^(-1)
  dm = 19 *exp(0.038*PAR_kd)
  sm = ((1+exp(-0.15*irr))/(1+exp(-3*Ig)))
  sm[sm < 1] = 1
  theta  = dm*sm
  sf_mod = (theta*chl)/bbp
  
  a = -16.8
  b= 1.57
  
  c <- ifelse(fe <= 0.99,30 ,47.027)
  d <- ifelse(fe <= 0.99,0.01 ,0.0125)
  
  sf <- ifelse(is.na(sf), sf_mod,sf)
  
  z    = c(0:200)
  uMax = 2
  
  IgFunc <- (1- exp(-5 * irr))
  
  lightmu <-  (I(1 / dm * a) + b) 
  nutmu <-  (I(1 / (dm*sm) * c) + d) 
  mu = (lightmu*nutmu) *  IgFunc
  
  carbon    = (dm*sm)*(chl)
  carbon = ifelse(fe <= 0.99, (dm*sm)*(chl*fe),carbon)
  
  # 3) Compute NPP through depth
  
  NPP     <- array(NA, c(length(chl), 201))
  mu.z   <- array(NA, c(length(chl),201))
  
  #Fill in NPP where data is present. Will overwrite NPP for  data beneath the mld
  #Start at shallowest mld  

  mld.shallowest <- ceiling(min(mld))
  
  for (k in 1:mld.shallowest){
    
    #Find cells below mld
   
     shallow <- which(mld>=z[k], arr.ind = T)
    
    irr.z  =  irr * exp(-kd_chl*z[k])
    
    IgFunc[shallow]      = (1 - exp(-5 * irr.z))
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm*sm) * c) + d) 
    
    mu[shallow] = (lightmu*nutmu) *  IgFunc[shallow]
    
    mu[shallow][ mu[shallow]>uMax] <- uMax
    
    NPP[shallow,k]  = mu[shallow] * carbon[shallow]
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  for (k in mld:201){
    
    #Find cells below mld

    deep <- which(mld<=z[k], arr.ind = T)
    
    irr.z     =  irr * exp(-kd_chl*z[k])
    
    IgFunc[deep]      = (1 - exp(-5 * irr.z))
    
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm) * c) + d) 
    
    mu[deep] = (lightmu*nutmu) *  IgFunc[deep]
    
    mu[deep][ mu[deep]>uMax] <- uMax
    
    NPP[deep,k]  = round(mu[deep] * carbon[deep],3)
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  return(NPP)
  
  
}

npp_func_var_c <- function(x){
  
  bbp   = x$bbp # change bbp for modelled carbon
  chl   = x$chl
  irr   = x$par
  mld   = x$mld
  yd = x$yd
  sf =  x$sf*12287
  fe = x$fe
  a490 = x$ap490 
  lat = x$lat
  
  mld[mld > 200] = 200

  kd_chl = 0.0166+0.0773*(chl^0.6715)
  PAR_kd = irr^0.45/kd_chl
  Ig = irr * exp(-0.5*kd_chl*mld)
  
  dm = 19 *exp(0.038*PAR_kd)
  sm = ((1+exp(-0.15*irr))/(1+exp(-3*Ig)))
  sm[sm < 1] = 1
  theta  = dm*sm
  sf_mod = (theta*chl)/bbp
  
  a = -16.8
  b= 1.57
  
  c <- ifelse(fe <= 0.99,30 ,47.027)
  d <- ifelse(fe <= 0.99,0.01 ,0.0125)
  
  sf <- ifelse(is.na(sf), sf_mod,sf)
  
  z    = c(0:200)
  uMax = 2
  
  IgFunc <- (1- exp(-5 * irr))
  
  lightmu <-  (I(1 / dm * a) + b) 
  nutmu <-  (I(1 / (dm*sm) * c) + d) 
  mu = (lightmu*nutmu) *  IgFunc
  
  #carbon    = (dm*sm)*(chl)
  carbon  = (bbp-0.00027)*(sf)

  # 3) Compute NPP through depth
  
  NPP     <- array(NA, c(length(chl), 201))
  mu.z   <- array(NA, c(length(chl),201))
  
  #Fill in NPP where data is present. Will overwrite NPP for  data beneath the mld
  #Start at shallowest mld  
  
  mld.shallowest <- ceiling(min(mld))
  
  for (k in 1:mld.shallowest){
    
    #Find cells below mld
    
    shallow <- which(mld>=z[k], arr.ind = T)
    
    irr.z  =  irr * exp(-kd_chl*z[k])
    
    IgFunc[shallow]      = (1 - exp(-5 * irr.z))
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm*sm) * c) + d) 
    
    mu[shallow] = (lightmu*nutmu) *  IgFunc[shallow]
    
    mu[shallow][ mu[shallow]>uMax] <- uMax
    
    NPP[shallow,k]  = mu[shallow] * carbon[shallow]
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  for (k in mld:201){
    
    #Find cells below mld
    
    deep <- which(mld<=z[k], arr.ind = T)
    
    irr.z     =  irr * exp(-kd_chl*z[k])
    
    
    IgFunc[deep]      = (1 - exp(-5 * irr.z))
    
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm) * c) + d) 
    
    mu[deep] = (lightmu*nutmu) *  IgFunc[deep]
    
    mu[deep][ mu[deep]>uMax] <- uMax
    
    NPP[deep,k]  = round(mu[deep] * carbon[deep],3)
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  return(NPP)
  
  
}




npp_func_bbp_c <- function(x){
  
  bbp   = x$bbp # change bbp for modelled carbon
  chl   = x$chl
  irr   = x$par
  kd  = x$kd_M07
  mld   = x$mld
  yd = x$yd
  sf =  x$sf*12287
  fe = x$fe
  a490 = x$ap490 
  lat = x$lat
  
  mld[mld > 200] = 200
  
  kd_chl = 0.0166+0.0773*(chl^0.6715)
  PAR_kd = irr^0.45/kd_chl
  Ig = irr * exp(-0.5*kd_chl*mld)
  
  dm = 19 *exp(0.038*PAR_kd)
  sm = ((1+exp(-0.15*irr))/(1+exp(-3*Ig)))
  sm[sm < 1] = 1
  theta  = dm*sm
  sf_mod = (theta*chl)/bbp
  
  a = -16.8
  b= 1.57
  
  c <- ifelse(fe <= 0.99,30 ,47.027)
  d <- ifelse(fe <= 0.99,0.01 ,0.0125)
  
  sf <- ifelse(is.na(sf), sf_mod,sf)
  
  z    = c(0:200)
  uMax = 2
  
  IgFunc <- (1- exp(-5 * irr))
  
  lightmu <-  (I(1 / dm * a) + b) 
  nutmu <-  (I(1 / (dm*sm) * c) + d) 
  mu = (lightmu*nutmu) *  IgFunc
  
  #carbon    = (dm*sm)*(chl)
  carbon  = (bbp)*(12128)
  
  # 3) Compute NPP through depth
  
  NPP     <- array(NA, c(length(chl), 201))
  mu.z   <- array(NA, c(length(chl),201))
  
  #Fill in NPP where data is present. Will overwrite NPP for  data beneath the mld
  #Start at shallowest mld  
  
  mld.shallowest <- ceiling(min(mld))
  
  for (k in 1:mld.shallowest){
    
    #Find cells below mld
    
    shallow <- which(mld>=z[k], arr.ind = T)
    
    irr.z  =  irr * exp(-kd_chl*z[k])
    
    IgFunc[shallow]      = (1 - exp(-5 * irr.z))
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm*sm) * c) + d) 
    
    mu[shallow] = (lightmu*nutmu) *  IgFunc[shallow]
    
    mu[shallow][ mu[shallow]>uMax] <- uMax
    
    NPP[shallow,k]  = mu[shallow] * carbon[shallow]
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  for (k in mld:201){
    
    #Find cells below mld
    
    deep <- which(mld<=z[k], arr.ind = T)
    

    irr.z     =  irr * exp(-kd_chl*z[k])
    
    
    IgFunc[deep]      = (1 - exp(-5 * irr.z))
    
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm) * c) + d) 
    
    mu[deep] = (lightmu*nutmu) *  IgFunc[deep]
    
    mu[deep][ mu[deep]>uMax] <- uMax
    
    NPP[deep,k]  = round(mu[deep] * carbon[deep],3)
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  return(NPP)
  
  
}

npp_func_var_c_naames <- function(x){
  
  bbp   = x$bbp # change bbp for modelled carbon
  chl   = x$chl
  irr   = x$par
  kd  = x$kd_M07
  mld   = x$mld
  yd = x$yd
  sf =  x$sf*12287
  fe = x$fe
  a490 = x$ap490 
  lat = x$lat
  
  mld[mld > 200] = 200
  a490 <- ifelse(is.na(a490), chl * 0.03,a490)
  
  kd_chl = 0.0166+0.0773*(chl^0.6715)
  PAR_kd = irr^0.45/kd_chl
  Ig = irr * exp(-0.5*kd_chl*mld)
  
  dm = 19 *exp(0.038*PAR_kd)
  sm = ((1+exp(-0.15*irr))/(1+exp(-3*Ig)))
  sm[sm < 1] = 1
  theta  = dm*sm
  sf_mod = (theta*chl)/bbp
  
  a = -16.8
  b= 1.57
  
  c <- ifelse(fe <= 0.99,30 ,47.027)
  d <- ifelse(fe <= 0.99,0.01 ,0.0125)
  
  sf <- ifelse(yd >= 315, sf_mod*1.75,sf)
  sf <- ifelse(is.na(sf), sf_mod,sf)
  
  z    = c(0:200)
  uMax = 2
  
  IgFunc <- (1- exp(-5 * irr))
  
  lightmu <-  (I(1 / dm * a) + b) 
  nutmu <-  (I(1 / (dm*sm) * c) + d) 
  mu = (lightmu*nutmu) *  IgFunc
  
  #carbon    = (dm*sm)*(chl)
  carbon  = (bbp-0.00027)*(sf)
  
  # 3) Compute NPP through depth
  
  NPP     <- array(NA, c(length(chl), 201))
  mu.z   <- array(NA, c(length(chl),201))
  
  #Fill in NPP where data is present. Will overwrite NPP for  data beneath the mld
  #Start at shallowest mld  
  
  mld.shallowest <- ceiling(min(mld))
  
  for (k in 1:mld.shallowest){
    
    #Find cells below mld
    
    shallow <- which(mld>=z[k], arr.ind = T)
    
    irr.z  =  irr * exp(-kd_chl*z[k])
    
    IgFunc[shallow]      = (1 - exp(-5 * irr.z))
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm*sm) * c) + d) 
    
    mu[shallow] = (lightmu*nutmu) *  IgFunc[shallow]
    
    mu[shallow][ mu[shallow]>uMax] <- uMax
    
    NPP[shallow,k]  = mu[shallow] * carbon[shallow]
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  for (k in mld:201){
    
    #Find cells below mld
    
    deep <- which(mld<=z[k], arr.ind = T)
    
    irr.z     =  irr * exp(-kd_chl*z[k])
    
    
    IgFunc[deep]      = (1 - exp(-5 * irr.z))
    
    
    lightmu <-  (I(1 / dm * a) + b) 
    nutmu <-  (I(1 / (dm) * c) + d) 
    
    mu[deep] = (lightmu*nutmu) *  IgFunc[deep]
    
    mu[deep][ mu[deep]>uMax] <- uMax
    
    NPP[deep,k]  = round(mu[deep] * carbon[deep],3)
    
    #Parameters dependent on level z-1
    
    print(k)
  }
  
  return(NPP)
  
  
}



