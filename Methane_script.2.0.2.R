


rm(list = ls())

start.time <- Sys.time()


##########################
#INSERT CODE DESCRIPTION##
##########################

#########
#MODEL
#########

#########
#DESCRIPTION OF INPUTS
#########

#########
#DESCIPTION OF OUTPUTS
#########


#########
#DEPENDENCIES
#########
library(limSolve) ##Linear solvers
library(pracma)  ##Linear algebra toolkit
library(ggplot2)  ##Plotting toolkit

#########
#OUTPUT FILE
#########

experName <-  "Temp_3c_b"
experName2 <- "Ts_b"

#########
#FIXED PARAMETERS
#########

Q10_moxid <- 2 #Q10 of methane oxidation
operch4 <- 4  #numbers of oxygen
R10 <- 1
rho <- 1000  #density of water (kg m-3)
g <- 9.81    #acceleration parameter (m s-2)
Patm <- 1e5   #atmospheric pressure (Pa)

tenc <- 283.15   # ten degrees C in Kelvin eqn 11 Loyd and Taylor
To <- 227   # Kelvin
Eo <- 308.56 
mmtom <- .001 #conversion of mm to meters

#inverse Henry constants for oxygen and methane (mol m-3)/Pa
iHenryO2_0 <-   1/769.23 * rho * 1e-5 * 8.13*298.15  
iHenryCH4_0 <-  1/714.29 * rho * 1e-5 * 8.13*298.15

mconsump_max <- 1.25e-5*86400*365 #from CLM but lower value
kch4 <- 5e-2 #5e-1 #(mol m-3) #10 times higher value
ko2 <- 2e-2  #(mol m-3)

gtmol <- 12

#########
#INPUTS
########

ranox <- 0.4
k10exu <- 13     
k10litter <- 0.1  
k10fast <- .03
k10slow <- .001
fexu <-  0.0
fatm <- 1.0
ffast <- 0.985
fslow <- 0.015
NPP <- (1.096300 *1000)/gtmol # Mol m-2 yr-1 ## 1 ) 1.096300 ### 2) 0.957200  ## 3) 0.954100 ## kg C m2  yr
ch4frac0 <- 0.25
bgNPP <- 0.5
LAI <- 3
ctiller <- 0.22
aeren_rad <- 3e-3
wilt <- 0.01
fsand <- .50
fom <- .50
fclay <- 0       
thetasatom <- 0.9 
psisatom <- -0.0103
bom <- 2.7   
wpos <-  0.5
eta <- 400
mresp <- 0.5
rturn <- 1.3
rgresp <- .33
ch4_bubbles <- 0.15
lambdaroot <- .2517
rL <- 3 
P <- 0.3
wtm <-   0
dm1 <- 3.33
dm2 <- 2
o2om <- 0.2

Tmean <-  293 # K 
A = 5 # amplitude of temp

## half of double 
lambdadryom <-  0.05 # Wm-1K-1dry thermal conductivity of OM Farouki,1981 in CLM
lambdaliq <-    0.57       #  Wm-1K-1 table 2.3 CLM
lambdasom <-  0.25 # W m-1 K-1 #heat exchange coefficient for organic matter


#cice <-   2.11727*10^3     # Specific heat capacity of ice J kg-1 K-1 Table 2.6 CLM
cliq <-  4.188*10^3  * rho # Specific heat capacity of water J kg-1 K-1  converted to J m-3 K-1 Table 2.6 CLM 
#csbedr <-  2 * 10^6 # heat capacity of bedrock J m-3 K-1
csom <- 2.5*10^6 # J m-3K-1 heat capacity of OM


###############END INPUT#############


# Ambient temmperature initial conditions??


Tair <- Tmean # Initial conditions?gerber
Tc <- Tair - 273.15  # initial Temp for diffusion coeffieciens in Celsius 



#########
#Iteration parameters
#########
duration <- 400

# Use formula to determine correct time step 
#https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis

#sg: commented this - we need to calculate everything else first to make timestep determination
#suggest to put in a fixed timestep, but manually control --

#timestepnew <- 0.5 * ((ci_h * ddepth^2 ) / lambdai)

timestep <- 5.e-6 # year  1.285039e-04 
nite <- duration/timestep 
nout <- 4800 # base 200 
ointerval <- nite %/% nout 
nout <- nout + 1



#INPUT and constant boundary conditions

#Air
p_O2 <- Patm*0.2
o2_air <- p_O2/(8.13*Tair)
pch4 <- Patm*1.8e-6
ch4_air <- pch4/(8.13*Tair)


#Grid setup
depth <- c(0.0071,0.0279,0.0623, 0.1189, 0.2122,0.3661,0.6198,1.0380,1.7276,2.8646, 4.7392, 7.8298, 12.9253, 21.32,35.1776 ) # gerber 
ndepth <- length(depth)
top <- c(0,(depth[1:(ndepth - 1)] + depth[2:ndepth])/2)  #location of upper boundary
bottom <- c(top[2:ndepth],depth[ndepth] + (depth[ndepth] - top[ndepth])) #location of lower boundary
depthborders <- c(top,bottom[ndepth]) # upper and lower boundary
ddepth <- depthborders[2:(ndepth + 1)] - depthborders[1:(ndepth)] #thickness of layer
dzm1 <- c(depth[1], depth[2:ndepth] - depth[1:(ndepth - 1)])#distance from midpoint to midpoint above
dzp1 <-  c(dzm1[2:ndepth],bottom[ndepth] - depth[ndepth]) #distance from midpoint to midpoint below

#####
#Soil layer boundary  calculations
####
#Centers
centers <- depth[1:(ndepth - 1)] + ((depth[2:ndepth] - depth[1:(ndepth - 1)])/2)

#Distances
Dist_up <- rep(NA,(ndepth - 2))
for (ind in 2:(ndepth - 1)) {
  Dist_up[ind - 1] <- depth[ind] - centers[ind - 1]  
}
Dist_down <- rep(NA,(ndepth - 2))
for (ind in 2:(ndepth - 1)) {
  Dist_down[ind - 1] <- centers[ind] - depth[ind]
}
#Weights
W_high <- rep(NA,(ndepth - 2))
for (ind in 1:(ndepth - 2)) {
  W_high[ind] <- Dist_down[ind]/(Dist_down[ind] + Dist_up[ind])
}
W_high <- c(W_high,0) #For the last level
W_low <- 1 - W_high



###################################################
#####Set up secondary variables from values above##
###################################################


## Soil Moisture calculations un changing 
psisatmin <- (-10.0 * 10^(1.88 - 0.0131*(fsand*100))) * mmtom   # water potential 
bmin <- 2.91 + 0.159*(fclay*100)   # water retention parameter
b <- (1 - fom)*bmin + fom*bom          #  water retention parameter
thetasatmin <- 0.489 - 0.00126 * (fsand*100)   #  saturation volume fraction
psisat <- (1 - fom)*psisatmin + fom*psisatom   # 
thetasat <- (1 - fom) * thetasatmin + fom*thetasatom  # old PSP 
bd = 2700*(1 - thetasat)       # bulk density ( kg m-3) eq 6.84 CLM c


#heat exchange parameterization, boundary conditions and initial conditions


lambdasmin <-  (8.80 * fsand + 2.92*fclay)/(fsand + fclay) # mineral soil solid therman conductivity ( unitless?)
lambdas <-  (1 - fom)*lambdasmin + fom*lambdasom   # W m-1 K-1  ?
lambdadrymin  <-  (0.135*bd + 64.7)/(2700 - 0.947*bd) # thermal conductivity if mineral soil ( W m-1 K-1) 6.84
lambdadryi <-  (1 - fom) * (lambdadrymin + fom * lambdadryom)  * (86400 * 365) #  (J yr-1 m-1 K-1)
lambdasat <-  lambdas^(1 - thetasat) * lambdaliq^(thetasat)  * (86400*365)   # (J yr-1 m-1 K-1)


csmin <- ((2.128*fsand + 2.385*fclay) / (fsand + fclay))*10^6  # heat capacity of mineral soil solids (J m-3 K-1)
cs <-  (1 - fom)*csmin + fom*csom  # heat capacity of soil solids ( J m-3 K-1)


#Araenchima area
A_aerench = bgNPP*NPP*LAI/ctiller*3.14*(aeren_rad)^2

#fraction of root with depth
froot <- exp(-top/lambdaroot) - exp(-bottom/lambdaroot) #.2517 lambda root # Wania  

##################
##################
###Containers for storing working values
##################
##################


#output variables
outfile = paste(experName,".dat",sep= "")
outfile2= paste(experName2,".dat", sep = "")
semdifch4 = 0
semaerenchch4 = 0
semch4 = 0
semebullch4 = 0
smineralizedc = 0
semo2 = 0
sco2_het = 0

#write header
cat("timestep","NPP", "wt","sco2_het","semch4","semo2","smineralizedc", "semdifch4 ","semaerenchch4 ","semebullch4","\n",file=outfile)  

cat("timestep","T0.0071","T0.0279","T0.0623", "T0.1189", "T0.2122","T0.3661","T0.6198","T1.0380","T1.7276","T2.8646", "T4.7392", "T7.8298", "T12.9253", "T21.32","T35.1776", "Tmean","\n", file=outfile2)




#litter<-matrix(0,ndepth)  #g/m2/yr?
#litter<-matrix(1,ndepth)*froot*NPP/(k10litter*Rtot)  #mol/m2 #estimate from equilibrium


## keep sg 
litter <- matrix(1,ndepth)*froot*NPP/(k10litter*0.2)  #mol/m2 #estimate from equilibrium
exu <- matrix(0,ndepth)
fast <- exu
slow <- exu

#co <- matrix(1,ndepth)*.9*o2_air  #mol m -3
#co <- matrix(1,ndepth)*1.*o2_air  #mol m -3
co <- matrix(1,ndepth) * 0    # keep sg
#  km <- (0.1875 + 0.0013*Ts) * 10^-4 * 86400 * 365  #methane exchange coeficcient m2 yr-1, currently not used
#cm <- matrix(1,ndepth)*1*ch4_air   # Concentration of methane UP initial in mol m-3
#cm <- matrix(1,ndepth)*1.*ch4_air   # Concentration of methane UP initial in mol m-3
cm <- matrix(1,ndepth) * 0 # keep sg
Ts  <- matrix(1,ndepth) * Tmean  # soil temp in K #keep sg
#Tc  <- matrix(1,ndepth) * Tmean - 273.15  # soil temp in C 
k_odif <- matrix(0,ndepth)    #oxygen diffusion coefficient in a staggered grid
k_mdif <- matrix(0,ndepth)   #methane diffusion coefficient in a staggered grid
#co2prod <- matrix(0,ndepth) #mol/m3/ yr
mprod <- matrix(0,ndepth) #methane production in soil (mol m-3 yr-1)
mconsump <- matrix(0,ndepth) #methane consumption in soil (mol m-3 yr-1)
maerench <- matrix(0,ndepth) #araenchyma transport (mol m-3 yr-1)
ebullition <- matrix(0,ndepth)
oresp <- matrix(0,ndepth) #oxygen consumption when respiring (mol m-3 yr-1)
oconsump <- matrix(0,ndepth) #oxygen consumption when respiring methane (mol m-3 yr-1)
oaerench <- matrix(0,ndepth) #oxygen araenchima transport (mol m-3 yr-1)
sch4 <-  matrix(0,ndepth)   #source of new methane
so2 <-  matrix(0,ndepth)   #source of new oxygen
fdifo2 <- matrix(0,ndepth)  #diffusion at each layer boundary (except lowest)
fdifch4 <- matrix(0,ndepth) #

Rm <- matrix(1,ndepth)      #modifier to calculate effective air concentrations
Ro <- matrix(1,ndepth)      
rm <- matrix(1,ndepth)      #modifier to calculate biological relevant concentration
ro <- matrix(1,ndepth)

thetaw <- matrix(0,ndepth)


###########################
#Data Collection for output
###########################
years <- 2 #Collect daily data for last "years" number of years
summary.df <- data.frame(semdifch4= rep(NA,nout) ,semaerenchch4 = rep(NA,nout), semebullch4 = rep(NA,nout))
output.df <- data.frame(emdifch4 = rep(NA,years*365),emaerenchch4 = rep(NA,years*365), emebull = rep(NA,years*365), wt = rep(NA,years*365))
output.mconsump <- data.frame(level = rep(1:length(depth),years*365), value = rep(NA, length(depth)*years*365), day = rep(NA, length(depth)*years*365))
output.mprod <- matrix(NA,length(depth),years*365)
k <- 1
j <- 1



################
#Run Main loop
################



#nite = 2 
for (i in 1:(nite - 1)) {

  if (i == nite %/% 2) {
   
    Tmean <-  Tmean + 0 # K 
     
  }
  # first step: update water table 
  wt = wtm + wpos*sin(2*pi*i*timestep)
  depthlesswt <- depth < wt
  sumdepth <- sum(depthlesswt) > 0
  thetaw <- rep(thetasat,length(depth))
  thetaw[depthlesswt] <- thetasat * ((psisat - (wt - depth))/(psisat))[depthlesswt]^(-1/b)
  thetaw <-  pmax(thetaw,wilt)  # Do not allow for thetaw to go bellow wilt point
  

  jwt <- which.min(abs(depth - wt))
  if (depth[jwt] < wt) jwt <- jwt + 1
  
  thetaa <- thetasat - thetaw
  
  
  ### Heat exchange Calculations ###
  
  
  Tair <- Tmean + A * sin(2*pi*i*timestep) # Ambient T Kelvin 
  
  # wliq= thetaw*rho*ddepth ## ddept is deltaz in CLM
  
  ci_h  <- cs * (1 - thetasat) + (thetaw) * cliq   # volumetric heat capacity (J m-3 K-1) 
  
  Sri <-  thetaw/thetasat # wetness of soil with respect to saturation 6.86 ( unitless)
  Kei <-  log(thetaw/thetasat) + 1  # unitless?
  lambdaweti <-    Kei *  lambdasat + (1 - Kei) * lambdadryi  #thermal conductivity iλ [W m-1 K-1]; Jyr-1m eq 6.28 CLM depends on Sr 
  lambdai <- ifelse(Sri <= 1.e-7,lambdadryi, lambdaweti)  #
  
  

  
  
  ###Crank Nicholson solution, and tridiagonal solver###
  ###Temperature###
  
  lambdai_m1 <- c( lambdai[1], 1/(W_low/lambdai[2:ndepth] + W_high/lambdai[1:(ndepth - 1)]))
  lambdai_p1 <- c(lambdai_m1[2:ndepth],0)
  
  # change to Ts (lambda)
  diam1_ts     <- -lambdai_m1/(2*ddepth*dzm1) #diagonal below center
  diap1_ts     <- -lambdai_p1/(2*ddepth*dzp1)  #diagonal above center
  dia_ts   <-  ci_h/timestep + (lambdai_m1/dzm1 + lambdai_p1/dzp1)/(2*ddepth) # center diagonal
  
  hts  <- c(Tair,Ts,Ts[ndepth]) #help variable, temp at the boundaries ##
  
  rhs_ts   <-  Ts*ci_h/timestep  + 
    1/(2*ddepth) *
    (lambdai_p1/dzp1*(hts[3:(ndepth + 2)] - hts[2:(ndepth + 1)]) - 
       lambdai_m1/dzm1*(hts[2:(ndepth + 1)] - hts[1:ndepth]))
  # ##### ( cm *ci_h * ddpepth) or wo ()  - Rm is ci_h * ddeppth
  
  
  #update righ hand side with boundary conditions, top and bottom
  #not nessecary, because it includes boundary conditions??
  rhs_ts[1] = rhs_ts[1] - diam1_ts[1]*hts[1]
  rhs_ts[ndepth] = rhs_ts[ndepth] + diap1_ts[ndepth]*hts[(ndepth + 2)]
  
  #Stefan additions
  T_old = Ts
  Ts <-  Solve.tridiag(diam1_ts[2:ndepth],dia_ts, diap1_ts[1:(ndepth - 1)],rhs_ts) #note the shorter vectors m1 and p1  
  deltaT = Ts - T_old
  
  # Calculations to check for conservation of Heat
  # htsnew   = c(Tair,Ts)
  # heatflux = (0.5*lambdai_m1*(htsnew[1:ndepth]-htsnew[2:(ndepth+1)]+hts[1:ndepth]-hts[2:(ndepth+1)]))/dzm1 
  # heatflux = c(heatflux,0)
  #sumheat = sum(ci_h * ddepth * (Ts - T_old))
  
  #heatflux[1]
  Tc <- Ts - 273.15   # Soil Temp needs to be in Celsius for Diffusion calculations
  
  
  #Soil properties
  #Diffusion 
  
  D_O2_g  <- (.1759 + .00117 * Tc )*1e-4 * 86400 * 365
  D_O2_aq <- (1.172 + 0.03443 * Tc + 0.0005048*Tc^2)*1e-9*86400*365
  D_CH4_g <-  (0.1875 + 0.0013*Tc)*1e-4*86400*365
  D_CH4_aq <- (0.9798 + 0.02986*Tc + 0.0004381 * Tc^2)*1e-9*86400*365
  Sc_CH4 <- 1898 - 110.1*Tc + 2.834*Tc^2 - 0.0279*Tc^3
  Sc_O2 <- 1800.6 - 120.1*Tc + 3.7818*Tc^2 - 0.0476*Tc^3

  
  #gas diffusion parameters
  Dm  <- rep(thetasat^dm2,ndepth)
  if (sumdepth) {
    Dm[depthlesswt] <- thetaa[depthlesswt]^(dm1)/(thetasat^dm2)
  }
  
  
  D_CH4 <- rep(D_CH4_aq*Dm,ndepth) 
  D_CH4[depthlesswt] <- D_CH4_g*Dm
  
  D_O2 <- rep(D_O2_aq*Dm,ndepth) 
  D_O2[depthlesswt] <- D_O2_g*Dm
  
  #Diffusion at the boundaries
  #kludge to deal with indundated surface
  if (wt > 0) D_CH4_S <- D_CH4_g[1] else D_CH4_S = D_CH4_aq[1]
  #D_CH4_S[wt <= 0] <- D_CH4_aq
  
  
  D_CH4_m1 <- c( 2/(1/D_CH4_S + 1/D_CH4[1]), 1/(W_low/D_CH4[2:ndepth] + W_high/D_CH4[1:(ndepth - 1)]))
  D_CH4_p1 <- c(D_CH4_m1[2:ndepth],0)
  D_O2_S  <- ifelse(wt <= 0,D_O2_aq,D_O2_g)
  D_O2_m1 <- c( 2/(1/D_O2_S + 1/D_O2[1]), 1/(W_low/D_O2[2:ndepth] + W_high/D_O2[1:(ndepth - 1)]))
  D_O2_p1 <- c(D_O2_m1[2:ndepth],0)
  
  ##############
  
  Rm <- thetaa + iHenryCH4_0*thetaw
  Ro <- thetaa + iHenryO2_0*thetaw
  
  
  #next timestep
  
  RT <- R10*exp(Eo*(1/(tenc - To) - 1 / (Ts - To)))  # Lloyd and Taylor ...describe more here carla 
  Rmoist <- rep(ranox, ndepth)
  Rmoist[thetaa > 0] <- 1
  Rtot <- RT*Rmoist  
  
  #rm <- ifelse(depth<wt,1,iHenryCH4_0)
  #ro <- ifelse(depth<wt,1,iHenryO2_0)
  rm <- Rm
  ro <- Ro
  
  ch4frac <- ch4frac0*1/(1 + eta * ro*co) # eq 19.4 clm
  co2frac <-  1 - ch4frac
  
  olitter <- litter*k10litter*Rtot
  oexu   <- exu*k10exu*Rtot
  ofast  <- fast*k10fast*Rtot
  oslow  <- slow*k10slow*Rtot	
  ilitter <- NPP*(1 - fexu)*froot   # root depth 
  iexu   <- fexu*NPP*froot
  islow  <- olitter*(1 - fatm)*fslow #do not input into som pools here here because... takes long time to reach equilibrium
  ifast  <- olitter*(1 - fatm)*ffast #to obtain mineralized c
  mineralizedc <- (olitter + oexu + ofast + oslow - ifast - islow)/ddepth
  
  #mineralizedc <- NPP*froot/ddepth
  
  co2prod <- mineralizedc*co2frac  
  mprod <- mineralizedc*ch4frac
  
  # Wania methanotrophy
  if (sumdepth) {mconsump <- pmin(cm,0.5*co) #replace ifelse()
  mconsump[!depthlesswt] <- 0}
  else{mconsump <- rep(0,ndepth)}    
  cm = cm - mconsump
  co = co - 2*mconsump
  mconsump = mconsump/timestep
  
  
  #maintenance and growth respiration estimated based on NPP 				
  #approximation of CLM base respiration using C/N = 60 and growth respiration is a 1/3 of NPP, and roots turn
  #over annually
  
  rrespiration <- (NPP*froot*mresp*rturn*RT + NPP*froot*rgresp)/ddepth    
  co2prod      <- co2prod + rrespiration   
  co2prod      <- co2prod
  oresp        <- co2prod  * (1 - o2om) #assume 20% of oxygen can be taken out of organic matter)
  

  maxrate <- mconsump_max/ko2*ro*co #oxygen limited rate for heterotrophs is the same as for methanotrophs
  scale <- maxrate/(oresp + oconsump) #remove ifelse
  scale[!((oresp + oconsump) > maxrate)] < -1
  scale1  <- 1 + ch4frac*(1 - scale)
  
  oexu <- oexu*scale1
  olitter <- olitter*scale1 
  ofast <- ofast*scale1
  oslow <- oslow*scale1
  islow <- islow*scale1  #do input into som pools here because...
  ifast <- ifast*scale1  #to obtain mineralized c
  
  rrespiration <- rrespiration*scale
  co2prod <- co2prod*scale
  mineralizedc <- mineralizedc*scale1    #scale back co2 but not methane
  
  oresp <- oresp*scale
  oconsump <- oconsump*scale
  mconsump <- mconsump*scale
  
  co2_het <- co2prod - rrespiration
  
  #update organic matter pools
  litter <- litter + (ilitter - olitter)*timestep
  
  #aerenchyma transport, formulated as a transport out (sink)
  #note: cm[j] already is gaseous concentration
  #division by ddepth to obtain a source term
  maerench <-  ((cm - ch4_air)*D_CH4_g*P*froot*A_aerench/(rL*depth) ) / ddepth
  
  oaerench <- (((co - o2_air)*D_O2_g*P*froot*A_aerench)/(rL*depth) )  / ddepth 
  
  #ebullition of methane -> ebullition going out
  ch4_crit <-  ch4_bubbles*(Patm + g*(depth - wt)*rho) / (8.31*Ts)
  #ebullition over a daily timescale, hence the factor 365
  ebullition <- 365*(cm - ch4_crit)  #Remove ifelse
  ebullition[!((cm > ch4_crit) & (depth > wt))] <- 0
  febull <- ebullition*ddepth*Rm # mol m-2 year-1
  sum_ebull  <- sum(febull)
  if (jwt == 1) { # ebullition into atmosphere
    emebull <- sum_ebull
  } 
  else  {
    ebullition[jwt - 1] <- -sum_ebull/ddepth[jwt - 1] 
    emebull <- 0.
    febull[jwt - 1]  <- -sum_ebull
  }
  
  #total sources and sinks
  sch4 <-  mprod - maerench - ebullition # Wania
  so2 <- -oresp - oaerench  # Wania 


  
  
  
  #Crank Nicholson solution, and tridiagonal solver
  #fist methane
  diam1_ch4 <- -D_CH4_m1/(2*ddepth*dzm1) #diagonal below center
  diap1_ch4 <- -D_CH4_p1/(2*ddepth*dzp1) #diagonal above center
  dia_ch4   <-  Rm/timestep + (D_CH4_m1/dzm1 + D_CH4_p1/dzp1)/(2*ddepth) #center diagonal
  hcm       <- c(ch4_air,cm,cm[ndepth]) #help variable, add concentration at the boundaries
  
  rhs_ch4   <-  cm*Rm/timestep  + 
    1/(2*ddepth) *
    (D_CH4_p1/dzp1*(hcm[3:(ndepth + 2)] - hcm[2:(ndepth + 1)]) - D_CH4_m1/dzm1*(hcm[2:(ndepth + 1)] - hcm[1:ndepth])) +
    sch4
  #update righ hand side with boundary conditions, top and bottom
  #not nessecary, because it includes boundary conditions??
  rhs_ch4[1] = rhs_ch4[1] - diam1_ch4[1]*hcm[1]
  rhs_ch4[ndepth] = rhs_ch4[ndepth] + diap1_ch4[ndepth]*hcm[(ndepth + 2)]
  
  cm = Solve.tridiag(diam1_ch4[2:ndepth],dia_ch4, diap1_ch4[1:(ndepth - 1)],rhs_ch4) #note the shorter vectors m1 and p1  
  
  #repeat the same for oxygen
  diam1_o2 <- -D_O2_m1/(2*ddepth*dzm1) #diagonal below center
  diap1_o2 <- -D_O2_p1/(2*ddepth*dzp1) #diagonal above center
  dia_o2   <-  Ro/timestep + (D_O2_m1/dzm1 + D_O2_p1/dzp1)/(2*ddepth) #center diagonal
  hco      <- c(o2_air,co,co[ndepth]) #help variable, add concentration at the boundaries
  
  rhs_o2   <-  co*Ro/timestep  + 
    1/(2*ddepth) *
    (D_O2_p1/dzp1*(hco[3:(ndepth + 2)] - hco[2:(ndepth + 1)]) - D_O2_m1/dzm1*(hco[2:(ndepth + 1)] - hco[1:ndepth])) +
    so2
  
  #update righ hand side with boundary conditions, top and bottom
  rhs_o2[1] = rhs_o2[1] - diam1_o2[1]*hco[1]
  rhs_o2[ndepth] = rhs_o2[ndepth] + diap1_o2[ndepth]*hco[(ndepth + 2)]
  
  co_old <- co
co = Solve.tridiag(diam1_o2[2:ndepth],dia_o2, diap1_o2[1:(ndepth - 1)],rhs_o2) #note the shorter vectors m1 and p1
  
  
  #diagnosis of diffusion, flux this time, i.e. mol m-2 yr-1
  fdifch4 = -diam1_ch4*ddepth*((hcm[1:ndepth] - hcm[2:(ndepth + 1)]) + (c(ch4_air,cm[1:(ndepth - 1)]) - cm) )   
  fdifo2 = -diam1_o2*ddepth*((hco[1:ndepth] - hco[2:(ndepth + 1)]) + (c(o2_air,co[1:(ndepth - 1)]) - co) ) 
  
  #diagnosis of aerenchima flux (from mol m-3 yr-1 to mol m-2 yr-1
  faerenchch4 = maerench*ddepth
  faerencho2  = oaerench*ddepth
  
  #ebullition has already been diagnosed
  
  
  #emission of methane (transport through surface)
  emdifch4 = -fdifch4[1]
  emaerenchch4 = sum(faerenchch4)
  emch4 = emdifch4 + emaerenchch4 + emebull
  
  # the same for oxygen
  emdifo2 = -fdifo2[1]
  emaerencho2 = sum(faerencho2)
  emo2 = emdifo2 + emaerencho2
  
  
  #cumulative emissions (sums)
  semdifch4 = semdifch4 + emdifch4/ointerval
  semaerenchch4 = semaerenchch4 + emaerenchch4/ointerval
  semch4 = semch4 + emch4/ointerval
  semebullch4 = semebullch4 + emebull/ointerval
  smineralizedc = smineralizedc + sum(mineralizedc*ddepth)/ointerval
  semo2 = semo2 + emo2/ointerval
  sco2_het = sco2_het + sum(co2_het*ddepth)/ointerval
  
  if (i %% ointerval == 0 || i == nite) {
    if (i*timestep > (duration - years)) {
    summary.df[k,] <- c(semdifch4 ,semaerenchch4, semebullch4) # sum over nout - creates average -- cannot use first value its 100myear avergage 
    #emch4.df[k,] <- c(semdifch4 ,semaerenchch4, semebullch4)
    }
    cat(i*timestep,NPP,wt,sco2_het, semch4, semo2, smineralizedc, semdifch4,semaerenchch4,semebullch4,"\n",file = outfile, append = TRUE)  
   
  
    y <- matrix(c(t(Ts),Tmean), nrow = 1, ncol = 16)
    cat(i*timestep,y, file = outfile2, append = T,"\n")
    
    
     semdifch4 = 0
    semaerenchch4 = 0
    semch4 = 0
    semebullch4 = 0
    smineralizedc = 0
    semo2 = 0
    sco2_het = 0
  }
  
  #cat("timestep","NPP", "wt","sco2_het","semch4","semo2","smineralizedc", "semdifch4 ","semaerenchch4 ","semebullch4","Ts","Tmean","\n",file=outfile)  
  
  
  
  #check for negative cm and co
  
  if (min(cm) < 0) {
    cat("negative cm, stop \n")
    cat("negative cm, stop",outfile = outfile)
    break
  }
  if (min(co) < 0) {
    cat("negative co, stop \n")
    cat("negative co, stop",outfile = outfile)
    break
  }
  if (i*timestep > (duration - years)) {
    if (j > (1/timestep)/365) {
      output.df[k,] <- c(emdifch4 ,emaerenchch4, emebull, wt)
      output.mconsump$value[((k - 1)*length(depth)):((k)*length(depth) - 1)] <- mconsump
      output.mconsump$day[((k - 1)*length(depth)):((k)*length(depth) - 1)] <- i
      output.mprod[,k] <- mprod
      j <- 0
       k <- k + 1
    }
    j <- j + 1
  }

  
 


  
  
 # cbind(i*timestep,Ts)
  #cat(as.matrix(c(i*timestep,Ts)))
} #end of loop over time
# cat("\t","timestep","NPP", "wt","sco2_het","semch4","semo2","smineralizedc", "semdifch4 ","semaerenchch4 ","semebullch4", "\n")  
# print(c(timestep,NPP,wt,emch4,emdifch4,emaerenchch4,emebull,sum(mineralizedc*ddepth),NPP - sum(mineralizedc*ddepth) - emch4))
# return(c(timestep,NPP,wt,emch4,emdifch4,emaerenchch4,emebull,sum(mineralizedc*ddepth),NPP - sum(mineralizedc*ddepth) - emch4))
 #print(semch4)
# 
# head(output.df2)
# output.df2 <- output.df[complete.cases(output.df),]
# output.mconsump2 <- output.mconsump[,1:(dim(output.mconsump)[2]-1)]
# output.mprod2 <- output.mprod[,1:(dim(output.mprod)[2]-1)]
# 

#cat(i*timestep, Ts,file = outfile2,append = T)
#cat(as.matrix(c(i*timestep,Ts)))
 
#save.image("~/Desktop/testscript.RData")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#source("Multiplot.R")
# p1 <- ggplot(data = output.df2, aes(x = 1:(dim(output.df2)[1]), y = output.df2$emdifch4)) + 
#   geom_line() + 
#   ggtitle("emdifch4")
# 
# p2 <- ggplot(data = output.df2, aes(x = 1:(dim(output.df2)[1]), y = output.df2$emaerenchch4)) +
#   geom_line() +
#   ggtitle("emaerenchch4 ")
# 
# p3 <- ggplot(data = output.df2, aes(x = 1:(dim(output.df2)[1]), y = output.df2$emebull)) +
#   geom_line() +
#   ggtitle("emebull ")
# 
# p4 <- ggplot(data = output.df2, aes(x = 1:(dim(output.df2)[1]), y = output.df2$wt)) +
#   geom_line() +
#   ggtitle("wt")
# 
# 
# multiplot(p1,p2,p3,p4,cols=2)
# output.mconsump2[,] <- 1:17
# 
# p5 <- ggplot(data = output.mconsump, aes(x = day , y = value, colour=as.factor(level))) 
# p5 +  geom_line()+
#   ggtitle("MConsump")
# 
# p6 <- ggplot(data = output.mconsump[output.mconsump$level==1,], aes(x = day, y = value)) 
# p6 +  geom_line() + 
#   ggtitle("Level 2")
# 
# p7 <- ggplot(data = output.mconsump[output.mconsump$level==2,], aes(x = day, y = value)) 
# p7 <- p7 +  geom_line() + 
#   ggtitle("Level 2")
# 
# p8 <- ggplot(data = output.mconsump[output.mconsump$level==3,], aes(x = day, y = value)) 
# p8 <- p8 +  geom_line() + 
#   ggtitle("Level 3")
# 
# multiplot(p5,p6,p7,p8,cols=2)


#Tnew = read.table("Ts.dat",sep='',skip=1)
# dim(Tnew)


#Tnew = read.table("Ts.dat",sep='')


                   