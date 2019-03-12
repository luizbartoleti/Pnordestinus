## ABC data simulation for P. nordestinus - Scenario 2 - LGM Refugia Hypothesis

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html


## set working directory
#setwd("~/elen")

## make directory for temporary files
#system(sprintf("mkdir temp"))

## number of simulations
numsim <- 100000


### variable declarations

  ## total sample sizes
  Nsam_ND2 = 168
  Nsam_16S = 151
  Nsam_SiaH = 155*2
  Nsam_Rhodo = 189*2

  ## sample sizes of group North
  PopN_ND2 = 131
  PopN_16S = 129
  PopN_SiaH = 130*2
  PopN_Rhodo = 146*2
  
  ## sample sizes of group South
  PopS_ND2 = 34
  PopS_16S = 18
  PopS_SiaH = 23*2
  PopS_Rhodo = 41*2

  ## sample sizes of group West
  PopW_ND2 = 3
  PopW_16S = 4
  PopW_SiaH = 2*2
  PopW_Rhodo = 2*2
    
  ## number of base pairs
  L_ND2 <- 813
  L_16S <- 1035
  L_SiaH <- 361
  L_Rhodo <- 378
  
  ## number of years per generation
  genlen <- 1


### create a variable to store all the parameters
scenario2_parameters <- data.frame()



### model4: single refugium in the South, with posterior colonization of North and West

model4 <- "model4"

for (i in 1:numsim) {

### Define parameters

  ## mutation rates
  mutrate_ND2 <- runif(1, 5*10^(-9), 1.5*10^(-8))  #95% HPD obtained in *BEAST analysis
  mutrate_16S <- runif(1, 2*10^(-9), 5*10^(-9))  #95% HPD obtained in *BEAST analysis
  mutrate_SiaH <- runif(1, 3*10^(-10), 1.5*10^(-9)) #95% HPD obtained in *BEAST analysis
  mutrate_Rhodo <- runif(1, 1*10^(-10), 1*10^(-9)) #95% HPD obtained in *BEAST analysis
 
  ## Ne distribution based on mitochondrial mutation rates and theta-W intervals calculated in dnasp (theta_ND2 = 21 - 52, theta_16S = 9 - 24)
  Ne <- runif(1, 2500000, 9000000)
  
  ## use a single Ne value to obtain theta values for all markers in each simulation
  theta_ND2 <- Ne*mutrate_ND2*L_ND2
  theta_16S <- Ne*mutrate_16S*L_16S
  theta_SiaH <- 4*Ne*mutrate_SiaH*L_SiaH
  theta_Rhodo <- 4*Ne*mutrate_Rhodo*L_Rhodo
  
  ## divergence time priors following uniform distributions
  DivTimeNS_0_LIG <- runif(1, 0, 120000) #split between N and S from present to LIG
  DivTimeSW_0_LIG <- runif(1, 0, 120000) #split between S and W from present to LIG
  DivTimeSW_0_LGM <- 0 #split between S and W from present to LGM (absent in model4)
  DivTimeNS_LGM_LIG <- 0 #split between N and S from LGM to LIG (absent in model4)
  DivTimeNW_0_LGM <- 0 #split between N and W from present to LGM (absent in model4)
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  coalDivTimeLIG <- 120000/(genlen*4*Ne)
  #coalDivTimeLGM <- 21000/(genlen*4*Ne)
  coalDivTimeNS_0_LIG <- DivTimeNS_0_LIG/(genlen*4*Ne)
  coalDivTimeSW_0_LIG <- DivTimeSW_0_LIG/(genlen*4*Ne)
  #coalDivTimeSW_0_LGM <- DivTimeSW_0_LGM/(genlen*4*Ne)
  #coalDivTimeNS_LGM_LIG <- DivTimeNS_LGM_LIG/(genlen*4*Ne)
  #coalDivTimeNW_0_LGM <- DivTimeNW_0_LGM/(genlen*4*Ne)
  
  # exponential growth rate in each group
  GrowthRatio <- runif(1, 0.25, 0.9) # ratio of Ne before expansion/Ne after expansion (expansion of 10% - 75%)
  GrowthRateN1 <- -(1/coalDivTimeNS_0_LIG)*log(GrowthRatio)
  GrowthRateS1 <- -(1/coalDivTimeLIG)*log(GrowthRatio)
  GrowthRateW1 <- -(1/coalDivTimeSW_0_LIG)*log(GrowthRatio)
  GrowthRateN2 <- 0
  GrowthRateS2 <- 0
  GrowthRateW2 <- 0
  GrowthRateW3 <- 0
  
  ## bottleneck intensity in each group, measured as the ratio of the Ne in the original population
  BotN <- runif(1, 0.01, 0.2)
  BotW <- runif(1, 0.01, 0.2)


### ms's command
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -ej %f 3 2 -ej %f 1 2 -eg %f 2 0 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, GrowthRateN1, GrowthRateS1, GrowthRateW1, coalDivTimeNS_0_LIG, BotN, coalDivTimeSW_0_LIG, BotW, coalDivTimeSW_0_LIG, coalDivTimeNS_0_LIG, coalDivTimeLIG, model4))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -ej %f 3 2 -ej %f 1 2 -eg %f 2 0 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, GrowthRateN1, GrowthRateS1, GrowthRateW1, coalDivTimeNS_0_LIG, BotN, coalDivTimeSW_0_LIG, BotW, coalDivTimeSW_0_LIG, coalDivTimeNS_0_LIG, coalDivTimeLIG, model4))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -ej %f 3 2 -ej %f 1 2 -eg %f 2 0 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, GrowthRateN1, GrowthRateS1, GrowthRateW1, coalDivTimeNS_0_LIG, BotN, coalDivTimeSW_0_LIG, BotW, coalDivTimeSW_0_LIG, coalDivTimeNS_0_LIG, coalDivTimeLIG, model4))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -ej %f 3 2 -ej %f 1 2 -eg %f 2 0 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, GrowthRateN1, GrowthRateS1, GrowthRateW1, coalDivTimeNS_0_LIG, BotN, coalDivTimeSW_0_LIG, BotW, coalDivTimeSW_0_LIG, coalDivTimeNS_0_LIG, coalDivTimeLIG, model4))

## save parameter values
scenario2_parameters <- rbind(scenario2_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTimeNS_0_LIG, DivTimeSW_0_LIG, DivTimeSW_0_LGM, DivTimeNS_LGM_LIG, DivTimeNW_0_LGM, GrowthRateN1, GrowthRateS1, GrowthRateW1, GrowthRateN2, GrowthRateS2, GrowthRateW2, GrowthRateW3, BotN, BotW))
}


### model5: North and South refugia, with posterior colonization of West by individuals from South

model5 <- "model5"

for (i in 1:numsim) {
  
  ## mutation rates
  mutrate_ND2 <- runif(1, 5*10^(-9), 1.5*10^(-8))  #95% HPD obtained in *BEAST analysis
  mutrate_16S <- runif(1, 2*10^(-9), 5*10^(-9))  #95% HPD obtained in *BEAST analysis
  mutrate_SiaH <- runif(1, 3*10^(-10), 1.5*10^(-9)) #95% HPD obtained in *BEAST analysis
  mutrate_Rhodo <- runif(1, 1*10^(-10), 1*10^(-9)) #95% HPD obtained in *BEAST analysis
  
  ## Ne distribution based on mitochondrial mutation rates and theta-W intervals calculated in dnasp (theta_ND2 = 21 - 52, theta_16S = 9 - 24)
  Ne <- runif(1, 2500000, 9000000)
  
  ## use a single Ne value to obtain theta values for all markers in each simulation
  theta_ND2 <- Ne*mutrate_ND2*L_ND2
  theta_16S <- Ne*mutrate_16S*L_16S
  theta_SiaH <- 4*Ne*mutrate_SiaH*L_SiaH
  theta_Rhodo <- 4*Ne*mutrate_Rhodo*L_Rhodo
  
  ## divergence time priors following uniform distributions
  DivTimeNS_0_LIG <- 0 #split between N and S from present to LIG (absent in model5)
  DivTimeSW_0_LIG <- 0 #split between S and W from present to LIG (absent in model5)
  DivTimeSW_0_LGM <- runif(1, 0, 21000) #split between S and W from present to LGM
  DivTimeNS_LGM_LIG <- runif(1, 21000, 120000) #split between N and S from LGM to LIG
  DivTimeNW_0_LGM <- 0 #split between N and W from present to LGM (absent in model5)
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  #coalDivTimeLIG <- 120000/(genlen*4*Ne)
  coalDivTimeLGM <- 21000/(genlen*4*Ne)
  #coalDivTimeNS_0_LIG <- DivTimeNS_0_LIG/(genlen*4*Ne)
  #coalDivTimeSW_0_LIG <- DivTimeSW_0_LIG/(genlen*4*Ne)
  coalDivTimeSW_0_LGM <- DivTimeSW_0_LGM/(genlen*4*Ne)
  coalDivTimeNS_LGM_LIG <- DivTimeNS_LGM_LIG/(genlen*4*Ne)
  #coalDivTimeNW_0_LGM <- DivTimeNW_0_LGM/(genlen*4*Ne)
  
  # exponential growth rate in each group
  GrowthRatio <- runif(1, 0.25, 0.9) # ratio of Ne before expansion/Ne after expansion (expansion of 10% - 75%)
  GrowthRateN1 <- 0
  GrowthRateW1 <- 0
  GrowthRateS1 <- 0
  GrowthRateN2 <- -(1/coalDivTimeLGM)*log(GrowthRatio)
  GrowthRateW2 <- -(1/coalDivTimeSW_0_LGM)*log(GrowthRatio)
  GrowthRateS2 <- -(1/coalDivTimeLGM)*log(GrowthRatio)
  GrowthRateW3 <- 0
 
  ## bottleneck intensity in each group, measured as the ratio of the Ne in the original population
  BotN <- 1
  BotW <- runif(1, 0.01, 0.2)
  
  
### ms's command
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 2 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, GrowthRateN2, GrowthRateS2, GrowthRateW2, coalDivTimeSW_0_LGM, BotW, coalDivTimeSW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model5))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 2 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, GrowthRateN2, GrowthRateS2, GrowthRateW2, coalDivTimeSW_0_LGM, BotW, coalDivTimeSW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model5))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 2 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, GrowthRateN2, GrowthRateS2, GrowthRateW2, coalDivTimeSW_0_LGM, BotW, coalDivTimeSW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model5))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 2 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, GrowthRateN2, GrowthRateS2, GrowthRateW2, coalDivTimeSW_0_LGM, BotW, coalDivTimeSW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model5))

## save parameter values
scenario2_parameters <- rbind(scenario2_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTimeNS_0_LIG, DivTimeSW_0_LIG, DivTimeSW_0_LGM, DivTimeNS_LGM_LIG, DivTimeNW_0_LGM, GrowthRateN1, GrowthRateS1, GrowthRateW1, GrowthRateN2, GrowthRateS2, GrowthRateW2, GrowthRateW3, BotN, BotW))
}


## model6: North and South refugia, with posterior colonization of West by individuals from North

model6 <- "model6"

for (i in 1:numsim) {
  
  ## mutation rates
  mutrate_ND2 <- runif(1, 5*10^(-9), 1.5*10^(-8))  #95% HPD obtained in *BEAST analysis
  mutrate_16S <- runif(1, 2*10^(-9), 5*10^(-9))  #95% HPD obtained in *BEAST analysis
  mutrate_SiaH <- runif(1, 3*10^(-10), 1.5*10^(-9)) #95% HPD obtained in *BEAST analysis
  mutrate_Rhodo <- runif(1, 1*10^(-10), 1*10^(-9)) #95% HPD obtained in *BEAST analysis
  
  ## Ne distribution based on mitochondrial mutation rates and theta-W intervals calculated in dnasp (theta_ND2 = 21 - 52, theta_16S = 9 - 24)
  Ne <- runif(1, 2500000, 9000000)
  
  ## use a single Ne value to obtain theta values for all markers in each simulation
  theta_ND2 <- Ne*mutrate_ND2*L_ND2
  theta_16S <- Ne*mutrate_16S*L_16S
  theta_SiaH <- 4*Ne*mutrate_SiaH*L_SiaH
  theta_Rhodo <- 4*Ne*mutrate_Rhodo*L_Rhodo
  
  ## divergence time priors following uniform distributions
  DivTimeNS_0_LIG <- 0 #split between N and S from present to LIG (absent in model6)
  DivTimeSW_0_LIG <- 0 #split between S and W from present to LIG (absent in model6)
  DivTimeSW_0_LGM <- 0 #split between S and W from present to LGM (absent in model6)
  DivTimeNS_LGM_LIG <- runif(1, 21000, 120000) #split between N and S from LGM to LIG
  DivTimeNW_0_LGM <- runif(1, 0, 21000) #split between N and W from present to LGM
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  #coalDivTimeLIG <- 120000/(genlen*4*Ne)
  coalDivTimeLGM <- 21000/(genlen*4*Ne)
  #coalDivTimeNS_0_LIG <- DivTimeNS_0_LIG/(genlen*4*Ne)
  #coalDivTimeSW_0_LIG <- DivTimeSW_0_LIG/(genlen*4*Ne)
  #coalDivTimeSW_0_LGM <- DivTimeSW_0_LGM/(genlen*4*Ne)
  coalDivTimeNS_LGM_LIG <- DivTimeNS_LGM_LIG/(genlen*4*Ne)
  coalDivTimeNW_0_LGM <- DivTimeNW_0_LGM/(genlen*4*Ne)
  
  # exponential growth rate in each group
  GrowthRatio <- runif(1, 0.25, 0.9) # ratio of Ne before expansion/Ne after expansion (expansion of 10% - 75%)
  GrowthRateN1 <- 0
  GrowthRateW1 <- 0
  GrowthRateS1 <- 0
  GrowthRateN2 <- -(1/coalDivTimeLGM)*log(GrowthRatio)
  GrowthRateW2 <- 0
  GrowthRateS2 <- -(1/coalDivTimeLGM)*log(GrowthRatio)
  GrowthRateW3 <- -(1/coalDivTimeNW_0_LGM)*log(GrowthRatio)
  
  ## bottleneck intensity in each group, measured as the ratio of the Ne in the original population
  BotN <- 1
  BotW <- runif(1, 0.01, 0.2)
  
  
### ms's command
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 1 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, GrowthRateN2, GrowthRateS2, GrowthRateW3, coalDivTimeNW_0_LGM, BotW, coalDivTimeNW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model6))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 1 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, GrowthRateN2, GrowthRateS2, GrowthRateW3, coalDivTimeNW_0_LGM, BotW, coalDivTimeNW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model6))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 1 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, GrowthRateN2, GrowthRateS2, GrowthRateW3, coalDivTimeNW_0_LGM, BotW, coalDivTimeNW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model6))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -ej %f 3 1 -eg %f 1 0 -eg %f 2 0 -ej %f 2 1 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, GrowthRateN2, GrowthRateS2, GrowthRateW3, coalDivTimeNW_0_LGM, BotW, coalDivTimeNW_0_LGM, coalDivTimeLGM, coalDivTimeLGM, coalDivTimeNS_LGM_LIG, model6))

## save parameter values
scenario2_parameters <- rbind(scenario2_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTimeNS_0_LIG, DivTimeSW_0_LIG, DivTimeSW_0_LGM, DivTimeNS_LGM_LIG, DivTimeNW_0_LGM, GrowthRateN1, GrowthRateS1, GrowthRateW1, GrowthRateN2, GrowthRateS2, GrowthRateW2, GrowthRateW3, BotN, BotW))
}



### export summary statistics (pi_total, ss, D, pi_N, pi_S, pi_W, pi_NS, pi_NW, pi_SW)
#setwd("/home/elen/Área de Trabalho/abc p_nordestinus/1_simulations/Scenario2/temp")
model4_ND2 <- read.table("model4_ND2.txt")
model4_ND2 <- data.frame(pi_total.ND2=(model4_ND2[,1]),
                        ss.ND2=(model4_ND2[,2]),
                        D.ND2=(model4_ND2[,3]),
                        pi_N.ND2=(model4_ND2[,6]),
                        pi_S.ND2=(model4_ND2[,7]),
                        pi_W.ND2=(model4_ND2[,8]),
                        pi_NS.ND2=(model4_ND2[,9]),
                        pi_NW.ND2=(model4_ND2[,10]),
                        pi_SW.ND2=(model4_ND2[,11]))
model4_16S <- read.table("model4_16S.txt")
model4_16S <- data.frame(pi_total.16S=(model4_16S[,1]),
                        ss.16S=(model4_16S[,2]),
                        D.16S=(model4_16S[,3]),
                        pi_N.16S=(model4_16S[,6]),
                        pi_S.16S=(model4_16S[,7]),
                        pi_W.16S=(model4_16S[,8]),
                        pi_NS.16S=(model4_16S[,9]),
                        pi_NW.16S=(model4_16S[,10]),
                        pi_SW.16S=(model4_16S[,11]))
model4_SiaH <- read.table("model4_SiaH.txt")
model4_SiaH <- data.frame(pi_total.SiaH=(model4_SiaH[,1]),
                         ss.SiaH=(model4_SiaH[,2]),
                         D.SiaH=(model4_SiaH[,3]),
                         pi_N.SiaH=(model4_SiaH[,6]),
                         pi_S.SiaH=(model4_SiaH[,7]),
                         pi_W.SiaH=(model4_SiaH[,8]),
                         pi_NS.SiaH=(model4_SiaH[,9]),
                         pi_NW.SiaH=(model4_SiaH[,10]),
                         pi_SW.SiaH=(model4_SiaH[,11]))
model4_Rhodo <- read.table("model4_Rhodo.txt")
model4_Rhodo <- data.frame(pi_total.Rhodo=(model4_Rhodo[,1]),
                          ss.Rhodo=(model4_Rhodo[,2]),
                          D.Rhodo=(model4_Rhodo[,3]),
                          pi_N.Rhodo=(model4_Rhodo[,6]),
                          pi_S.Rhodo=(model4_Rhodo[,7]),
                          pi_W.Rhodo=(model4_Rhodo[,8]),
                          pi_NS.Rhodo=(model4_Rhodo[,9]),
                          pi_NW.Rhodo=(model4_Rhodo[,10]),
                          pi_SW.Rhodo=(model4_Rhodo[,11]))

model5_ND2 <- read.table("model5_ND2.txt")
model5_ND2 <- data.frame(pi_total.ND2=(model5_ND2[,1]),
                         ss.ND2=(model5_ND2[,2]),
                         D.ND2=(model5_ND2[,3]),
                         pi_N.ND2=(model5_ND2[,6]),
                         pi_S.ND2=(model5_ND2[,7]),
                         pi_W.ND2=(model5_ND2[,8]),
                         pi_NS.ND2=(model5_ND2[,9]),
                         pi_NW.ND2=(model5_ND2[,10]),
                         pi_SW.ND2=(model5_ND2[,11]))
model5_16S <- read.table("model5_16S.txt")
model5_16S <- data.frame(pi_total.16S=(model5_16S[,1]),
                         ss.16S=(model5_16S[,2]),
                         D.16S=(model5_16S[,3]),
                         pi_N.16S=(model5_16S[,6]),
                         pi_S.16S=(model5_16S[,7]),
                         pi_W.16S=(model5_16S[,8]),
                         pi_NS.16S=(model5_16S[,9]),
                         pi_NW.16S=(model5_16S[,10]),
                         pi_SW.16S=(model5_16S[,11]))
model5_SiaH <- read.table("model5_SiaH.txt")
model5_SiaH <- data.frame(pi_total.SiaH=(model5_SiaH[,1]),
                          ss.SiaH=(model5_SiaH[,2]),
                          D.SiaH=(model5_SiaH[,3]),
                          pi_N.SiaH=(model5_SiaH[,6]),
                          pi_S.SiaH=(model5_SiaH[,7]),
                          pi_W.SiaH=(model5_SiaH[,8]),
                          pi_NS.SiaH=(model5_SiaH[,9]),
                          pi_NW.SiaH=(model5_SiaH[,10]),
                          pi_SW.SiaH=(model5_SiaH[,11]))
model5_Rhodo <- read.table("model5_Rhodo.txt")
model5_Rhodo <- data.frame(pi_total.Rhodo=(model5_Rhodo[,1]),
                           ss.Rhodo=(model5_Rhodo[,2]),
                           D.Rhodo=(model5_Rhodo[,3]),
                           pi_N.Rhodo=(model5_Rhodo[,6]),
                           pi_S.Rhodo=(model5_Rhodo[,7]),
                           pi_W.Rhodo=(model5_Rhodo[,8]),
                           pi_NS.Rhodo=(model5_Rhodo[,9]),
                           pi_NW.Rhodo=(model5_Rhodo[,10]),
                           pi_SW.Rhodo=(model5_Rhodo[,11]))

model6_ND2 <- read.table("model6_ND2.txt")
model6_ND2 <- data.frame(pi_total.ND2=(model6_ND2[,1]),
                         ss.ND2=(model6_ND2[,2]),
                         D.ND2=(model6_ND2[,3]),
                         pi_N.ND2=(model6_ND2[,6]),
                         pi_S.ND2=(model6_ND2[,7]),
                         pi_W.ND2=(model6_ND2[,8]),
                         pi_NS.ND2=(model6_ND2[,9]),
                         pi_NW.ND2=(model6_ND2[,10]),
                         pi_SW.ND2=(model6_ND2[,11]))
model6_16S <- read.table("model6_16S.txt")
model6_16S <- data.frame(pi_total.16S=(model6_16S[,1]),
                         ss.16S=(model6_16S[,2]),
                         D.16S=(model6_16S[,3]),
                         pi_N.16S=(model6_16S[,6]),
                         pi_S.16S=(model6_16S[,7]),
                         pi_W.16S=(model6_16S[,8]),
                         pi_NS.16S=(model6_16S[,9]),
                         pi_NW.16S=(model6_16S[,10]),
                         pi_SW.16S=(model6_16S[,11]))
model6_SiaH <- read.table("model6_SiaH.txt")
model6_SiaH <- data.frame(pi_total.SiaH=(model6_SiaH[,1]),
                          ss.SiaH=(model6_SiaH[,2]),
                          D.SiaH=(model6_SiaH[,3]),
                          pi_N.SiaH=(model6_SiaH[,6]),
                          pi_S.SiaH=(model6_SiaH[,7]),
                          pi_W.SiaH=(model6_SiaH[,8]),
                          pi_NS.SiaH=(model6_SiaH[,9]),
                          pi_NW.SiaH=(model6_SiaH[,10]),
                          pi_SW.SiaH=(model6_SiaH[,11]))
model6_Rhodo <- read.table("model6_Rhodo.txt")
model6_Rhodo <- data.frame(pi_total.Rhodo=(model6_Rhodo[,1]),
                           ss.Rhodo=(model6_Rhodo[,2]),
                           D.Rhodo=(model6_Rhodo[,3]),
                           pi_N.Rhodo=(model6_Rhodo[,6]),
                           pi_S.Rhodo=(model6_Rhodo[,7]),
                           pi_W.Rhodo=(model6_Rhodo[,8]),
                           pi_NS.Rhodo=(model6_Rhodo[,9]),
                           pi_NW.Rhodo=(model6_Rhodo[,10]),
                           pi_SW.Rhodo=(model6_Rhodo[,11]))


scenario2_models <- rep(c("model4", "model5", "model6"), each=numsim)

## join data
model4<-cbind(model4_ND2, model4_16S, model4_SiaH, model4_Rhodo)
model5<-cbind(model5_ND2, model5_16S, model5_SiaH, model5_Rhodo)
model6<-cbind(model6_ND2, model6_16S, model6_SiaH, model6_Rhodo)
scenario2_sust <- rbind(model4, model5, model6)
names(scenario2_sust) <- c("pi_total.ND2", "ss.ND2", "D.ND2", "pi_N.ND2", "pi_S.ND2", "pi_W.ND2", "pi_NS.ND2", "pi_NW.ND2", "pi_SW.ND2", 
                           "pi_total.16S", "ss.16S", "D.16S", "pi_N.16S", "pi_S.16S", "pi_W.16S", "pi_NS.16S", "pi_NW.16S", "pi_SW.16S", 
                           "pi_total.SiaH", "ss.SiaH", "D.SiaH", "pi_N.SiaH", "pi_S.SiaH", "pi_W.SiaH", "pi_NS.SiaH", "pi_NW.SiaH", "pi_SW.SiaH", 
                           "pi_total.Rhodo", "ss.Rhodo", "D.Rhodo", "pi_N.Rhodo", "pi_S.Rhodo", "pi_W.Rhodo", "pi_NS.Rhodo", "pi_NW.Rhodo", "pi_SW.Rhodo")

##optional lines to save the files in the folder, not only in the R environment
write.table(scenario2_models, file="scenario2_models.txt", quote=F, row.names=F, col.names=F)
write.table(scenario2_parameters, file="scenario2_parameters.txt", quote=F, row.names=F, col.names=F)
write.table(scenario2_sust, file="scenario2_sust.txt", quote=F, row.names=F, col.names=F)

## erase temporary files
system(sprintf("rm -rf model4_ND2.txt model4_16S.txt model4_SiaH.txt model4_Rhodo.txt model5_ND2.txt model5_16S.txt model5_SiaH.txt model5_Rhodo.txt model6_ND2.txt model6_16S.txt model6_SiaH.txt model6_Rhodo.txt"))
