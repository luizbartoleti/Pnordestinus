## Script for data simulation - ABC - P. nordestinus - Scenario 1 - Riverine Barrier Hypothesis

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html


## set working directory
#setwd("/home/elen/Área de Trabalho/abc p_nordestinus/1_simulations/Scenario1")

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
scenario1_parameters <- data.frame()



### model1: vicariance among groups North, South and West

scenario1_model1 <- "scenario1_model1"

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
  
  ## divergence time prior following an uniform distribution from 0.12 to 5.33 million years ago
  DivTime <- runif(1, 120000, 5333000)
  ## divergence time prior following an uniform distribution from present to DivTime (posterior colonization event, absent in model 1)
  DivTime2 <- 0
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  coalDivTime <- DivTime/(genlen*4*Ne)
  coalDivTime2 <- DivTime2/(genlen*4*Ne)
  
  ## bottleneck intensity in group West, measured as the ratio of the Ne in the original population (no bottleneck in model 1)
  Bot <- 1


### ms's command
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, coalDivTime, coalDivTime, scenario1_model1))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, coalDivTime, coalDivTime, scenario1_model1))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, coalDivTime, coalDivTime, scenario1_model1))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, coalDivTime, coalDivTime, scenario1_model1))

## save parameter values
scenario1_parameters <- rbind(scenario1_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTime, DivTime2, Bot))
}



### model2: vicariance between groups North and South, colonization of West by individuals from South

scenario1_model2 <- "scenario1_model2"

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
  
  ## divergence time prior following an uniform distribution from 0.12 to 5.33 million years ago
  DivTime <- runif(1, 120000, 5333000)
  ## divergence time prior following an uniform distribution from present to DivTime (posterior colonization event)
  DivTime2 <- runif(1, 0, DivTime)
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  coalDivTime <- DivTime/(genlen*4*Ne)
  coalDivTime2 <- DivTime2/(genlen*4*Ne)
  
  ## bottleneck intensity in group West, measured as the ratio of the Ne in the original population (prior following an uniform distribution from 1 to 20% of the Ne)
  Bot <- runif(1, 0.01, 0.2)
  
  
### ms's command
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model2))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model2))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model2))
system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model2))
  
## save parameter values
scenario1_parameters <- rbind(scenario1_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTime, DivTime2, Bot))
}



## model3: vicariance between groups North and South, colonization of West by individuals from North

scenario1_model3 <- "scenario1_model3"

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
  
  ## divergence time prior following an uniform distribution from 0.12 to 5.33 million years ago
  DivTime <- runif(1, 120000, 5333000)
  ## divergence time prior following an uniform distribution from present to DivTime (posterior colonization event)
  DivTime2 <- runif(1, 0, DivTime)
  
  ## use the DivTime in years to calculte divergence time in colaescent units (required by ms)
  coalDivTime <- DivTime/(genlen*4*Ne)
  coalDivTime2 <- DivTime2/(genlen*4*Ne)
  
  ## bottleneck intensity in group West, measured as the ratio of the Ne in the original population (prior following an uniform distribution from 1 to 20% of the Ne)
  Bot <- runif(1, 0.01, 0.2)
  
  
  ### ms's command
  system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> %s_ND2.txt", Nsam_ND2, theta_ND2, PopN_ND2, PopS_ND2, PopW_ND2, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model3))
  system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> %s_16S.txt", Nsam_16S, theta_16S, PopN_16S, PopS_16S, PopW_16S, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model3))
  system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> %s_SiaH.txt", Nsam_SiaH, theta_SiaH, PopN_SiaH, PopS_SiaH, PopW_SiaH, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model3))
  system(sprintf("./ms %d 1 -t %f -I 3 %d %d %d -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> %s_Rhodo.txt", Nsam_Rhodo, theta_Rhodo, PopN_Rhodo, PopS_Rhodo, PopW_Rhodo, coalDivTime2, Bot, coalDivTime2, coalDivTime, scenario1_model3))
  
  ## save parameter values
  scenario1_parameters <- rbind(scenario1_parameters, data.frame(Ne, theta_ND2, theta_16S, theta_SiaH, theta_Rhodo, DivTime, DivTime2, Bot))
}



### export summary statistics (pi_total, ss, D, pi_N, pi_S, pi_W, pi_NS, pi_NW, pi_SW)
#setwd("/home/elen/Área de Trabalho/abc p_nordestinus/1_simulations/Scenario1/temp")
scenario1_model1_ND2 <- read.table("scenario1_model1_ND2.txt")
scenario1_model1_ND2 <- data.frame(pi_total.ND2=(scenario1_model1_ND2[,1]),
                        ss.ND2=(scenario1_model1_ND2[,2]),
                        D.ND2=(scenario1_model1_ND2[,3]),
                        pi_N.ND2=(scenario1_model1_ND2[,6]),
                        pi_S.ND2=(scenario1_model1_ND2[,7]),
                        pi_W.ND2=(scenario1_model1_ND2[,8]),
                        pi_NS.ND2=(scenario1_model1_ND2[,9]),
                        pi_NW.ND2=(scenario1_model1_ND2[,10]),
                        pi_SW.ND2=(scenario1_model1_ND2[,11]))
scenario1_model1_16S <- read.table("scenario1_model1_16S.txt")
scenario1_model1_16S <- data.frame(pi_total.16S=(scenario1_model1_16S[,1]),
                        ss.16S=(scenario1_model1_16S[,2]),
                        D.16S=(scenario1_model1_16S[,3]),
                        pi_N.16S=(scenario1_model1_16S[,6]),
                        pi_S.16S=(scenario1_model1_16S[,7]),
                        pi_W.16S=(scenario1_model1_16S[,8]),
                        pi_NS.16S=(scenario1_model1_16S[,9]),
                        pi_NW.16S=(scenario1_model1_16S[,10]),
                        pi_SW.16S=(scenario1_model1_16S[,11]))
scenario1_model1_SiaH <- read.table("scenario1_model1_SiaH.txt")
scenario1_model1_SiaH <- data.frame(pi_total.SiaH=(scenario1_model1_SiaH[,1]),
                         ss.SiaH=(scenario1_model1_SiaH[,2]),
                         D.SiaH=(scenario1_model1_SiaH[,3]),
                         pi_N.SiaH=(scenario1_model1_SiaH[,6]),
                         pi_S.SiaH=(scenario1_model1_SiaH[,7]),
                         pi_W.SiaH=(scenario1_model1_SiaH[,8]),
                         pi_NS.SiaH=(scenario1_model1_SiaH[,9]),
                         pi_NW.SiaH=(scenario1_model1_SiaH[,10]),
                         pi_SW.SiaH=(scenario1_model1_SiaH[,11]))
scenario1_model1_Rhodo <- read.table("scenario1_model1_Rhodo.txt")
scenario1_model1_Rhodo <- data.frame(pi_total.Rhodo=(scenario1_model1_Rhodo[,1]),
                          ss.Rhodo=(scenario1_model1_Rhodo[,2]),
                          D.Rhodo=(scenario1_model1_Rhodo[,3]),
                          pi_N.Rhodo=(scenario1_model1_Rhodo[,6]),
                          pi_S.Rhodo=(scenario1_model1_Rhodo[,7]),
                          pi_W.Rhodo=(scenario1_model1_Rhodo[,8]),
                          pi_NS.Rhodo=(scenario1_model1_Rhodo[,9]),
                          pi_NW.Rhodo=(scenario1_model1_Rhodo[,10]),
                          pi_SW.Rhodo=(scenario1_model1_Rhodo[,11]))

scenario1_model2_ND2 <- read.table("scenario1_model2_ND2.txt")
scenario1_model2_ND2 <- data.frame(pi_total.ND2=(scenario1_model2_ND2[,1]),
                         ss.ND2=(scenario1_model2_ND2[,2]),
                         D.ND2=(scenario1_model2_ND2[,3]),
                         pi_N.ND2=(scenario1_model2_ND2[,6]),
                         pi_S.ND2=(scenario1_model2_ND2[,7]),
                         pi_W.ND2=(scenario1_model2_ND2[,8]),
                         pi_NS.ND2=(scenario1_model2_ND2[,9]),
                         pi_NW.ND2=(scenario1_model2_ND2[,10]),
                         pi_SW.ND2=(scenario1_model2_ND2[,11]))
scenario1_model2_16S <- read.table("scenario1_model2_16S.txt")
scenario1_model2_16S <- data.frame(pi_total.16S=(scenario1_model2_16S[,1]),
                         ss.16S=(scenario1_model2_16S[,2]),
                         D.16S=(scenario1_model2_16S[,3]),
                         pi_N.16S=(scenario1_model2_16S[,6]),
                         pi_S.16S=(scenario1_model2_16S[,7]),
                         pi_W.16S=(scenario1_model2_16S[,8]),
                         pi_NS.16S=(scenario1_model2_16S[,9]),
                         pi_NW.16S=(scenario1_model2_16S[,10]),
                         pi_SW.16S=(scenario1_model2_16S[,11]))
scenario1_model2_SiaH <- read.table("scenario1_model2_SiaH.txt")
scenario1_model2_SiaH <- data.frame(pi_total.SiaH=(scenario1_model2_SiaH[,1]),
                          ss.SiaH=(scenario1_model2_SiaH[,2]),
                          D.SiaH=(scenario1_model2_SiaH[,3]),
                          pi_N.SiaH=(scenario1_model2_SiaH[,6]),
                          pi_S.SiaH=(scenario1_model2_SiaH[,7]),
                          pi_W.SiaH=(scenario1_model2_SiaH[,8]),
                          pi_NS.SiaH=(scenario1_model2_SiaH[,9]),
                          pi_NW.SiaH=(scenario1_model2_SiaH[,10]),
                          pi_SW.SiaH=(scenario1_model2_SiaH[,11]))
scenario1_model2_Rhodo <- read.table("scenario1_model2_Rhodo.txt")
scenario1_model2_Rhodo <- data.frame(pi_total.Rhodo=(scenario1_model2_Rhodo[,1]),
                           ss.Rhodo=(scenario1_model2_Rhodo[,2]),
                           D.Rhodo=(scenario1_model2_Rhodo[,3]),
                           pi_N.Rhodo=(scenario1_model2_Rhodo[,6]),
                           pi_S.Rhodo=(scenario1_model2_Rhodo[,7]),
                           pi_W.Rhodo=(scenario1_model2_Rhodo[,8]),
                           pi_NS.Rhodo=(scenario1_model2_Rhodo[,9]),
                           pi_NW.Rhodo=(scenario1_model2_Rhodo[,10]),
                           pi_SW.Rhodo=(scenario1_model2_Rhodo[,11]))

scenario1_model3_ND2 <- read.table("scenario1_model3_ND2.txt")
scenario1_model3_ND2 <- data.frame(pi_total.ND2=(scenario1_model3_ND2[,1]),
                         ss.ND2=(scenario1_model3_ND2[,2]),
                         D.ND2=(scenario1_model3_ND2[,3]),
                         pi_N.ND2=(scenario1_model3_ND2[,6]),
                         pi_S.ND2=(scenario1_model3_ND2[,7]),
                         pi_W.ND2=(scenario1_model3_ND2[,8]),
                         pi_NS.ND2=(scenario1_model3_ND2[,9]),
                         pi_NW.ND2=(scenario1_model3_ND2[,10]),
                         pi_SW.ND2=(scenario1_model3_ND2[,11]))
scenario1_model3_16S <- read.table("scenario1_model3_16S.txt")
scenario1_model3_16S <- data.frame(pi_total.16S=(scenario1_model3_16S[,1]),
                         ss.16S=(scenario1_model3_16S[,2]),
                         D.16S=(scenario1_model3_16S[,3]),
                         pi_N.16S=(scenario1_model3_16S[,6]),
                         pi_S.16S=(scenario1_model3_16S[,7]),
                         pi_W.16S=(scenario1_model3_16S[,8]),
                         pi_NS.16S=(scenario1_model3_16S[,9]),
                         pi_NW.16S=(scenario1_model3_16S[,10]),
                         pi_SW.16S=(scenario1_model3_16S[,11]))
scenario1_model3_SiaH <- read.table("scenario1_model3_SiaH.txt")
scenario1_model3_SiaH <- data.frame(pi_total.SiaH=(scenario1_model3_SiaH[,1]),
                          ss.SiaH=(scenario1_model3_SiaH[,2]),
                          D.SiaH=(scenario1_model3_SiaH[,3]),
                          pi_N.SiaH=(scenario1_model3_SiaH[,6]),
                          pi_S.SiaH=(scenario1_model3_SiaH[,7]),
                          pi_W.SiaH=(scenario1_model3_SiaH[,8]),
                          pi_NS.SiaH=(scenario1_model3_SiaH[,9]),
                          pi_NW.SiaH=(scenario1_model3_SiaH[,10]),
                          pi_SW.SiaH=(scenario1_model3_SiaH[,11]))
scenario1_model3_Rhodo <- read.table("scenario1_model3_Rhodo.txt")
scenario1_model3_Rhodo <- data.frame(pi_total.Rhodo=(scenario1_model3_Rhodo[,1]),
                           ss.Rhodo=(scenario1_model3_Rhodo[,2]),
                           D.Rhodo=(scenario1_model3_Rhodo[,3]),
                           pi_N.Rhodo=(scenario1_model3_Rhodo[,6]),
                           pi_S.Rhodo=(scenario1_model3_Rhodo[,7]),
                           pi_W.Rhodo=(scenario1_model3_Rhodo[,8]),
                           pi_NS.Rhodo=(scenario1_model3_Rhodo[,9]),
                           pi_NW.Rhodo=(scenario1_model3_Rhodo[,10]),
                           pi_SW.Rhodo=(scenario1_model3_Rhodo[,11]))


scenario1_models <- rep(c("scenario1_model1", "scenario1_model2", "scenario1_model3"), each=numsim)

## join data
scenario1_model1<-cbind(scenario1_model1_ND2, scenario1_model1_16S, scenario1_model1_SiaH, scenario1_model1_Rhodo)
scenario1_model2<-cbind(scenario1_model2_ND2, scenario1_model2_16S, scenario1_model2_SiaH, scenario1_model2_Rhodo)
scenario1_model3<-cbind(scenario1_model3_ND2, scenario1_model3_16S, scenario1_model3_SiaH, scenario1_model3_Rhodo)
scenario1_sust <- rbind(scenario1_model1, scenario1_model2, scenario1_model3)
names(scenario1_sust) <- c("pi_total.ND2", "ss.ND2", "D.ND2", "pi_N.ND2", "pi_S.ND2", "pi_W.ND2", "pi_NS.ND2", "pi_NW.ND2", "pi_SW.ND2", 
                           "pi_total.16S", "ss.16S", "D.16S", "pi_N.16S", "pi_S.16S", "pi_W.16S", "pi_NS.16S", "pi_NW.16S", "pi_SW.16S", 
                           "pi_total.SiaH", "ss.SiaH", "D.SiaH", "pi_N.SiaH", "pi_S.SiaH", "pi_W.SiaH", "pi_NS.SiaH", "pi_NW.SiaH", "pi_SW.SiaH", 
                           "pi_total.Rhodo", "ss.Rhodo", "D.Rhodo", "pi_N.Rhodo", "pi_S.Rhodo", "pi_W.Rhodo", "pi_NS.Rhodo", "pi_NW.Rhodo", "pi_SW.Rhodo")

##optional lines to save the files in the folder, not only in the R environment
write.table(scenario1_models, file="scenario1_models.txt", quote=F, row.names=F, col.names=F)
write.table(scenario1_parameters, file="scenario1_parameters.txt", quote=F, row.names=F, col.names=F)
write.table(scenario1_sust, file="scenario1_sust.txt", quote=F, row.names=F, col.names=F)

## erase temporary files
system(sprintf("rm -rf scenario1_model1_ND2.txt scenario1_model1_16S.txt scenario1_model1_SiaH.txt scenario1_model1_Rhodo.txt scenario1_model2_ND2.txt scenario1_model2_16S.txt scenario1_model2_SiaH.txt scenario1_model2_Rhodo.txt scenario1_model3_ND2.txt scenario1_model3_16S.txt scenario1_model3_SiaH.txt scenario1_model3_Rhodo.txt"))
