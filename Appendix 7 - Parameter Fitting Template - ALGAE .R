#This script executes Monod-type microbial kinetics using deSolve 
library(deSolve)
library(FME)



Period=10 #number of days
N=100

Microbial <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { 
    
    #model parameters
    
    ualg <- 1.6 #d^-1
    kresp <- 0.1 #d^-1
    kdeath <- 0.1 #d^-1
    KC <- 0.00432 #gC m^-3
    ICO2 <- 120 #gC m^-3
    KN <- 0.1 #gC m^-3
    KO2 <- 0.2 #gO2 m^-3
    Kpr <- 0.01 
    Tau <- 2 
    SsatO2 <- 7.1904 #gO2 m^-3
    Topt <- 25 #degrees celsius
    s <- 13
    alpha <- 0.001935 #(uE m^-2)^-1
    beta <- 0.00000057848 #(uE m^-2)^-1
    gamma <- 0.1460 #s^-1
    delta <- 0.0004796 #s^-1
    Ef <- 1.74 #uE J^-1
    N <- 0.74
    rho <- 1353 #W m^-2
    KaO2 <- 4 #d^-1
    KaCO2 <- 0.7 #d^-1
    KaNH3 <- 0.7 #d^-1
    Tact <- 13.5 
    Keq1 <- 1 #10^(17.843-(3404.71/(273.15+Tact))-0.032786*(273.15+Tact))
    Keq2 <- 1 #10^(9.494-(2902.39/(273.15+Tact))-0.02379*(273.15+Tact))
    Keq3 <- 1 #10^(2.891-(2727/(273.15+Tact)))
    Keqw <- 1 #10^(-(4470.99/(273.15+Tact))+12.0875-0.01706*(273.15+Tact))
    keq1 <- .1 #d^-1
    keq2 <- .1 #d^-1
    keq3 <- .1 #d^-1
    keqw <- .1 #gm^-1 d^-1
    iC <- 0.387 #gC gCOD^-1
    iH <- 0.075 #gH gCOD^-1
    iO <- 0.538 #gO2 gCOD^-1
    iN <- 0.065 #gN gCOD^-1
    I <- 200 
    SwatO2 <- 10.5 
    SwatCO2 <- 220 
    SwatNH3 <- 16 
    #SO2 <- 8
    
    
    #stoichiometric coefficients
    v11a = -iN   #gN gCOD^-1
    v41a = (8/3)*iC+(8*iH)-iO-(12/7)*iN #gO2 gCOD^-1
    v51a = -iC  #gC gCOD^-1
    v81a = iN/14 #gH gCOD^-1
    v101a = 1 #gCOD gCOD^-1
    v31b = -iN #gN gCOD^-1
    v41b =(8/3)*iC+(8*iH)-iO-(20/7)*iN #gO2 gCOD^-1
    v51b = -iC  #gC gCOD^-1
    v81b = -iN/14 #gH gCOD^-1
    v101b = 1 #gCOD gCOD^-1
    v12 = iN #gN gCOD^-1
    v42 = (iO)-8*(iH)-(8/3)*(iC)+(12/7)*(iN) #gO2 gCOD^-1
    v52 = iC #gC gCOD^-1
    v82 = -(1/14)*(iN) #gH gCOD^-1
    v102 = -1 #gCOD gCOD^-1
    v13 = iN #gN gCOD^-1
    v43 = (iO)-8*(iH)-(8/3)*(iC)+(12/7)*(iN) #gO2 gCOD^-1
    v53 = iC #gC gCOD^-1
    v83 = -(1/14)*(iN) #gH gCOD^-1
    v103 = -1 #gCOD gCOD^-1
    v54 = -1  #WAS NEGATIVE #gC gC^-1
    v64 = 1 #gC gC^-1
    v84 = 1/12 #gH gC^-1
    v65 = -1 #gC gC^-1
    v75 = 1 #gC gC^-1
    v85 = 1/12 #gH gC^-1
    v16 = -1 #WAS NEGATIVE #gN gN^-1
    v26 = 1 #gN gN^-1
    v86 = 1/14 #gH gN^-1
    v87 = 1 #gH gH^-1
    v97 = 1 #gH gH^-1
    v4O2 = 1
    v5CO2 = 1
    v2NH3 = 1
    
    
    #Process rates
    fL <- (alpha*delta*I)/((alpha*beta*((I)^2))+((alpha+beta)*delta*I)+(gamma*delta))
    #fL = 1
    #fPR <- if (SO2>(Tau*SsatO2)) {
    #0} else {
    #1-tanh((Kpr*(SO2/(Tau*SsatO2)))/(1-(SO2/(Tau*SsatO2))))}
    fPR = 0.5
    nps <- fL*fPR
    #nps = 1
    fT <- exp(-((Tact-Topt)/s)^2)
    #fT = 1
    
    
    p1a = ualg*fT*nps*((SCO2+SHCO3)/(KC+SCO2+SHCO3+(((SCO2)^2)/ICO2)))*((SNH3+SNH4)/(KN+SNH3+SNH4))*Xalg #Growth on ammonia
    p1b = ualg*fT*nps*((SCO2+SHCO3)/(KC+SCO2+SHCO3+(((SCO2)^2)/ICO2)))*(SNO3/(KN+SNO3))*(KN/(KN+SNH3+SNH4))*Xalg #Growth on nitrate
    p2 = kresp*fT*(SO2/(KO2+SO2))*Xalg #Endogenous respiration
    p3 = kdeath*fT*Xalg #Inactivation
    p4 = keq1*(SCO2-((SH*SHCO3)/Keq1)) #CO2<->HCO3-
    p5 = keq2*(SHCO3-((SH*SCO3)/Keq2)) #HCO#-<->CO32-
    p6 = keq3*(SNH4-((SH*SNH3)/Keq3)) #NH4+<->NH3
    p7 = keqw*(1-((SH*SOH)/Keqw)) #H+<->OH-
    pO2 = KaO2*(SwatO2-SO2)
    pCO2 = KaCO2*(SwatCO2-SCO2)
    pnh3 = KaNH3*(SwatNH3-SNH3)
    
    #System of DAEs
    dSNH4 = p1a*v11a+p2*v12+p3*v13+p6*v16
    dSNH3 = p6*v26+pnh3*v2NH3
    dSNO3 = p1b*v31b
    dSO2 = p1a*v41a+p1b*v41b+p2*v42+p3*v43+pO2*v4O2
    dSCO2 = p1a*v51a+p1b*v51b+p2*v52+p3*v53+p4*v54+pCO2*v5CO2
    dSHCO3 = p4*v64+p5*v65
    dSCO3 = p7*v75
    dSH = p1a*v81a+p1b*v81b+p2*v82+p3*v83+p4*v84+p5*v85+p6*v86+p7*v87
    dSOH = p7*v97
    dXalg = p1a*v101a+p1b*v101b+p2*v102+p3*v103
    
    
    return(list(c(dSNH4, dSNH3, dSNO3, dSO2, dSCO2, dSHCO3, dSCO3, dSH, dSOH, dXalg)))
    
  })
}
dat <- data.frame (
  time = seq(from = 0, to = 9, by = 1),
  SNH4 = c(6.885, 5.988344, 6.4486636, 7.879068, 9.8168387, 11.99050, 14.27054, 1.598758, 18.94811, 21.30462),
  SNH3 = c(0.582, 7.341017, 9.9305313, 10.62804, 10.665681, 10.54405, 10.42312, 10.33779, 10.28618, 10.25955),
  SNO3 = c(9.665, 9.655609, 9.6509340, 9.647483, 9.6446015, 9.642094, 9.639859, 9.637853, 9.636025, 9.634351),
  SO2 = c(6.8000, 10.58095, 9.7564105, 9.365003, 9.1747385, 9.071812, 9.011224, 8.972813, 8.946718, 8.927814),
  SCO2 = c(0.680, 91.50242, 137.98161, 161.9230, 174.16874, 180.3239, 183.3213, 184.6996, 185.2657, 185.4400), 
  SHCO3 = c(85.0, 77.78401, 71.242809, 64.95223, 59.166452, 63.58839, 54.05013, 45.88668, 42.72116, 40.05481),
  SCO3 = c(0.9945, 1.07643, 1.1432874, 1.189526, 1.2168286, 1.229480, 1.231899, 1.227760, 1.219819, 1.210008),
  SH = c(0.000002686, 0.888924, 1.5916114, 2.1030037, 2.4885756, 2.798343, 3.062307, 3.297562, 3.513766, 3.71634),
  SOH = c(0.00241, 0.08433, 0.1511926, 0.197432, 0.2247337, 0.237386, 0.239805, 0.235666, .2277244, 0.217913),
  Xalg = c(85.000, 93.7326, 97.796835, 99.88826, 101.05973, 101.7453, 101.1505, 102.3812, 102.4952, 102.5259)
)
parms <- c(ualg=1.6, kresp=0.1, kdeath=0.1, KC=0.00432, ICO2=120, KN=0.1, KO2=0.2, Kpr=0.01, Tau=2, SsatO2=7.1904,  
           Topt=25, s=13, alpha=0.001935, beta=0.00000057848, gamma=0.1460, delta=0.0004796, Ef=1.74, N=0.74, rho=1353, 
           KaO2=4, KaCO2=0.7, KaNH3=0.7, Tact=13.5, Keq1=1, Keq2=1, Keq3=1, Keqw=1, keq1=.1, keq2=.1, keq3=.1, keqw=.1, 
           iC=0.387, iH=0.075, iO=0.538, iN=0.065, I=200, SwatO2=10.5, SwatCO2=220, SwatNH3=16)
state <- c(SNH4=8.1, SNH3=0.685, SNO3=11.37,SO2=8, SCO2=0.8, SHCO3=100, SCO3=1.17, SH=0.00000316, SOH=0.00283, Xalg=100)
times <- seq(from = 0, to = 9, by = 0.01)
out<- ode(y = state, times = times, func = Microbial, parms = parms, method = "rk4")
#return(out)
#tiff(file="Algae_model.tiff", width = 8, height = 10, units = "in", res = 300)
par(mar=c(1, 1, 1, 1))
#plot(out, lwd=2, xlab= "Time (days)", ylab= "Concentration (gm^-3)", col= "red", fg="gray", xlim=c(0,9))#cex.axis=2, cex.lab=1.45, cex.main=3.5, xlim=c(0,9))
plot(out, obs = dat)
#define cost function (least squares)
Cost <- function(P, dat) {
  P <- parms
  dat <- dat
  parms["ualg"]<-P[1]
  parms["kresp"]<-P[2]
  parms["kdeath"]<-P[3]
  parms["KC"]<-P[4]
  parms["ICO2"]<-P[5]
  parms["KN"]<-P[6]
  parms["KO2"]<-P[7]
  parms["Kpr"]<-P[8]
  parms["Tau"]<-P[9]
  parms["SsatO2"]<-P[10]
  parms["Topt"]<-P[11]
  parms["s"]<-P[12]
  parms["alpha"]<-P[13]
  parms["beta"]<-P[14]
  parms["gamma"]<-P[15]
  parms["delta"]<-P[16]
  parms["Ef"]<-P[17]
  parms["N"]<-P[18]
  parms["rho"]<-P[19]
  parms["KaO2"]<-P[20]
  parms["KaCO2"]<-P[21]
  parms["KaNH3"]<-P[22]
  parms["Tact"]<-P[23]
  parms["Keq1"]<-P[24]
  parms["Keq2"]<-P[25]
  parms["Keq3"]<-P[26]
  parms["Keqw"]<-P[27]
  parms["keq1"]<-P[28]
  parms["keq2"]<-P[29]
  parms["keq3"]<-P[30]
  parms["keqw"]<-P[31]
  parms["iC"]<-P[32]
  parms["iH"]<-P[33]
  parms["iO"]<-P[34]
  parms["iN"]<-P[35]
  parms["I"]<-P[36]
  parms["SwatO2"]<-P[37]
  parms["SwatCO2"]<-P[38]
  parms["SwatNH3"]<-P[39]
  out1 <- ode(y = state, times = times, func = Microbial, parms = P)
  modCost(model = out1, obs = dat, weight = "std")
}


#fit the model 
parms <- c(ualg=1.6, kresp=0.1, kdeath=0.1, KC=0.00432, ICO2=120, KN=0.1, KO2=0.2, Kpr=0.01, Tau=2, SsatO2=7.1904,  
           Topt=25, s=13, alpha=0.001935, beta=0.00000057848, gamma=0.1460, delta=0.0004796, Ef=1.74, N=0.74, rho=1353, 
           KaO2=4, KaCO2=0.7, KaNH3=0.7, Tact=13.5, Keq1=1, Keq2=1, Keq3=1, Keqw=1, keq1=.1, keq2=.1, keq3=.1, keqw=.1, 
           iC=0.387, iH=0.075, iO=0.538, iN=0.065, I=200, SwatO2=10.5, SwatCO2=220, SwatNH3=16)
fit <- modFit(f = Cost, p = parms)
summary(fit)

#compare output with data
out2 <- ode(state, times, Microbial, coef(fit))
plot(out, out2, obs = dat)

#fit parameters and initial values
Cost <- function(P, dat){
  yy <- P[c("SNH4", "SNH3", "SNO3","SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg")]
  pp <- P[c("ualg", "kresp", "kdeath", "KC", "ICO2", "KN", "KO2", "Kpr", "Tau", "SsatO2",  
            "Topt", "s", "alpha", "beta", "gamma", "delta", "Ef", "N", "rho", 
            "KaO2", "KaCO2", "KaNH3", "Tact", "Keq1", "Keq2", "Keq3", "Keqw", "keq1", "keq2", "keq3", "keqw", 
            "iC", "iH", "iO", "iN", "I", "SwatO2", "SwatCO2", "SwatNH3")]
  out1 <- ode(y=yy, times=times, func=Microbial, parms=pp)
  modCost(model=out1, obs=dat, weght="std")
}

state <- coef(fit)[c("SNH4", "SNH3", "SNO3","SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg")]
pp <- coef(fit)[c("ualg", "kresp", "kdeath", "KC", "ICO2", "KN", "KO2", "Kpr", "Tau", "SsatO2",  
                  "Topt", "s", "alpha", "beta", "gamma", "delta", "Ef", "N", "rho", 
                  "KaO2", "KaCO2", "KaNH3", "Tact", "Keq1", "Keq2", "Keq3", "Keqw", "keq1", "keq2", "keq3", "keqw", 
                  "iC", "iH", "iO", "iN", "I", "SwatO2", "SwatCO2", "SwatNH3")]
out3 <- ode(y=state, times=times, func=Microbial, parms=pp)
plot(out1, out2, out3, obs=dat)
summary(fit)

