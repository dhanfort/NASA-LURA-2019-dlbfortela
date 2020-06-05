##----------------------------------------------------------------
# This R-script was adopted from Fortela et al. (2019) Water Environment Research, Vol 91(9), P.865-876
# DOI: https://doi.org/10.1002/wer.1127
# and revised accordingly by A. DeLattre and D.L.B. Fortela (2019-2020)
# for the LURA 2019 Project sponsored by LaSPACE
# Sub-award No.: PO-0000105697
# Primary NASA Agreement No.: NNX15AH82H
# LaSPACE - Louisiana Space Consortium: https://laspace.lsu.edu/about-us/
#
# This LURA project was conducted at the University of Louisiana at Lafayette
# Department of Chemical Engineering: https://chemical.louisiana.edu/
#
# This script executes GSA-fPCA calculations on the Activated Sludge Model No.3 (ASM3)
# Title: Activated Sludge Models: ASM1, ASM2, ASM2d and ASM3 (2002)
# Authors: Henze, M., Gujer, W., Mino, T., & van Loosdrecht, M.
# Publisher: London, UK: IWA Publishing. 
# DOI: https://doi. org/10.2166/9781780402369.
#
# This must be executed in the R-statictical software, an open-source computing software
# R installer may be downloaded from the R-project website: https://www.r-project.org/
# We executed R via RStudio, which is an integrated development environment (IDE) for R
# RStudio may be downloaded for free from RStudio website: https://rstudio.com/
# We suggest the user to run R via RStudio for more organized workflow
##----------------------------------------------------------------


# Load the R-packages
library(deSolve); # To solve the ODE's
library(fda); # To perform PCA
library(sensitivity); # To perform GSA

# Set the working directory
directory = "..."; #paste the path to your working director in "..."; must use forward slash "/"
setwd(directory);

# Create a dummy csv file in the directory folder to check validity of path
write.csv(c(seq(0,3,by=1), format(Sys.time(), "%a %b %d %X %Y")), quote = FALSE, "testoutput.csv");
Run_code="Period10Band30_BSplines_Basis_r100";

r = 50; # Setting the number of EE's for the Morris GSA; good value is r=50; use r=2 for testing
RNG = 3; # Number of Random Numer Generator (RNG) seeds
k = 43; # Number of ASM3 model parameters under GSA study

# Initialize the dataframes for outputs
mu_RNG_ind_ALL <- data.frame(); # Results data frame for EE's mu
mu.star_RNG_ind_ALL <- data.frame(); # Results data frame for EE's mu*
sigma_RNG_ind_ALL <- data.frame(); # Results data frame for EE's sigma
PC1to3var_RNG_ind_ALL <- data.frame(); # Results data frame for PC1, PC2 & PC3 variances


##-----Loop through random number sets RNG
for(i in 1:RNG) {
  
  print(paste("Simulating the DAEs for i =", i, sep = " ")); # Print on screen to track progress
  
    ##-------------
    ## Sample parameters values using Morris method
    
    set.seed(i); # Set Randomized sampling seed for reproducibility of results
    
    mo <- morris(model = NULL, factors = k, r = r,
               design = list(type = "oat", levels = 20, grid.jump = 10),
               
               #  Set lower bounds and upper bounds of model parameter sampling range
               binf = c(0.00575, 0.0115, 0.0115, 0.0115, 0.0345, 0.000, 0.115, 0.7475, 0.575, 0.575, 0.23, 0.23, 0.23, 0.115, 0.023, 4.600, 4.60, 1.15, 0.920, 0.575, 0.115, 0.115, 0.0575, 0.115, 0.115, 0.115, 0.2875, 0.4025, 0.115, 0.23, 0.0575, 0.0575, 5.75, 0.00115, 0.115, 0.115, 0.0115, 0.0115, 0.23, 0.23, 0.0575, 0.0575, 0.115),
               bsup = c(0.0085 , 0.0425, 0.0595, 0.0425, 0.0935, 0.017, 0.255, 0.8075, 0.680, 0.680, 0.85, 0.85, 0.85, 0.255, 0.085, 12.75, 17.0, 5.10, 0.935, 0.850, 0.425, 0.425, 0.2125, 0.289, 0.850, 0.850, 0.8500, 0.8500, 0.850, 1.70, 0.4250, 4.2500, 17.0, 0.08500, 8.500, 8.500, 0.8500, 0.8500, 4.25, 4.25, 0.4250, 0.4250, 0.850),
               #        iNSI ,iNSS ,iNXI ,iNXS , iNBM, fSI , fXI, yH02,yHNO3,yHNO2,ySTOO2,YSTONO3,YSTONO2,YAOB,YNOB, kH,kSTO,uH ,uAOB,uNOB,bHO2,bSTOO2,bAOB ,bNOB ,nHNO3,nHNO2,nHendNO3,nHendNO2,nNend,KX  ,KHO2 ,KHO2inh,KHSS,KHNH4 ,KHNO3,KHNO2,KHALK,KHSTO,KAOBO2,KNOBO2,KAOBNH4,KNOBNO2,KNALK
               scale = TRUE 
               );
    
    N = 200 # no. of simulation samples to take
    Period = 10 # Simulated time in [days]
    
    ##--- Defining the ASM3 ----
    
  Microbial <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), { 
      
      
      # Model parameters
      
      iNSI = parameters[1]
      iNSS = parameters[2]
      iNXI = parameters[3]
      iNXS = parameters[4]
      iNBM = parameters[5]
      fSI  = parameters[6]
      fXI  = parameters[7]
      
      YHO2    = parameters[8]
      YHNO3   = parameters[9]
      YHNO2   = parameters[10]
      YSTOO2  = parameters[11]
      YSTONO3 = parameters[12]
      YSTONO2 = parameters[13]
      YAOB    = parameters[14]
      YNOB    = parameters[15]
      
      kH   = parameters[16]
      kSTO = parameters[17]
      muH   = parameters[18]
      muAOB = parameters[19]
      muNOB = parameters[20]
      
      bHO2     = parameters[21]
      bSTOO2   = parameters[22]
      bAOB     = parameters[23]
      bNOB     = parameters[24]
      nHNO3    = parameters[25]
      nHNO2    = parameters[26]
      nHendNO3 = parameters[27]
      nHendNO2 = parameters[28]
      nNend    = parameters[29]
      
      KX      = parameters[30]
      KHO2    = parameters[31]
      KHO2inh = parameters[32]
      KHSS    = parameters[33]
      KHNH4   = parameters[34]
      KHNO3   = parameters[35]
      KHNO2   = parameters[36]
      KHALK   = parameters[37]
      KHSTO   = parameters[38]
      KAOBO2  = parameters[39]
      KNOBO2  = parameters[40]
      KAOBNH4 = parameters[41]
      KNOBNO2 = parameters[42]
      KNALK   = parameters[43]
      
      
      # Stoichiometric coefficients
      
      u12 = 1-fSI
      u17 = (iNXS-fSI*iNSI-(1-fSI)*iNSS)/14
      u18 = fSI
      u111 = -1
      
      u21 = YSTOO2-1
      u22 = -1
      u23 = iNSS
      u27 = iNSS/14
      u212 = YSTOO2
      u32 = -1
      u33 = iNSS
      u34 = (1-YSTONO3)/1.14
      u35 = -(1-YSTONO3)/1.14
      u37 = iNSS/14
      u312 = YSTONO3
      
      u42 = -1
      u43 = iNSS
      u44 = -(1-YSTONO2)/1.72
      u46 = (1-YSTONO2)/1.72
      u47 = (iNSS-(YSTONO2-1)/1.72)/14
      u412 = YSTONO2
      
      u51 = (1-(1/YHO2))
      u53 = -iNBM

      u57 = -iNBM/14
      u510 = 1
      u512 = -1/YHO2
      
      u63 = -iNBM
      u64 = (1-YHNO3)/(1.14*YHNO3)
      u65 = -(1-YHNO3)/(1.14*YHNO3)
      u67 = -iNBM/14
      u610 = 1
      u612 = -1/YHNO3
      
      u73 = -iNBM
      u74 = -(1-YHNO2)/(1.72*YHNO2)
      u76 = (1-YHNO2)/(1.72*YHNO2)
      u77 = ((YHNO2-1)/(1.72*YHNO2)-iNBM)/14
      u710 = 1
      u712 = -1/YHNO2
      
      u81 = -(1-fXI)
      u83 = -fXI*iNXI + iNBM
      u87 = (-fXI*iNXI + iNBM)/14
      u89 = fXI
      u810 = -1
      
      u93 = -fXI*iNXI + iNBM
      u94 = (1-fXI)/1.14
      u95 = -(1-fXI)/1.14
      u97 = (-fXI*iNXI + iNBM)/14
      u99 = fXI
      u910 = -1
      
      u103 = -fXI*iNXI + iNBM
      u104 = -(1-fXI)/1.72
      u106 = (1-fXI)/1.72
      u107 = (iNBM - fXI*iNXI-(fXI-1)/1.72)/14
      u109 = fXI
      u1010 = -1
      
      u111 = -1
      u1112 = -1
      
      u124 = 1/1.14
      u125 = -1/1.14
      u1211 = -1
      u1212 = -1
      
      u134 = -1/1.72 
      u136 = 1/1.72
      u137 = 1/24.08
      u1312 = -1
      
      u141 = 1-(38/(14*YAOB))
      u143 = -(1/YAOB) - iNBM
      u144 = 1/YAOB
      u147 = (-1/14)*((2/YAOB)+iNBM)
      u1413 = 1
      
      u151 = -(1-fXI)
      u153 = -fXI*iNXI + iNBM
      u157 = (-fXI*iNXI + iNBM)/14
      u159 = fXI
      u1513 = -1
      
      u163 = -fXI*iNXI + iNBM
      u165 = -(1-fXI)/2.86
      u166 = (1-fXI)/2.86
      u167 = (iNBM - fXI*iNXI-(fXI-1)/2.86)/14
      u169 = fXI
      u1613 = -1
      
      u171 = 1-(16/(14*YAOB))
      u173 = -iNBM
      u174 = -1/YNOB
      u175 = 1/YNOB
      u177 = -iNBM/14
      u1714 = 1
      
      u181 = -(1-fXI)
      u183 = -fXI*iNXI + iNBM 
      u187 = (-fXI*iNXI + iNBM)/14
      u189 = fXI
      u1814 = -1
      
      u193 = -fXI*iNXI + iNBM 
      u195 = -(1-fXI)/2.86
      u196 = (1-fXI)/2.86
      u197 = (iNBM - fXI*iNXI - (fXI - 1)/2.86)/14
      u199 = fXI
      u1914 = -1
      
      # Aeration
      SO2sat = 8 # Saturation concentration of oxygen in gO2/m3
      kLA_O2=3*24 # kLA=3/hr for oxygen (3*24/day) in typical sewage aeration
      
      # Reaction rates
      
      j1  = kH*((XS/XH)/((KX+XS)/XH))*XH # Hydrolysis
      j2  = kSTO*((SO2)/(KHO2+SO2))*((SS)/(KHSS+SS))*XH # Aerobic Storage of Ss
      j3a = kSTO*nHNO3*((KHO2inh)/(KHO2inh+SO2))*((SS)/(KHSS+SS))*((SNO3)/(KHNO3+SNO3))*XH # Anoxic Storage of Ss NO3-NO2
      j3b = kSTO*nHNO2*((KHO2inh)/(KHO2inh+SO2))*((SS)/(KHSS+SS))*((SNO2)/(KHNO2+SNO2))*XH # Anoxic Storage pf Ss NO2-N2
      j4  = muH*((SO2)/(KHO2+SO2))*((SNH4)/(KHNH4+SNH4))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*XH # Aerobic Growth of Xh
      j5a = muH*nHNO3*((KHO2inh)/(KHO2inh+SO2))*((SNH4)/(KHNH4+SNH4))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*((SNO3)/(KHO2inh+SNO3))*XH # Anoxic Growth of Xh NO3-NO2
      j5b = muH*nHNO2*((KHO2inh)/(KHO2inh+SO2))*((SNH4)/(KHNH4+SNH4))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*((SNO2)/(KHO2inh+SNO2))*XH # Anoxic Growth of Xh NO2-N2
      j6  = bHO2*((SO2)/(KHO2+SO2))*XH # Aerobic Endog. Resp. of Xh
      j7a = bHO2*nHendNO3*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XH # Anoxic Endog. Resp. of Xh NO3-NO2
      j7b = bHO2*nHendNO2*((KHO2inh)/(KHO2inh+SO2))*((SNO2)/(KHNO2+SNO2))*XH # Anoxic Endog. Resp. of Xh NO2-N2
      j8  = bSTOO2*((SO2)/(KHO2+SO2))*XSTO # Aerobic Endog. Resp. of Xsto
      j9a = bSTOO2*nHendNO3*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XSTO # Anoxic Endog. Resp. of Xsto NO3-NO2
      j9b = bSTOO2*nHendNO2*((KHO2inh)/(KHO2inh+SO2))*((SNO2)/(KHNO2+SNO2))*XSTO # Anoxic Endog. Resp. of Xsto NO2-N2
      j10a = muAOB*((SO2)/(KAOBO2+SO2))*((SNH4)/(KAOBNH4+SNH4))*((SALK)/(KNALK+SALK))*XAOB # Aerobic Growth of Xaob Nitritation
      j10b = muNOB*((SO2)/(KNOBO2+SO2))*((SNH4)/(KHNH4+SNH4))*((SNO2)/(KNOBNO2+SNO2))*((SALK)/(KNALK+SALK))*XNOB # Aerobic Growth of Xaob Nitritation
      j11a = bAOB*((SO2)/(KHO2+SO2))*XAOB # Aerobic Endog. Resp. of Xaob
      j11b = bNOB*((SO2)/(KHO2+SO2))*XNOB # Aerobic Endog. Resp. of Xnob
      j12a = bAOB*nNend*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XAOB # Anoxic Endog. Resp. of Xaob
      j12b = bNOB*nNend*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XNOB # Anoxic Endog. Resp. of Xnob
      jO2 = kLA_O2*(SO2sat-SO2) # Aeration rate with sewage aeration parameters levels
      
      print(paste("r=", r, "i=", i, "Done with the reaction rates", sep = " "));
      
      
      # System of DAEs

      dSO2 = u21*j2+u51*j4+u81*j6+u111*j8+u141*j10a+u151*j11a+u171*j10b+u181*j11b+jO2
      dSS = u12*j1+u22*j2+u32*j3a+u42*j3b
      dSNH4 = u23*j2+u33*j3a+u43*j3b+u53*j4+u63*j5a+u73*j5b+u83*j6+u93*j7a+u103*j7b+u143*j10a+u153*j11a+u163*j12a+u173*j10b+
        u183*j11b+u193*j12b
      dSNO2 = u34*j3a+u44*j3b+u64*j5a+u74*j5b+u94*j7a+u104*j7b+u124*j9a+u134*j9b+u144*j10a+u174*j10b
      dSNO3 = u35*j3a+u65*j5a+u95*j7a+u125*j9a+u165*j12a+u175*j10b+u195*j12b
      dSN2 = u46*j3b+u76*j5b+u106*j7b+u136*j9b+u166*j12a+u196*j12b
      dSALK = u17*j1+u27*j2+u37*j3a+u47*j3b+u57*j4+u67*j5a+u77*j5b+u87*j6+u97*j7a+u107*j7b+u137*j9b+u147*j10a+u157*j11a+
        u167*j12a+u177*j10b+u187*j11b+u197*j12b
      dSI = u18*j1
      dXI = u89*j6+u99*j7a+u109*j7b+u159*j11a+u169*j12a+u189*j11b+u199*j12b
      dXH = u510*j4+u610*j5a+u710*j5b+u810*j6+u910*j7a+u1010*j7b
      dXS = u111*j1
      dXSTO = u212*j2+u312*j3a+u412*j3b+u512*j4+u612*j5a+u712*j5b+u1112*j8+u1212*j9a+u1312*j9b
      dXAOB = u1413*j10a+u1513*j11a+u1613*j12a
      dXNOB = u1714*j10b+u1814*j11b+u1914*j12b
      
      
      print(paste("r=", r, "i=", i, "Done with the ODEs", sep = " "));
      
      return(list(c(dSO2, dSS, dSNH4, dSNO2, dSNO3, dSN2, dSALK, dSI, dXI, dXH, dXS, dXSTO, dXAOB, dXNOB)))
      

    })
  }
    
    
    ##--- ASM3 Simulation ----
    ## Simulate ASM3 using ode() to generate data
    PCA_data <- function(doe) {
    
    # Specifying the constants
    #parameters = doe[]
    parameters = c(
      
      iNSI = doe[1],
      iNSS = doe[2],
      iNXI = doe[3],
      iNXS = doe[4],
      iNBM = doe[5],
      fSI  = doe[6],
      fXI  = doe[7],
      
      YHO2    = doe[8],
      YHNO3   = doe[9],
      YHNO2   = doe[10],
      YSTOO2  = doe[11],
      YSTONO3 = doe[12],
      YSTONO2 = doe[13],
      YAOB    = doe[14],
      YNOB    = doe[15],
      
      kH   = doe[16],
      kSTO = doe[17],
      muH   = doe[18],
      muAOB = doe[19],
      muNOB = doe[20],
      
      bHO2     = doe[21],
      bSTOO2   = doe[22],
      bAOB     = doe[23],
      bNOB     = doe[24],
      nHNO3    = doe[25],
      nHNO2    = doe[26],
      nHendNO3 = doe[27],
      nHendNO2 = doe[28],
      nNend    = doe[29],
      
      KX      = doe[30],
      KHO2    = doe[31],
      KHO2inh = doe[32],
      KHSS    = doe[33],
      KHNH4   = doe[34],
      KHNO3   = doe[35],
      KHNO2   = doe[36],
      KHALK   = doe[37],
      KHSTO   = doe[38],
      KAOBO2  = doe[39],
      KNOBO2  = doe[40],
      KAOBNH4 = doe[41],
      KNOBNO2 = doe[42],
      KNALK   = doe[43]
      
    )
    
    state     = c(SO2=8, SS=60, SNH4=16, SNO2=0, SNO3=0, SN2=0, SALK=5, SI=30, XI=25, XH=30, XS=115, XSTO=0.25, XAOB=0.25, XNOB=0.25)

    times      <- seq(0, Period, length = N)
    
    # Solving the ASM3
    out <- ode(y = state, times = times, func = Microbial, parms = parameters, method = "bdf")
    
    print(paste("r=", r, "i=", i, "Done evaluating ASM3 using ode()", sep = " "));
    
    return(out)
    
    
    }
    
    
    ##--- PCA -----
    
    # Generate the output for each row of the input design
    y <-apply(mo$X,1, PCA_data)
    
    print(paste("r=", r, "i=", i, "Done simulating the DAEs", sep = " "));
    
    
    for(j in 1:14){
    
      component <- c('SO2[gO2/m^3]', 'SS[gCOD/m^3]', 'SNH4[gN/m^3]', 'SNO2[gN/m^3]', 'SNO3[gN/m^3]', 'SN2[gN/m^3]', 'SALK[moleHCO3/m^3]', 'SI[gCOD/m^3]', 
                     'XI[gCOD/m^3]', 'XH[gCOD/m^3]', 'XS[gCOD/m^3]', 'XSTO[gCOD/m^3]', 'XAOB[gCOD/m^3]', 'XNOB[gCOD/m^3]');
      y_range_start <- (j*N) + 1;
      y_range_end <- (j+1)*N;
    
      y_j <- y[y_range_start:y_range_end,];
      
      
      print(paste("Performing fPCA on i =", i, sep = " "));
      print(paste("i = ", i, "; j = ", j, component[j], sep = " ")); # Print on screen to track progress
      
        
      setEPS();
      mypath <- file.path(directory,paste("ASM3_", "r", r, "i", i, "j", j, ".eps", sep=""));
      postscript(file=mypath, width = 5, height = 5);
      matplot_title <- paste(component[j], "Simulations", sep=" ");
      matplot(x=seq(0, Period, length = N), y=y_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], " ", sep = " "));
      dev.off();
      
      
      # Turn the vectors of model values into functional objects
      basisobjBSplines <- create.bspline.basis(rangeval = c(0,Period))
      
      yf <- Data2fd(seq(0,Period,length = N), y_j, basisobj = basisobjBSplines)
      
      
      # Then calculate the first 3 PCs using the fda
      ypca <- pca.fd(yf,3)
      
      # Plot the PCAs to see the variation they describe
      setEPS();
      mypath <- file.path(directory,paste("ASM3_PCsPlots_", "r", r, "i", i, "j", j, ".eps", sep=""));
      postscript(file=mypath, width = 6, height = 2.5);
      op <- par(mfrow=c(1,3), oma = c(0,0,0,0));
      plot.pca.fd(ypca);
      par(op);
      dev.off();
      
      
      # Here only the look at the first 3 PCs as they describe pretty much all the variance
      100*ypca$varprop
      
      # Now pass the coefficients (or scores) of the PCA back to the Morris method
      # just consider the first 3 as they account for all the variance
      y1 <- tell(mo,ypca$scores[,1:3])
      
      
      ##--------------------
      ## Retrieving the Morris SA parameters
      
      mu <- apply(y1$ee, 3, function(M){
        apply(M, 2, mean)
      })
      mu.star <- apply(abs(y1$ee), 3, function(M){
        apply(M, 2, mean)
      })
      sigma <- apply(y1$ee, 3, function(M){
        apply(M, 2, sd)
      })
      
      # Create codenames for the model parameters
      X = c(1:k); 
      Parameters = paste('P', 1:k, sep = "");
      
      # Organize the data row entries
      mu_RNG_ind <- cbind(mu,i,j,r,X,Parameters);
      mu.star_RNG_ind <- cbind(mu.star,i,j,r,X,Parameters);
      sigma_RNG_ind <- cbind(sigma,i,j,r,X,Parameters);
      PC1to3var_RNG_ind <- cbind(ypca$varprop[1],ypca$varprop[2],ypca$varprop[3],i,j,r);
      
      # Append the rows to create the  data frames for outputting
      mu_RNG_ind_ALL <- rbind(mu_RNG_ind_ALL, mu_RNG_ind);
      mu.star_RNG_ind_ALL <- rbind(mu.star_RNG_ind_ALL, mu.star_RNG_ind);
      sigma_RNG_ind_ALL <- rbind(sigma_RNG_ind_ALL, sigma_RNG_ind);
      PC1to3var_RNG_ind_ALL <- rbind(PC1to3var_RNG_ind_ALL, PC1to3var_RNG_ind);
      
      
    };



};


# File names for csv outputs
outfile_mu <- paste(Run_code,"-r", r, "-ASM3-GSA-Expt1-mu_RNG_ind_ALL.csv", sep = "");
outfile_mustar <- paste (Run_code,"-r", r, "-ASM3-GSA-Expt1-mu.star_RNG_ind_ALL.csv", sep = "");
outfile_sigma <- paste(Run_code, "-r", r, "-ASM3-GSA-Expt1-sigma_RNG_ind_ALL.csv", sep = "");
outfile_PC1to3var <- paste(Run_code, "-r", r, "-ASM3-GSA-Expt1-PC1to3var_RNG_ind_ALLL.csv", sep = "");

# Save output data files as csv format  
write.csv(mu_RNG_ind_ALL, quote = FALSE, outfile_mu);
write.csv(mu.star_RNG_ind_ALL,quote = FALSE, outfile_mustar);
write.csv(sigma_RNG_ind_ALL, quote = FALSE, outfile_sigma);
write.csv(PC1to3var_RNG_ind_ALL, quote = FALSE, outfile_PC1to3var)
