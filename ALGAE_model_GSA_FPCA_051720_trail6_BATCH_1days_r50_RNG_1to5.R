# ALGAE Model
# Batch mode, check line 47-48

library(deSolve)
library(fda)
library(sensitivity)

# Set the working directory where outputs will be saved
directory = "/Users/dhan-lordfortela/Documents/EIL/Proposals/LURA 2019/LURA_working/trial6_batch_1day_r50_RNG_1to5";
setwd(directory);

write.csv(seq(0,5,by=1), quote = FALSE, "testoutput.csv"); # Test directory is okay before the long computation stage


r = 50; # Setting the number of EE's for the Morris GSA
RNG = 5; # Number of Random Numer Generator (RNG) seeds
k = 20;

mu_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu
mu.star_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu*
sigma_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's sigma
PC1to3var_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances




#Period=10 #number of days
#N=100
##-----Loop through random number sets RNG
for(i in 1:RNG) {
  
  print(paste("Simulating the DAEs for i =", i, sep = " ")); # Print on screen to track progress
  
  ##--------Sample parameters using Morris method-----
  set.seed(i)
  
  doe3 <- morris(model = NULL, factors = k, r = r,
                 design = list(type = "oat", levels = 6, grid.jump = 3),
                 binf = c(1.36, 0.085, 0.085, 0.003672, 102, 0.085, 0.17, 0.0085, 1.7, 21.25, 11.05, 0.00164475, 0.000000491708, 0.1241, 0.00040766,  3.4, 0.595, 0.595, 20,     170),
                 bsup = c(1.84, 0.115, 0.115, 0.004968, 138, 0.115, 0.23, 0.0115, 2.3, 28.75, 14.95, 0.00222525, 0.000000665252, 0.1679, 0.00055154,  4.6, 0.805, 0.805, 40,     230),
                 #        ualg, kresp, kdeath, KC,     ICO2,  KN,   KO2,  Kpr,    Tau, Topt,   s,      alpha,     beta,          gamma,     delta,    KaO2, KaCO2, KaNH3, Tact,  I
                 #        1       2     3       4       5       6     7     8       9  10     11        12        13              14        15        16    17     18    19      20 
                 scale = TRUE)
  
  Period = 1 #Simulation time period in Days  
  N = 100 # no. of simulation samples to take for GSA-fPCA
  
  ##--------Setup simulation conditions------
  #Q = 0.075; # Continuous mode (CSTR), m^3/day
  Q = 0.0; # Batch mode, m^3/day
  Vliq=0.450; #m^3
  Vgas=0.060; #m^3
  tau=Q/Vliq; #/day
  
  
  #######---------Defining the Algae Model---------######
  Microbial <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), { 
      
      #model parameters
      
      doe3 <- parameters
      
      ualg = doe3[1] #1.6 #d^-1
      kresp  = doe3[2] #0.1 #d^-1
      kdeath = doe3[3] #0.1 #d^-1
      KC = doe3[4] #0.00432 #gC m^-3
      ICO2 = doe3[5] #120 #gC m^-3
      KN = doe3[6] #0.1 #gC m^-3
      KO2 = doe3[7] #0.2 #gO2 m^-3
      Kpr = doe3[8] #0.01 
      Tau = doe3[9] #2 #originally 4
      Topt = doe3[10] #25 #degrees celsius
      s = doe3[11] #13
      alpha = doe3[12] #0.001935 #(uE m^-2)^-1
      beta = doe3[13] #0.00000057848 #(uE m^-2)^-1
      gamma = doe3[14] #0.1460 #s^-1
      delta = doe3[15] #0.0004796 #s^-1
      #Ef = doe3[16] #1.74 #uE J^-1
      #N = doe3[17] #0.74
      KaO2 = doe3[16] #4 #d^-1
      KaCO2 = doe3[17] #0.7 #d^-1
      KaNH3 = doe3[18] #0.7 #d^-1
      Tact = doe3[19] #13.5 #PLACEHOLDER, good range accdg to the paper is 15C to 25C
      I = doe3[20] #200 #PLACEHOLDER, range is 3.25 to 665 uE/(m^2-s)
      
      SsatO2 = 7.1904 #gO2 m^-3
      rho = 1353 #W m^-2
      Keq1 = 10^(17.843-(3404.71/(273.15+Tact))-0.032786*(273.15+Tact))
      Keq2 = 10^(9.494-(2902.39/(273.15+Tact))-0.02379*(273.15+Tact))
      Keq3 = 10^(2.891-(2727/(273.15+Tact)))
      Keqw = 10^(-(4470.99/(273.15+Tact))+12.0875-0.01706*(273.15+Tact))
      keq1 = 10000 #d^-1
      keq2 = 1000 #d^-1
      keq3 = 1000 #d^-1
      keqw = 1000 #gm^-1 d^-1
      iC = 0.387 #gC gCOD^-1
      iH = 0.075 #gH gCOD^-1
      iO = 0.538 #gO2 gCOD^-1
      iN = 0.065 #gN gCOD^-1
      
      #SwatO2 = doe3[37] #10.5 #PLACEHOLDER
      SwatO2 = 7.19 #7.19 g/m^3 #PLACEHOLDER
      #SwatCO2 = doe3[38] #1440 g/m^3
      SwatCO2 = 220 #1440 g/m^3
      #SwatNH3 = doe3[39] #180000 g/m^3
      SwatNH3 = 16 #180000 g/m^3
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
      
      
      p1a = ualg*fT*nps*((SCO2+SHCO3)/(KC+SCO2+SHCO3+(((SCO2)^2)/ICO2)))*((SNH3+SNH4_alg)/(KN+SNH3+SNH4_alg))*Xalg #Growth on ammonia
      p1b = ualg*fT*nps*((SCO2+SHCO3)/(KC+SCO2+SHCO3+(((SCO2)^2)/ICO2)))*(SNO3/(KN+SNO3))*(KN/(KN+SNH3+SNH4_alg))*Xalg #Growth on nitrate
      p2 = kresp*fT*(SO2/(KO2+SO2))*Xalg #Endogenous respiration
      p3 = kdeath*fT*Xalg #Inactivation
      p4 = keq1*(SCO2-((SH*SHCO3)/Keq1)) #CO2<->HCO3-
      p5 = keq2*(SHCO3-((SH*SCO3)/Keq2)) #HCO#-<->CO32-
      p6 = keq3*(SNH4_alg-((SH*SNH3)/Keq3)) #NH4+<->NH3
      p7 = keqw*(1-((SH*SOH)/Keqw)) #H+<->OH-
      pO2 = KaO2*(SwatO2-SO2)
      pCO2 = KaCO2*(SwatCO2-SCO2)
      pnh3 = KaNH3*(SwatNH3-SNH3)
      
      # Influent comepoenets levels - also equal to initial conditions
      SNH4_alg_in=1.50 #g/m^3
      SNH3_in=0.685 #g/m^3
      SNO3_in=9 #g/m^3
      SO2_in=1 #g/m^3
      SCO2_in=0.8 #g/m^3
      SHCO3_in=100 #g/m^3
      SCO3_in=1.17#g/m^3
      SH_in=0.00000316#g/m^3
      SOH_in=0.00283 #g/m^3
      Xalg_in=80 #gCOD/m^3; COD/TSS=0.80
      
      
      
      #System of DAEs
      dSNH4_alg = tau*SNH4_alg_in-tau*SNH4_alg+p1a*v11a+p2*v12+p3*v13+p6*v16
      dSNH3 = tau*SNH3_in-tau*SNH3+p6*v26+pnh3*v2NH3
      dSNO3 = tau*SNO3_in-tau*SNO3+p1b*v31b
      dSO2 = tau*SO2_in-tau*SO2+p1a*v41a+p1b*v41b+p2*v42+p3*v43+pO2*v4O2
      dSCO2 = tau*SCO2_in-tau*SCO2+p1a*v51a+p1b*v51b+p2*v52+p3*v53+p4*v54+pCO2*v5CO2
      dSHCO3 = tau*SHCO3_in-tau*SHCO3+p4*v64+p5*v65
      dSCO3 = tau*SCO3_in-tau*SCO3+p7*v75
      dSH = tau*SH_in-tau*SH+p1a*v81a+p1b*v81b+p2*v82+p3*v83+p4*v84+p5*v85+p6*v86+p7*v87
      dSOH = tau*SOH_in-tau*SOH+p7*v97
      dXalg = tau*Xalg_in-tau*Xalg+p1a*v101a+p1b*v101b+p2*v102+p3*v103
      
      
      return(list(c(dSNH4_alg, dSNH3, dSNO3, dSO2, dSCO2, dSHCO3, dSCO3, dSH, dSOH, dXalg)))
      
    })
  }
  #SNH4_alg_sum = unname(V36)
  #SNO3_sum = unname(V38)
  #SO2_sum = unname(V34)
  
  ##--------Define the simulation function w/ ode()------
  PCA_data3 <- function(doe3) {
    
    
    # Specifying the constants
    #GSA_parameters = params
    ##states initial condition, liquid within the digester, not the input
    
    state <- c(SNH4_alg=1.50, SNH3=0.685, SNO3=9, SO2=1, SCO2=0.8, SHCO3=100, SCO3=1.17, SH=0.00000316, SOH=0.00283, Xalg=80)
    
    parameters=c(ualg = doe3[1], kresp  = doe3[2], kdeath = doe3[3], KC = doe3[4], ICO2 = doe3[5], KN = doe3[6], KO2 = doe3[7], Kpr = doe3[8], Tau = doe3[9], 
                 Topt = doe3[10], s = doe3[11], alpha = doe3[12], beta = doe3[13], gamma = doe3[14], delta = doe3[15], 
                 KaO2 = doe3[16], KaCO2 = doe3[17], KaNH3 = doe3[18], Tact = doe3[19], I = doe3[20])
    times <- seq(0, Period, length = N)
    #######-----Solving the algae model--------#####      
    out3<- ode(y = state, times = times, func = Microbial, parms = parameters, method = "bdf")
    
    return(out3)
    
  }
  
  #print(paste("r=", r, "i=", i, "Done evaluating algae using ode()", sep = " "));
  
  ##--------Pass GSA doe parameters & generate the algae output(s)-----
  yAlg <-apply(doe3$X,1, PCA_data3);
  
  ##------Run fPCA for each DAE simulation output-------
  for(j in 1:10){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
    
    component <- c("SNH4_alg", "SNH3", "SNO3", "SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg");
    y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
    y_range_end <- (j+1)*N;
    
    #
    
    y_j <- yAlg[y_range_start:y_range_end,];
    
    y_j_gsa <- yAlg[y_range_end,];
    
    # Plot all runs for j within i and save plots file to directory
    #plot(y_j[,1])
    setEPS(); # format as postscript for high-quality graphics
    mypath <- file.path(directory,paste("ALGAE_model_Responses_","r", r, "i", i, "j", j, ".eps", sep=""));
    postscript(file=mypath, width = 5, height = 5);
    matplot_title <- paste(component[j], "Simulations for", "r =", r, "&","RNG seed =", i, sep = " ");
    matplot(x=seq(0, Period, length = N), y=y_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], "(g/m^3)", sep = " "));
    dev.off();
    
    print(paste("Performing GSA calcs for i =", i,", j =", j, sep = " "));
    ##--------Run fPCA using fda package & perform Morris GSA on fPCA scores------
      # Turn the vectors of model values into functional objects
      yf <- Data2fd(y_j, seq(0, Period, length = N));
      #plot(yf)
      
      # Then calculate the first 3 PCs using the fda
      ypca <- pca.fd(yf,3);
      
      # Plot the first 2 PCAs to see the variation they describe
      #plot.pca.fd(ypca)
      
      
      
      setEPS(); # format as postscript for high-quality graphics
      mypath <- file.path(directory,paste("ALGAE_model_PCsPlots_","r", r, "i", i, "j", j, ".eps", sep=""));
      postscript(file=mypath, width = 10, height = 3);
      op <- par(mfrow=c(1,3), oma = c(0,0,0,0));
      plot.pca.fd(ypca);
      par(op);
      #title("My Title");
      dev.off();
      
      
      # Here only the look at the first 3 PCs as they describe pretty much all the variance
      #100*ypca$varprop;
      # Now pass the coefficients (or scores) of the PCA back to the Morris method
      # just consider the first 3 as they account for all the variance
      y1Alg <- tell(doe3,ypca$scores[,1:3]);
      
      
      ##--------Calculate the Morris GSA indices------
      
      mu <- apply(y1Alg$ee, 3, function(M){
        apply(M, 2, mean)
      });
      mu.star <- apply(abs(y1Alg$ee), 3, function(M){
        apply(M, 2, mean)
      });
      sigma <- apply(y1Alg$ee, 3, function(M){
        apply(M, 2, sd)
      });
      
      
      ##--------Plot Morris S.D. vs mu*-------
      #plot(x = mu.star[,1], y = sigma[,1], main = "PC1 Morris - r=40 RNGset.set.seed(1) - Try1.12/01/17", xlab = "mean, mu*", ylab = "S.D.", col = "yellow", cex = 1)
      #text(mu.star[,1], sigma[,1], labels = row.names(mu.star), cex = 0.7, offset = 15)
      
      
      X = c(1:k);
      Parameters = c("ualg", "kresp", "kdeath", "KC", "ICO2", "KN", "KO2", "Kpr", "Tau", "Topt", "s", "alpha", "beta", "gamma", "delta", "KaO2", "KaCO2", "KaNH3", "Tact",  "I");
      mu_RNG_ind <- cbind(mu,i,j,r,X,Parameters);
      mu.star_RNG_ind <- cbind(mu.star,i,j,r,X,Parameters);
      sigma_RNG_ind <- cbind(sigma,i,j,r,X,Parameters);
      PC1to3var_RNG_ind <- cbind(ypca$varprop[1],ypca$varprop[2],ypca$varprop[3],i,j,r);
      
      
      mu_RNG_ind_ALL <- rbind(mu_RNG_ind_ALL, mu_RNG_ind);
      mu.star_RNG_ind_ALL <- rbind(mu.star_RNG_ind_ALL, mu.star_RNG_ind);
      sigma_RNG_ind_ALL <- rbind(sigma_RNG_ind_ALL, sigma_RNG_ind);
      PC1to3var_RNG_ind_ALL <- rbind(PC1to3var_RNG_ind_ALL, PC1to3var_RNG_ind);
      
    
    };
    
  };
  

# File names for csv outputs
outfile_mu <- paste("r", r, "-ALGAE-GSA-FPCA_batch1days_05192020_trial6_mu_ALL.csv", sep = "");
outfile_mustar <- paste ("r", r, "-ALGAE-GSA-FPCA_batch1days_05192020_trial6_mu_star_ALL.csv", sep = "");
outfile_sigma <- paste("r", r, "-ALGAE-GSA-FPCA_batch1days_05192020_trial6_sigma_ALL.csv", sep = "");
outfile_PC1to3var <- paste("r", r, "-ALGAE-GSA-FPCA_batch1days_05192020_trial6-PC1to3var_ALL.csv", sep = "");

# Save output data files as csv format  
write.csv(mu_RNG_ind_ALL, quote = FALSE, outfile_mu);
write.csv(mu.star_RNG_ind_ALL,quote = FALSE, outfile_mustar);
write.csv(sigma_RNG_ind_ALL, quote = FALSE, outfile_sigma);
write.csv(PC1to3var_RNG_ind_ALL, quote = FALSE, outfile_PC1to3var)

