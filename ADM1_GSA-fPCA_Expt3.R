##----------------------------------------------------------------
# This R-script was adopted from Fortela et al. (2019) Biochemical Engineering Journal, Vol 143, P.212-223
# DOI: https://doi.org/10.1016/j.bej.2019.01.001
# and revised accordingly by A. DeLattre and D.L.B. Fortela (2019-2020)
# for the LURA 2019 Project sponsored by LaSPACE
#
# This script executes GSA-fPCA calculations on the Anaerobic Digestion Model No.1 (ADM1)
# based on the model Handbook published by the Internaitonal Water Association (IWA)
# Title: Anaerobic Digestion Model No.1 (ADM1), 2002
# Authored by D.J. Batstone; J. Keller;I. Angelidaki;S.V. Kalyuzhnyi;
# S.G. Pavlostathis;A. Rozzi;W.T.M. Sanders;H. Siegrist;V.A. Vavilin
# ISBN13: 9781900222785 
# The constants are based on the nominal values suggested in the ADM1 handbook
#
# This must be executed in the R-statictical software, an open-source computing software
# R installer may be downloaded from the R-project website: https://www.r-project.org/
# We executed R via RStudio, which is an integrated development environment (IDE) for R
# RStudio may be downloaded for free from RStudio website: https://rstudio.com/
# We suggest the user to run R via RStudio for more organized workflow
##----------------------------------------------------------------


##----------------------------------------------------------------
library(deSolve) #attach the 'deSolve' package; functions for solving the ODE's
library(fda) #attach the 'fda' package; functions for fPCA
library(sensitivity) #attach the 'sensitivity' package; functions for GSA (Morris GSA)


# Set the working directory where outputs will be saved
directory = "..."; #paste the path to your working directory in "..."
setwd(directory); #setting the working directory

write.csv(seq(0,5,by=1), quote = FALSE, "testoutput.csv"); # Test directory is okay before the long computation stage

# Computational Runs Index Codes:
# j is the index of the ODE output from the ode() solver
# i is the integer agrument to the set.seed(i) for RNG
# r is the number of elementary effects (EE's) for Morris GSa

r = 2; # Setting the number of EE's for the Morris GSA
RNG = 1; # Number of Random Numer Generator (RNG) seeds
k = 22; # Number of model parameters undergoing GSA analysis
  
mu_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu
mu.star_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu*
sigma_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's sigma
PC1to3var_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances

  ##-----Loop through random number sets RNG
  for(i in 1:RNG) {
    
    print(paste("Simulating the DAEs for i =", i, sep = " ")); # Print on screen to track progress
    
    ##--------Sample parameters using Morris method-----
    set.seed(i)
    
    doe <- morris(model = NULL, factors = k, r = r,
                  design = list(type = "oat", levels = 6, grid.jump = 3),
                  binf = c(0.91, 0.16, 0.09, 0.18, 0.18, 0.03, 0.28, 0.28, 0.13, 0.04, 0.18, 0.07, 0.14, 0.14, 0.18, 0.06, 0.07, 0.04, 0.04, 0.02, 0.03, 0.04),
                  bsup = c(0.98, 0.30, 0.17, 0.34, 0.35, 0.07, 0.53, 0.53, 0.25, 0.08, 0.32, 0.13, 0.26, 0.26, 0.32, 0.10, 0.13, 0.08, 0.08, 0.05, 0.06, 0.08),
                  scale = TRUE)
    # params: Ffa_li, Fva_aa, Fbu_su, Fbu_aa, Fpro_su, Fpro_aa, Fac_su, Fac_aa, Fh2_su, Fh2_aa, FxI_xc, FsI_xc, Fpr_xc, Fch_xc, Fli_xc, Yaa, Ysu, Yc4, Yfa, Ypro, Yac, Yh2      
    
    Period = 30 # Simulation time period in Days  
    N = 100 # no. of simulation samples to take for GSA-fPCA
    
    ##--------Setup simulation conditions------
    Q = 0.0090; # Continuous mode (CSTR)
    #Q = 0.0; # Batch mode
    Vliq=0.054;
    Vgas=0.006;
    tau=Q/Vliq;
    
    
    ##--------Define the ADM1 ode function-------
    ADM1_A <- function(t,state,parameters){ 
      with(as.list(c(state,parameters)), {
        
        ## Constants
        # GSA doe
        doe <- parameters
        
        Ffa_li=doe[1];
        Fva_aa=doe[2];
        Fbu_su=doe[3];
        Fbu_aa=doe[4];
        Fpro_su=doe[5];
        Fpro_aa=doe[6];
        Fac_su=doe[7];
        Fac_aa=doe[8];
        Fh2_su=doe[9];
        Fh2_aa=doe[10];
        FxI_xc=doe[11];
        FsI_xc=doe[12];
        Fpr_xc=doe[13];
        Fch_xc=doe[14];
        Fli_xc=doe[15];
        
        # Yield coefficients
        Yaa=doe[16];
        Ysu=doe[17];
        Yc4=doe[18];
        Yfa=doe[19];
        Ypro=doe[20];
        Yac=doe[21];
        Yh2=doe[22];
        
        # Yaa=0.08;
        # Ysu=0.1;
        # Yc4=0.06;
        # Yfa=0.06;
        # Ypro=0.04;
        # Yac=0.05;
        # Yh2=0.06;
        # 
        # Nitrogen coeff
        Nbac=0.08/14;
        Naa =0.007; 
        Nxc=0.0376/14;
        Ni=0.06/14;
        
        # Rate equation parameters
        Kdis=0.5;
        Khyd_ch=13;
        Khyd_pr=10;
        Khyd_li=10.5;
        Km_su=30;
        Ks_su=0.5;
        Km_aa=50;
        Ks_aa=0.3;
        Km_fa=6;
        Ks_fa=0.4;
        Km_c4=20;
        Ks_c4=0.2;
        Km_pro=13;
        Ks_pro=0.1;
        Km_ac=8;
        Ks_ac=0.15;
        Km_h2=35;
        Ks_h2=7e-6;
        Kdec_xsu= Kdec_xaa= Kdec_xfa=Kdec_xc4= Kdec_xpro= Kdec_xac= Kdec_xh2=0.02;
        Cxc=0.02786;
        CsI=0.03;
        Cch=0.0313;
        Cpr=0.03;
        Cli=0.022;
        CxI=0.03;
        Csu=0.0313;
        Caa=0.03;
        Cbu=0.025;
        Cpro=0.0268;
        Cac=0.0313;
        Cbac =0.0313;
        Cva=0.024;
        Cfa=0.0217;
        Cch4=0.0156;
        pHuL_aa=5.5;
        pHlL_aa=4;
        pHuL_ac=7;
        pHlL_ac=6;
        pHuL_h2=6;
        pHlL_h2=5;
        Ks_in=1e-4;
        Ki_h2_fa=5e-6;
        Ki_h2_c4=1e-5;
        Ki_h2_pro=3.5e-6;
        Ki_nh3=0.0018;
        Ka_bva= Ka_bbu= Ka_bpro= Ka_bac= Ka_bco2= Ka_bin=1e10;
        Ka_va=10^(-4.86); 
        Ka_bu=10^(-4.82);
        Ka_pro=10^(-4.88);
        Ka_ac=10^(-4.76);
        KL=200;
        R=0.083145;
        Tbase=298.15;
        Top=308.15;
        Pbar= 1.013;
        Patm = 1;
        Kw=(exp(55900/(R*100)*(1/Tbase-1/Top)))*(10^(-14))
        Ka_co2=10^(-6.35)*exp(7646/(R*100)*(1/Tbase-1/Top)) 
        Ka_in=10^(-9.25)*exp(51965/(R*100)*(1/Tbase-1/Top))
        Kh_h2=(7.8e-4)*exp(-4180/(R*100)*(1/Tbase-1/Top)) 
        Kh_ch4=0.0014*exp(-14240/(R*100)*(1/Tbase-1/Top)) 
        Kh_co2=0.035*exp(-19410/(R*100)*(1/Tbase-1/Top)) 
        Pgas_h2o=0.0313*exp(5290*(1/Tbase-1/Top))
        
        print(paste("r=", r, "i=", i, "Done with model parameters", sep = " "));
        
        
        ##input values
        Ssu_in=0.01; 
        Saa_in=0.001*0.1; # nominal(0.001) x 0.1
        Sfa_in=0.001;
        Sva_in=0.001;
        Sbu_in=0.001;
        Spro_in=0.001;
        Sac_in=0.001;
        Sh2_in=1e-8;
        Sch4_in=1e-5;
        Sic_in=0.01;
        Sin_in=0.02;
        Si_in=1.2;
        Xc_in=5.2;
        Xch_in=24.50;
        Xpr_in=7.9;
        Xli_in=13.90;
        Xsu_in=0.00;
        Xaa_in=0.01;
        Xfa_in=0.001;
        Xc4_in=0.01;
        Xpro_in=0.01;
        Xac_in=0.01;
        Xh2_in=0.01;
        Xi_in=16;
        Scation_in=0.04;
        Sanion_in=0.02;
        
        print(paste("r=", r, "i=", i, "Done with initial states", sep = " "));
        
        # Algebraic equ
        Snh4=Sin-Snh3
        Sco2=Sic-Shco3_m 
        Z=Scation+Snh4-Shco3_m-Sac_m/64-Spro_m/112-Sbu_m/160-Sva_m/208-Sanion 
        Kw=(exp(55900/(R*100)*(1/Tbase-1/Top)))*(10^(-14)) 
        Sh=-Z*.5+.5*sqrt(Z^2+4*Kw)
        
        print(paste("r=", r, "i=", i, "Done with algebraic equations", sep = " "));
        
        #pH
        pH = -log10(Sh)
        
        
        # Inhibition factors
        IpH_aa<- if ( pH<pHuL_aa) exp(-3*((pH-pHuL_aa)/(pHuL_aa-pHlL_aa))^2) else 1
        IpH_ac<- if ( pH<pHuL_ac) exp(-3*((pH-pHuL_ac)/(pHuL_ac-pHlL_ac))^2) else 1
        IpH_h2<- if ( pH<pHuL_h2) exp(-3*((pH-pHuL_h2)/(pHuL_h2-pHlL_h2))^2) else 1
        Iin_lim = 1/(1+Ks_in/Sin)
        Ih2_fa = 1/(1+Sh2/Ki_h2_fa)
        Ih2_c4 = 1/(1+Sh2/Ki_h2_c4)
        Ih2_pro = 1/(1+Sh2/Ki_h2_pro)
        Inh3 = 1/(1+Snh3/Ki_nh3)
        I5=I6=IpH_aa*Iin_lim
        I7=IpH_aa*Iin_lim*Ih2_fa
        I8=I9=IpH_aa*Iin_lim*Ih2_c4
        I10=IpH_aa*Iin_lim*Ih2_pro
        I11=IpH_ac*Iin_lim*Inh3
        I12=IpH_h2*Iin_lim
        
        print(paste("r=", r, "i=", i, "Done with inhibition effects", sep = " "));
        
        
        # Process rate equations
        P1=Kdis*Xc
        P2=Khyd_ch*Xch
        P3=Khyd_pr*Xpr
        P4=Khyd_li*Xli
        P5=Km_su*Ssu/(Ks_su+Ssu)*Xsu*I5
        P6=Km_aa*Saa/(Ks_aa+Saa)*Xaa*I6
        P7=Km_fa*Sfa/(Ks_fa+Sfa)*Xfa*I7
        P8=Km_c4*Sva/(Ks_c4+Sva)*Xc4*Sva/(Sbu+Sva+1e-6)*I8
        P9=Km_c4*Sbu/(Ks_c4+Sbu)*Xc4*Sbu/(Sva+Sbu+1e-6)*I9
        P10=Km_pro*Spro/(Ks_pro+Spro)*Xpro*I10
        P11=Km_ac*Sac/(Ks_ac+Sac)*Xac*I11
        P12=Km_h2*Sh2/(Ks_h2+Sh2)*Xh2*I12
        P13=Kdec_xsu*Xsu
        P14=Kdec_xaa*Xaa
        P15=Kdec_xfa*Xfa
        P16=Kdec_xc4*Xc4
        P17=Kdec_xpro*Xpro
        P18=Kdec_xac*Xac
        P19=Kdec_xh2*Xh2
        
        print(paste("r=", r, "i=", i, "Done with rxn rates", sep = " "));
        
        
        # Inorganic carbon rate coefficents
        S1=-Cxc+FsI_xc*CsI+Fch_xc*Cch+Fpr_xc*Cpr+Fli_xc*Cli+FxI_xc*CxI 
        S2=-Cch+Csu
        S3=-Cpr+Caa
        S4=-Cli+(1-Ffa_li)*Csu+Ffa_li*Cfa 
        S5=-Csu+(1-Ysu)*(Fbu_su*Cbu+Fpro_su*Cpro+Fac_su*Cac)+Ysu*Cbac 
        S6=-Caa+(1-Yaa)*(Fva_aa*Cva+Fbu_aa*Cbu+Fpro_aa*Cpro+Fac_aa*Cac)+Yaa*Cbac 
        S7=-Cfa+(1-Yfa)*0.7*Cac+Yfa*Cbac 
        S8=-Cva+(1-Yc4)*.54*Cpro+(1-Yc4)*.31*Cac+Yc4*Cbac 
        S9=-Cbu+(1-Yc4)*.8*Cac+Yc4*Cbac
        S10=-Cpro+(1-Ypro)*.57*Cac+ Ypro*Cbac
        S11=-Cac+(1-Yac)*Cch4+ Yac*Cbac
        S12=(1-Yh2)*Cch4+ Yh2*Cbac
        S13=-Cbac+Cxc
        
        print(paste("r=", r, "i=", i, "Done with inorganic carbon rate coefficient", sep = " "));
        
        #acid-base rates:
        Pa_4=Ka_bva*(Sva_m*(Ka_va+Sh)-Ka_va*Sva) 
        Pa_5=Ka_bbu*(Sbu_m*(Ka_bu+Sh)-Ka_bu*Sbu) 
        Pa_6=Ka_bpro*(Spro_m*(Ka_pro+Sh)-Ka_pro*Spro) 
        Pa_7=Ka_bac*(Sac_m*(Ka_ac+Sh)-Ka_ac*Sac) 
        Ka_co2=10^(-6.35)*exp(7646/(R*100)*(1/Tbase-1/Top)) 
        Pa_10=Ka_bco2*(Shco3_m*(Ka_co2+Sh)-Ka_co2*Sic) 
        Ka_in=10^(-9.25)*exp(51965/(R*100)*(1/Tbase-1/Top)) 
        Pa_11=Ka_bin*(Snh3*(Ka_in+Sh)-Ka_in*Sin)
        
        print(paste("r=", r, "i=", i, "Done with acid-base rates", sep = " "));
        
        
        #gas transfer equ&as transfer rates
        Pgas_h2=Sgas_h2*R*Top/16 
        Pgas_ch4=Sgas_ch4*R*Top/64 
        Pgas_co2=Sgas_co2*R*Top 
        Kh_h2=(7.8e-4)*exp(-4180/(R*100)*(1/Tbase-1/Top)) 
        Pt_8=KL*(Sh2-16*Kh_h2*Pgas_h2) 
        Kh_ch4=0.0014*exp(-14240/(R*100)*(1/Tbase-1/Top)) 
        Pt_9=KL*(Sch4-64*Kh_ch4*Pgas_ch4) 
        Kh_co2=0.035*exp(-19410/(R*100)*(1/Tbase-1/Top)) 
        Pt_10=KL*(Sco2-Kh_co2*Pgas_co2) 
        Pgas_h2o=0.0313*exp(5290*(1/Tbase-1/Top))
        Qgas = R*Top/(Patm-Pgas_h2o)*Vliq*(Pt_8/16+Pt_9/64+Pt_10)
        
        print(paste("r=", r, "i=", i, "Done with gas transfer rates", sep = " "));
        
        
        #Components dff equ.
        dSsu = tau*Ssu_in-tau*Ssu+(P2+(1-Ffa_li)*P4-P5) #C1 components 
        dSaa = tau*Saa_in-tau*Saa+(P3-P6)#C2
        dSfa = tau*Sfa_in-tau*Sfa+(Ffa_li*P4-P7)#C3
        dSva = tau*Sva_in-tau*Sva+((1-Yaa)*Fva_aa*P6-P8)#C4
        dSbu = tau*Sbu_in-tau*Sbu+((1-Ysu)*Fbu_su*P5+(1-Yaa)*Fbu_aa*P6-P9)#C5
        dSpro= tau*Spro_in-tau*Spro+((1-Ysu)*Fpro_su*P5+(1-Yaa)*Fpro_aa*P6+(1- Yc4)*.54*P8-P10)#C6
        dSac = tau*Sac_in-tau*Sac+((1-Ysu)*Fac_su*P5+(1-Yaa)*Fac_aa*P6+.50*(1-Yfa)*P7 +.20*(1-Yc4)*P8+.4*(1-Yc4)*P9+.50*(1-Ypro)*P10-P11)#C7
        dSh2 = tau*Sh2_in-tau*Sh2+((1-Ysu)*Fh2_su*P5+(1-Yaa)*Fh2_aa*P6+.50*(1-Yfa)*P7+.26*(1-Yc4)*P8+.6*(1-Yc4)*P9+.50*(1- Ypro)*P10-P12-Pt_8)#C8
        dSch4 = tau*Sch4_in-tau*Sch4+((1-Yac)*P11+(1-Yh2)*P12-Pt_9)#C9
        dSic = tau*Sic_in-tau*Sic- (sum(S1*P1,S2*P2,S3*P3,S4*P4,S5*P5,S6*P6,S7*P7,S8*P8,S9*P9,S10*P10,S11*P11, S12*P12)+S13*(P13+P14+P15+P16+P17+P18+P19))-Pt_10 #C10
        dSin = tau*Sin_in-tau*Sin-Ysu*Nbac*P5+(Naa-Yaa*Nbac)*P6-Yfa*Nbac* P7- Yc4*Nbac*P8-Yc4*Nbac*P9-Ypro*Nbac*P10-Yac*Nbac*P11-Yh2*Nbac*P12+(Nbac-Nxc)*sum(P13,P14,P15,P16,P17,P18,P19)+(Nxc-FxI_xc*Ni-FsI_xc*Ni- Fpr_xc*Naa)*P1 #C11
        dSi = tau*Si_in-tau*Si+FsI_xc*P1 #C12
        dXc = tau*Xc_in-tau*Xc +(-P1+sum(P13,P14,P15,P16,P17,P18,P19)) #C13 
        dXch = tau*Xch_in-tau*Xch +(Fch_xc*P1-P2) #C14
        dXpr = tau*Xpr_in-tau*Xpr +(Fpr_xc*P1-P3) #C15
        dXli = tau*Xli_in-tau*Xli +(Fli_xc*P1-P4) #C16
        dXsu = tau*Xsu_in-tau*Xsu +(Ysu*P5-P13) #C17
        dXaa = tau*Xaa_in-tau*Xaa +(Yaa*P6-P14) #C18
        dXfa =tau*Xfa_in-tau*Xfa +(Yfa*P7-P15) #C19
        dXc4 =tau*Xc4_in-tau*Xc4 +(Yc4*P8+Yc4*P9-P16) #C20
        dXpro =tau*Xpro_in-tau*Xpro +(Ypro*P10-P17) #C21
        dXac =tau*Xac_in-tau*Xac +(Yac*P11-P18) #C22
        dXh2 =tau*Xh2_in-tau*Xh2 +(Yh2*P12-P19) #C23
        dXi =tau*Xi_in-tau*Xi +(FxI_xc*P1) #C24
        dScation =tau*Scation_in-tau*Scation #C25 cations and anions
        dSanion =tau*Sanion_in-tau*Sanion #C26
        dSva_m = -Pa_4 #C27 ion states
        dSbu_m = -Pa_5 #C28
        dSpro_m = -Pa_6 #C29
        dSac_m = -Pa_7 #C30
        dShco3_m = -Pa_10 #C31
        dSnh3 = -Pa_11 #C32
        dSgas_h2 =-Sgas_h2*Qgas/Vgas+Pt_8*Vliq/Vgas #33
        dSgas_ch4 =-Sgas_ch4*Qgas/Vgas+Pt_9*Vliq/Vgas #34
        dSgas_co2 =-Sgas_co2*Qgas/Vgas+Pt_10*Vliq/Vgas #35
        
        print(paste("r=", r, "i=", i, "Done with the ODEs", sep = " "));
        
        
        return(list(c(dSsu,dSaa,dSfa,dSva,dSbu,dSpro,dSac,dSh2,dSch4,dSic,dSin,dSi,dXc,dXch,dXpr,dXli,dXsu,dXaa,dXfa,dXc4,dXpro,dXac,dXh2,dXi,dScation,dSanion,dSva_m,dSbu_m,dSpro_m,dSac_m,dShco3_m,dSnh3,dSgas_h2,dSgas_ch4,dSgas_co2)))
        
        
        
      })
 
      
    }
    
    # End of the ADM1 DAE system of equations function
    
    ##--------Define the simulation function w/ ode()------
    PCA_data <- function(doe) {
      
      
      # Specifying the constants
      #GSA_parameters = params
      ##states initial condition, liquid within the digester, not the input
      state=c(Ssu=0.011,Saa=0.005,Sfa=0.093, Sva=0.013, Sbu=0.013, Spro=0.0153, Sac=0.193, Sh2=2.3e-7,
              Sch4=0.055, Sic=0.04, Sin=0.01, Si=0.02, Xc=0.3, Xch=0.026, Xpr=0.1, Xli=0.03, Xsu=0.4, Xaa=1.17,
              Xfa=0.20, Xc4=0.41, Xpro=0.137, Xac=0.7, Xh2=0.317, Xi=5, Scation=0.04, Sanion=0.02,
              Sva_m=0.0601, Sbu_m=0.0905,Spro_m=0.13, Sac_m=0.159, Shco3_m=0.0090,
              Snh3=0.0165, Sgas_h2=0.03, Sgas_ch4=0.029, Sgas_co2=0.0378)
      #state      <- c(Cc = 1, Cs = 10, Cp = 0.01)
      
      ##parameters
      parameters=c(Ffa_li=doe[1], Fva_aa=doe[2], Fbu_su=doe[3], Fbu_aa=doe[4], Fpro_su=doe[5], Fpro_aa=doe[6],
                   Fac_su=doe[7], Fac_aa=doe[8], Fh2_su=doe[9], Fh2_aa=doe[10], FxI_xc=doe[11], FsI_xc=doe[12],
                   Fpr_xc=doe[13], Fch_xc=doe[14], Fli_xc=doe[15],
                    Yaa=doe[16], Ysu=doe[17], Yc4=doe[18], Yfa=doe[19],Ypro=doe[20], Yac=doe[21], Yh2=doe[22])
      #              Nbac=Nbac, Naa=Naa, Nxc=Nxc, Ni=Ni,
      #              Kdis= Kdis,Khyd_ch= Khyd_ch,
      #              Khyd_pr= Khyd_pr,Khyd_li= Khyd_li,Km_su= Km_su,Ks_su= Ks_su,Km_aa= Km_aa,Ks_aa= Ks_aa,Km_fa= Km_fa,
      #              Ks_fa= Ks_fa,Km_c4= Km_c4,Ks_c4= Ks_c4,Km_pro= Km_pro,Ks_pro= Ks_pro,Km_ac= Km_ac,Ks_ac= Ks_ac,
      #              Km_h2= Km_h2,Ks_h2= Ks_h2,Kdec_xsu=Kdec_xsu, Kdec_xaa= Kdec_xaa, Kdec_xfa= Kdec_xfa,
      #              Kdec_xc4= Kdec_xc4, Kdec_xpro= Kdec_xpro, Kdec_xac=Kdec_xac, Kdec_xh2= Kdec_xh2,Cxc= Cxc,
      #              CsI= CsI,Cch= Cch,Cpr= Cpr,Cli= Cli,CxI= CxI,Csu= Csu,Caa= Caa,Cbu= Cbu,
      #              Cpro= Cpro,Cac= Cac,Cbac = Cbac,Cva= Cva,Cfa= Cfa,Cch4= Cch4,pHuL_aa= pHuL_aa,pHlL_aa= pHlL_aa,
      #              pHuL_ac= pHuL_ac,pHlL_ac= pHlL_ac, pHuL_h2= pHuL_h2,pHlL_h2= pHlL_h2,Ks_in= Ks_in,Ki_h2_fa= Ki_h2_fa,
      #              Ki_h2_c4= Ki_h2_c4,Ki_h2_pro= Ki_h2_pro,Ki_nh3= Ki_nh3, Ka_bva= Ka_bva, Ka_bbu= Ka_bbu,
      #              Ka_bpro= Ka_bpro, Ka_bac= Ka_bac, Ka_bco2= Ka_bco2, Ka_bin= Ka_bin, Ka_va= Ka_va,
      #              Ka_bu= Ka_bu,Ka_pro= Ka_pro,Ka_ac= Ka_ac,KL= KL,R= R, Tbase=Tbase,Top= Top,Patm= Patm,
      #              Kh_h2= Kh_h2, Kh_ch4= Kh_ch4,Kh_co2= Kh_co2,Ka_in= Ka_in,Pgas_h2o=Pgas_h2o)
      
      times <- seq(0, Period, length = N);
      
      # Solving the DAEs
      out <- ode(y = state, times = times, func = ADM1_A, parms = parameters, method = "bdf", maxsteps = 20000)
      
      print(paste("r=", r, "i=", i, "Done evaluating the ADM1 sys of equations using ode()", sep = " "));
      
      return(out) # Return C1 component values
      
    }
    
    
    ##--------Pass GSA doe parameters & generate the ADM1 output(s)-----
    y <-apply(doe$X,1, PCA_data);
    
    print(paste("r=", r, "i=", i, "Done simulating the DAEs", sep = " "));
    
    ##------Run fPCA for each DAE simulation output-------
    for(j in 1:35){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
      
      component <- c('Ssu','Saa','Sfa','Sva','Sbu','Spro','Sac','Sh2','Sch4','Sic','Sin','Si','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','Xi','Scation','Sanion','Sva_m','Sbu_m','Spro_m','Sac_m','Shco3_m','Snh3','Sgas_h2','Sgas_ch4','Sgas_co2');
      y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
      y_range_end <- (j+1)*N;
      
      y_j <- y[y_range_start:y_range_end,];
      
      print(paste("Performing fPCA on i =", i, sep = " "));
      print(paste("i = ", i, "; j = ", j, component[j], sep = " ")); # Print on screen to track progress
      
      # Plot all runs for j within i and save plots file to directory
      #plot(y_j[,1])
      setEPS(); # format as postscript for high-quality graphics
      mypath <- file.path(directory,paste("AMD1_MSThesis_Responses_","r", r, "i", i, "j", j, ".eps", sep=""));
      postscript(file=mypath, width = 5, height = 5);
      matplot_title <- paste(component[j], "Simulations for", "r =", r, "&","RNG seed =", i, sep = " ");
      matplot(x=seq(0, Period, length = N), y=y_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], "(kgCOD/m^3)", sep = " "));
      dev.off();
      
      
      ##--------Run fPCA using fda package & perform Morris GSA on fPCA scores------
      # Turn the vectors of model values into functional objects
      yf <- Data2fd(y_j, seq(0, Period, length = N));
      #plot(yf)
      
      # Then calculate the first 3 PCs using the fda
      ypca <- pca.fd(yf,3);
      
      # Plot the first 2 PCAs to see the variation they describe
      #plot.pca.fd(ypca)
      
      
      
      setEPS(); # format as postscript for high-quality graphics
      mypath <- file.path(directory,paste("ADM1_MSThesis_PCsPlots_","r", r, "i", i, "j", j, ".eps", sep=""));
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
      y1 <- tell(doe,ypca$scores[,1:3]);
      
      
      ##--------Calculate the Morris GSA indices------
      
      mu <- apply(y1$ee, 3, function(M){
        apply(M, 2, mean)
      });
      mu.star <- apply(abs(y1$ee), 3, function(M){
        apply(M, 2, mean)
      });
      sigma <- apply(y1$ee, 3, function(M){
        apply(M, 2, sd)
      });
      
      
      ##--------Plot Morris S.D. vs mu*-------
      #plot(x = mu.star[,1], y = sigma[,1], main = "PC1 Morris - r=40 RNGset.set.seed(1) - Try1.12/01/17", xlab = "mean, mu*", ylab = "S.D.", col = "yellow", cex = 1)
      #text(mu.star[,1], sigma[,1], labels = row.names(mu.star), cex = 0.7, offset = 15)
      
  
      X = c(1:k);
      Parameters = c("fa_li", "va_aa", "bu_su", "bu_aa", "pro_su", "pro_aa", "ac_su", "ac_aa", "h2_su", "h2_aa", "xI_xc", "sI_xc", "pr_xc", "ch_xc", "li_xc", "aa", "su", "c4", "fa", "pro", "ac", "h2");
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
outfile_mu <- paste("r", r, "-ADM1-GSA-Expt3-mu_RNG_ind_ALL.csv", sep = "");
outfile_mustar <- paste ("r", r, "-ADM1-GSA-Expt3-mu.star_RNG_ind_ALL.csv", sep = "");
outfile_sigma <- paste("r", r, "-ADM1-GSA-Expt3-sigma_RNG_ind_ALL.csv", sep = "");
outfile_PC1to3var <- paste("r", r, "-ADM1-GSA-Expt3-PC1to3var_RNG_ind_ALLL.csv", sep = "");

# Save output data files as csv format  
write.csv(mu_RNG_ind_ALL, quote = FALSE, outfile_mu);
write.csv(mu.star_RNG_ind_ALL,quote = FALSE, outfile_mustar);
write.csv(sigma_RNG_ind_ALL, quote = FALSE, outfile_sigma);
write.csv(PC1to3var_RNG_ind_ALL, quote = FALSE, outfile_PC1to3var)
  

