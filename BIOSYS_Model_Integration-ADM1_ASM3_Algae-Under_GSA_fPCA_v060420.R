#GSA on INTEGRATED MODEL FOR ADM1, ASM3, AND ALGAE MODELS



#This script executes ADM1 using deSolve
library(deSolve)
library(fda)
library(sensitivity)

# Set the working directory where outputs will be saved
directory = "/Users/dhan-lordfortela/Documents/EIL/Proposals/LURA 2019/LURA_working/BIOSYS_Bioprocesses_INTEGRATED/GSA Model Integration-ADM1_ASM3_Algae/BIOSYS_Model_Integration-ADM1_ASM3_Algae-Under_GSA_fPCA_v060420"
setwd(directory);

write.csv(seq(0,5,by=1), quote = FALSE, "testoutput.csv"); # Test directory is okay before the long computation stage

# Computational Runs Index Codes:
# j is the index of the ODE output from the ode() solver
# i is the integer agrument to the set.seed(i) for RNG
# r is the number of elementary effects (EE's) for Morris GSa

r = 5; # Setting the number of EE's for the Morris GSA
RNG = 1; # Number of Random Numer Generator (RNG) seeds
k = 62; # Number of model parameters under GSA-fPCA calcs

mu_RNG_ind_ALL_ADM1 <- data.frame(); # Initialize the results data frame for EE's mu
mu.star_RNG_ind_ALL_ADM1 <- data.frame(); # Initialize the results data frame for EE's mu*
sigma_RNG_ind_ALL_ADM1 <- data.frame(); # Initialize the results data frame for EE's sigma
PC1to3var_RNG_ind_ALL_ADM1 <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances

mu_RNG_ind_ALL_ASM3 <- data.frame(); # Initialize the results data frame for EE's mu
mu.star_RNG_ind_ALL_ASM3 <- data.frame(); # Initialize the results data frame for EE's mu*
sigma_RNG_ind_ALL_ASM3 <- data.frame(); # Initialize the results data frame for EE's sigma
PC1to3var_RNG_ind_ALL_ASM3 <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances

mu_RNG_ind_ALL_Algae <- data.frame(); # Initialize the results data frame for EE's mu
mu.star_RNG_ind_ALL_Algae <- data.frame(); # Initialize the results data frame for EE's mu*
sigma_RNG_ind_ALL_Algae <- data.frame(); # Initialize the results data frame for EE's sigma
PC1to3var_RNG_ind_ALL_Algae <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances


i=1

doe <- morris(model = NULL, factors = k, r = r,
              design = list(type = "oat", levels = 6, grid.jump = 3),
              binf = c(0.91, 0.16, 0.09, 0.18, 0.18, 0.03, 0.28, 0.28, 0.13, 0.04, 0.18, 0.07, 0.14, 0.14, 0.18, 0.06, 0.07, 0.04, 0.04, 0.02, 0.03, 0.04,
              #ADM1 params: Ffa_li, Fva_aa, Fbu_su, Fbu_aa, Fpro_su, Fpro_aa, Fac_su, Fac_aa, Fh2_su, Fh2_aa, FxI_xc, FsI_xc, Fpr_xc, Fch_xc, Fli_xc, Yaa, Ysu, Yc4, Yfa, Ypro, Yac, Yh2
                       0.00575, 0.0115, 0.0115, 0.0115, 0.0345, 0.000, 0.115, 0.7475, 0.575, 0.575, 0.23, 0.23, 0.23, 0.115, 0.023, 4.600, 4.60, 1.15, 0.920, 0.575,
              #ASM3 params: iNSI ,iNSS ,iNXI ,iNXS , iNBM, fSI , fXI, yH02,yHNO3,yHNO2,ySTOO2,YSTONO3,YSTONO2,YAOB,YNOB, kH,kSTO,uH ,uAOB,uNOB
                       1.36, 0.085, 0.085, 0.003672, 102, 0.085, 0.17, 0.0085, 1.7, 21.25, 11.05, 0.00164475, 0.000000491708, 0.1241, 0.00040766,  3.4, 0.595, 0.595, 20, 170
              #Algae params: ualg, kresp, kdeath, KC, ICO2, KN, KO2,  Kpr, Tau, Topt, s, alpha, beta, gamma, delta, KaO2, KaCO2, KaNH3, Tact, I              
                       ),
              
              bsup = c(0.98, 0.30, 0.17, 0.34, 0.35, 0.07, 0.53, 0.53, 0.25, 0.08, 0.32, 0.13, 0.26, 0.26, 0.32, 0.10, 0.13, 0.08, 0.08, 0.05, 0.06, 0.08,
              #ADM1 params: Ffa_li, Fva_aa, Fbu_su, Fbu_aa, Fpro_su, Fpro_aa, Fac_su, Fac_aa, Fh2_su, Fh2_aa, FxI_xc, FsI_xc, Fpr_xc, Fch_xc, Fli_xc, Yaa, Ysu, Yc4, Yfa, Ypro, Yac, Yh2
                       0.0085 , 0.0425, 0.0595, 0.0425, 0.0935, 0.017, 0.255, 0.8075, 0.680, 0.680, 0.85, 0.85, 0.85, 0.255, 0.085, 12.75, 17.0, 5.10, 0.935, 0.850,
              #ASM3 params: iNSI ,iNSS ,iNXI ,iNXS , iNBM, fSI , fXI, yH02,yHNO3,yHNO2,ySTOO2,YSTONO3,YSTONO2,YAOB,YNOB, kH,kSTO,uH ,uAOB,uNOB
                       1.84, 0.115, 0.115, 0.004968, 138, 0.115, 0.23, 0.0115, 2.3, 28.75, 14.95, 0.00222525, 0.000000665252, 0.1679, 0.00055154,  4.6, 0.805, 0.805, 40, 230
              #Algae params: ualg, kresp, kdeath, KC, ICO2, KN, KO2,  Kpr, Tau, Topt, s, alpha, beta, gamma, delta, KaO2, KaCO2, KaNH3, Tact, I              
                       ),
              scale = TRUE)


Period = 1 # Simulation time period in Days  
N = 100 # no. of simulation samples to take for GSA-fPCA

##--------Setup simulation conditions------
#Q = 0.0090; # Continuous mode (CSTR)
Q = 0.0; # Batch mode
Vliq=0.054;
Vgas=0.006;
tau=Q/Vliq;


##--------Define the ADM1 ode function-------  
ADM1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { 
    
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
    
    print(paste("ADM1: r=", r, "Done with model parameters", sep = " "));
    
    
    ##input values
    Ssu_in=0.01; 
    Saa_in=0.001; # nominal(0.001)
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
    
    #print(paste("r=", r, "i=", i, "Done with initial states", sep = " "));
    
    # Algebraic equ
    Snh4=Sin-Snh3
    Sco2=Sic-Shco3_m 
    Z=Scation+Snh4-Shco3_m-Sac_m/64-Spro_m/112-Sbu_m/160-Sva_m/208-Sanion 
    Kw=(exp(55900/(R*100)*(1/Tbase-1/Top)))*(10^(-14)) 
    Sh=-Z*.5+.5*sqrt(Z^2+4*Kw)
    
    #print(paste("r=", r, "i=", i, "Done with algebraic equations", sep = " "));
    
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
    
    print(paste("ADM1: r=", r, "Done with inhibition effects", sep = " "));
    
    
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
    
    #print(paste("r=", r, "i=", i, "Done with rxn rates", sep = " "));
    
    
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
    
    #print(paste("r=", r, "i=", i, "Done with inorganic carbon rate coefficient", sep = " "));
    
    #acid-base rates:
    Pa_4=Ka_bva*(Sva_m*(Ka_va+Sh)-Ka_va*Sva) 
    Pa_5=Ka_bbu*(Sbu_m*(Ka_bu+Sh)-Ka_bu*Sbu) 
    Pa_6=Ka_bpro*(Spro_m*(Ka_pro+Sh)-Ka_pro*Spro) 
    Pa_7=Ka_bac*(Sac_m*(Ka_ac+Sh)-Ka_ac*Sac) 
    Ka_co2=10^(-6.35)*exp(7646/(R*100)*(1/Tbase-1/Top)) 
    Pa_10=Ka_bco2*(Shco3_m*(Ka_co2+Sh)-Ka_co2*Sic) 
    Ka_in=10^(-9.25)*exp(51965/(R*100)*(1/Tbase-1/Top)) 
    Pa_11=Ka_bin*(Snh3*(Ka_in+Sh)-Ka_in*Sin)
    
    #print(paste("r=", r, "i=", i, "Done with acid-base rates", sep = " "));
    
    
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
    
    #print(paste("r=", r, "i=", i, "Done with gas transfer rates", sep = " "));
    
    
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
    
    #print(paste("r=", r, "i=", i, "Done with the ODEs", sep = " "));
    
    
    return(list(c(dSsu,dSaa,dSfa,dSva,dSbu,dSpro,dSac,dSh2,dSch4,dSic,dSin,dSi,dXc,dXch,dXpr,dXli,dXsu,dXaa,dXfa,dXc4,dXpro,dXac,dXh2,dXi,dScation,dSanion,dSva_m,dSbu_m,dSpro_m,dSac_m,dShco3_m,dSnh3,dSgas_h2,dSgas_ch4,dSgas_co2)))
    
#    return(list(c(dSsu, dSaa, dSfa, dSva, dSbu, dSpro, dSac, dSh2, dSch4, dSic, dSin, dSi, dXc, dXch, dXpr, dXli, dXsu, dXaa, dXfa, dXc4, dXpro, dXac, dXh2, dXi, dSva_m, dSbu_m, dSpro_m, dSac_m, dShco3_m, dSnh3_m, dSgas_h2, dSgas_ch4, dSgas_co2)))
    
  })
}
# End of the ADM1 DAE system of equations function

##--------Define the simulation function w/ ode()------
PCA_data1 <- function(doe) {
  
  
  # Specifying the constants
  #GSA_parameters = params
  ##states initial condition, liquid within the digester, not the input
  state1=c(Ssu=0.011,Saa=0.005,Sfa=0.093, Sva=0.013, Sbu=0.013, Spro=0.0153, Sac=0.193, Sh2=2.3e-7,
          Sch4=0.055, Sic=0.04, Sin=0.01, Si=0.02, Xc=0.3, Xch=0.026, Xpr=0.1, Xli=0.03, Xsu=0.4, Xaa=1.17,
          Xfa=0.20, Xc4=0.41, Xpro=0.137, Xac=0.7, Xh2=0.317, Xi=5, Scation=0.04, Sanion=0.02,
          Sva_m=0.0601, Sbu_m=0.0905,Spro_m=0.13, Sac_m=0.159, Shco3_m=0.0090,
          Snh3=0.0165, Sgas_h2=0.03, Sgas_ch4=0.029, Sgas_co2=0.0378)
  ##parameters
  parameters1=c(ffa_li=doe[1], fva_aa=doe[2], fbu_su=doe[3], fbu_aa=doe[4], fpro_su=doe[5], fpro_aa=doe[6],
               fac_su=doe[7], fac_aa=doe[8], fh2_su=doe[9], fh2_aa=doe[10], fxi_xc=doe[11], FsI_xc=doe[12], 
               fpr_xc=doe[13], fch_xc=doe[14], fli_xc=doe[15], Yaa=doe[16], Ysu=doe[17], Yc4=doe[18], Yfa=doe[19],
               Ypro=doe[20], Yac=doe[21], Yh2=doe[22])
  
  times1 <- seq(0, Period, length=N)
  
  # Solving the DAEs
  out1<- ode(y = state1, times = times1, func = ADM1, parms = parameters1, method = "bdf", maxsteps=5000)
  
  return(out1)
  
}
##--------Pass GSA doe parameters & generate the ADM1 output(s)-----
yADM1 <-apply(doe$X[,1:22],1, PCA_data1);

 ##------Run fPCA for each DAE simulation output-------
 for(j in 1:35){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
   
   component <- c('Ssu','Saa','Sfa','Sva','Sbu','Spro','Sac','Sh2','Sch4','Sic','Sin','Si','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','Xi','Scation','Sanion','Sva_m','Sbu_m','Spro_m','Sac_m','Shco3_m','Snh3','Sgas_h2','Sgas_ch4','Sgas_co2');
   y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
   y_range_end <- (j+1)*N;
   
   #yADM1 <- as.data.frame(yADM1)
   
   y1_j <- yADM1[y_range_start:y_range_end,];
   
   # Plot all runs for j within i and save plots file to directory
   #plot(y_j[,1])
   setEPS(); # format as postscript for high-quality graphics
   mypath <- file.path(directory,paste("AMD1_Responses_","r", r, "i", i, "j", j, ".eps", sep=""));
   postscript(file=mypath, width = 5, height = 5);
   matplot_title <- paste(component[j], "Simulations for", "r =", r, "&","RNG seed =", i, sep = " ");
   matplot(x=seq(0, Period, length = N), y=y1_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], "(kgCOD/m^3)", sep = " "));
   dev.off();
   
   
   ##--------Run fPCA using fda package & perform Morris GSA on fPCA scores------
   # Turn the vectors of model values into functional objects
   yf_1j <- Data2fd(y1_j, seq(0, Period, length = N));
   #plot(yf)
   
   # Then calculate the first 3 PCs using the fda
   ypca_1j <- pca.fd(yf_1j,3);
   

   # Plot the first 3 PCAs to see the variation they describe
   #plot.pca.fd(ypca)
   
   setEPS(); # format as postscript for high-quality graphics
   mypath <- file.path(directory,paste("ADM1_PCsPlots_","r", r, "i", i, "j", j, ".eps", sep=""));
   postscript(file=mypath, width = 10, height = 3);
   op <- par(mfrow=c(1,3), oma = c(0,0,0,0));
   plot.pca.fd(ypca_1j);
   par(op);
   #title("My Title");
   dev.off();
 

    # Here only the look at the first 3 PCs as they describe pretty much all the variance
    #100*ypca$varprop;
    # Now pass the coefficients (or scores) of the PCA back to the Morris method
    # just consider the first 3 as they account for all the variance

   
    y1ADM1_scores <- tell(doe,ypca_1j$scores[,1:3]);
   
   ##--------Calculate the Morris GSA indices------
   
    mu_y1ADM1 <- apply(y1ADM1_scores$ee, 3, function(M){
      apply(M, 2, mean)
    });
    mu.star_y1ADM1 <- apply(abs(y1ADM1_scores$ee), 3, function(M){
      apply(M, 2, mean)
    });
    sigma_y1ADM1 <- apply(y1ADM1_scores$ee, 3, function(M){
      apply(M, 2, sd)
    });
   
   
##-----Summarize the GSA indices per parameter and PC
    X = c(1:k);
    Parameters = paste('Parm', 1:k, sep = "");
    
    mu_RNG_ind_ADM1 <- cbind(mu_y1ADM1,i,j,r,X,Parameters);
    mu.star_RNG_ind_ADM1 <- cbind(mu.star_y1ADM1,i,j,r,X,Parameters);
    sigma_RNG_ind_ADM1 <- cbind(sigma_y1ADM1,i,j,r,X,Parameters);
    PC1to3var_RNG_ind_ADM1 <- cbind(ypca_1j$varprop[1],ypca_1j$varprop[2],ypca_1j$varprop[3],i,j,r);
   
   
    mu_RNG_ind_ALL_ADM1 <- rbind(mu_RNG_ind_ALL_ADM1, mu_RNG_ind_ADM1);
    mu.star_RNG_ind_ALL_ADM1 <- rbind(mu.star_RNG_ind_ALL_ADM1, mu.star_RNG_ind_ADM1);
    sigma_RNG_ind_ALL_ADM1 <- rbind(sigma_RNG_ind_ALL_ADM1, sigma_RNG_ind_ADM1);
    PC1to3var_RNG_ind_ALL_ADM1 <- rbind(PC1to3var_RNG_ind_ALL_ADM1, PC1to3var_RNG_ind_ADM1);
   
 };


#Save the GSA-fPCA results for the ADM1 section

# File names for csv outputs
outfile_mu_ADM1 <- paste("r", r, "-BIOSYS-GSAfPCA-ADM1section-mu_RNG_ind_ALL.csv", sep = "");
outfile_mustar_ADM1 <- paste ("r", r, "-BIOSYS-GSAfPCA-ADM1section-mu_star_RNG_ind_ALL.csv", sep = "");
outfile_sigma_ADM1 <- paste("r", r, "-BIOSYS-GSAfPCA-ADM1section-sigma_RNG_ind_ALL.csv", sep = "");
outfile_PC1to3var_ADM1 <- paste("r", r, "-BIOSYS-GSAfPCA-ADM1section-PC1to3var_RNG_ind_ALL.csv", sep = "");

# Save output data files as csv format  
write.csv(mu_RNG_ind_ALL_ADM1, quote = FALSE, outfile_mu_ADM1);
write.csv(mu.star_RNG_ind_ALL_ADM1,quote = FALSE, outfile_mustar_ADM1);
write.csv(sigma_RNG_ind_ALL_ADM1, quote = FALSE, outfile_sigma_ADM1);
write.csv(PC1to3var_RNG_ind_ALL_ADM1, quote = FALSE, outfile_PC1to3var_ADM1)







#Start of ADM1 to ASM3
   y1_j_all <- data.frame();
   
   for(j1 in 1:35){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
     
     component <- c('Ssu','Saa','Sfa','Sva','Sbu','Spro','Sac','Sh2','Sch4','Sic','Sin','Si','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','Xi','Scation','Sanion','Sva_m','Sbu_m','Spro_m','Sac_m','Shco3_m','Snh3','Sgas_h2','Sgas_ch4','Sgas_co2');
     y_range_start <- (j1*N) + 1; # the first j*N are the Time values, so skip them
     y_range_end <- (j1+1)*N;
     
     #yADM1 <- as.data.frame(yADM1)
     y1_j <- yADM1[y_range_end,]
     y1_j_all <- rbind(y1_j_all, y1_j)
     
     
   }
   
   NumCol1<-ncol(y1_j_all)
   colnames(y1_j_all) <- 1:NumCol1
   rownames(y1_j_all) <- c('Ssu','Saa','Sfa','Sva','Sbu','Spro','Sac','Sh2','Sch4','Sic','Sin','Si','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','Xi','Sva_m','Sbu_m','Spro_m','Sac_m','Shco3_m','Snh3','Sgas_h2','Sgas_ch4','Sgas_co2','Scation','Sanion');
   
   nth=1
   V1=y1_j_all['Ssu', nth]
   V2=y1_j_all['Saa', nth]
   V3=y1_j_all['Sfa', nth]
   V4=y1_j_all['Sva', nth]
   V5=y1_j_all['Sbu', nth]
   V6=y1_j_all['Spro', nth]
   V7=y1_j_all['Sac', nth]
   V8=y1_j_all['Sh2', nth]
   V9=y1_j_all['Sch4', nth]
   V10=y1_j_all['Sic', nth]
   V11=y1_j_all['Sin', nth]
   V12=y1_j_all['Si', nth]
   V13=y1_j_all['Xc', nth]
   V14=y1_j_all['Xch', nth]
   V15=y1_j_all['Xpr', nth]
   V16=y1_j_all['Xli', nth]
   V17=y1_j_all['Xsu', nth]
   V18=y1_j_all['Xaa', nth]
   V19=y1_j_all['Xfa', nth]
   V20=y1_j_all['Xc4', nth]
   V21=y1_j_all['Xpro', nth]
   V22=y1_j_all['Xac', nth]
   V23=y1_j_all['Xh2', nth]
   V24=y1_j_all['Xi', nth]
   V25=y1_j_all['Sva_m', nth]
   V26=y1_j_all['Sbu_m', nth]
   V27=y1_j_all['Spro_m', nth]
   V28=y1_j_all['Sac_m', nth]
   V29=y1_j_all['Shco3_m', nth]
   V30=y1_j_all['Snh3_m', nth]
   V31=y1_j_all['Sgas_h2', nth]
   V32=y1_j_all['Sgas_ch4', nth]
   V33=y1_j_all['Sgas_co2', nth]
   




###
# This script executes GSA-fPCA on the Activated Sludge Model No. 3 (ASM3)
# Code by Fortela DLB, et al. (2018)
# The Energy Institute, University of Louisiana, Lafayette, LA 70504 USA
###

# Load the R-packages
#library(deSolve); # To solve the ODE's
#
#directory = "/Users/Alyssa/Documents/LURA/ASM3";
#setwd(directory);
#
##-----Loop through random number sets RNG
#for(i in 1:RNG) {
#  
#  print(paste("Simulating the DAEs for i =", i, sep = " ")); # Print on screen to track progress
#  
#  ##--------Sample parameters using Morris method-----
#  set.seed(i)
#  
#  doe2 <- morris(model = NULL, factors = k, r = r,
#                 design = list(type = "oat", levels = 6, grid.jump = 3),
#                 binf = c(0.006375, 0.0255, 0.034, 0.0255, 0.0595, 0.0085, 0.17, 0.68, 0.5525, 0.5525, 0.51, 0.51, 0.17, 0.051, 8.0750, 10.2, 2.975, 0.8075, 0.6375, 0.255, 0.255, 0.1275, 0.187, 0.4675, 0.4675, 0.53125, 0.57375, 0.4675, 0.935, 2.14625, 2.14625, 10.625, 0.042925, 4.2925, 4.2925, 0.42925, 0.42925, 2.21, 2.21, 0.23375, 0.23375, 0.4675),
#                 bsup = c(0.008625, 0.0345, 0.046, 0.0345, 0.0805, 0.0115, 0.23, 0.92, 0.7475, 0.7475, 0.69, 0.69, 0.23, 0.069, 10.925, 13.8, 4.025, 1.0925, 0.8625, 0.345, 0.345, 0.1725, 0.253, 0.6325, 0.6325, 0.71875, 0.77625, 0.6325, 1.265, 2.77750, 2.77750, 14.375, 0.058075, 5.8075, 5.8075, 0.58075, 0.58075, 2.99, 2.99, 0.31625, 0.31625, 0.6325),
#                 scale = TRUE)
#  #parameters <- c(0.0075, 0.03, 0.04, 0.03, 0.07, 0.01, 0.2, 0.80, 0.65, 0.65 , 0.6, 0.6, 0.6, 0.2, 0.06, 9.5, 12, 3.5, 0.95, 0.75, 0.3, 0.3, 0.15, 0.22, 0.55, 0.55, 0.625, 0.675, 0.55, 1.1, 2.525, 2.525, 12.5, 0.0505, 5.05, 5.05, 0.505, 0.505, 2.6, 2.6, 0.275, 0.275, 0.55)
#  
#  #        iNSI ,iNSS ,iNXI ,iNXS , iNBM, fSI , fXI, yH02,yHNO3,yHNO2,ySTOO2,YSTONO3,YSTONO2,YAOB,YNOB, kH,kSTO,uH ,uAOB,uNOB,bHO2,bSTOO2,bAOB ,bNOB ,nHNO3,nHNO2,nHendNO3,nHendNO2,nNend,KX  ,KHO2 ,KHO2inh,KHSS,KHNH4 ,KHNO3,KHNO2,KHALK,KHSTO,KAOBO2,KNOBO2,KAOBNH4,KNOBNO2,KNALK
  
  

  #Conversion of ADM1 response values 
  #Need to convert kgCOD/m3 to gCOD/m3
   
   XS_sum = unname(1000*(V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23))
   SS_sum = unname(1000*(V1+V2+V3+V4+V5+V6+V7))
   SALK_sum = unname(1000*((V7+V6+V5+V4+V11+V10)/1)) #negative 1 (-1) denominator?
  
  ##--- Defining the ASM3 ----
  

  ASM3 <- function(t, state, parameters2) {
    with(as.list(c(state, parameters2)), { 
      
      
      # Model parameters
      doe <- parameters2
      
      iNSI = doe[1]
      iNSS = doe[2]
      iNXI = doe[3]
      iNXS = doe[4]
      iNBM = doe[5]
      fSI  = doe[6]
      fXI  = doe[7]
      
      YHO2    = doe[8]
      YHNO3   = doe[9]
      YHNO2   = doe[10]
      YSTOO2  = doe[11]
      YSTONO3 = doe[12]
      YSTONO2 = doe[13]
      YAOB    = doe[14]
      YNOB    = doe[15]
      
      kH   = doe[16]
      kSTO = doe[17]
      muH   = doe[18]
      muAOB = doe[19]
      muNOB = doe[20]

      bHO2=0.3
      bSTOO2=0.3
      bAOB=0.15
      bNOB=0.22
      nHNO3=0.5
      nHNO2=0.5
      nHendNO3=0.5
      nHendNO2=0.5
      nNend=0.5
      KX=1.0
      KHO2=0.2
      KHO2inh=0.2
      KHSS=10
      KHNH4=0.01
      KHNO3=0.5
      KHNO2=0.5
      KHALK=0.1
      KHSTO=0.1
      KAOBO2=0.8
      KNOBO2=0.8
      KAOBNH4=0.14
      KNOBNO2=0.28
      KNALK=0.5
    
      print(paste("ASM3: r=", r, "Done with model parameters", sep = " "));
      
#      bHO2     = doe2[21]
#      bSTOO2   = doe2[22]
#      bAOB     = doe2[23]
#      bNOB     = doe2[24]
#      nHNO3    = doe2[25]
#      nHNO2    = doe2[26]
#      nHendNO3 = doe2[27]
#      nHendNO2 = doe2[28]
#      nNend    = doe2[29]
      
#      KX      = doe2[30]
#      KHO2    = doe2[31]
#      KHO2inh = doe2[32]
#      KHSS    = doe2[33]
#      KHNH4   = doe2[34]
#      KHNO3   = doe2[35]
#      KHNO2   = doe2[36]
#      KHALK   = doe2[37]
#      KHSTO   = doe2[38]
#      KAOBO2  = doe2[39]
#      KNOBO2  = doe2[40]
#      KAOBNH4 = doe2[41]
#      KNOBNO2 = doe2[42]
#      KNALK   = doe2[43]
      
      
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
      j4  = muH*((SO2)/(KHO2+SO2))*((SNH4_asm3)/(KHNH4+SNH4_asm3))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*XH # Aerobic Growth of Xh
      j5a = muH*nHNO3*((KHO2inh)/(KHO2inh+SO2))*((SNH4_asm3)/(KHNH4+SNH4_asm3))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*((SNO3)/(KHO2inh+SNO3))*XH # Anoxic Growth of Xh NO3-NO2
      j5b = muH*nHNO2*((KHO2inh)/(KHO2inh+SO2))*((SNH4_asm3)/(KHNH4+SNH4_asm3))*((SALK)/(KHALK+SALK))*((XSTO/XH)/(KHSTO+XSTO/XH))*((SNO2)/(KHO2inh+SNO2))*XH # Anoxic Growth of Xh NO2-N2
      j6  = bHO2*((SO2)/(KHO2+SO2))*XH # Aerobic Endog. Resp. of Xh
      j7a = bHO2*nHendNO3*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XH # Anoxic Endog. Resp. of Xh NO3-NO2
      j7b = bHO2*nHendNO2*((KHO2inh)/(KHO2inh+SO2))*((SNO2)/(KHNO2+SNO2))*XH # Anoxic Endog. Resp. of Xh NO2-N2
      j8  = bSTOO2*((SO2)/(KHO2+SO2))*XSTO # Aerobic Endog. Resp. of Xsto
      j9a = bSTOO2*nHendNO3*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XSTO # Anoxic Endog. Resp. of Xsto NO3-NO2
      j9b = bSTOO2*nHendNO2*((KHO2inh)/(KHO2inh+SO2))*((SNO2)/(KHNO2+SNO2))*XSTO # Anoxic Endog. Resp. of Xsto NO2-N2
      j10a = muAOB*((SO2)/(KAOBO2+SO2))*((SNH4_asm3)/(KAOBNH4+SNH4_asm3))*((SALK)/(KNALK+SALK))*XAOB # Aerobic Growth of Xaob Nitritation
      j10b = muNOB*((SO2)/(KNOBO2+SO2))*((SNH4_asm3)/(KHNH4+SNH4_asm3))*((SNO2)/(KNOBNO2+SNO2))*((SALK)/(KNALK+SALK))*XNOB # Aerobic Growth of Xaob Nitritation
      j11a = bAOB*((SO2)/(KHO2+SO2))*XAOB # Aerobic Endog. Resp. of Xaob
      j11b = bNOB*((SO2)/(KHO2+SO2))*XNOB # Aerobic Endog. Resp. of Xnob
      j12a = bAOB*nNend*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XAOB # Anoxic Endog. Resp. of Xaob
      j12b = bNOB*nNend*((KHO2inh)/(KHO2inh+SO2))*((SNO3)/(KHNO3+SNO3))*XNOB # Anoxic Endog. Resp. of Xnob
      jO2 = kLA_O2*(SO2sat-SO2) # Aeration rate with sewage aeration parameters levels
      
      #print(paste("r=", r, "i=", i, "Done with the reaction rates", sep = " "));
      
      print(paste("ASM3: r=", r, "Evaluating the DAE's", sep = " "));
      
      # System of DAEs
      
      dSO2 = u21*j2+u51*j4+u81*j6+u111*j8+u141*j10a+u151*j11a+u171*j10b+u181*j11b+jO2
      dSS = u12*j1+u22*j2+u32*j3a+u42*j3b
      dSNH4_asm3 = u23*j2+u33*j3a+u43*j3b+u53*j4+u63*j5a+u73*j5b+u83*j6+u93*j7a+u103*j7b+u143*j10a+u153*j11a+u163*j12a+u173*j10b+
        u183*j11b+u193*j12b
      dSNO2 = u34*j3a+u44*j3b+u64*j5a+u74*j5b+u94*j7a+u104*j7b+u124*j9a+u134*j9b+u144*j10a+u174*j10b
      dSNO3 = u35*j3a+u65*j5a+u95*j7a+u125*j9a+u165*j12a+u175*j10b+u195*j12b
      dSN2 = u46*j3b+u76*j5b+u106*j7b+u136*j9b+u166*j12a+u196*j12b
      dSALK = u17*j1+u27*j2+u37*j3a+u47*j3b+u57*j4+u67*j5a+u77*j5b+u87*j6+u97*j7a+u107*j7b+u137*j9b+u147*j10a+u157*j11a+
        u167*j12a+u177*j10b+u187*j11b+u197*j12b
      dSI_asm3 = u18*j1
      dXI_asm3 = u89*j6+u99*j7a+u109*j7b+u159*j11a+u169*j12a+u189*j11b+u199*j12b
      dXH = u510*j4+u610*j5a+u710*j5b+u810*j6+u910*j7a+u1010*j7b
      dXS = u111*j1
      dXSTO = u212*j2+u312*j3a+u412*j3b+u512*j4+u612*j5a+u712*j5b+u1112*j8+u1212*j9a+u1312*j9b
      dXAOB = u1413*j10a+u1513*j11a+u1613*j12a
      dXNOB = u1714*j10b+u1814*j11b+u1914*j12b
      
      
      #print(paste("r=", r, "i=", i, "Done with the ODEs", sep = " "));
      
      return(list(c(dSO2, dSS, dSNH4_asm3, dSNO2, dSNO3, dSN2, dSALK, dSI_asm3, dXI_asm3, dXH, dXS, dXSTO, dXAOB, dXNOB)))
      
    })
  } 
  

  ##--------Define the simulation function w/ ode()------
  PCA_data2 <- function(doe) {
    
    
    # Specifying the constants
    #GSA_parameters = params
    ##states initial condition, liquid within the digester, not the input
    state2 = c(SO2=8, SS=SS_sum, SNH4_asm3=16, SNO2=0, SNO3=0, SN2=0, SALK=SALK_sum, SI_asm3=V12, XI_asm3=V24, XH=30, XS=XS_sum, XSTO=0.25, XAOB=0.25, XNOB=0.25)
    #state2 = c(SO2=8, SS=60, SNH4_asm3=16, SNO2=0, SNO3=0, SN2=0, SALK=5, SI_asm3=30, XI_asm3=25, XH=30, XS=115, XSTO=0.25, XAOB=0.25, XNOB=0.25)
    #
    parameters2=c(iNSI=doe[1], iNSS=doe[2], iNXI=doe[3], iNXS=doe[4], iNBM=doe[5], fSI=doe[6], fXI=doe[7], YHO2=doe[8],
                 YHNO3=doe[9], YHNO2=doe[10], YSTOO2=doe[11], YSTONO3=doe[12], YSTONO2=doe[13], YAOB=doe[14], YNOB=doe[15],
                 kH=doe[16], kSTO=doe[17], muH=doe[18], muAOB=doe[19], muNOB=doe[20])
    
    times2 <- seq(0, Period, length = N)
    
    # Solving the ASM3
    out2 <- ode(y = state2, times = times2, func = ASM3, parms = parameters2, method = "bdf", maxsteps=5000)
    #colnames(out2)[9]<- "SI_asm3"
    #colnames(out2)[10]<- "XI_asm3"
    return(out2)
    
  }
    #print(paste("r=", r, "i=", i, "Done evaluating ASM3 using ode()", sep = " "));
    
    ##--------Pass GSA doe parameters & generate the ASM3 output(s)-----
    yASM3 <-apply(doe$X[,23:42],1, PCA_data2);
    
     ##------Run fPCA for each DAE simulation output-------
     for(j in 1:14){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
       
       component <- c('SO2','SS','SNH4_asm3','SNO2','SNO3','SN2','SALK','SI_asm3','XI_asm3','XH','XS','XSTO','XAOB','XNOB');
       y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
       y_range_end <- (j+1)*N;
       
       #yASM3 <- as.data.frame(yASM3)
       
       y2_j <- yASM3[y_range_start:y_range_end,];
       
       # Plot all runs for j within i and save plots file to directory
       #plot(y_j[,1])
       setEPS(); # format as postscript for high-quality graphics
       mypath <- file.path(directory,paste("ASM3_Responses_","r", r, "i", i, "j", j, ".eps", sep=""));
       postscript(file=mypath, width = 5, height = 5);
       matplot_title <- paste(component[j], "Simulations for", "r =", r, "&","RNG seed =", i, sep = " ");
       matplot(x=seq(0, Period, length = N), y=y2_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], "(kgCOD/m^3)", sep = " "));
       dev.off();
       
       ##--------Run fPCA using fda package & perform Morris GSA on fPCA scores------
       # Turn the vectors of model values into functional objects
       yf_2j <- Data2fd(y2_j, seq(0, Period, length = N));
       #plot(yf)
       
       # Then calculate the first 3 PCs using the fda
       ypca_2j <- pca.fd(yf_2j,3);
    
           
  
      # Plot the first 3 PCAs to see the variation they describe
      #plot.pca.fd(ypca)

        setEPS(); # format as postscript for high-quality graphics
        mypath <- file.path(directory,paste("ASM3_PCsPlots_","r", r, "i", i, "j", j, ".eps", sep=""));
        postscript(file=mypath, width = 10, height = 3);
        op <- par(mfrow=c(1,3), oma = c(0,0,0,0));
        plot.pca.fd(ypca_2j);
        par(op);
        #title("My Title");
        dev.off();
   
       
       # Here only the look at the first 3 PCs as they describe pretty much all the variance
       #100*ypca$varprop;
       # Now pass the coefficients (or scores) of the PCA back to the Morris method
       # just consider the first 3 as they account for all the variance

        yASM3_scores <- tell(doe,ypca_2j$scores[,1:3]);
       
       
       ##--------Calculate the Morris GSA indices------
       
       mu_yASM3 <- apply(yASM3_scores$ee, 3, function(M){
         apply(M, 2, mean)
       });
       mu.star_yASM3 <- apply(abs(yASM3_scores$ee), 3, function(M){
         apply(M, 2, mean)
       });
       sigma_yASM3 <- apply(yASM3_scores$ee, 3, function(M){
         apply(M, 2, sd)
       });
       
       
       ##-----Summarize the GSA indices per parameter and PC
       X = c(1:k);
       Parameters = paste('Parm', 1:k, sep = "");
       
       mu_RNG_ind_ASM3 <- cbind(mu_yASM3,i,j,r,X,Parameters);
       mu.star_RNG_ind_ASM3 <- cbind(mu.star_yASM3,i,j,r,X,Parameters);
       sigma_RNG_ind_ASM3 <- cbind(sigma_yASM3,i,j,r,X,Parameters);
       PC1to3var_RNG_ind_ASM3 <- cbind(ypca_2j$varprop[1],ypca_2j$varprop[2],ypca_2j$varprop[3],i,j,r);
       
       
       mu_RNG_ind_ALL_ASM3 <- rbind(mu_RNG_ind_ALL_ASM3, mu_RNG_ind_ASM3);
       mu.star_RNG_ind_ALL_ASM3 <- rbind(mu.star_RNG_ind_ALL_ASM3, mu.star_RNG_ind_ASM3);
       sigma_RNG_ind_ALL_ASM3 <- rbind(sigma_RNG_ind_ALL_ASM3, sigma_RNG_ind_ASM3);
       PC1to3var_RNG_ind_ALL_ASM3 <- rbind(PC1to3var_RNG_ind_ALL_ASM3, PC1to3var_RNG_ind_ASM3);
       
     };


#Save the GSA-fPCA results for the ADM1 section

# File names for csv outputs
outfile_mu_ASM3 <- paste("r", r, "-BIOSYS-GSAfPCA-ASM3section-mu_RNG_ind_ALL.csv", sep = "");
outfile_mustar_ASM3 <- paste ("r", r, "-BIOSYS-GSAfPCA-ASM3section-mu_star_RNG_ind_ALL.csv", sep = "");
outfile_sigma_ASM3 <- paste("r", r, "-BIOSYS-GSAfPCA-ASM3section-sigma_RNG_ind_ALL.csv", sep = "");
outfile_PC1to3var_ASM3 <- paste("r", r, "-BIOSYS-GSAfPCA-ASM3section-PC1to3var_RNG_ind_ALL.csv", sep = "");

# Save output data files as csv format  
write.csv(mu_RNG_ind_ALL_ASM3, quote = FALSE, outfile_mu_ASM3);
write.csv(mu.star_RNG_ind_ALL_ASM3,quote = FALSE, outfile_mustar_ASM3);
write.csv(sigma_RNG_ind_ALL_ASM3, quote = FALSE, outfile_sigma_ASM3);
write.csv(PC1to3var_RNG_ind_ALL_ASM3, quote = FALSE, outfile_PC1to3var_ASM3)

    #   
    #   
    #   ##--------Plot Morris S.D. vs mu*-------
    #   #plot(x = mu.star[,1], y = sigma[,1], main = "PC1 Morris - r=40 RNGset.set.seed(1) - Try1.12/01/17", xlab = "mean, mu*", ylab = "S.D.", col = "yellow", cex = 1)
    #   #text(mu.star[,1], sigma[,1], labels = row.names(mu.star), cex = 0.7, offset = 15)
    #   
    #   
    #   X = c(1:k);
    #   Parameters = c("iNSI", "iNSS", "iNXI", "iNXS", "iNBM", "fSI", "fXI", "YHO2", "YHNO3", "YHNO2", "YSTOO2", "YSTONO3", "YSTONO2", "YAOB", "YNOB", "kH", "kSTO", "muH", "muAOB", "muNOB", "bHO2", "bSTOO2", "bAOB", "bNOB", "nHNO3", "nHNO2", "nHendNO3", "nHendNO2", "nNend", "KX", "KHO2", "KHO2inh", "KHSS", "KHNH4", "KHNO3", "KHNO2", "KHALK", "KHSTO", "KAOBO2","KNOBO2", "KAOBNH4", "KNOBNO2", "KNALK");
    #   mu_RNG_ind <- cbind(mu,i,j,r,X,Parameters);
    #   mu.star_RNG_ind <- cbind(mu.star,i,j,r,X,Parameters);
    #   sigma_RNG_ind <- cbind(sigma,i,j,r,X,Parameters);
    #   PC1to3var_RNG_ind <- cbind(ypca$varprop[1],ypca$varprop[2],ypca$varprop[3],i,j,r);
    #   
    #   
    #   mu_RNG_ind_ALL <- rbind(mu_RNG_ind_ALL, mu_RNG_ind);
    #   mu.star_RNG_ind_ALL <- rbind(mu.star_RNG_ind_ALL, mu.star_RNG_ind);
    #   sigma_RNG_ind_ALL <- rbind(sigma_RNG_ind_ALL, sigma_RNG_ind);
    #   PC1to3var_RNG_ind_ALL <- rbind(PC1to3var_RNG_ind_ALL, PC1to3var_RNG_ind);
    #   

##-------------Passing ASM3 outputs to Algae------------##       

y2_j_all <- data.frame();

for(j in 1:14){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
  
  component <- c('SO2','SS','SNH4_asm3','SNO2','SNO3','SN2','SALK','SI_asm3','XI_asm3','XH','XS','XSTO','XAOB','XNOB');
  y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
  y_range_end <- (j+1)*N;
  
  #yADM1 <- as.data.frame(yADM1)
  y2_j <- yASM3[y_range_end,]
  y2_j_all <- rbind(y2_j_all, y2_j)
  
  
}

NumCol2<-ncol(y2_j_all)
colnames(y2_j_all) <- 1:NumCol2
rownames(y2_j_all) <- c('SO2','SS','SNH4_asm3','SNO2','SNO3','SN2','SALK','SI_asm3','XI_asm3','XH','XS','XSTO','XAOB','XNOB');


      

      mth= 1
      V34=y2_j_all['SO2',mth]
      V35=y2_j_all['SS',mth]
      V36=y2_j_all['SNH4_asm3',mth]
      V37=y2_j_all['SNO2',mth]
      V38=y2_j_all['SNO3',mth]
      V39=y2_j_all['SN2',mth]
      V40=y2_j_all['SALK',mth]
      V41=y2_j_all['SI_asm3',mth]
      V42=y2_j_all['XI_asm3',mth]
      V43=y2_j_all['XH',mth]
      V44=y2_j_all['XS',mth]
      V45=y2_j_all['XSTO',mth]
      V46=y2_j_all['XAOB',mth]
      V47=y2_j_all['XNOB',mth]
      
      
      
#      # ALGAE Model
#      # Batch mode, check line 47-48
#      
#      library(deSolve)
#      library(fda)
#      library(sensitivity)
#      
#      # Set the working directory where outputs will be saved
#      directory = "/Users/dhan-lordfortela/Documents/EIL/Proposals/LURA 2019/LURA_working/trial6_batch_1day_r50_RNG_1to5";
#      setwd(directory);
#      
#      write.csv(seq(0,5,by=1), quote = FALSE, "testoutput.csv"); # Test directory is okay before the long computation stage
#      
#      
#      r = 50; # Setting the number of EE's for the Morris GSA
#      RNG = 5; # Number of Random Numer Generator (RNG) seeds
#      k = 20;
#      
#      mu_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu
#      mu.star_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's mu*
#      sigma_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for EE's sigma
#      PC1to3var_RNG_ind_ALL <- data.frame(); # Initialize the results data frame for PC1, PC2 & PC3 variances
#      
#      
#      
#      
#      #Period=10 #number of days
#      #N=100
#      ##-----Loop through random number sets RNG
#      for(i in 1:RNG) {
#        
#        print(paste("Simulating the DAEs for i =", i, sep = " ")); # Print on screen to track progress
#        
#        ##--------Sample parameters using Morris method-----
#        set.seed(i)
#        
#        doe3 <- morris(model = NULL, factors = k, r = r,
#                       design = list(type = "oat", levels = 6, grid.jump = 3),
#                       binf = c(1.36, 0.085, 0.085, 0.003672, 102, 0.085, 0.17, 0.0085, 1.7, 21.25, 11.05, 0.00164475, 0.000000491708, 0.1241, 0.00040766,  3.4, 0.595, 0.595, 20,     170),
#                       bsup = c(1.84, 0.115, 0.115, 0.004968, 138, 0.115, 0.23, 0.0115, 2.3, 28.75, 14.95, 0.00222525, 0.000000665252, 0.1679, 0.00055154,  4.6, 0.805, 0.805, 40,     230),
#                       #        ualg, kresp, kdeath, KC,     ICO2,  KN,   KO2,  Kpr,    Tau, Topt,   s,      alpha,     beta,          gamma,     delta,    KaO2, KaCO2, KaNH3, Tact,  I
#                       #        1       2     3       4       5       6     7     8       9  10     11        12        13              14        15        16    17     18    19      20 
#                       scale = TRUE)
#        
#        Period = 1 #Simulation time period in Days  
#        N = 100 # no. of simulation samples to take for GSA-fPCA
#        
#        ##--------Setup simulation conditions------
#        #Q = 0.075; # Continuous mode (CSTR), m^3/day
#        Q = 0.0; # Batch mode, m^3/day
#        Vliq=0.450; #m^3
#        Vgas=0.060; #m^3
#        tau=Q/Vliq; #/day
        
        
        #######---------Defining the Algae Model---------######
        ALGAE_model <- function(t, state, parameters) {
          with(as.list(c(state, parameters)), { 
            
            #model parameters
            
            doe <- parameters
            
            ualg = doe[1] #1.6 #d^-1
            kresp  = doe[2] #0.1 #d^-1
            kdeath = doe[3] #0.1 #d^-1
            KC = doe[4] #0.00432 #gC m^-3
            ICO2 = doe[5] #120 #gC m^-3
            KN = doe[6] #0.1 #gC m^-3
            KO2 = doe[7] #0.2 #gO2 m^-3
            Kpr = doe[8] #0.01 
            Tau = doe[9] #2 #originally 4
            Topt = doe[10] #25 #degrees celsius
            s = doe[11] #13
            alpha = doe[12] #0.001935 #(uE m^-2)^-1
            beta = doe[13] #0.00000057848 #(uE m^-2)^-1
            gamma = doe[14] #0.1460 #s^-1
            delta = doe[15] #0.0004796 #s^-1
            KaO2 = doe[16] #4 #d^-1
            KaCO2 = doe[17] #0.7 #d^-1
            KaNH3 = doe[18] #0.7 #d^-1
            Tact = doe[19] #13.5 #PLACEHOLDER, good range accdg to the paper is 15C to 25C
            I = doe[20] #200 #PLACEHOLDER, range is 3.25 to 665 uE/(m^2-s)
            
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
            
            print(paste("ALGAE: r=", r, "Done with model parameters", sep = " "));
            
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
            
            
            print(paste("ALGAE: r=", r, "Evaluating the DAE's", sep = " "));
            
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

      
      SNH4_alg_sum = unname(V36)
      SNO3_sum = unname(V38)
      SO2_sum = unname(V34)        
      
      
        ##--------Define the simulation function w/ ode()------
        PCA_data3 <- function(doe) {
          
          
          # Specifying the constants
          #GSA_parameters = params
          ##states initial condition, liquid within the digester, not the input
          
          #state3 <- c(SNH4_alg=SNH4_alg_sum, SNH3=0.685, SNO3=SNO3_sum,SO2=SO2_sum, SCO2=0.8, SHCO3=100, SCO3=1.17, SH=0.00000316, SOH=0.00283, Xalg=80)
          state3 <- c(SNH4_alg=1.50, SNH3=0.685, SNO3=9, SO2=1, SCO2=0.8, SHCO3=100, SCO3=1.17, SH=0.00000316, SOH=0.00283, Xalg=80)
          parameters3=c(ualg = doe[1], kresp  = doe[2], kdeath = doe[3], KC = doe[4], ICO2 = doe[5], KN = doe[6], KO2 = doe[7], Kpr = doe[8], Tau = doe[9], 
                       Topt = doe[10], s = doe[11], alpha = doe[12], beta = doe[13], gamma = doe[14], delta = doe[15], 
                       KaO2 = doe[16], KaCO2 = doe[17], KaNH3 = doe[18], Tact = doe[19], I = doe[20])
          times3 <- seq(0, Period, length = N)
          #######-----Solving the algae model--------#####
          out3<- ode(y = state3, times = times3, func = ALGAE_model, parms = parameters3, method = "bdf")
          
          return(out3)
          
        }
        
        #print(paste("r=", r, "i=", i, "Done evaluating algae using ode()", sep = " "));
        
        ##--------Pass GSA doe parameters & generate the algae output(s)-----
        yAlg <-apply(doe$X[,43:62],1, PCA_data3);
        
      
        ##------Run fPCA for each DAE simulation output-------
        for(j in 1:10){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
           
           component <- c("SNH4_alg", "SNH3", "SNO3", "SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg");
           y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
           y_range_end <- (j+1)*N;
           
           #yAlg <- as.data.frame(yAlg)
           
           y3_j <- yAlg[y_range_start:y_range_end,];
           
          # y_j_gsa <- yAlg[y_range_end,];
          # 
           # Plot all runs for j within i and save plots file to directory
           #plot(y_j[,1])
           setEPS(); # format as postscript for high-quality graphics
           mypath <- file.path(directory,paste("Algae_Responses_","r", r, "i", i, "j", j, ".eps", sep=""));
           postscript(file=mypath, width = 5, height = 5);
           matplot_title <- paste(component[j], "Simulations for", "r =", r, "&","RNG seed =", i, sep = " ");
           matplot(x=seq(0, Period, length = N), y=y3_j, type="l", main = matplot_title, xlab = "t (days)", ylab = paste(component[j], "(kgCOD/m^3)", sep = " "));
           dev.off();
           
          # print(paste("Performing GSA calcs for i =", i,", j =", j, sep = " "));
           ##--------Run fPCA using fda package & perform Morris GSA on fPCA scores------
           # Turn the vectors of model values into functional objects
          
           
           yf_3j <- Data2fd(y3_j, seq(0, Period, length = N));
           #plot(yf)
           
           # Then calculate the first 3 PCs using the fda
           ypca_3j <- pca.fd(yf_3j,3);
          
            
          # # Plot the first 2 PCAs to see the variation they describe
          # #plot.pca.fd(ypca)
          # 
          # 
          # 
            setEPS(); # format as postscript for high-quality graphics
            mypath <- file.path(directory,paste("Algae_PCsPlots_","r", r, "i", i, "j", j, ".eps", sep=""));
            postscript(file=mypath, width = 10, height = 3);
            op <- par(mfrow=c(1,3), oma = c(0,0,0,0));
            plot.pca.fd(ypca_3j);
            par(op);
            #title("My Title");
            dev.off();
           
           
          # 
          # # Here only the look at the first 3 PCs as they describe pretty much all the variance
          # #100*ypca$varprop;
          # # Now pass the coefficients (or scores) of the PCA back to the Morris method
          # # just consider the first 3 as they account for all the variance
           yAlg_scores <- tell(doe,ypca_3j$scores[,1:3]);
           
           
           ##--------Calculate the Morris GSA indices------
           
           mu_yAlg <- apply(yAlg_scores$ee, 3, function(M){
             apply(M, 2, mean)
           });
           mu.star_yAlg <- apply(abs(yAlg_scores$ee), 3, function(M){
             apply(M, 2, mean)
           });
           sigma_yAlg <- apply(yAlg_scores$ee, 3, function(M){
             apply(M, 2, sd)
           });
           
           
           ##-----Summarize the GSA indices per parameter and PC
           X = c(1:k);
           Parameters = paste('Parm', 1:k, sep = "");
           
           mu_RNG_ind_Algae <- cbind(mu_yAlg,i,j,r,X,Parameters);
           mu.star_RNG_ind_Algae <- cbind(mu.star_yAlg,i,j,r,X,Parameters);
           sigma_RNG_ind_Algae <- cbind(sigma_yAlg,i,j,r,X,Parameters);
           PC1to3var_RNG_ind_Algae <- cbind(ypca_3j$varprop[1],ypca_3j$varprop[2],ypca_3j$varprop[3],i,j,r);
           
           
           mu_RNG_ind_ALL_Algae <- rbind(mu_RNG_ind_ALL_Algae, mu_RNG_ind_Algae);
           mu.star_RNG_ind_ALL_Algae <- rbind(mu.star_RNG_ind_ALL_Algae, mu.star_RNG_ind_Algae);
           sigma_RNG_ind_ALL_Algae <- rbind(sigma_RNG_ind_ALL_Algae, sigma_RNG_ind_Algae);
           PC1to3var_RNG_ind_ALL_Algae <- rbind(PC1to3var_RNG_ind_ALL_Algae, PC1to3var_RNG_ind_Algae);
           
        };
      
      
      #Save the GSA-fPCA results for the ADM1 section
      
      # File names for csv outputs
      outfile_mu_Algae <- paste("r", r, "-BIOSYS-GSAfPCA-ALGAEsection-mu_RNG_ind_ALL.csv", sep = "");
      outfile_mustar_Algae <- paste ("r", r, "-BIOSYS-GSAfPCA-ALGAEsection-mu_star_RNG_ind_ALL.csv", sep = "");
      outfile_sigma_Algae <- paste("r", r, "-BIOSYS-GSAfPCA-ALGAEsection-sigma_RNG_ind_ALL.csv", sep = "");
      outfile_PC1to3var_Algae <- paste("r", r, "-BIOSYS-GSAfPCA-ALGAEsection-PC1to3var_RNG_ind_ALL.csv", sep = "");
      
      # Save output data files as csv format  
      write.csv(mu_RNG_ind_ALL_Algae, quote = FALSE, outfile_mu_Algae);
      write.csv(mu.star_RNG_ind_ALL_Algae,quote = FALSE, outfile_mustar_Algae);
      write.csv(sigma_RNG_ind_ALL_Algae, quote = FALSE, outfile_sigma_Algae);
      write.csv(PC1to3var_RNG_ind_ALL_Algae, quote = FALSE, outfile_PC1to3var_Algae)
        
      
           y3_j_all <- data.frame();
           
           for(j in 1:10){ # j = 1 to 24 are the 12 soluble components and 12 particulate components
             
             component <- c("SNH4_alg", "SNH3", "SNO3", "SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg");
             y_range_start <- (j*N) + 1; # the first j*N are the Time values, so skip them
             y_range_end <- (j+1)*N;
             
             #yADM1 <- as.data.frame(yADM1)
             y3_j <- yASM3[y_range_end,]
             y3_j_all <- rbind(y3_j_all, y3_j)
             
             
           }
           
           NumCol3<-ncol(y3_j_all)
           colnames(y3_j_all) <- 1:NumCol3
           rownames(y3_j_all) <- c("SNH4_alg", "SNH3", "SNO3", "SO2", "SCO2", "SHCO3", "SCO3", "SH", "SOH", "Xalg")
           
      
      #return(out3)
      #jth <- nrow(out3)
      #V48=yAlg[jth, 'SNH4_alg']
      #V49=yAlg[jth, 'SNH3']
      #V50=yAlg[jth, 'SNO3']
      #V51=yAlg[jth, 'SO2']
      #V52=yAlg[jth, 'SCO2']
      #V53=yAlg[jth, 'SHCO3']
      #V54=yAlg[jth, 'SCO3']
      #V55=yAlg[jth, 'SH']
      #V56=yAlg[jth, 'SOH']
      #V57=yAlg[jth, 'Xalg']
      
# # Combine all GSA-fPCA results from ADM1-ASM3-Algae sections
#            
#            mu_RNG_ind_ALL = mu_RNG_ind_ALL_Algae #Copy one full set for mu, e.g. Algae
#            mu.star_RNG_ind_ALL = mu.star_RNG_ind_ALL_Algae #Copy one full set for mu_star, e.g. Algae
#            sigma_RNG_ind_ALL = sigma_RNG_ind_ALL_Algae #Copy one full set for sigma, e.g. Algae
#            PC1to3var_RNG_ind_ALL = PC1to3var_RNG_ind_ALL_Algae #Copy one full set for PC's variations, e.g. Algae
#            
#            
#            mu_ADM1 <- mu_RNG_ind_ALL_ADM1[,c('V1','V2','V3')];
#            mu_ASM3 <- mu_RNG_ind_ALL_ASM3[,c('V1','V2','V3')];
#            mu_Algae <- mu_RNG_ind_ALL_Algae[,c('V1','V2','V3')];
#            
#            mu_ALL <- mu_ADM1 + mu_ASM3 + mu_Algae;
#            
#            
#            # File names for csv outputs
#            outfile_mu <- paste("r", r, "-BIOSYS-GSAfPCA-ALLsections-mu_RNG_ind_ALL.csv", sep = "");
#            outfile_mustar <- paste ("r", r, "-BIOSYS-GSAfPCA-ALLsections-mu_star_RNG_ind_ALL.csv", sep = "");
#            outfile_sigma <- paste("r", r, "-BIOSYS-GSAfPCA-ALLsections-sigma_RNG_ind_ALL.csv", sep = "");
#            outfile_PC1to3var <- paste("r", r, "-BIOSYS-GSAfPCA-ALLsections-PC1to3var_RNG_ind_ALL.csv", sep = "");
#            
#            # Save output data files as csv format  
#            write.csv(mu_RNG_ind_ALL, quote = FALSE, outfile_mu);
#            write.csv(mu.star_RNG_ind_ALL,quote = FALSE, outfile_mustar);
#            write.csv(sigma_RNG_ind_ALL, quote = FALSE, outfile_sigma);
#            write.csv(PC1to3var_RNG_ind_ALL, quote = FALSE, outfile_PC1to3var)
# 
