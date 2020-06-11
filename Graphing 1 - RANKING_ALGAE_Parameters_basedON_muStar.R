# This is coded for the Algae model, but you may modify it for AMD1 and ASM3 data analysis
# Algae-GSA_indices: Script for the avrages of mu_star values and calc of their rankings

# Load the data from directory
setwd("..."); # paste the path to your working directory into "..."
mydata <- read.csv("r50-ALGAE-GSA-FPCA_batch1days_05192020_trial6_mu_star_ALL.csv"); # paste the filename of the mu_Star GSA
mydata <- mydata[,-1];

# Set analysis conditions/givens
Run = "Nominal";
Subs_Var = "Nominal";
Level = "Nominal";
r=50;
k=20;
RNG=5;
Total_runs = r*(k+1);

# Create the columns for PC calcs
mydata["PC1_TOT"] <- mydata$V1*Total_runs; #consider only the PC1
#mydata["PC2_TOT"] <- mydata$V2*Total_runs; #consider only the PC2
#mydata["PC3_TOT"] <- mydata$V3*Total_runs; #consider only the PC3

# Initialize the PC1 calc variables
PC1_mustar_AVE_perParam <- c()
PC1_mustar_AVE_perj <- c()
PC1_mustar_AVE <- c()

# # Initialize the PC2 calc variables
# PC2_mustar_AVE_perParam <- c()
# PC2_mustar_AVE_perj <- c()
# PC2_mustar_AVE <- c()

# # Initialize the PC3 calc variables
# PC3_mustar_AVE_perParam <- c()
# PC3_mustar_AVE_perj <- c()
PC3_mustar_AVE <- c()

Parameters = paste('P', 1:k, sep = "");

# Process the PC1 mu_stars
for (l in 1:10) {
  for (m in 1:k) {
    PC1_mustar_AVE_perParam[m] <- sum(mydata[which(mydata$X==m & mydata$j==l), "PC1_TOT"])/(Total_runs*RNG);
    PC1_mustar_AVE_perj <- as.data.frame(cbind("j" = l, "PC_AVE_mustar" = PC1_mustar_AVE_perParam, "PC"=1));
  }
  PC1_mustar_AVE_perj["Rank"] <- rank(-PC1_mustar_AVE_perj$PC_AVE, ties.method = "random");
  PC1_mustar_AVE_perj_named <- cbind("X" = c(1:m), PC1_mustar_AVE_perj, Subs_Var, Level, Run, Parameters);
  PC1_mustar_AVE <-rbind(PC1_mustar_AVE, PC1_mustar_AVE_perj_named);
}


# # Process the PC2 mu_stars
# for (l in 1:10) {
#   for (m in 1:k) {
#     PC2_mustar_AVE_perParam[m] <- sum(mydata[which(mydata$X==m & mydata$j==l), "PC2_TOT"])/(Total_runs*RNG);
#     PC2_mustar_AVE_perj <- as.data.frame(cbind("j" = l, "PC_AVE_mustar" = PC2_mustar_AVE_perParam, "PC"=2));
#   }
#   PC2_mustar_AVE_perj["Rank"] <- rank(-PC2_mustar_AVE_perj$PC_AVE, ties.method = "random");
#   PC2_mustar_AVE_perj_named <- cbind("X" = c(1:m), PC2_mustar_AVE_perj, Subs_Var, Level, Run, Parameters);
#   PC2_mustar_AVE <-rbind(PC2_mustar_AVE, PC2_mustar_AVE_perj_named);
# }
# 
# 
# # Process the PC2 mu_stars
# for (l in 1:10) {
#   for (m in 1:k) {
#     PC3_mustar_AVE_perParam[m] <- sum(mydata[which(mydata$X==m & mydata$j==l), "PC3_TOT"])/(Total_runs*RNG);
#     PC3_mustar_AVE_perj <- as.data.frame(cbind("j" = l, "PC_AVE_mustar" = PC3_mustar_AVE_perParam, "PC"=3));
#   }
#   PC3_mustar_AVE_perj["Rank"] <- rank(-PC3_mustar_AVE_perj$PC_AVE, ties.method = "random");
#   PC3_mustar_AVE_perj_named <- cbind("X" = c(1:m), PC3_mustar_AVE_perj, Subs_Var, Level, Run, Parameters);
#   PC3_mustar_AVE <-rbind(PC3_mustar_AVE, PC3_mustar_AVE_perj_named);
# }
# 
# 
# # Append all PC1, PC2 & PC3 into single data frame
# PC12_mustar_AVE <- rbind(PC1_mustar_AVE, PC2_mustar_AVE);
# PC123_mustar_AVE <- rbind(PC1_mustar_AVE, PC2_mustar_AVE, PC3_mustar_AVE);


# File names for csv outputs
PC1_outfile_name <- paste(Run, "-muStarRanks","-r50-ALGAE-GSA-FPCA_05172020_trial6-PC1.csv", sep = "");
# PC2_outfile_name <- paste(Run, "-muStarRanks","-r50-ALGAE-GSA-FPCA_05172020_trial4-PC2.csv", sep = "");
# PC3_outfile_name <- paste(Run, "-muStarRanks","-r50-ALGAE-GSA-FPCA_05172020_trial4-PC3.csv", sep = "");
# PC12_outfile_name <- paste(Run, "-muStarRanks","-r50-ALGAE-GSA-FPCA_05172020_trial4-PC12.csv", sep = "");
# PC123_outfile_name <- paste(Run, "-muStarRanks","-r50-ALGAE-GSA-FPCA_05172020_trial4-PC123.csv", sep = "");


# Save output data files as csv format  
write.csv(PC1_mustar_AVE, quote = FALSE, PC1_outfile_name)
# write.csv(PC2_mustar_AVE, quote = FALSE, PC2_outfile_name)
# write.csv(PC3_mustar_AVE, quote = FALSE, PC3_outfile_name)
# write.csv(PC12_mustar_AVE, quote = FALSE, PC12_outfile_name)
# write.csv(PC123_mustar_AVE, quote = FALSE, PC123_outfile_name)

