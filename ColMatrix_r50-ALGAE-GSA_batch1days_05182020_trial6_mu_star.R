

library(ggplot2)
library(RColorBrewer)
library(stats)
library(grid);
library(cowplot);
library(grDevices);
library(ggthemes);


setwd("/Users/dhan-lordfortela/Documents/EIL/Proposals/LURA 2019/LURA_working/Processed/trial6_batch_1day_r50_RNG_1to5")
mydata <- read.csv("Nominal-SigmaRanks-r50-ALGAE-GSA-FPCA_05172020_trial6-PC1.csv")
#j_names <- c('SO2', 'SS', 'SNH4', 'SNO2', 'SNO3', 'SN2', 'SALK', 'SI', 'XI', 'XH', 'XS', 'XSTO', 'XAOB', 'XNOB');
j_names <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10');
mydata$y1 <- j_names[mydata$j];


# levels(mydata$y1) <- list("SO2"="SO2", "SS"="SS", "SNH4"="SNH4", "SNO2"="SNO2", "SNO3"="SNO3", "SN2"="SN2", "SALK"="SALK", "SI"="SI", 
#                        "XI"="XI", "XH"="XH", "XS"="XS", "XSTO"="XSTO", "XAOB"="XAOB", "XNOB"="XNOB");

levels(mydata$y1) <- list("1"="01", "2"="02", "3"="03", "4"="04", "5"="05", "6"="06", "7"="07", "8"="08", 
                          "9"="09", "10"="10");

levels(mydata$Parameters) <- list("p20"="P20", "p19"="P19", "p18"="P18", "p17"="P17", "p16"="P16", "p15"="P15", "p14"="P14", "p13"="P13", "p12"="P12",
  "p11"="P11", "p10"="P10", "p9"="P9","p8"="P8","p7"="P7", "p6"="P6", "p5"="P5", "p4"="P4",  "p3"="P3","p2"="P2", "p1"="P1");

# Set the output graph format and name to be saved to the directory

tiff(file="r50-ALGAE-GSA-FPCA_batch1day_05182020_trial6_Sigma-PC1_purpleyellow.tiff", width = 3.5, height = 4.5, units = 'in', res = 300);



# There is no specific heatmap plotting function in ggplot2, but combining geom_tile with 
# a smooth gradient fill does the job very well.

(p <- ggplot(mydata, aes(x=Parameters, y=y1)) + 
  geom_tile(aes(fill = Rank), colour = "white") + 
  scale_fill_gradient(low = "dark blue", high = "orange", name = "Rank of p\nper S or X\nbased on\n\n", limits = c(1,20), breaks = c(1:20), guide = guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE)))
  # scale_fill_hue(c=45, l=80))
  # scale_fill_gradientn(colours = brewer.pal(22, "Paired")))
# A few finishing touches to the formatting, and the heatmap plot is ready for presentation.
base_size <- 9
p + theme_bw(base_size = base_size) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, ncol = 1, keywidth = 0.7, keyheight = 0.7)) +
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0), breaks = c(paste('p', 1:20, sep = "")), drop=FALSE) +
  #scale_y_discrete(expand = c(0, 0), breaks = c("SO2", "SS", "SNH4", "SNO2", "SNO3", "SN2", "SALK", "SI", 
  #                                              "XI", "XH", "XS", "XSTO", "XAOB", "XNOB"), drop=FALSE) +
  #scale_y_discrete(expand = c(0, 0), breaks = c('SO2', 'SS', 'SNH4', 'SNO2', 'SNO3', 'SN2', 'SALK', 'SI', 'XI', 'XH', 'XS', 'XSTO', 'XAOB', 'XNOB'), drop=FALSE) +
  scale_y_discrete(expand = c(0, 0), breaks = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10'), drop=FALSE) +
  geom_text(aes(label = mydata$Rank), size =2.5, col = "white") +
  
  coord_flip() +
  theme(plot.margin = unit(c(1,0.1,1.4,1.2), "cm"), 

        #axis.text.x = element_text(size=8.3, hjust=1, vjust =0.5, angle = 90),
        axis.text.x = element_blank(),
        #axis.text.y = element_text(size=8.3),
        axis.text.y = element_blank(),
        panel.spacing.y = unit(0.05,"lines"),
        strip.text.x = element_blank(),
        legend.position = "right") +

  facet_wrap(~PC, scales = "free_x", nrow = 1)


###
grid.text(expression(sigma), x=unit(0.81, "npc"), y=unit(0.785, "npc"), hjust =0, rot = 0, gp = gpar(fontsize = 20, fontface = "bold.italic", col = "red"));
#grid.text(expression('*'), x=unit(0.845, "npc"), y=unit(0.79, "npc"), hjust = 0, rot = 0, gp = gpar(fontsize = 20, fontface = "bold.italic", col = "red"));


#---- Annotate the y-labels (parameter names) ----
grid.text(expression(mu['alg']), x=unit(0.06, "npc"), y=unit(0.89, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p3
grid.text(expression('k'['resp']), x=unit(0.05, "npc"), y=unit(0.855, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p5
grid.text(expression('k'['death']), x=unit(0.05, "npc"), y=unit(0.82, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p6
grid.text(expression('K'['C']), x=unit(0.08, "npc"), y=unit(0.775, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p8
grid.text(expression('I'['CO2']), x=unit(0.06, "npc"), y=unit(0.735, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p10
grid.text(expression('K'['N']), x=unit(0.08, "npc"), y=unit(0.695, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p12
grid.text(expression('K'['O2']), x=unit(0.07, "npc"), y=unit(0.655, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p14
grid.text(expression('K'['pr']), x=unit(0.08, "npc"), y=unit(0.615, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p16
grid.text(expression(tau), x=unit(0.10, "npc"), y=unit(0.578, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9.5, fontface = "bold.italic")); # p17
grid.text(expression('T'['opt']), x=unit(0.06, "npc"), y=unit(0.54, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p19
grid.text(expression('s'), x=unit(0.10, "npc"), y=unit(0.50, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p21
grid.text(expression(alpha), x=unit(0.10, "npc"), y=unit(0.465, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9.5, fontface = "bold.italic")); # p23
grid.text(expression(beta), x=unit(0.10, "npc"), y=unit(0.425, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9.5, fontface = "bold.italic")); # p25
grid.text(expression(gamma), x=unit(0.10, "npc"), y=unit(0.385, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9.5, fontface = "bold.italic")); # p27
grid.text(expression(delta), x=unit(0.10, "npc"), y=unit(0.345, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9.5, fontface = "bold.italic")); # p29
grid.text(expression('Ka'['O2']), x=unit(0.05, "npc"), y=unit(0.305, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p34
grid.text(expression('Ka'['C02']), x=unit(0.04, "npc"), y=unit(0.265, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p36
grid.text(expression('Ka'['NH3']), x=unit(0.04, "npc"), y=unit(0.225, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p38
grid.text(expression('T'['act']), x=unit(0.07, "npc"), y=unit(0.19, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # p40
grid.text(expression('I'), x=unit(0.11, "npc"), y=unit(0.15, "npc"), rot = 0, hjust = 0, gp = gpar(fontsize = 9, fontface = "bold.italic")); # p42


#---- Annotate Variables names ----

#j_names <- c('(a) SNH4_alg','(b) SNH3','(c) SNO3','(d) SO2','(e) SCO2','(f) SHCO3','(g) SCO3','(h) SH','(i) SOH','(j) Xalg');
grid.text(expression('S'['NH4']), x=unit(0.165, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SO2
grid.text(expression('S'['NH3']), x=unit(0.225, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNH4
grid.text(expression('S'['NO3']), x=unit(0.29, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO2
grid.text(expression('S'['O2']), x=unit(0.355, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SN2
grid.text(expression('S'['CO2']), x=unit(0.415, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SALK
grid.text(expression('S'['HCO3']), x=unit(0.475, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SI
grid.text(expression('S'['CO3']), x=unit(0.53, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XH
grid.text(expression('S'['H']), x=unit(0.60, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XS
grid.text(expression('S'['OH']), x=unit(0.66, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XAOB
grid.text(expression('X'['alg']), x=unit(0.72, "npc"), y=unit(0.11, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XNOB

# grid.text(expression('S'['NH4']), x=unit(0.49, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SO2
# grid.text(expression('S'['NH3']), x=unit(0.525, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNH4
# grid.text(expression('S'['NO3']), x=unit(0.565, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO2
# grid.text(expression('S'['O2']), x=unit(0.60, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SN2
# grid.text(expression('S'['CO2']), x=unit(0.64, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SALK
# grid.text(expression('S'['HCO3']), x=unit(0.675, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SI
# grid.text(expression('S'['CO3']), x=unit(0.715, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XH
# grid.text(expression('S'['H']), x=unit(0.75, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XS
# grid.text(expression('X'['AOB']), x=unit(0.79, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XAOB
# grid.text(expression('X'['alg']), x=unit(0.825, "npc"), y=unit(0.06, "npc"), rot = 90, hjust = 1, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XNOB


#---- PCs projection percentage ----
# PC1
#grid.text(expression('PC1'), x=unit(0.30, "npc"), y=unit(0.93, "npc"), rot = 00, gp = gpar(fontsize = 12.5, fontface = "bold.italic"));
#grid.lines(x=unit(c(0.1,0.45), "npc"), y=unit(c(0.975,0.975), "npc"));
#j_names <- c('SO2', 'SS', 'SNH4', 'SNO2', 'SNO3', 'SN2', 'SALK', 'SI', 'XI', 'XH', 'XS', 'XSTO', 'XAOB', 'XNOB');
# grid.text(expression('65.8'['%']), x=unit(0.105, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SO2
# grid.text(expression('80.7'['%']), x=unit(0.133, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SS
# grid.text(expression('87.3'['%']), x=unit(0.16, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNH4
# grid.text(expression('98.7'['%']), x=unit(0.186, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO2
# grid.text(expression('98.8'['%']), x=unit(0.213, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO3
# grid.text(expression('97.8'['%']), x=unit(0.24, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SN2
# grid.text(expression('85.7'['%']), x=unit(0.265, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SALK
# grid.text(expression('98.3'['%']), x=unit(0.29, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SI
# grid.text(expression('99.1'['%']), x=unit(0.316, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XI
# grid.text(expression('84.0'['%']), x=unit(0.345, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XH
# grid.text(expression('88.2'['%']), x=unit(0.37, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XS
# grid.text(expression('72.1'['%']), x=unit(0.394, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XST0
# grid.text(expression('99.7'['%']), x=unit(0.422, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XAOB
# grid.text(expression('99.0'['%']), x=unit(0.448, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XNOB

# PC2
#grid.text(expression('PC2'), x=unit(0.67, "npc"), y=unit(0.93, "npc"), rot = 00, gp = gpar(fontsize = 12.5, fontface = "bold.italic"));
#grid.lines(x=unit(c(0.48,0.82), "npc"), y=unit(c(0.975,0.975), "npc"));
#j_names <- c('SO2', 'SS', 'SNH4', 'SNO2', 'SNO3', 'SN2', 'SALK', 'SI', 'XI', 'XH', 'XS', 'XSTO', 'XAOB', 'XNOB');
# grid.text(expression('15.6'['%']), x=unit(0.485, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SO2
# grid.text(expression('13.5'['%']), x=unit(0.513, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SS
# grid.text(expression('  9.9'['%']), x=unit(0.537, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNH4
# grid.text(expression('  1.1'['%']), x=unit(0.563, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO2
# grid.text(expression('  1.0'['%']), x=unit(0.59, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SNO3
# grid.text(expression('  1.8'['%']), x=unit(0.615, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SN2
# grid.text(expression('11.9'['%']), x=unit(0.642, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SALK
# grid.text(expression('  1.5'['%']), x=unit(0.665, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # SI
# grid.text(expression('  0.8'['%']), x=unit(0.693, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XI
# grid.text(expression('12.9'['%']), x=unit(0.72, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XH
# grid.text(expression('  8.7'['%']), x=unit(0.745, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XS
# grid.text(expression('19.8'['%']), x=unit(0.771, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XST0
# grid.text(expression('  0.2'['%']), x=unit(0.80, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XAOB
# grid.text(expression('  0.9'['%']), x=unit(0.825, "npc"), y=unit(0.93, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # XNOB


##----- Arrow ----

grid.lines(x=unit(c(0.93,0.93), "npc"), y=unit(c(0.68,0.73), "npc"), gp = gpar(fill = "black"), arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed"));
grid.text(expression('Increasing rank (increasing GSA index)'), x=unit(0.93, "npc"), y=unit(0.21, "npc"), rot = 90, hjust = 0, gp = gpar(fontsize = 8.5, fontface = "bold.italic")); # Arrow note
grid.lines(x=unit(c(0.93,0.93), "npc"), y=unit(c(0.13,0.20), "npc") );

dev.off()