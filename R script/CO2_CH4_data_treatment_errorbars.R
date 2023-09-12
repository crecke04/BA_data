#Plots for Substrates Charlotte Recke
library(ggplot2)
library(reshape2)
library(hues)
library(RColorBrewer)
library(tidyr)

setwd("H:/Desktop/PM4_qPCR/Plot")

display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 8, name = "Dark2")

##### H2 
H2 <- read.table("h2.txt", sep="\t", header=T, check.names = FALSE, dec=",")

names(H2) <- c("Sample", "Replicat", "Subsample", "Day", "h2", "2h","H2")

str(H2)

H2$Sample <- factor(H2$Sample, levels=c("C", "B", "LB01",
                                         "LB05",
                                        "BHI01",
                                        "BHI05",
                                        "BHI05LB05"),
                          labels = c("Widdel medium Control",
                                     "Control with Sample",
                                     "LB 0.1%", "LB 0.5%",
                                     "BHI 0.1%",
                                     "BHI 0.5%",
                                     "BHI 0.5% + LB 0.5%"), ordered=T)



cols=c("#1D3557","#E5383B","#457B9D", "#FDBF6F")
cols2=c("#E63946","#457B9D")

ggplot(H2, aes(x=Day, y=H2, col=Replicat, shape=Replicat)) +
  geom_line(aes(group=Subsample),size=.75 ) +
  geom_point(size=2) +
  scale_color_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0,size = 12,face="bold"),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0,size = 12),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size=23, face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title= element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = NA, fill = NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 80)) +
  labs(title ="H2 concentration in liquid",
       x = "Incubation time [d]",
       y = "Concentration [µM]",
       colour="Replicate",
       shape="Replicate") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA")),shape = guide_legend(override.aes = list(color = c("#1D3557","#E5383B","#457B9D")))) +
  facet_wrap(~ Sample, scales = "free", nrow=4) +
  theme(strip.background = element_rect(fill = "White", colour = "black"),
        strip.text = element_text(angle = 0, face = "bold", size = 12),
        panel.spacing = unit(0.2,"lines") ) 


##### SO4
SO4<- read.table("SO4.txt", sep="\t", header=T, check.names = FALSE, dec=",")

names(SO4) <- c("Value", "Replicat","Subsample", "Day", "Sample")

str(SO4)

SO4$Sample <- factor(SO4$Sample, levels=c("C",
                                          "B",
                                          "LB 01",
                                          "LB05",
                                          "BHI01",
                                          "BHI05",
                                          "BHI05LB05"),
                     labels = c("Widdel medium Control",
                                "Control with Sample",
                                " LB 0.1%", "LB 0.5 %",
                                "BHI 0.1%",
                                "BHI 0.5%",
                                "BHI 0.5 % + LB 0.5%"), ordered=T)

ggplot(SO4, aes(x=Day, y=Value, col=Replicat, shape=Replicat)) +
  geom_smooth(aes(group=Subsample),size=.75 , method= lm, formula= y~x) +
  geom_point(size=2) +
  scale_color_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 0, hjust=0.0, vjust=0,size = 12,face="bold"),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0,size = 12),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size=23, face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title= element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = NA, fill = NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 80)) +
  labs(title ="SO4 concentration",
       x = "Incubation time [d]",
       y = "Concentration [mM]",
       colour = "Replicates",
       shape="Replicates") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA"))) +
  facet_wrap(~ Sample, scales = "free", nrow=4) +
  theme(strip.background = element_rect(fill = "White", colour = "black"),
        strip.text = element_text(angle = 0, face = "bold", size = 12),
        panel.spacing = unit(0.2,"lines") ) 
#### gen 5 ch4
CH4_5 <- read.table("ch4_gen5.txt", sep="\t", header=T, check.names = FALSE, dec=",")
str(CH4_5)
names(CH4_5) <- c("Sample1", "Gas","Value1")
CH4_5$Sample1 <- factor(CH4_5$Sample1, levels=c("C1","C2","C3",
                                            "B1","B2","B3",
                                            "LB01A", "LB01B","LB01C",
                                           "LB05A", "LB05B","LB05C",
                                            "BHI01A","BHI01B","BHI01C",
                                            "BHI05A","BHI05B", "BHI05C",
                                            "LB05BHI05"),
                      labels = c("Widdel medium Control A","Widdel medium Control B","Widdel medium Control C",
                                 "Control with Sample A","Control with Sample B","Control with Sample c",
                                 " LB 0.1% A"," LB 0.1% B"," LB 0.1% c", "LB 0.5 % A","LB 0.5 % B", "LB 0.5 % C",
                                 "BHI 0.1% A","BHI 0.1% B","BHI 0.1% C",
                                 "BHI 0.5% A","BHI 0.5% B", "BHI 0.5% C",
                                 "BHI 0.5 % + LB 0.5%"), ordered=T)

ggplot(CH4_5, aes(x=Sample1, y=Value1, col=Gas, shape=Gas)) +
  geom_point(size=1.5) +
  scale_color_manual(values=cols2)+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.0,size = 9,face="bold"),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17, face="bold",hjust=0.0),
        legend.title=element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 14000)) +
  labs(title ="CH4 and CO2 measurements",
       x = "Sample",
       y = "Concentration [µM]") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA")))


#### Ch4 and Co2
CH4 <- read.table("CH4.txt", sep="\t", header=T, check.names = FALSE, dec=",")
CH4_2 <- read.table("CH4_2.txt", sep="\t", header=T, check.names = FALSE, dec=",")

str(CH4)
str(CH4_2)
names(CH4) <- c("Sample1", "Gas","Value1")
names(CH4_2) <- c("Sample2", "Gas","Value2")

CH4$Sample1 <- factor(CH4$Sample1, levels=c("WM C1",
                                        "WM C2",
                                        "WM C3",
                                        "1A","1B","1C",
                                        "2A","2B","2C",
                                        "3A","3B","3C",
                                        "4A","4B","4C",
                                        "5A","5B","5C",
                                        "6A","6B","6C"),
                     labels = c("C1",
                               "C2",
                               "C3",
                               "1A","1B","1C",
                               "2A","2B","2C",
                               "3A","3B","3C",
                               "4A","4B","4C",
                               "5A","5B","5C",
                               "6A","6B","6C"), ordered=T)

CH4_2$Sample2 <- factor(CH4_2$Sample2, levels=c("WM1","WM2", "WM3",
                                            "E3-1A","E3-1B","E3-1C",
                                            "E3-2A","E3-2B","E3-2C"),
                                  labels = c("C1","C2", "C3",
                                             "E3-1A","E3-1B","E3-1C",
                                             "E3-2A","E3-2B","E3-2C"), ordered=T)

ggplot(CH4, aes(x=Sample1, y=Value1, color=Gas, shape=Gas)) +
  geom_bar(aes(color = Gas, fill = Gas),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  labs(title ="CH4 and CO2 measurements generation 4",
       x = "Sample",
       y = "Concentration [µM]",
       color="Gas",
       shape="Gas") +
  scale_color_manual(values=cols2)+
  scale_fill_manual("Gas", values=cols2) +
 theme(strip.background = element_rect(fill="white"),
       axis.text.x = element_text(size=12),
       axis.text.y =element_text(size=12),
       axis.ticks.x = element_blank(),
       axis.title.x = element_text(),
       axis.title.y = element_text(),
       plot.title = element_text(size=17, face="bold"),
       panel.background = element_rect(fill = "white"),
       panel.border = element_rect(colour = "black", fill = "NA"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position = "right",
       legend.title = element_text(size=12, face="bold"),
       legend.text = element_text(size=12),
       legend.key = element_rect(colour = NA, fill = NA),
       legend.key.size = unit(0.35, "cm"))+
  scale_y_continuous(expand=c(0,0),limits = c(0, 12000))


ggplot(CH4, aes(x=Sample1, y=Value1, col=Gas, shape=Gas)) +
  geom_point(size=3) +
  scale_color_manual(values=cols2)+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0,size = 9,face="bold"),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17, face="bold"),
        legend.title=element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12000)) +
  labs(title ="CH4 and CO2 measurements generation 4",
       x = "Sample",
       y = "Concentration [µM]") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA")))


ggplot(CH4_2, aes(x=Sample2, y=Value2, col=Gas, shape=Gas)) +
  geom_point(size=3) +
  scale_color_manual(values=cols2)+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0,size = 9,face="bold"),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17, face="bold"),
        legend.title=element_text(face="bold"),
        legend.position = "right",
        legend.text = element_text(size=9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10000)) +
  labs(title ="CH4 and CO2 measurements generation 3",
       x = "Sample",
       y = "Concentration [µM]") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA")))


ggplot(CH4_2, aes(x=Sample2, y=Value2, color=Gas)) +
  geom_bar(aes(color = Gas, fill = Gas),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  labs(title ="CH4 and CO2 measurements generation 3",
       x = "Sample",
       y = "Concentration [µM]",
       color="Gas") +
  scale_color_manual(values=cols2)+
  scale_fill_manual("Gas", values=cols2) +
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text( face="bold"),
        axis.title.y = element_text( face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  scale_y_continuous(expand=c(0,0),limits = c(0, 10000))

##### Protein
Protein<- read.table("Protein.txt", sep="\t", header=T, check.names = FALSE, dec=",")

names(Protein) <- c("Subsample", "Day","value", "SD","Sample", "Replicate")

str(Protein)

Protein$Sample <- factor(Protein$Sample, levels=c("C",
                                          "B",
                                          "LB01",
                                          "LB05",
                                          "BHI01",
                                          "BHI05",
                                          "BHI05LB05"),
                     labels = c("Widdel medium Control",
                                "Control with Sample",
                                " LB 0.1%", "LB 0.5 %",
                                "BHI 0.1%",
                                "BHI 0.5%",
                                "BHI 0.5 % + LB 0.5%"), ordered=T)

ggplot(Protein, aes(x=Day, y=value, col=Replicate, shape=Replicate)) +
  geom_line(aes(group=Subsample),size=.75 ) +
  geom_point(size=1.5) +  
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD), lty="solid", width=.75)+
  scale_color_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0,size = 12,face="bold"),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0,size = 12),
        axis.title = element_text(face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size=23, face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title= element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = NA, fill = NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,7)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 80)) +
  labs(title ="Protein concentration",
       x = "Incubation time [d]",
       y = "Concentration [mg/ml]",
       colour = "Replicates",
       shape="Replicates") +
  guides(fill=guide_legend(ncol=1), color=guide_legend(override.aes=list(fill="NA"))) +
  facet_wrap(~ Sample, scales = "free", nrow=4) +
  theme(strip.background = element_rect(fill = "White", colour = "black"),
        strip.text = element_text(angle = 0, face = "bold", size = 12),
        panel.spacing = unit(0.2,"lines") )
