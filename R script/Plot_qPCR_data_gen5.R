library(ggplot2)
library(reshape2)
library(hues)
library(RColorBrewer)
library(tidyr)
library(patchwork)
library(gridExtra)
library(grid)

setwd("H:/Desktop/PM4_qPCR/Plots")

qPCR <- read.table("Gen5_qpcr.txt", sep="\t", header=T, check.names = FALSE, dec=",")
names(qPCR) <- c("Taxa","Sample","Day","value", "SD","Treatment", "Replicate")

str(qPCR)
qPCR$Treatment <- factor(qPCR$Treatment, levels=c( "B", "LB 01",
                                            "LB05",
                                            "BHI01",
                                            "BHI05",
                                            "BHI05LB05"),
                      labels = c("Control with Sample",
                                 "LB 0.1%", "LB 0.5%",
                                 "BHI 0.1%",
                                 "BHI 0.5%",
                                 "BHI 0.5 % + LB 0.5%"), ordered=T)

qPCR$Taxa <- factor(qPCR$Taxa, levels=c("Loki2b","Archaea","Izemo","Bacteria", "DSR"), ordered=T)

qPCR$Day <- factor(qPCR$Day, levels=c("14",
                                      "28",
                                      "42",
                                      "74"), ordered= T)

cols=c("#457B9D","#1D3557","#E5383B","#660708","#FDBF6F")

p1 <-  ggplot(qPCR[qPCR$Treatment == "Control with Sample",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
             lty="solid", width=0.2)+
  ggtitle("Control with Sample")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  facet_grid(~Replicate)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)

p2 <-  ggplot(qPCR[qPCR$Treatment == "LB 0.1%",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid", width=0.2)+  ggtitle("LB 0.1%")+
  labs( y = "Gene copies / mL Slurry (log10 scale)") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  facet_grid(~Replicate)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)

p3 <-ggplot(qPCR[qPCR$Treatment == "LB 0.5%",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid", width=0.2)+  ggtitle("LB 0.5%")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  facet_grid(~Replicate)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)

p4 <-ggplot(qPCR[qPCR$Treatment == "BHI 0.1%",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid", width=0.2)+  ggtitle("BHI 0.1%")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  facet_grid(~Replicate)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)

p5 <- ggplot(qPCR[qPCR$Treatment == "BHI 0.5%",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid", width=0.2)+  ggtitle("BHI 0.5%")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  facet_grid(~Replicate)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)

p6 <-  ggplot(qPCR[qPCR$Treatment == "BHI 0.5 % + LB 0.5%",], 
             aes(x = Day, y = value, color = Taxa, shape= Taxa))+ 
  geom_line(aes(group=Taxa),size=.5 ) +
  geom_point(size=2)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid", width=0.2)+
  labs( y = "Gene copies / mL Slurry (log10 scale)") +
  ggtitle("BHI 0.5% + LB 0.5%")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(title = element_text(size=10),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=9),
        axis.text.y =element_text(size=9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.35, "cm"))+
  scale_fill_manual("Taxa", values=cols)+ 
  scale_y_continuous(expand=c(0,0),limits = c(1,100000000000000), trans="log10",breaks=c(10,1e03,1e05,1e07,1e09,1e11,1e13), labels = scales::scientific)


(p1 + p2 + p3) / (p4 + p5 +( p6 +plot_spacer()+ plot_spacer()))+plot_layout(guides="collect")

title=textGrob("Rt-qPCR results from the 5th generation", gp=gpar(fontface="bold", fontsize=23))
y=textGrob("Gene copies / mL Slurry (log10 scale)",rot = 90, gp=gpar(fontface="bold"))
x=textGrob("Incubation day", gp=gpar(fontface="bold"), )
grid.arrange(patchworkGrob((p1 + p2 + p3) / (p4 + p5 +( p6 +plot_spacer()+ plot_spacer()))+plot_layout(guides="collect")), left = y, bottom= x, top= title, padding = unit(1, "line"), widths=1)
