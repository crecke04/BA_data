library(ggplot2)
library(reshape2)
library(hues)
library(RColorBrewer)
library(tidyr)
library(patchwork)
library(gridExtra)
library(grid)

setwd("H:/Desktop/PM4_qPCR/Plots")
 
qPCR <- read.table("Gen4_qpcr.txt", sep="\t", header=T, check.names = FALSE)
names(qPCR) <- c("Taxa","Sample","Day","value", "SD","Treatment")

str(qPCR)

qPCR$Sample <- factor(qPCR$Sample, levels=c("1A","1B","1C",
                                            "2A","2B","2C",
                                            "3A","3B","3C",
                                            "4A","4B","4C",
                                            "5A","5B","5C",
                                            "6A","6B","6C"), ordered=T)
                      

qPCR$Day <- factor(qPCR$Day, levels=c("0",
                                      "21",
                                      "35"), ordered= T)

qPCR$Treatment <- factor(qPCR$Treatment, levels=c("1","2","3","4","5","6"), ordered= T)

cols=c("#1D3557","#E5383B","#457B9D")

p1 <- ggplot(qPCR[qPCR$Treatment == "1",], 
       aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 1")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)


p2 <- ggplot(qPCR[qPCR$Treatment == "2",], 
       aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 2")+
  labs( y = "Gene copies / mL Slurry (log10 scale)") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)


p3 <- ggplot(qPCR[qPCR$Treatment == "3",], 
             aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 3")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)


p4 <- ggplot(qPCR[qPCR$Treatment == "4",], 
             aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 4")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)



p5 <- ggplot(qPCR[qPCR$Treatment == "5",], 
             aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 5")+
  labs( y = "Gene copies / mL Slurry (log10 scale)") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)


p6 <- ggplot(qPCR[qPCR$Treatment == "6",], 
             aes(x = Day, y = value, color = Taxa))+ 
  geom_bar(aes(color = Taxa, fill = Taxa),
           stat = "identity", position = position_dodge(0.9),  width = 0.85)+
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=Taxa), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4)+
  ggtitle("Sample 6")+
  labs(color="Taxa") +
  scale_color_manual(values=cols)+
  theme(strip.background = element_rect(fill="white"),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
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
  facet_grid(~Sample)+
  scale_fill_manual("Taxa", values=cols) +
  scale_y_continuous(expand=c(0,0),limits = c(1, 200000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07), labels = scales::scientific)


p1/ p2 / p3 / p4 / p5 /p6 + plot_layout(guides= "collect" )

title=textGrob("Rt-qPCR results from the 4th generation", gp=gpar(fontface="bold", fontsize=25))
y=textGrob("Gene copies / mL Slurry (log10 scale)",rot = 90, gp=gpar(fontface="bold"))
x=textGrob("Incubation day", gp=gpar(fontface="bold"), )
grid.arrange(patchworkGrob(p1/ p2 / p3 / p4 / p5 /p6 + plot_layout(guides= "collect")), left = y, bottom= x, top= title, padding = unit(1, "line"))

             