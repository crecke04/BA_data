## script for subsetting a taxmap object to a final long format table with taxa above a certain threshold

# input:
# taxmap object (created in first steps here)
# undesired taxa sorted out
# asv counts, no relative abundance calculated yet
# mdata associated with asv table

# output:
# 1. list of taxa above certain threshold (can be set) in at least one sample for each rank (genus, family, etc.) 
# 2. long format table with samples and individual ASVs with additional columns for each rank renamed taxa name if ASV was below threshold
# 3. long format table with samples and taxa summed up per taxa from table 2


# packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(phyloseq)
library(hues)
library(cowplot)
library(devEMF)


# import data ####
setwd("H:/Desktop/PM4_qPCR/Sequencing")

## ASV table raw counts
asv <- read_tsv("asv_table_tax_Lib80_bac.txt")

## metadata, need to have same sample names as in taxmap
mdata <- read_tsv("meta_arc_bac.txt")

# create taxmap object and do some quality control ####

## get sample names
sample.names <- mdata$Sample

## make taxmap object ####
obj <- parse_tax_data(
  tax_data = asv,            # which table to use
  class_cols = c("Kingdom", "Phylum", "Class", "Order", "Family","Genus"),   # which columns have the taxonomy
  named_by_rank = T          # keep the column names as rank names in command 'taxon_ranks'
)
obj
# the object has some additional columns now with the ASV ID and the taxonomy

### remove taxonomy columns from ASV count table, only keep ASV_ID and taxon_ID column
obj$data$tax_data <- obj$data$tax_data[c("taxon_id", "ASV_ID", sample.names)]

### rename table to asv_table
names(obj$data) <- "asv_table"

obj
# 385 taxa (including subtaxa)
# 202 ASVs (=number of rows of tax_data table)
ini_reads <- sum(obj$data$asv_table[, sample.names])
ini_reads
# 1630323 reads in original table


## modify taxmap object ####
### remove any ASVs which are not Bacteria
obj <- metacoder::filter_taxa(obj,
                              taxon_names == "Bacteria",
                              subtaxa = T)
sum(obj$data$asv_table[, sample.names])/ini_reads
# there were no none Bacteria reads, nothing was removed

### remove doubletons and singletons
obj$data$asv_table <- metacoder::zero_low_counts(obj, "asv_table",
                                                 min_count = 3,
                                                 use_total = T,
                                                 other_cols = T)
# Zeroing 20 of 202 rows with total counts less than 3

### check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 20 empty ASVs

# remove empty ASVs
obj <- metacoder::filter_obs(obj, "asv_table",
                             ! no_reads,
                             drop_taxa = T)
obj
# taxa reduced to 349
# ASVs reduced to 182
sum(obj$data$asv_table[, sample.names])/ini_reads
# 99.99% of the reads kept

# check taxonomy table if everything looks ok
print(taxonomy_table(obj), n = 300)

## calculate further tables ####
## Calculate relative abundance
# relative abundance per ASV
obj$data$rel_abd <- calc_obs_props(obj, "asv_table", other_cols = T)
# relative abundance per taxon
obj$data$rel_tax_abd <- calc_taxon_abund(obj, "rel_abd")
print(obj)

## save or load created taxmap object ####
save(obj, file = "taxmap-bac.RData")

load("taxmap-bac.RData")
# here named 'obj'


# create phyloseq object ####
# maybe need to change names of otu_table, otu_id_col and sample names (= column name of column with sample names)
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "ASV_ID", sample_data = mdata, sample_id_col = "Sample")

ntaxa(objphy)                  # how many taxa (including subtaxa)
nsamples(objphy)               # how many samples
sample_names(objphy)           # check sample names
sample_variables(objphy)       # sample variables from mdata
otu_table(objphy)              # ASV count table
tax_table(objphy)              # taxonomy for ASVs
taxa_names(objphy)             # ASV IDs

rank_names(objphy)             # check rank names
## if not already named like this, change rank names
colnames(tax_table(objphy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")


## calculate relative abundance
opr <- transform_sample_counts(objphy, function(x) x/sum(x))

## melt phyloseq object into dataframe
phylo_melt <- psmelt(opr)

### save object phylo_melt used to perform all following operations and remove not needed objects
save(phylo_melt, file = "bac_rawtable_phylomelt.RData")
rm(obj, objphy, opr, asv)

load("bac_rawtable_phylomelt.RData")



# rename taxa depending on max abundance ####
## use function 'sort_abundant_taxa'
source("filter_taxa_above_threshold.R")

## taxa above 1% in at least one samples
# enter threshold in 'abundance_threshold' in percent
# this function creates a list as output
# for each rank (Phylum, Class, etc.) it contains a vector with all taxa which are in at least one sample
#above the set threshold
# it contains an abundance table similar to the input, with additional modified taxon names
#in which, for each rank, taxa which were below the threshold were renamed as e.g. 'other_o_Aminicenantales_<1%'
# but be careful, it does this via making it NA, so any before unclassified ASV on that rank, will be renamed
#by its supertaxon! but you can always recheck in the not modified taxon name

# ta1 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 1, Abundance = "Abundance")
# str(ta1)
# 2 genera, 2 families, 2 orders, 2 classes, 2 phyla

## taxa above 0.1% in at least one samples
# ta01 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.1, Abundance = "Abundance")
# str(ta01)
# 12 genera, 14 families, 12 orders, 7 classes, 7 phyla


################## If you want to change something you need everything starting here

load("bac_long_rel_abd_table_bel05perc.RData")

## taxa above 0.5% in at least one samples
ta05 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.5, Abundance = "Abundance")
str(ta05)
# 3 genera, 3 families, 3 orders, 3 classes, 3 phyla


# use 0.5% threshold
# rm(ta01, ta1)


# order additional metadata columns for later plotting ####
ta05_original <- ta05
ta05$ASV_table_taxa_abv_0.5_perc$nucleic_acids <- factor(ta05$ASV_table_taxa_abv_0.5_perc$nucleic_acids, levels=c("DNA","RNA"), ordered=T)

ta05$ASV_table_taxa_abv_0.5_perc$Sample <- factor(ta05$ASV_table_taxa_abv_0.5_perc$Sample, levels=c("bac.2C","bac.2C.Gen4","bac.B1", "bac.B2", "bac.B3",
                                                                                              "bac.BHI01A", "bac.BHI01B", "bac.BHI01C",
                                                                                              "bac.BHI05A", "bac.BHI05B", "bac.BHI05C",
                                                                                              "bac.LB01A", "bac.LB01B", "bac.LB01C",
                                                                                              "bac.LB05A", "bac.LB05B", "bac.LB05C",
                                                                                              "bac.BHI05LB05","RNA.BHI05A", "RNA.BHI05B", "RNA.BHI05C",
                                                                                              "RNA.LB05A", "RNA.LB05B", "RNA.LB05C",
                                                                                              "RNA.BHI05LB05"),
                                               labels = c("2C","2C Gen4","Control A","Control B","Control C","BHI 0.1% A","BHI 0.1% B","BHI 0.1% C","BHI 0.5% A","BHI 0.5% B","BHI 0.5% C",
                                                          "LB 0.1% A","LB 0.1% B","LB 0.1% C", "LB 0.5% A","LB 0.5% B","LB 0.5% C", "BHI 0.5% LB 0.5%",
                                                          "BHI 0.5% A","BHI 0.5% B","BHI 0.5% C","LB 0.5% A","LB 0.5% B","LB 0.5% C", "BHI 0.5% LB 0.5%"), ordered = T)

ta05$ASV_table_taxa_abv_0.5_perc$treatment <- factor(ta05$ASV_table_taxa_abv_0.5_perc$treatment, levels=c("WM",
                                                                                                    "LB", "BHI", "BHI_LB"), 
                                                  labels = c("WM", "WM + LB", "WM + BHI", "WM + LB + BHI"),ordered=T)


ta05$ASV_table_taxa_abv_0.5_perc$concentration <- factor(ta05$ASV_table_taxa_abv_0.5_perc$concentration, levels=c("1", "0.1", "0.5"), ordered = T)

save(ta05, file = "bac_long_rel_abd_table_bel05perc.RData")

#######################

load("bac_long_rel_abd_table_bel05perc.RData")

# plots ####
## plot on class level, not sorted ####

### which classes will be plotted
unique(ta05$ASV_table_taxa_abv_0.5_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta05_class <- select(ta05$ASV_table_taxa_abv_0.5_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, nucleic_acids, treatment, concentration) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("Desulfovibrionia", "other_p_Desulfobacterota_<0.5%","Bacilli","other_p_Firmicutes_<0.5%",
                                              "Gammaproteobacteria","other_p_Proteobacteria_<0.5%","other_Bacteria_<0.5%"), ordered = T), 
         Phylum_mod=factor(Phylum_mod)) %>% 
  group_by(Class_mod)


### sum ASVs of same taxonomy
ta05_class_s <- ta05_class %>% 
  group_by(Sample, Class_mod, Phylum_mod, nucleic_acids, treatment, concentration) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup()

### get automatic colors
# function to define multiple colour pallete based on subgroups from here
# https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar?rq=1
#  ColourPalleteMulti <- function(df, group, subgroup){
#    
#    # Find how many colour categories to create and the number of colours in each
#    categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
#    category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
#    category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom
# 
#    # Build Colour pallette
#    colours <- unlist(lapply(1:nrow(categories),
#                             function(i){
#                               colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
#    return(colours)
#  }
# 
# ## create colour gradient with color pallet function
#  colours <- ColourPalleteMulti(ta05_class_s, "Phylum_mod", "Class_mod")

cols2 <- c("#457B9D","#e0e0e0","#E5383B","#e0e0e0","#FDBF6F","#e0e0e0","#969696")


### plot 
p1 <- ggplot(ta05_class_s[ta05_class_s$nucleic_acids == "DNA" & ta05_class_s$treatment %in% c("WM","WM + LB", "WM + BHI",  "WM + LB + BHI") ,], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.2) + 
  scale_fill_manual("Taxa", values = cols2) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0.5),
        axis.text.x=element_text(size=7, angle=90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  #ggtitle("16S rRNA gene abundance archaea (class level)") + 
 labs(x="", y = "DNA")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid( ~treatment  ,scales = "free")

p2 <- ggplot(ta05_class_s[ta05_class_s$nucleic_acids == "RNA" & ta05_class_s$treatment %in% c("WM","WM + LB", "WM + BHI", "WM + LB + BHI") ,], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.2) + 
  scale_fill_manual("Taxa", values = cols2) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0.5),
        axis.text.x=element_text(size=7, angle=90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Treatment", y = "RNA")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid( ~treatment  ,scales = "free")

library(patchwork)
library(grid)
library(gridExtra)

p1/ (p2+plot_spacer()) +plot_layout(guides="collect")&
  theme(legend.position='bottom')

title=textGrob("16S rRNA gene abundance bacteria (class level)", gp=gpar(fontface="bold", fontsize=20))
y=textGrob("Relative abundance",rot = 90, gp=gpar(fontface="bold"))
x=textGrob("", gp=gpar(fontface="bold"), )
grid.arrange(patchworkGrob(p1/ (p2+plot_spacer()) +plot_layout(guides="collect")), left = y, bottom= x, top= title, padding = unit(1, "line"), widths=1)


ggsave("bac_relabundance_05.png", width = 12, height = 8, units = "cm")

