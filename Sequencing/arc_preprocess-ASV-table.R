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
asv <- read_tsv("asv_table_tax_Lib80_arc.txt")

## metadata, need to have same sample names as in taxmap
sessmdata <- read_tsv("meta_arc.txt")

# create taxmap object and do some quality control ####

## get sample names
sample.names <- sessmdata$Sample

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
# 51 taxa (including subtaxa)
# 38 ASVs (=number of rows of tax_data table)
ini_reads <- sum(obj$data$asv_table[, sample.names])
ini_reads
# 254977 reads in original table


## modify taxmap object ####
### remove any ASVs which are not Bacteria
obj <- metacoder::filter_taxa(obj,
                              taxon_names == "Archaea",
                              subtaxa = T)
sum(obj$data$asv_table[, sample.names])/ini_reads
# there were no none Bacteria reads, nothing was removed

### remove doubletons and singletons
obj$data$asv_table <- metacoder::zero_low_counts(obj, "asv_table",
                                                 min_count = 3,
                                                 use_total = T,
                                                 other_cols = T)
# Zeroing 2 of 38 rows with total counts less than 3

### check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 2 empty ASVs

# remove empty ASVs
obj <- metacoder::filter_obs(obj, "asv_table",
                             ! no_reads,
                             drop_taxa = T)
obj
# taxa reduced to 46
# ASVs reduced to 36
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
save(obj, file = "taxmap-arc.RData")

load("taxmap-arc.RData")
# here named 'obj'


# create phyloseq object ####
# maybe need to change names of otu_table, otu_id_col and sample names (= column name of column with sample names)
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "ASV_ID", sample_data = sessmdata, sample_id_col = "Sample")

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
save(phylo_melt, file = "arc_rawtable_phylomelt.RData")
rm(obj, objphy, opr, asv)

load("arc_rawtable_phylomelt.RData")



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

## taxa above 1% in at least one samples
 # ta1 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 1, Abundance = "Abundance")
 # str(ta1)
# 0 genera, 0 families, 1 orders, 2 classes, 3 phyla

## taxa above 0.1% in at least one samples
 # ta01 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.1, Abundance = "Abundance")
 # str(ta01)
# 2 genera, 2 families, 3 orders, 3 classes, 4 phyla

################################

load("arc_long_rel_abd_table_bel05perc.RData")
 
## taxa above 0.5% in at least one samples
ta05 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.5, Abundance = "Abundance")
# str(ta05)
# 1 genera, 1 families, 2 orders, 3 classes, 4 phyla


# use 0.5% threshold
# rm(ta01, ta1)


# order additional metadata columns for later plotting ####
ta05_original <- ta05
ta05$ASV_table_taxa_abv_0.5_perc$nucleic_acids <- factor(ta05$ASV_table_taxa_abv_0.5_perc$nucleic_acids, levels=c("DNA","RNA"), ordered=T)

ta05$ASV_table_taxa_abv_0.5_perc$Sample <- factor(ta05$ASV_table_taxa_abv_0.5_perc$Sample, levels=c("arc.2C","arc.2C.Gen4","arc.B1", "arc.B2", "arc.B3",
                                                                                                    "arc.BHI01A", "arc.BHI01B", "arc.BHI01C",
                                                                                                    "arc.BHI05A", "arc.BHI05B", "arc.BHI05C",
                                                                                                    "arc.LB01A", "arc.LB01B", "arc.LB01C",
                                                                                                    "arc.LB05A", "arc.LB05B", "arc.LB05C", "arc.BHI05LB05"),
                                                  labels = c("2C","2C Gen4","Control A","Control B","Control C","BHI 0.1% A","BHI 0.1% B","BHI 0.1% C","BHI 0.5% A","BHI 0.5% B","BHI 0.5% C",
                                                             "LB 0.1% A","LB 0.1% B","LB 0.1% C", "LB 0.5% A","LB 0.5% B","LB 0.5% C", "BHI 0.5% LB 0.5%"), ordered = T)

ta05$ASV_table_taxa_abv_0.5_perc$treatment <- factor(ta05$ASV_table_taxa_abv_0.5_perc$treatment, levels=c("WM",
                                                                                                          "LB", "BHI", "BHI_LB"), 
                                                     labels = c("WM", "WM + LB", "WM + BHI", "WM + LB + BHI"),ordered=T)


ta05$ASV_table_taxa_abv_0.5_perc$concentration <- factor(ta05$ASV_table_taxa_abv_0.5_perc$concentration, levels=c("1", "0.1", "0.5"), ordered = T)

save(ta05, file = "arc_long_rel_abd_table_bel05perc.RData")

#######################

load("arc_long_rel_abd_table_bel05perc.RData")

# plots ####
## plot on class level, not sorted ####

### which classes will be plotted
unique(ta05$ASV_table_taxa_abv_0.5_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta05_class <- select(ta05$ASV_table_taxa_abv_0.5_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, nucleic_acids, treatment, concentration) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("other_p_Asgardarchaeota_<0.5%", "Bathyarchaeia","Methanosarcinia","other_p_Halobacterota_<0.5%",
                                              "Thermoplasmata","other_Archaea_<0.5%"), 
                          ordered = T), 
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

cols2 <- c("#457B9D","#E5383B","#660708","#e0e0e0","#FDBF6F","#969696")

### plot 
ggplot(ta05_class_s, aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.2) + 
  scale_fill_manual("Taxa", values = cols2) + theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.5),
        axis.text.x=element_text(size=7, angle=90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("16S rRNA gene abundance archaea (class level)") + 
  labs(x="Treatment", y = "Relative abundance (DNA)")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~treatment  ,scales = "free")

ggsave("arc_relabundance_05.png", width = 20, height = 20, units = "cm")

