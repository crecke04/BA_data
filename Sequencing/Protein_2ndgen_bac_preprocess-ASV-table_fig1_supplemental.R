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
setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Enrichments_I_2020/2nd_generation_Nadja/")

## ASV table raw counts
asv <- read_tsv("bac_ASV_protein.txt")

## metadata, need to have same sample names as in taxmap
mdata <- read_tsv("meta_protein_2ndgen.txt")

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
obj$data$tax_data <- obj$data$tax_data[c("taxon_id", "ASV", sample.names)]

### rename table to asv_table
names(obj$data) <- "asv_table"

obj
# 2128 taxa (including subtaxa)
# 8120 ASVs (=number of rows of tax_data table)
ini_reads <- sum(obj$data$asv_table[, sample.names])
ini_reads
# 453749 reads in original table


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
# Zeroing 1058 of 8120 rows with total counts less than 3

### check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 5291 empty ASVs

# remove empty ASVs
obj <- metacoder::filter_obs(obj, "asv_table",
                             ! no_reads,
                             drop_taxa = T)
obj
# taxa reduced to 1179
# ASVs reduced to 2829
sum(obj$data$asv_table[, sample.names])/ini_reads
# 99.66% of the reads kept

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
save(obj, file = "2ndgen_Protein_taxmap-bac.RData")

load("2ndgen_Protein_taxmap-bac.RData")
# here named 'obj'


# create phyloseq object ####
# maybe need to change names of otu_table, otu_id_col and sample names (= column name of column with sample names)
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "ASV", sample_data = mdata, sample_id_col = "Sample")

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
save(phylo_melt, file = "2ndgen_Protein_bac_rawtable_phylomelt.RData")
rm(obj, objphy, opr, asv)
load("2ndgen_Protein_bac_rawtable_phylomelt.RData")



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

ta1 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 1, Abundance = "Abundance")
str(ta1)
# 6 genera, 7 families, 7 orders, 9 classes, 8 phyla

## taxa above 0.1% in at least one samples
ta01 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.1, Abundance = "Abundance")
str(ta01)
# 5 genera, 7 families, 9 orders, 12 classes, 8 phyla

## taxa above 2% in at least one samples
ta05 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.5, Abundance = "Abundance")
str(ta05)
# 3 genera, 2 families, 6 orders, 9 classes, 6 phyla


# use 1% threshold
rm(ta01, ta05)


# order additional metadata columns for later plotting ####
ta1_original <- ta1
ta1$ASV_table_taxa_abv_1_perc$Incubation_time <- factor(ta1$ASV_table_taxa_abv_1_perc$Incubation_time, levels=c("98","157"), ordered=T)

ta1$ASV_table_taxa_abv_1_perc$Sample <- factor(ta1$ASV_table_taxa_abv_1_perc$Sample, levels=c("D98_Control_1", "D98_Control_2",
                                                                                              "D98_Protein_1", "D98_Protein_2",
                                                                                              "D157_Control_1", "D157_Control_2",
                                                                                              "D157_Protein_1", "D157_Protein_2"),
                                               labels = c("98","*98","98","*98",
                                                          "157","*157","157","*157"), ordered = T)

ta1$ASV_table_taxa_abv_1_perc$Treatment <- factor(ta1$ASV_table_taxa_abv_1_perc$Treatment, levels=c("Control",
                                                                                                    "Protein"), 
                                                  labels = c("Control + AB", "Protein + AB"),ordered=T)


save(ta1, file = "2ndgen_protein_bac_long_rel_abd_table_bel1perc.RData")

load("2ndgen_protein_bac_long_rel_abd_table_bel1perc.RData")

# plots ####
## plot on class level, not sorted ####

### which classes will be plotted
unique(ta1$ASV_table_taxa_abv_1_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta1_class <- select(ta1$ASV_table_taxa_abv_1_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, Incubation_time, Treatment) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("Aminicenantia", "other_p_Acidobacteriota_<1%","JS1","Campylobacteria","Desulfobacteria","Desulfovibrionia","other_p_Desulfobacterota_<1%",
                                              "47209","other_p_Elusimicrobiota_<1%" ,"Bacilli","other_p_Firmicutes_<1%","Latescibacteria","other_p_Latescibacterota_<1%","Gammaproteobacteria",
                                              "other_p_Proteobacteria_<1%","other_Bacteria_<1%"), ordered = T), 
         Phylum_mod=factor(Phylum_mod)) %>% 
  group_by(Class_mod)

# levels=c("Aenigmarchaeia", "other_p_Aenigmarchaeota_<0.5%","Heimdallarchaeia","Lokiarchaeia","Odinarchaeia","other_p_Asgardarchaeota_<0.5%",
#          "Bathyarchaeia","other_p_Crenarchaeota_<0.5%","ANME-1","Methanosarcinia","other_p_Halobacterota_<0.5%","Micrarchaeia","Thermoplasmata","other_p_Thermoplasmatota_<0.5%",
#          "other_Archaea_<0.5%"), ordered = T)

### sum ASVs of same taxonomy
ta1_class_s <- ta1_class %>% 
  group_by(Sample, Class_mod, Phylum_mod, Incubation_time, Treatment) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup()

### get automatic colors
# function to define multiple colour pallete based on subgroups from here
# https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar?rq=1
 ColourPalleteMulti <- function(df, group, subgroup){
   
   # Find how many colour categories to create and the number of colours in each
   categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
   category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
   category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom

   # Build Colour pallette
   colours <- unlist(lapply(1:nrow(categories),
                            function(i){
                              colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
   return(colours)
 }

## create colour gradient with color pallet function
 colours <- ColourPalleteMulti(ta1_class_s, "Phylum_mod", "Class_mod")

# cols <- c("#a6cee3","#e7e1ef","#1f78b4","#e7e1ef","#b2df8a","#e7e1ef","#33a02c","#fb9a99","#e7e1ef","#e31a1c","#e7e1ef","#fdbf6f","#ff7f00","#ffff99","#e7e1ef","#cab2d6","#6a3d9a","#e7e1ef","#e7e1ef","#878787","#e0e0e0")
# 
 cols2 <- c("#BA0098","#e0e0e0","#4393c3","#053061","#abd9e9","#e0e0e0","#a6dba0","#e0e0e0","#6a3d9a","#cab2d6","#e0e0e0","#1b7837","#FFFB5E","#e0e0e0","#e0e0e0")

cols2 <- c("#dfc27d","#e0e0e0","#016c59","#cab2d6","#abd9e9","#053061","#e0e0e0","#BA0098","#e0e0e0","#FFFB5E","#e0e0e0","#a6dba0","#e0e0e0","#4393c3","#e0e0e0","#969696")


### plot 
ggplot(ta1_class_s, aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.2) + 
  scale_fill_manual("Taxa", values = cols2) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0.5),
        axis.text.x=element_text(size=7, angle=90, hjust = 1, vjust = 0.5),
        axis.text.y=element_text(size=7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  #ggtitle("16S rRNA gene abundance archaea (class level)") + 
  labs(x="Day", y = "Relative abundance")+
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(~Treatment  ,scales = "free")

ggsave("../../Manuscript_Thermo_ENA/Plots/Relative_abundance_bac_class_ctrl_protein.png", width = 12, height = 8, units = "cm")

