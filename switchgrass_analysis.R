library(phyloseq)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)
library(grid)
library(decontam)
library(dplyr)
library(vegan)

otu = read.csv("otu_table.csv", sep=",", row.names=1)
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("metadata.csv", sep=",", row.names=1)##come back to these steps after adding negative controls to otu table, remove contaminants, then go through all steps again from start
##or do decontam on original otutable (exp plus neg controls) and merge with the current otu table for active taxa. if overlap is found then remove those taxa/contaminants
#metadata = read.csv("sample_metadata_dna_cdna.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

sample_names(phyloseq_merged)

##check column names of taxonomy table
colnames(tax_table(phyloseq_merged))

##remove mito, chloroplast and archaea, eukarya
phyloseq_merged_clean <- phyloseq_merged %>%
  subset_taxa(
    Domain == "Bacteria" 
    #      Family   != " Chloroplast" &
    #      Order   != "Mitochondria" &
    #      Class  != "Chloroplast" &
    #      # OTU     != "0a3cf58d4ca062c13d42c9db4ebcbc53" &
    #      # OTU     != "22f1fa0bdcc19746dee080bcc12a1840"
    #      Family  != " Mitochondria"
  )
phyloseq_merged_clean

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

##identify contaminants by prevalence. in this method, the distribution 
#of the frequency of each sequence feature as a function of the 
#prevalence is used to identify contaminants. In this method, 
#the prevalence (presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence in negative controls 
#to identify contaminants.
sample_data(phyloseq_merged_clean)$is.neg <- sample_data(phyloseq_merged_clean)$Sample_or_Control == "Control"
contamdf.prev0.5 <- isContaminant(phyloseq_merged_clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev0.5$contaminant)

#  FALSE   TRUE
#   8660     42   SWITCHGRASS
#   9913      9   BEAN

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
write.csv(df.pa, "contaminant-table-dna-rna.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5-dna-rna.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-prev-0.5-dna-rna.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam

##export otu table out of phyloseq object 

OTU1 = as(otu_table(phyloseq_merged_clean_decontam), "matrix")
#if(taxa_are_rows(phyloseq_merged_clean_decontam)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTUdf, "OTU_clean.csv")

otu = read.csv("OTU_clean.csv", sep=",", row.names=1)#with no contaminants, mito, chloroplast, neg control 
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("metadata.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied.csv")

library(vegan)
rarecurve(t(otu_table(phyloseq_merged)), step=50)

?rarecurve


phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 10000, rngseed = TRUE, trimOTUs=FALSE)

sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like df1-df2, as it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied.csv")
df2 <- read.csv("samples_rarefied.csv")
df= df1 %>% anti_join(df2,by="sampleid")

View(df)

write.csv(df, "samples-lost-rarefaction.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")
#if(taxa_are_rows(phyloseq_merged_clean_decontam)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_clean_noneg_rarefied15k.csv")
# MANUALLY READD DRY SAMPLES FOR SW

otu_final <- read.csv("OTU_clean_noneg_rarefied15k.csv")


to_remove <- read.csv("samples-lost-equiv.csv")
to_remove[,1]


otu_final <- select(otu_final, -c(to_remove[,1]))

write.csv(otu_final, "OTU_clean_noneg_rarefied15k_removedsamples.csv")

##################################
#divide into dna and rna dataset
library(tidyverse)
require(dplyr) 
require(tibble)
dna_rna<- read.csv("OTU_clean_noneg_rarefied15k_removedsamples.csv", row.names=1) ##new rarefied and cleaned OTU table, with 15k reads per sample, and exact match of DNA RNA samples
dna_rna = dna_rna[rowSums(dna_rna[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(dna_rna)
dna_rna ##otus reduced to 7374
dna_rna1 <-as_tibble(dna_rna, rownames="OTUID") 
dna_rna1
has_rownames(dna_rna1)
rownames(dna_rna1)
dna<- dna_rna1 %>% select(OTUID, starts_with("DNA"))
View(dna)
rna <- dna_rna1 %>% select(OTUID, starts_with("rna"))
View(rna)
write.csv(dna, "DNA.csv")
write.csv(rna, "RNA.csv")

##already removed all negative samples from DNA and rna dataset
##check that all rowsums are >1, remove any rowsums=0 , rarefied reads should still be the same.

DNA<- read.csv("DNA.csv", row.names=1)
View(DNA)##7374 otus
dnanozero = DNA[rowSums(DNA[])>0, ,drop=FALSE]
dnanozero
View(dnanozero)##5578 otus
write.csv(dnanozero, "dna_noneg_nozero.csv")

dnanozero["Total"] <- rowSums(dnanozero)
dnanozero
dnanozero1 = dnanozero %>% filter(dnanozero$Total== 0)

View(dnanozero1)

##subset rna dataframe to include only those OTU that appear in at least one sample
rna<- read.csv("RNA.csv", row.names=1)
View(rna)##7374 otus
rnanozero = rna[rowSums(rna[])>0, ,drop=FALSE]
rnanozero
View(rnanozero)##5481 otus
write.csv(rnanozero, "rna_noneg_nozero.csv")

##merge to see overlap in OTUs between dna and cdna datasets

A = read.csv("dna_noneg_nozero.csv")

B = read.csv("rna_noneg_nozero.csv")

C<- inner_join(A, B) ##3685 taxa across 279 columns, join keeps both sets of data from A and B

View(C)


###combine two matrices one for dna one for cdna -- remove headers and rownames , use files DNA_initial_noneg_edited 
##and cDNA_initial_noneg_edited
library(tidyverse)

dna<- read.csv("DNA_mat_noheader.csv", header=FALSE)
rna<- read.csv("RNA_mat_noheader.csv", header=FALSE)
View(dna)
View(rna)

##OTU table where all DNA=0 is replaced with 1

dna.mat<- as.matrix(dna)
dna.mat[dna.mat == 0] <- 1

View(dna.mat)

rna.mat <- as.matrix(rna)

View(rna.mat)
merge.mat<- rna.mat/dna.mat


View(merge.mat)
write.csv(merge.mat, "rna-dna-ratio-method1.csv")

##selecting ratios greater than equal to 1

merge.mat[merge.mat <= 1] <- 0 ##second try modified <1 to <= 1, to only keep values of ratio greater than one
View(merge.mat)
#merge.mat[is.na(merge.mat)] <- 0

merge.mat[!is.finite(merge.mat)] <- 0
View(merge.mat)

write.csv(merge.mat, "rna-dna-allfinite-greaterthanone.csv")

##replacing ratios greater than 1 with abund values in rna table
##change this part : replacing ratios greater than 1 with abund values in **dna table**
##this will help reduce biases pertaining to copy number issue in RNA reads, for example: there may be DNA for which there are over-expressed copies of RNA 
##so to avoid that bias we will filter DNA abundance to active
##next step should also be to make binary (presence absence dataset) to calculate activity dynamics. this is to look at piles of similarly responsive OTUs
##and look at key transitions. for example: are there specific OTUs that are consistently detected across treatments (specifically drought planted ##treatments)
rna_ratio<- read.csv("rna-dna-allfinite-greaterthanone.csv", row.names = 1)
View(rna_ratio)
rna_ratio_mat<- as.matrix(rna_ratio)
View(rna_ratio_mat)



#dna_abun <- read.csv("DNA_mat.csv", row.names = 1) #replace above code with DNA_mat.csv
dna_abun_mat <- dna.mat #as.matrix(dna_abun)
View(dna_abun_mat)

##use this
#dna abundance of active OTUs
mask <- rna_ratio_mat > 0
mask

rna_ratio_mat[mask] <- dna_abun_mat[mask]
rna_ratio_mat
View(rna_ratio_mat)
write.csv(rna_ratio_mat, "finalactiverna_method2_redo.csv")

View(dna.mat)
View(merge.mat)

rna_active<- read.csv("finalactiverna_method2_redo.csv", row.names=1)
activefinalnozero = rna_active[rowSums(rna_active[])>0, ,drop=FALSE]


View(activefinalnozero) 

write.csv(activefinalnozero, "activefinalnozero.csv")


#dna abundance of inactive OTUs
rna_ratio<- read.csv("rna-dna-allfinite-greaterthanone.csv", row.names = 1)
rna_ratio_mat<- as.matrix(rna_ratio)
view(rna_ratio_mat)

rna_ratio_mat[rna_ratio_mat > 1] <- 10000
view(rna_ratio_mat)

mask2 <- rna_ratio_mat < 1
mask2

rna_ratio_mat[mask2] <- dna_abun_mat[mask2]
view(rna_ratio_mat)

rna_ratio_mat[rna_ratio_mat > 9999] <- 0
view(rna_ratio_mat)

write.csv(rna_ratio_mat, "finalinactiverna_method2_redo.csv")

View(dna.mat)
View(merge.mat)

rna_inactive<- read.csv("finalinactiverna_method2_redo.csv", row.names=1)
inactivefinalnozero = rna_inactive[rowSums(rna_active[])>0, ,drop=FALSE]


View(inactivefinalnozero) 

write.csv(inactivefinalnozero, "inactivefinalnozero.csv")


otu = read.csv("activefinalnozero.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("metadata_RNA.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged


phyloseq_phylum <- phyloseq_merged %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at family/phylum/class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.00001) %>%                         # Filter out low abundance taxa
  arrange(Phylum) # Sort data frame alphabetically by phylum/family etc

phyloseq_phylum$Phylum[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

sum(phyloseq_phylum$Abundance)

unique(phyloseq_phylum$Phylum)


library(randomcoloR)
n <- 36
palette <- distinctColorPalette(n)

#make the plot!

plot1 <- ggplot(phyloseq_phylum, aes(x = Treatment_Timepoint, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU") + 
  scale_fill_manual(values = palette)

plot1

ggsave(filename = "SW_Active_2.png", plot = plot1,
       width = 50,
       height = 30, units = c("cm"),
       dpi = 300)


##richness, diversity, etc

library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)

library(tibble)
df <- tibble::rownames_to_column(df, "VALUE")

otu = read.csv("activefinalnozero_flip.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = sample_df
OTU = otu_table(otu, taxa_are_rows = FALSE)
TAX = tax_table(tax)
meta = sample_data(Final_new)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged


data_richness <- estimateR(otu)
view(data_richness)

data_evenness <- diversity(otu) / log(specnumber(otu))
view(data_evenness)

data_shannon <- diversity(otu, index = "shannon")
view(data_shannon)

data_alphadiv <- cbind(Final_new, t(data_richness), data_shannon, data_evenness)     

data_alphadiv_tidy <- 
  data_alphadiv %>%
  mutate(sample_id = rownames(data_alphadiv)) %>%
  gather(key   = alphadiv_index,
         value = obs_values,
         -sample_id, -Treatment, -Timepoint, -Treatment_Timepoint, -Cat, -PercActive, -PercPhantom)

head(data_alphadiv_tidy)







richness <- ggplot(data_alphadiv, aes(x=Treatment_Timepoint, y=S.obs)) +
  geom_boxplot(fill=c("chocolate4","chocolate3","lightsteelblue4","lightsteelblue1","darkorchid2","darkorchid3","darkorchid4","firebrick2","firebrick3","firebrick4","seagreen2","seagreen3","seagreen4","steelblue2","steelblue3","steelblue4")) +
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"))+ 
  labs(title= 'Richness', x= ' ', y= '') +
  geom_point() 
P1



ggsave(filename = "SW_Richness.png", plot = P1,
       width = 70,
       height = 30, units = c("cm"),
       dpi = 300)







P2 <- ggplot(data_alphadiv, aes(x=Treatment_Timepoint, y=S.chao1)) +
  geom_boxplot(fill=c("black","white","black","white","black","white","black","white","black","white","black","white","black","white","black","white")) +
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"))+ 
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point()
P2

P3 <- ggplot(data_alphadiv, aes(x=Treatment_Timepoint, y=data_evenness)) +
  geom_boxplot(fill=c("black","white","black","white","black","white","black","white","black","white","black","white","black","white","black","white")) +
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"))+ 
  labs(title= 'Eveness', x= ' ', y= '', tag = "C") +
  geom_point()
P3

P4 <- ggplot(data_alphadiv, aes(x=Treatment_Timepoint, y=data_shannon)) +
  geom_boxplot(fill=c("black","white","black","white","black","white","black","white","black","white","black","white","black","white","black","white")) +
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"))+ 
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point()
P4

# all plots together using the patchwork package
library(patchwork)
(P1 | P2) / (P3 | P4)


#########PERCENT ACTIVE AND PERCENT PHANTOM


#read in metadata
metadata = read.csv("metadata_RNA.csv", sep=",", row.names=1)

#Make OTU table containing only 16sRNA
oturna = read.csv("RNA.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
OTUsRNA = otu_table(oturna, taxa_are_rows = TRUE)

#Make OTU table containing only 16sDNA
otudna = read.csv("DNA.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
OTUsDNA = otu_table(otudna, taxa_are_rows = TRUE)


#Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
#First, add the RNA and DNA datasets together for each sample.
TotalOTUs <- data.frame(OTUsRNA+OTUsDNA)

#Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
TotalOTUs <- colSums(TotalOTUs !=0)
view(TotalOTUs)

view(OTUsDNA)
#add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
OTUsDNA<-replace(OTUsDNA, OTUsDNA == 0, 1)


#Make table of 16s ratios (RNA/DNA).
OTUs_16sRatio <- data.frame(OTUsRNA/OTUsDNA)
view(OTUs_16sRatio)


#Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
threshold <- 1
OTUs_16sRatio[OTUs_16sRatio<threshold] = 0
OTUs_16sRatio[OTUs_16sRatio>=threshold] = 1

view(OTUs_16sRatio)


#Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
#This gives us the Percent of OTUs that are active.

OTUs_Active <- colSums(OTUs_16sRatio == 1)
view(OTUs_Active)



PercActive <- colSums(OTUs_16sRatio, na.rm=TRUE)/TotalOTUs*100
view(PercActive)


PercActive <- as.data.frame(PercActive)



#Make OTU table containing only 16sRNA
oturna = read.csv("RNA.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
OTUsRNA = otu_table(oturna, taxa_are_rows = TRUE)
view(OTUsRNA)

#Make OTU table containing only 16sDNA
otudna = read.csv("DNA.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
OTUsDNA = otu_table(otudna, taxa_are_rows = TRUE)
view(OTUsDNA)

#Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
#First, add the RNA and DNA datasets together for each sample.
TotalOTUs <- data.frame(OTUsRNA+OTUsDNA)

#Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
TotalOTUs <- colSums(TotalOTUs !=0)
view(TotalOTUs)


#in RNA table, replace RNA>=1 as 1, and keep all zeros as zero.
OTUsRNA[OTUsRNA>=1] = 1
OTUsRNA[OTUsRNA<1] = 0
view(OTUsRNA)

#In DNA table, replace all 0's with 1, and everything else with 0.
#First, make things that are currently 1 into a higher number.
OTUsDNA[OTUsDNA>=1] = 2

#then, make everything that's currently 0 into a 1.
OTUsDNA[OTUsDNA<1] = 1

#then, make everything greater than 1 into 0.
OTUsDNA[OTUsDNA>1] = 0

view(OTUsDNA)


#sum the two tables
Phantom <- data.frame(OTUsRNA+OTUsDNA)

#Count the number of 2's for each sample (RNA=1 coded as 1 and DNA=0 coded as 1)
Phantom <- colSums(Phantom ==2)

view(Phantom)

#calculate percent for each sample
PercPhantom <- Phantom/TotalOTUs*100
PercPhantom <- as.data.frame(PercPhantom)
view(PercPhantom)

#read in metadata
metadata = read.csv("metadata_RNA.csv", sep=",", row.names=1)
Final <- cbind(PercActive, PercPhantom, metadata)

Final_new <- Final

Final_new$Cat <- factor(Final_new$Cat,
                        levels=c("Pre-hormone treatment", "Abscisic acid", "Salicylic acid", "Methanol control", "Water control"))

Final_new$Timepoint <- factor(Final_new$Timepoint,
                              levels=c("Field soil", "Predry", "Postdry", "Water acclimated", "Day 1","Day 7","Day 14"))

Final_new$Treatment <- factor(Final_new$Treatment,
                              levels=c("Field soil", "Predry", "Postdry", "Water acclimated", "Abscisic acid","Salicylic acid","Methanol control","Water control"))

Final_new$Timepoint2 <- factor(Final_new$Timepoint2,
                               levels=c("Pre-hormone treatment", "Day 1","Day 7","Day 14"))


library(ggplot2)

# Basic violin plot
PercActivePlot <- ggplot(metadata, aes(x=Timepoint, y=PercActive,color=Cat), show.legend=FALSE) + 
  geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.07)) +
  scale_color_manual(values = c("goldenrod4", "darkolivegreen3", "darkorchid2", "firebrick3", "deepskyblue3")) +
  facet_wrap(~Cat, scale="free_x", nrow = 5) +
  xlab("") +
  ylab("%") +
  labs(title="(A) Percent Active OTUs") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 10),
        strip.text = element_text(size = 11, face = "bold"))

PercActivePlot

ggsave(filename = "Percent_active_final.png", plot = PercActivePlot,
       width = 10,
       height = 35, units = c("cm"),
       dpi = 300)




# Basic violin plot
PercPhantomPlot <- ggplot(metadata, aes(x=Timepoint, y=PercPhantom,color=Cat), show.legend=FALSE) + 
  geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.07)) +
  scale_color_manual(values = c("goldenrod4", "darkolivegreen3", "darkorchid2", "firebrick3", "deepskyblue3")) +
  facet_wrap(~Cat, scale="free_x", nrow = 5) +
  xlab("") +
  ylab("") +
  labs(title="(B) Percent Phantom OTUs") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 9, colour = "black", vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 10),
        strip.text = element_text(size = 11, face = "bold"))

PercPhantomPlot

ggsave(filename = "Percent_phantom_final.png", plot = PercPhantomPlot,
       width = 10,
       height = 35, units = c("cm"),
       dpi = 300)


library(patchwork)
Perc <- (PercActivePlot | PercPhantomPlot)

ggsave(filename = "Percent_active_phantom_final.png", plot = Perc,
       width = 17,
       height = 35, units = c("cm"),
       dpi = 300)


#PCOA 
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
#meta_ABA = sample_data(Final_new_ABA)
#meta_SA = sample_data(Final_new_SA)
#meta_MC = sample_data(Final_new_MC)
#meta_WC = sample_data(Final_new_WC)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge

phyloseq = phyloseq(OTU, TAX, meta)
#phyloseq_ABA = phyloseq(OTU, TAX, meta_ABA)
#phyloseq_SA = phyloseq(OTU, TAX, meta_SA)
#phyloseq_MC = phyloseq(OTU, TAX, meta_MC)
#phyloseq_WC = phyloseq(OTU, TAX, meta_WC)

phyloseq_pcoa

#PCOA all
phyloseq_pcoa <- ordinate(
  physeq = phyloseq, 
  method = "PCoA", 
  distance = "bray")

# Plot 
pcoa_all<-plot_ordination(
  physeq = phyloseq,
  ordination = phyloseq_pcoa,
  color = "Cat",
  shape = "Timepoint2",
  title = "") +
  scale_color_manual(values = c("darkorange", "darkolivegreen4","darkorchid3","firebrick3","deepskyblue3")) +
  scale_shape_manual(values = c(18, 20, 20, 20)) +
  geom_point(aes(size=Timepoint2), alpha=0.4) +
  scale_size_manual (values= c(4,6,8,10)) +
  geom_point(colour = "grey90")

pcoa_all

ggsave(filename = "pcoa_all.png", plot = pcoa_all,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)


phyloseq_1 <- subset_samples(phyloseq, Cat=="Methanol control"|Cat=="Water control")

phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
sample_df <- data.frame(sample_data(phyloseq_merged))

adonis2(phyloseq_bray ~ Treatment*Timepoint, data = sample_df, permutations=999)

phyloseq <- subset_samples(phyloseq_original, Cat==4|B==14)
phyloseq_bray <- phyloseq::distance(phyloseq, method = "bray")
sample_df <- data.frame(sample_data(phyloseq_original))
adonisvalues <- adonis(phyloseq_bray ~ Treatment_Timepoint, data = sample_df, permutations=999)
adonisvalues[1,4]
adonisvalues[1,5]



library(dplyr)
library(stats)


wilcox <- pairwise.wilcox.test(metadata$Richness, metadata$Treatment_Timepoint,
                               p.adjust.method = "bonferroni")

view(wilcox$p.value)

write.csv(wilcox$p.value,"richnesswilcox.csv")



#DESeq

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)
library(grid)
library(decontam)
library(dplyr)
library(vegan)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

otu = read.csv("activefinalnozero.csv", sep=",", row.names=1)
tax = read.csv("taxonomy.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("metadata_RNA.csv", sep=",", row.names=1)

metadata$Cat <- factor(metadata$Cat,
                       levels=c("Pre-hormone treatment", "Abscisic acid", "Salicylic acid", "Methanol control", "Water control"))

metadata$Timepoint <- factor(metadata$Timepoint,
                             levels=c("Field soil", "Predry", "Postdry", "Water acclimated", "Day 1","Day 7","Day 14"))

metadata$Treatment <- factor(metadata$Treatment,
                             levels=c("Field soil", "Predry", "Postdry", "Water acclimated", "Abscisic acid","Salicylic acid","Methanol control","Water control"))

metadata$Timepoint2 <- factor(metadata$Timepoint2,
                              levels=c("Pre-hormone treatment", "Day 1","Day 7","Day 14"))

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
##rownames(meta) <- meta$sampleid (this is not needed to merge)
##merge
phyloseq_original = phyloseq(OTU, TAX, meta)
phyloseq_original

sample_data(phyloseq_original)


phyloseq_1 <- subset_samples(phyloseq_original, Cat=="Salicylic acid"|Cat=="Abscisic acid")


head(sample_data(phyloseq_1))


diagdds = phyloseq_to_deseq2(phyloseq_1, ~ Cat)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_1)[rownames(sigtab), ], "matrix"))
head(sigtab)

write.csv(sigtab,"SA_vs_WC.csv")

sigtab <- read.csv("ABA_vs_SA.csv", row.names = 1)


library(randomcoloR)
n <- 10
palette <- distinctColorPalette(n)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))

h1 <- ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Class)) + 
  geom_point(size=6,alpha=0.5,shape=18,position = position_jitter(seed = 1, width = 0.2)) + 
  theme(axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust=0.5))

h1

