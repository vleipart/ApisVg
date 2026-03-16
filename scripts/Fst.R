#Fst Calculations
library(pegas)
library(adegenet)
library(hierfstat)
library(ape)
library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyverse)
library(here)


#Importing sequences including and excluding putative recombinant sequences
alignment_norecomb <- read.dna(here("..", "results", "vgseq_combined_msa_exons_norecomb.fasta"), format = "fasta")
alignment_all <- read.dna(here("..", "results", "vgseq_combined_msa_exons.fasta"), format = "fasta")

# Extract sequence names
sequence_names_norecomb <- rownames(alignment_norecomb)
sequence_names_all <- rownames(alignment_all)

# Extract population codes from sequence names
populations_norecomb <- sub(".*_(\\w{2})\\d+_.*", "\\1", sequence_names_norecomb)
populations_all <- sub(".*_(\\w{2})\\d+_.*", "\\1", sequence_names_all)

# Convert to factor
populations_norecomb <- factor(populations_norecomb)
populations_all <- factor(populations_all)

genind_obj_norecomb <- DNAbin2genind(alignment_norecomb, pop = populations_norecomb)
genind_obj_all <- DNAbin2genind(alignment_all, pop = populations_all)

hierfstat_obj_norecomb <- genind2hierfstat(genind_obj_norecomb)
hierfstat_obj_all <- genind2hierfstat(genind_obj_all)
fst_results_norecomb <- wc(hierfstat_obj_norecomb)
fst_results_all <- wc(hierfstat_obj_all)
fst_results_norecomb #Global Fst 
fst_results_all #Global Fst 

pairwise_fst_norecomb <- pairwise.WCfst(hierfstat_obj_norecomb)
pairwise_fst_all <- pairwise.WCfst(hierfstat_obj_all)

# View the pairwise Fst matrix
# Convert the matrix to a long format data frame
fst_melt_norecomb <- melt(as.matrix(pairwise_fst_norecomb), varnames = c('Population1', 'Population2'), value.name = "Fst")
fst_melt_all <- melt(as.matrix(pairwise_fst_all), varnames = c('Population1', 'Population2'), value.name = "Fst")

pop_order <- c("FN", "RN", "SW", "DK", "TX", "IR", "CS", "FR", 
               "PL", "PO", "SL", "IT", "MK", "TR", "MT", 
               "CA", "AZ", "NC", "MD", "IL", "MN")

fst_melt_norecomb$Population1 <- factor(fst_melt_norecomb$Population1, levels = pop_order)
fst_melt_norecomb$Population2 <- factor(fst_melt_norecomb$Population2, levels = pop_order)

fst_melt_all$Population1 <- factor(fst_melt_all$Population1, levels = pop_order)
fst_melt_all$Population2 <- factor(fst_melt_all$Population2, levels = pop_order)

fst_melt_norecomb <- fst_melt_norecomb[order(fst_melt_norecomb$Population1, fst_melt_norecomb$Population2), ]
fst_melt_all <- fst_melt_all[order(fst_melt_all$Population1, fst_melt_all$Population2), ]

# Define the population groups
brownbee <- c("CS", "DK", "FN", "FR", "IR", "PL", "RN", "SW", "TX")
southEurope <- c("IT", "MK", "MT", "PO", "SL", "TR")
usa <- c("AZ", "CA", "IL", "MD", "MN", "NC")

individual_pops_norecomb <- pop(genind_obj_norecomb)
individual_pops_all <- pop(genind_obj_all)

# Assign group labels based on population
group_labels_norecomb <- ifelse(individual_pops_norecomb %in% brownbee, "European Dark\n honey bee samples",
                       ifelse(individual_pops_norecomb %in% southEurope, "Other European\n samples",
                              ifelse(individual_pops_norecomb %in% usa, "USA samples", NA)))

group_labels_all <- ifelse(individual_pops_all %in% brownbee, "European Dark\n honey bee samples",
                       ifelse(individual_pops_all %in% southEurope, "Other European\n samples",
                              ifelse(individual_pops_all %in% usa, "USA samples", NA)))

# Convert to factor
group_labels_norecomb <- factor(group_labels_norecomb)
group_labels_all <- factor(group_labels_all)

# Double-check length and names match genind object
length(group_labels_norecomb) == nInd(genind_obj_norecomb)  # Should be TRUE
length(group_labels_all) == nInd(genind_obj_all)  # Should be TRUE

# Subset genind object by group and calculate global Fst
calc_group_fst_norecomb <- function(group) {
  keep_inds_norecomb <- which(group_labels_norecomb == group)
  group_genind_norecomb <- genind_obj_norecomb[keep_inds_norecomb, ]
  group_hierfstat_norecomb <- genind2hierfstat(group_genind_norecomb)
  wc(group_hierfstat_norecomb)$FST
}

calc_group_fst_all <- function(group) {
  keep_inds_all <- which(group_labels_all == group)
  group_genind_all <- genind_obj_all[keep_inds_all, ]
  group_hierfstat_all <- genind2hierfstat(group_genind_all)
  wc(group_hierfstat_all)$FST
}

# Calculate Fst for each group
fst_brownbee_norecomb <- calc_group_fst_norecomb("European Dark\n honey bee samples")
fst_southEurope_norecomb <- calc_group_fst_norecomb("Other European\n samples")
fst_usa_norecomb <- calc_group_fst_norecomb("USA samples")

fst_brownbee_all <- calc_group_fst_all("European Dark\n honey bee samples")
fst_southEurope_all <- calc_group_fst_all("Other European\n samples")
fst_usa_all <- calc_group_fst_all("USA samples")

fst_brownbee_norecomb  # Global Fst for Brownbee group
fst_southEurope_norecomb  # Global Fst for SouthEurope group
fst_usa_norecomb  # Global Fst for USA group

fst_brownbee_all  # Global Fst for Brownbee group
fst_southEurope_all  # Global Fst for SouthEurope group
fst_usa_all  # Global Fst for USA group

# Reassign the populations in genind_obj to the new group labels
genind_grouped_norecomb <- genind_obj_norecomb
genind_grouped_all <- genind_obj_all
pop(genind_grouped_norecomb) <- group_labels_norecomb
pop(genind_grouped_all) <- group_labels_all

# Convert the grouped genind object to a hierfstat object
hierfstat_grouped_norecomb <- genind2hierfstat(genind_grouped_norecomb)
hierfstat_grouped_all <- genind2hierfstat(genind_grouped_all)

# Calculate pairwise Fst between the three groups
pairwise_fst_groups_norecomb <- pairwise.WCfst(hierfstat_grouped_norecomb)
pairwise_fst_groups_all <- pairwise.WCfst(hierfstat_grouped_all)

# View the pairwise Fst matrix
pairwise_fst_groups_norecomb
pairwise_fst_groups_all

fst_melt_groups_norecomb <- melt(as.matrix(pairwise_fst_groups_norecomb), varnames = c('Population1', 'Population2'), value.name = "Fst")
fst_melt_groups_all <- melt(as.matrix(pairwise_fst_groups_all), varnames = c('Population1', 'Population2'), value.name = "Fst")


#Creating dataframe with Global Fst
fst <- c(fst_results_norecomb$FST, fst_brownbee_norecomb, fst_southEurope_norecomb, fst_usa_norecomb, 
         fst_results_all$FST, fst_brownbee_all,fst_southEurope_all, fst_usa_all)
Recomb <- c("Excluding recombinant alleles", "Excluding recombinant alleles", 
            "Excluding recombinant alleles", "Excluding recombinant alleles", 
            "Including recombinant alleles", "Including recombinant alleles", 
            "Including recombinant alleles", "Including recombinant alleles")
pop <- c("all", "European Dark honey bee samples","Other European samples", "USA samples", 
         "all", "European Dark honey bee samples","Other European samples", "USA samples")
size <- c(nrow(genind_grouped_norecomb$tab),as.integer(table(genind_grouped_norecomb@pop)), nrow(genind_grouped_all$tab), as.integer(table(genind_grouped_all@pop)))
Fstdf <- data.frame(fst,Recomb,pop,size)

#Plotting Global Fst for all groups (Supplement Figure S4)
Fst_plot <- ggplot(Fstdf, aes(x=fct_inorder(pop), y=fst, fill= Recomb, group = Recomb)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title = "", x="", y=expression("Global F"[ST])) +
  geom_text(label=size, position = position_dodge(width=0.9), vjust = 2, color = "black", size=5)  +
  scale_fill_manual(values = c("Excluding recombinant alleles" = "#9fc8c8", 
                               "Including recombinant alleles" = "#f0b077"), name = "") +
  theme_classic() +
  theme(legend.title = element_text(size = 10),legend.position = "top")

ggsave(filename = here("..", "figures", "FigureS4.png"), plot = Fst_plot, width = 8, height = 6, dpi = 300)


#Plotting pairwise Fst (Supplement Figure S7)
european_dark <- c("FN", "RN", "SW", "DK", "TX", "IR", "CS", "FR", "PL")
other_euro <- c("PO", "SL", "IT", "MK", "TR", "MT")
usa_mixed <- c("CA", "AZ", "NC", "MD", "IL", "MN")
pop_levels <- levels(fst_melt_all$Population1)
get_bounds <- function(group, levels) {
  pos <- match(group, levels)
  c(xmin = min(pos) - 0.5,
    xmax = max(pos) + 0.5,
    ymin = min(pos) - 0.5,
    ymax = max(pos) + 0.5)
}

bounds_dark <- get_bounds(european_dark, pop_levels)
bounds_euro <- get_bounds(other_euro, pop_levels)
bounds_usa <- get_bounds(usa_mixed, pop_levels)

hm_recomb <- ggplot(fst_melt_norecomb, aes(Population1, Population2, fill = Fst)) +
  geom_tile() +
  annotate("rect", xmin = bounds_dark["xmin"], xmax = bounds_dark["xmax"],
           ymin = bounds_dark["ymin"], ymax = bounds_dark["ymax"],
           color = "#D91A1A", fill = NA, size = 1.2) +
  annotate("rect", xmin = bounds_euro["xmin"], xmax = bounds_euro["xmax"],
           ymin = bounds_euro["ymin"], ymax = bounds_euro["ymax"],
           color = "#377FB8", fill = NA, size = 1.2) +
  annotate("rect", xmin = bounds_usa["xmin"], xmax = bounds_usa["xmax"],
           ymin = bounds_usa["ymin"], ymax = bounds_usa["ymax"],
           color = "#5DB45A", fill = NA, size = 1.2) +
  scale_fill_gradient(low = "#BF7F30", high = "yellow", na.value = "white", limits = c(0, 0.15)) +
  theme_minimal() +
  labs(title = "D) Apiary (excluding\nrecombinant haplotypes)", fill = "Fst", x="", y="") 

hm_all <- ggplot(fst_melt_all, aes(Population1, Population2, fill = Fst)) +
  geom_tile() +
  annotate("rect", xmin = bounds_dark["xmin"], xmax = bounds_dark["xmax"],
           ymin = bounds_dark["ymin"], ymax = bounds_dark["ymax"],
           color = "#D91A1A", fill = NA, size = 1.2) +
  annotate("rect", xmin = bounds_euro["xmin"], xmax = bounds_euro["xmax"],
           ymin = bounds_euro["ymin"], ymax = bounds_euro["ymax"],
           color = "#377FB8", fill = NA, size = 1.2) +
  annotate("rect", xmin = bounds_usa["xmin"], xmax = bounds_usa["xmax"],
           ymin = bounds_usa["ymin"], ymax = bounds_usa["ymax"],
           color = "#5DB45A", fill = NA, size = 1.2) +
  scale_fill_gradient(low = "#BF7F30", high = "yellow", na.value = "white", limits = c(0, 0.15)) +
  theme_minimal() +
  labs(title = "C) Apiary (including\nrecombinant haplotypes)", fill = "Fst", x="", y="") 

hm_groups <- ggplot(fst_melt_groups_norecomb, aes(Population1, Population2, fill = Fst)) +
  geom_tile() +
  scale_fill_gradient(low = "#BF7F30", high = "yellow", na.value = "white", limits = c(0, 0.15)) +
  theme_minimal() +
  labs(title = "B) Sampling categories\n(excluding recombinant haplotypes)", fill = "Fst", x="", y="") + 
  geom_text(aes(label = ifelse(Fst > 0, round(Fst, 3), "")), color = "black")

hm_groups_all <- ggplot(fst_melt_groups_all, aes(Population1, Population2, fill = Fst)) +
  geom_tile() +
  scale_fill_gradient(low = "#BF7F30", high = "yellow", na.value = "white", limits = c(0, 0.15)) +
  theme_minimal() +
  labs(title = "A) Sampling categories\n(including recombinant haplotypes)", fill = "Fst", x="", y="") + 
  geom_text(aes(label = ifelse(Fst > 0, round(Fst, 3), "")), color = "black")


layout <- "
AAABBB
CCCDDD"

FigS7 <- hm_groups_all +hm_groups + hm_all +  hm_recomb + plot_layout(design = layout) 

ggsave(filename = here("..", "figures", "FigureS7.png"), plot = FigS7, width = 13, height = 8, dpi = 300)

fst_pairs_all <- fst_melt_all %>% filter(Population1 != Population2) %>% mutate(p1=as.character(Population1), p2=as.character(Population2), lo=pmin(p1,p2), hi=pmax(p1,p2)) %>% distinct(lo, hi, .keep_all=TRUE) %>% transmute(Population1=lo, Population2=hi, Fst=Fst)
fst_pairs_all$Population1 <- factor(fst_pairs_all$Population1, levels = pop_order); fst_pairs_all$Population2 <- factor(fst_pairs_all$Population2, levels = pop_order)
fst_pairs_all <- fst_pairs_all %>% arrange(Population1, Population2)
write.csv(fst_pairs_all, here("..","results","pairwise_Fst_apiaries_all.csv"), row.names = FALSE)
fst_pairs_norecomb <- fst_melt_norecomb %>% filter(Population1 != Population2) %>% mutate(p1=as.character(Population1), p2=as.character(Population2), lo=pmin(p1,p2), hi=pmax(p1,p2)) %>% distinct(lo, hi, .keep_all=TRUE) %>% transmute(Population1=lo, Population2=hi, Fst=Fst); fst_pairs_norecomb$Population1 <- factor(fst_pairs_norecomb$Population1, levels = pop_order); fst_pairs_norecomb$Population2 <- factor(fst_pairs_norecomb$Population2, levels = pop_order); fst_pairs_norecomb <- fst_pairs_norecomb %>% arrange(Population1, Population2); write.csv(fst_pairs_norecomb, here("..","results","pairwise_Fst_apiaries_norecomb.csv"), row.names = FALSE)

fst_pairs_groups_all <- fst_melt_groups_all %>% filter(Population1 != Population2) %>% mutate(p1=as.character(Population1), p2=as.character(Population2), lo=pmin(p1,p2), hi=pmax(p1,p2)) %>% distinct(lo, hi, .keep_all=TRUE) %>% transmute(Population1=lo, Population2=hi, Fst=Fst)
grp_order <- c("European Dark\n honey bee samples","Other European\n samples","USA samples"); fst_pairs_groups_all$Population1 <- factor(fst_pairs_groups_all$Population1, levels=grp_order); fst_pairs_groups_all$Population2 <- factor(fst_pairs_groups_all$Population2, levels=grp_order)
fst_pairs_groups_all <- fst_pairs_groups_all %>% arrange(Population1, Population2)
write.csv(fst_pairs_groups_all, here("..","results","pairwise_Fst_all.csv"), row.names = FALSE)
fst_pairs_groups_norecomb <- fst_melt_groups_norecomb %>% filter(Population1 != Population2) %>% mutate(p1=as.character(Population1), p2=as.character(Population2), lo=pmin(p1,p2), hi=pmax(p1,p2)) %>% distinct(lo, hi, .keep_all=TRUE) %>% transmute(Population1=lo, Population2=hi, Fst=Fst); fst_pairs_groups_norecomb$Population1 <- factor(fst_pairs_groups_norecomb$Population1, levels=grp_order); fst_pairs_groups_norecomb$Population2 <- factor(fst_pairs_groups_norecomb$Population2, levels=grp_order); fst_pairs_groups_norecomb <- fst_pairs_groups_norecomb %>% arrange(Population1, Population2); write.csv(fst_pairs_groups_norecomb, here("..","results","pairwise_Fst_norecomb.csv"), row.names = FALSE)

# Fst Soni et al. 2022 test
# importing alignment with samples in only haplogroup 1 (hg1) or haplogroup 2 (hg2) 
alignment_hgs <- read.dna(here("..", "results", "vgseq_combined_msa_exons_hg1and2.fasta"), format = "fasta")

# Extract sequence names
sequence_names_hgs <- rownames(alignment_hgs)

# Extract population codes from sequence names
hgs <- sub(".*(hg[1-9]).*", "\\1", sequence_names_hgs)

# Convert to factor
hgs <- factor(hgs)

genind_obj_hgs <- DNAbin2genind(alignment_hgs, pop = hgs)
summary(genind_obj_hgs)

hierfstat_obj_hgs <- genind2hierfstat(genind_obj_hgs)

fst_results_hgs <- wc(hierfstat_obj_hgs)
fst_results_hgs #Global Fst all hgs

pairwise_fst_hgs <- pairwise.WCfst(hierfstat_obj_hgs)
pairwise_fst_hgs #Pairwise between hg1 and hg2 #0.15 

#per site
fst_per_locus <- fst_results_hgs$per.loc
fst_site <- fst_per_locus[, "FST"]

loc_pos <- as.integer(locNames(genind_obj_hgs))  
fst_df  <- data.frame(
  position = loc_pos,
  Fst      = fst_site
)

#Plotting Pairwise Fst per nucletoide in Vg
png(here("..", "figures", "FigureS9.png"), width = 13, height = 4, units = "in", res = 300)
par(mar = c(5, 4, 4, 2) + 0.5)
plot(fst_df$position, fst_df$Fst, pch = 16,
     xlab = "Vg gene position (bp)",
     ylab = "FST (hg1 vs hg2)",
     main = "Pairwise FST per nucleotide")
model <- lm(Fst ~ position, data = fst_df)
abline(model, col = "blue", lwd = 2)
mtext(label, side = 3, line = 0.5, adj = 1, col = "blue", cex = 0.9)
dev.off()

