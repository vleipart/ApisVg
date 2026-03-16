library("devtools")
library("ggplot2")
library(iMKT)
library(readxl)
library(tidyverse)
library(patchwork)
library(ggrepel)

daf_all_full <- read.delim(here("..", "results", "Daf_all_full.txt")) #no recombinated sequences
div_all_full <- read.delim(here("..", "results", "Div_all_full.txt")) #no recombinated seqeunces
daf_all_full_recomb <- read.delim(here("..", "results", "Daf_all_full_recomb.txt"))
div_all_full_recomb <- read.delim(here("..", "results", "Div_all_full_recomb.txt"))


#StandardMKT analysis
sMKTall <- standardMKT(daf = daf_all_full, divergence = div_all_full)
sMKTall_recomb <- standardMKT(daf = daf_all_full_recomb, divergence = div_all_full_recomb)

#FWW correction 
FWWall <- FWW(daf = daf_all_full, divergence = div_all_full, listCutoffs=c(0.10), plot=TRUE)
FWWall_recomb <- FWW(daf = daf_all_full_recomb, divergence = div_all_full_recomb, listCutoffs=c(0.10), plot=FALSE)

#Extended MKT
DGRPall <- DGRP(daf = daf_all_full, divergence = div_all_full, listCutoffs=c(0.10), plot=TRUE)
DGRPall_recomb <- DGRP(daf = daf_all_full_recomb, divergence = div_all_full_recomb, listCutoffs=c(0.10), plot=FALSE)

#Extracting results
FWWallResults <- FWWall$Results
FWWallResults_recomb <- FWWall_recomb$Results
DGRPallResults <- DGRPall$Results
DGRPallResults_recomb <- DGRPall_recomb$Results

#Making dataframe with results
a <- c(sMKTall$alpha.symbol, FWWallResults$alpha.symbol, DGRPallResults$alpha.symbol, sMKTall_recomb$alpha.symbol, FWWallResults_recomb$alpha.symbol, DGRPallResults_recomb$alpha.symbol)

p <- c(sMKTall$`Fishers exact test P-value`, FWWallResults$`Fishers exact test P-value`, DGRPallResults$`Fishers exact test P-value`, sMKTall_recomb$`Fishers exact test P-value`, FWWallResults_recomb$`Fishers exact test P-value`, DGRPallResults_recomb$`Fishers exact test P-value`)

Method <- c("sMKT", "FWW", "DGRP", "sMKT", "FWW", "DGRP")
Recomb <- c("Excluding recombinant alleles", "Excluding recombinant alleles", "Excluding recombinant alleles", 
            "Including recombinant alleles", "Including recombinant alleles", "Including recombinant alleles")
Diversity <- data.frame(Method,a,p)

#Plotting Results from all methods (Supplement Figure S10)
pall <- ggplot(Diversity, aes(y=a, x=fct_inorder(Method), fill = Recomb)) +
  geom_bar(stat="identity", position = "dodge") +
  labs(title = "", x="Methods", y="\u03b1")+
  scale_fill_manual(values = c("Excluding recombinant alleles" = "#9fc8c8", "Including recombinant alleles" = "#f0b077"),
                    name = "") +
  theme_classic() +
  theme(legend.position="top") +
  geom_text(aes(label=round(p,4)), position = position_dodge(width=0.9), vjust = -0.5, color = "black", size=4) +
  coord_cartesian(ylim=c(0, 0.8))

ggsave(filename = here("..", "figures", "FigureS10.png"), plot = pall, width = 8, height = 6, dpi = 300)

#Only DGRP results, but with del, neutral and weakly del
DGRP_a <- c(DGRPallResults$alpha.symbol, DGRPallResults_recomb$alpha.symbol)
DGRP_p <- c(DGRPallResults$`Fishers exact test P-value`, DGRPallResults_recomb$`Fishers exact test P-value`)
DiversityallGDPR <- data.frame(DGRP_a,DGRP_p)
colnames(DiversityallGDPR) <- c("\u03b1", "pvalue")

DGRP_Divergencemetrices <- DGRPall$`Divergence metrics`
DGRP_Divergencemetrices_recomb <- DGRPall_recomb$`Divergence metrics`

DGRP_estimates <- DGRP_Divergencemetrices$`Estimates by cutoff`
DGRP_estimates_recomb <- DGRP_Divergencemetrices_recomb$`Estimates by cutoff`

DiversityallGDPR$omegaA <- c(DGRP_estimates$omegaA.symbol, DGRP_estimates_recomb$omegaA.symbol)
DiversityallGDPR$omegaNA <- c(DGRP_estimates$omegaD.symbol, DGRP_estimates_recomb$omegaD.symbol)

DiversityallGDPR$type <- c("Vg")

df_Vg <- melt(DiversityallGDPR, id.vars = "type", measure.vars = c("\u03b1", "omegaA", "omegaNA"))
df_Vg$p <- c(DiversityallGDPR$p, rep(NA, nrow(df_Vg) - nrow(DiversityallGDPR)))
df_Vg$Recomb <- c("Excluding recombinant alleles", "Including recombinant alleles", "Excluding recombinant alleles", 
                  "Including recombinant alleles", "Excluding recombinant alleles", "Including recombinant alleles")

DGRP_fractions <- DGRPall$Fractions
DGRP_fractions_recomb <- DGRPall_recomb$Fractions
DGRP_fractions$type <- c("strongly deleterious (d)", "neutral (f)", "weakly deleterious (b)")
DGRP_fractions_recomb$type <- c("strongly deleterious (d)", "neutral (f)", "weakly deleterious (b)")
DGRP_fractions$Recomb <- "Excluding recombinant alleles"
DGRP_fractions_recomb$Recomb <- "Including recombinant alleles"
colnames(DGRP_fractions) <- c('fractions', 'type', 'Recomb')
colnames(DGRP_fractions_recomb) <- c('fractions', 'type', 'Recomb')

df_Vg_fractions <- rbind(DGRP_fractions, DGRP_fractions_recomb)
df_Vg_fractions$type <- factor(df_Vg_fractions$type,
                               levels = c("strongly deleterious (d)", "weakly deleterious (b)", "neutral (f)")
)
df_Vg_fractions <- df_Vg_fractions[order(df_Vg_fractions$Recomb, df_Vg_fractions$type), ]
df_Vg_fractions$fractions <- abs(df_Vg_fractions$fractions)

#Estimates for Table 1
df_Vg
df_Vg_fractions

#Estimates for Table S2
DGRPall$`MKT tables`
DGRPall_recomb$`MKT tables`
