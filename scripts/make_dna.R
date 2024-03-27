# Truncate the msa to start and stop columns of A. mellifera reference.

library(tidyverse)
library(ape)

# Hard code inputs and outputs because we will have several outputs
input_fasta <- "data/apis_vg.msa.fasta"
input_part <- "data/raw/Apis_mellifera_406088.part.csv"
output_fasta <- "results/apis_vg.dna.fasta"
output_part <- "results/apis_vg.dna.part.txt"

#### MAIN ######################################################################

# read inputs
aln <- read.dna(input_fasta, format="fasta", as.character = TRUE)
part <- read_csv(input_part, skip=1, show_col_types = FALSE)

# map columns to residues in the reference sequence
ref <- aln[1,]
pos <- cumsum(ref != "-")

# trim the stop codon
part$end[nrow(part)] <- part$end[nrow(part)] - 3

# map columns to intervals
starts <- sort(c(part$start,part$end+1))
int <- findInterval(pos, starts)

# truncate the first and last intervals
b <- int != max(int) & int != min(int)
aln <- aln[,b]
int <- int[b]

# remove frameshift deletions identified via the ref
exons <- which(int %% 2 == 1)
exon_aln <- aln[, exons]
exon_ref <- str_flatten(exon_aln[1, ])

gaps <- str_locate_all(exon_ref, "-+")[[1]]
gaps <- gaps[(gaps[,2]-gaps[,1]) %% 3 != 2, ]
gaps <- unlist(map2(gaps[,1], gaps[,2], seq))

aln <- aln[, -exons[gaps]]
int <- int[-exons[gaps]]

# save alignment
write.dna(aln, output_fasta, format="fasta", nbcol = 8, colsep = "")

# save partition
write_lines(int, output_part)
