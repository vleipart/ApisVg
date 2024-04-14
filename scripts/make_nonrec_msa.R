# Create a fasta file that does not contain any sequences
# that we have identified as recombinant.

library(tidyverse)
library(ape)

`%notin%` <- Negate(`%in%`)

# Hard code inputs and outputs
input_fasta <- here::here("results/apis_vg.dna.fasta")
input_tab <- here::here("results/full_haplotypes.csv.gz")
output_fasta <- here::here("results/apis_vg.nonrec.fasta")

# load data
aln <- read.dna(input_fasta, format="fasta", as.character = TRUE)
tab <- read_csv(input_tab)

# identify recombinant sequences
id <- tab |> filter(is_recomb == TRUE) |> pull(id)

# remove these sequences
aln <- aln[rownames(aln) %notin% id, ]

# save results
write.dna(aln, output_fasta, format="fasta", nbcol = 8, colsep = "")
