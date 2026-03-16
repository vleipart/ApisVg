# ApisVg
Analysis of the Vitellogenin gene in Apis mellifera

## Data Cleaning

Ns were added to the Apis florea sequences to correct abiological frameshifts
in two locations. This matches XM_003689645.3.

## Species Tree

Reference sequences were aligned via mafft and neighbor joining plus midpoint
rooting inferred the following tree.

(((dorsata,laboriosa),(mellifera,cerana)),florea);

## Population genetic study of vitellogenin 

This repository contains the scripts used to generate the VCF file, perform the PCoA, MKT and FST analyses, and produce the figures for our population genetic study of honey bee vitellogenin, see (DOI). 

The script  (`scripts/`) uses input files from (`results/`), which includes:

- altaa.csv: A csv file for all biallelic polymorphisms and their cDNA position (POS) and if this polymorphisms leads to an amino acid change (altaa)
- apis_mellifera_vg.dna.vcf: VCF file
- Daf and Div files, outputs from iMKT webserver. One set with all samples (_recomb), one set where putative recombinant haplotypes. 
- Three .fasta files. vg_seq_combined_msa_exons.fasta is the CDS regions for all samples. The same is found in vg_seq_combined_msa_exons_norecomb.fasta, but putative recombinant haplotypes are removed. The vg_seq_combined_msa_exons_hg1and2.fasta is filtered to only contain samples considered to be in haplogroup 1 or 2 and 'hg1' or 'hg2' is added to the fasta header (>). 

The (`results/`) folder also includes four output files from FST analysis showing the pairs compared:

 - pairwise_Fst_all.csv: Comparing three sampling categories
 - pairwise_Fst_all_norecomb.csv: Comparing three sampling categories, removed putative recombinant haplotypes.
 - pairwise_Fst_apiaries.csv: Comparing all apiaries. 
 - pairwise_Fst_apiaries_norecomb.csv: Comparing all apiaries, but removed putative recombinant haplotypes. 

The (`figures/`) folder have all the generated main and supplement figures from the scripts. 
