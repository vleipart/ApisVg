default: all

all:

.PHONY: default all

data/apis_vg.fasta: data/raw/Apis_mellifera_406088.fasta data/raw/Apis_cerana_108000069.fasta \
	data/raw/Apis_dorsata_102673109.fasta data/raw/Apis_laboriosa_122712558.fasta \
	data/raw/Apis_florea_100870965.fasta data/raw/Apis_mellifera_Vg_sequences.fasta
	cat $^ > $@

data/apis_vg.msa.fasta: data/apis_vg.fasta
	bash scripts/run.bash mafft --maxiterate 1000 --thread 8 $< > $@

# Truncate the msa to start and stop columns of A. mellifera reference.
# Assumes that first sequence is reference
data/apis_vg.cds.fasta data/apis_vg.cds.part.txt: data/apis_vg.msa.fasta data/raw/Apis_mellifera_406088.part.csv
	bash scripts/run.bash Rscript scripts/make_cds.R

