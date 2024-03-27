# Create a VCF of variation in Apis melifera with AA information derived from
# sister species

# Assumptions:
# - First sequence is the reference
# - All outgroup sequences end in _ref
# - All HP1 and HP2 alternate and refer to the same sample

library(tidyverse)
library(ape)

# Hard code inputs and outputs
input_fasta <- "data/apis_vg.cds.fasta"
input_part <- "data/apis_vg.cds.part.txt"
output_vcf <- "data/apis_mellifera_vg.cds.vcf"

#### HELPERS ###################################################################

`%notin%` <- Negate(`%in%`)

alleles_to_indexes <- function(x, outgroups) {
    names(x) <- NULL
    outgroups <- outgroups[outgroups != 1]
    x <- str_to_upper(x)
    ref <- x[1]
    out <- setNames(x[outgroups], names(outgroups))
    x <- x[-c(1, outgroups)]
    alleles <- sort(table(x), decreasing=TRUE)
    alt <- names(alleles)
    refalt <- c(ref, setdiff(alt, ref))
    nums <- match(x, refalt)
    if(length(refalt) == 1) {
        # No polymophism in mellifera
        aa <- refalt
    } else if(n_distinct(out) == 1 ) {
        # All outgroups have same allele
        aa <- out[1]
    } else if(n_distinct(out[-length(out)]) == 1) {
        # All outgroups but the farthest one have same allele
        aa <- out[1]
    } else if(alleles[[1]]/sum(alleles) >= 0.99) {
        # One allele is nearly fixed
        aa <- names(alleles)[1]
    } else {
        # Find the closest outgroup that matches.
        # This is only approximate because dorsata and laboriosa are sisters
        o <- match(out, refalt)
        o <- o[!is.na(o)]
        if(length(o) == 0) {
            # no outgroup matches
            aa <- NA_character_
        } else {
            aa <- refalt[o[1]]
        }
    }

    ret <- list(refalt = refalt,
                aa = aa,
                alleles = nums)
    ret <- c(ret, out)
    ret
}

#### LOAD DATA #################################################################

aln <- read.dna(input_fasta, format="fasta", as.character = TRUE)
rownames(aln) <- str_extract(rownames(aln), "\\S+")
part <- scan(input_part, quiet = TRUE)

## Split reference and non reference sequences
is_out <- str_detect(rownames(aln), "_ref$")
out <- aln[is_out,]
mel <- aln[!is_out, ]

outgroups <- str_replace(rownames(aln)[is_out], "^Apis_(.+)_ref$", "\\1")
outgroups <- setNames(which(is_out), outgroups)

#### IDENTIFY VARIANTS #########################################################

## Split alignment into segments based on gap patterns
m <- mel != "-"
no_gaps <- apply(m, 2, all)
groups <- cumsum(no_gaps)
groups <- split(seq_along(groups), groups)

## Convert to strings and extract variants
str <- apply(aln, 1, str_flatten)
start <- map_int(groups, 1L)
end <- map_int(groups, -1L)
variants <- str_sub_all(str, start, end)

vmat <- simplify(variants) |> str_remove_all(fixed("-")) |>
    matrix(nrow=length(groups))
vmat[vmat == ""] <- NA_character_
colnames(vmat) <- rownames(aln)
variants <- split(vmat, seq.int(nrow(vmat)))

pos <- cumsum(aln[1,] != "-")[start]
names(variants) <- as.character(pos)

# identify variants at each site
sites <- map(variants, alleles_to_indexes, outgroups)

#### CONSTRUCT VCF VARIANTS ####################################################

fix <- sites |> map(function(x) {
    y <- list(
        Ref = x$refalt[1],
        Alt = str_flatten_comma(x$refalt[-1]),
        N = str_flatten_comma(tabulate(x$alleles)),
        AA = x$aa,
        Cerana = x$cerana,
        Dorsata = x$dorsata,
        Laboriosa = x$laboriosa,
        Florea = x$florea
    )
    as_tibble_row(y)
}) |> list_rbind()

#### WRITE OUTPUT ##############################################################

header <- '##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=make_vcf.R
##contig=<ID=Vg,length=6094>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=Cerana,Number=1,Type=String,Description="Apis cerana Allele">
##INFO=<ID=Dorsata,Number=1,Type=String,Description="Apis dorsata Allele">
##INFO=<ID=Laboriosa,Number=1,Type=String,Description="Apis laboriosa Allele">
##INFO=<ID=Florea,Number=1,Type=String,Description="Apis florea Allele">
##INFO=<ID=Exon,Number=1,Type=Integer,Description="Which exon the variant is in">
##INFO=<ID=Intron,Number=1,Type=Integer,Description="Which intron the variant is in">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">
##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">
##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">
##INFO=<ID=HWE,Number=A,Type=Float,Description="HWE test (PMID:15789306); 1=good, 0=bad">
##INFO=<ID=ExcHet,Number=A,Type=Float,Description="Test excess heterozygosity; 1=good, 0=bad">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">
##FORMAT=<ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">
'

outfile <- "combined_msa_trunc.vcf"
write_file(header, outfile)
write_tsv(tab, outfile, na = ".", append=TRUE, col_names=TRUE)