# Create a VCF of variation in Apis melifera with AA information derived from
# sister species

# Assumptions:
# - First sequence is the reference
# - All outgroup sequences end in _ref
# - All HP1 and HP2 alternate and refer to the same sample

library(tidyverse)
library(ape)

# Hard code inputs and outputs
input_fasta <- "results/apis_vg.dna.fasta"
input_part <- "results/apis_vg.dna.part.txt"
output_vcf <- "results/apis_mellifera_vg.dna.vcf"

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
part <- part[start]

names(variants) <- as.character(pos)

# identify variants at each site
sites <- map(variants, alleles_to_indexes, outgroups)

#### CONSTRUCT VCF DATA ########################################################

tab <- sites |> map(function(x) {
    alt <- if(length(x$refalt) < 2) NA_character_ else 
        str_flatten(x$refalt[-1], ",")

    bees <- seq(1, length(x$alleles), 2)

    k <- tabulate(x$alleles)
    f <- sort(k, decreasing = TRUE)/sum(k)
    maf <- coalesce(f[2], 0)

    y <- list(
        ref = x$refalt[1],
        alt = alt,
        AA = x$aa,
        Cerana = x$cerana,
        Dorsata = x$dorsata,
        Laboriosa = x$laboriosa,
        Florea = x$florea,
        AN = sum(!is.na(x$alleles)),
        NS = sum(!is.na(x$alleles[bees]) & !is.na(x$alleles[bees+1])),
        MAF = maf,
        AC = str_flatten(k[-1], ",")
    )
    as_tibble_row(y)
}) |> list_rbind()

info <- tab |> select(AA:MAF) |> replace_na(list(AA = ".",
    Cerana = ".", Dorsata = ".", Laboriosa = ".", Florea = ".")) |>
    mutate(MAF = round(MAF, 4))

str <- str_c(names(info), "={", names(info), "}", collapse=";")

info <- str_glue_data(info, str)

ac <- if_else(tab$AC == "", "", str_c(";AC=", tab$AC))

pos <- as.integer(names(sites))
pos_part <- if_else(part %% 2 == 1,
    paste0("Exon=", (part+1) %/% 2),
    paste0("Intron=",(part+1) %/% 2) )

info <- str_glue("{pos_part};{info}{ac}")

fix <- tab |> transmute(`#CHROM` = "Vg",
    POS = pos,
    ID = NA_character_,
    REF = ref,
    ALT = alt,
    QUAL = NA_character_,
    FILTER = "PASS",
    INFO = info,
    FORMAT = "GT")

# gt data
gt <- sites |> map(pluck, "alleles") |>
    list_c() |> (`+`)(-1) |> as.character() |>
    coalesce(".") |>
    matrix(ncol = 2, byrow = TRUE)

# merge haplotypes into diplotypes
gt <- str_glue("{gt[,1]}|{gt[,2]}") |>
    matrix(nrow = length(sites), byrow = TRUE)

bees <- seq(1,nrow(mel),2)

colnames(gt) <- rownames(mel)[bees] |> str_replace("_HP1","")

df <- bind_cols(fix, gt)

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

write_file(header, output_vcf)
write_tsv(df, output_vcf, na = ".", append=TRUE, col_names=TRUE)
