library(stringdist)
library(tidyverse)
library(vcfR)
library(here)
library(ggrepel)
library(patchwork)
library(ggforce)
library(ggVennDiagram)

# Reading VCF created by make_vcf.R based on 543 diploid samples. Investigating biallelic sites
raw_vcf <- read.vcfR(here("..", "results", "apis_mellifera_vg.dna.vcf"),verbose = FALSE)
raw_vcf_fix <- vcfR2tidy(raw_vcf, info_only = TRUE)$fix

vcf <- raw_vcf[is.polymorphic(raw_vcf) &
                 is.biallelic(raw_vcf) & 
                 raw_vcf_fix$AA != "."]
vcf_fix <- vcfR2tidy(vcf, info_only = TRUE)$fix

#Analyzing biallelic sites
vg_positions <- data.frame(position = seq_len(nrow(vcf_fix)))
vg_positions_df <- cbind(vg_positions, POS = vcf_fix$POS, MAF = vcf_fix$MAF, Exon = vcf_fix$Exon, Intron = vcf_fix$Intron)
vg_positions_df$region <- ifelse(is.na(vg_positions_df$Exon), "Intron", "Exon")

polymorphic_sites_count <- sum(vg_positions_df$MAF > 0, na.rm = TRUE)

#Plotting all polymorphisms along the Vg gene (Supplement Figure 1A)
S1A <- ggplot(vg_positions_df, aes(x = POS, y = MAF, color = region)) +  geom_point() +
  annotate("rect", xmin = 0, xmax = 31, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 369, xmax = 1087, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 1236, xmax = 2135, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 2213, xmax = 3260, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 3356, xmax = 4981, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 5047, xmax = 5477, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 5540, xmax = 6094, ymin = 0, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("text", x = mean(c(369, 1087)), y = 0.57, label = "Exon 2", size = 3, color = "black") +
  annotate("text", x = mean(c(1236, 2135)), y = 0.57, label = "Exon 3", size = 3, color = "black") +
  annotate("text", x = mean(c(2213, 3260)), y = 0.57, label = "Exon 4", size = 3, color = "black") +
  annotate("text", x = mean(c(3356, 4981)), y = 0.57, label = "Exon 5", size = 3, color = "black") +
  annotate("text", x = mean(c(5047, 5477)), y = 0.57, label = "Exon 6", size = 3, color = "black") +
  annotate("text", x = mean(c(5540, 6094)), y = 0.57, label = "Exon 7", size = 3, color = "black")
S1A <- S1A + scale_color_manual(values = c("Exon" = "purple", "Intron" = "orange"))
S1A <- S1A + labs(title = "A) Polymorphisms MAF", x = "Vg Position", y = "MAF", color = "Region")
S1A <- S1A + theme_minimal() + theme(legend.position = "top")

# Extract phased haplotypes and encode variants as 0 = ancestral and 1 = derived
full_haplotypes_mat <- (extract.haps(vcf, verbose = FALSE) != vcf_fix$AA)*1L
colnames(full_haplotypes_mat) <- colnames(full_haplotypes_mat) |>
  str_replace_all(c("_0$" = "_HP1", "_1$" = "_HP2"))
full_haplotypes_str <- apply(full_haplotypes_mat, 2, str_flatten)
full_reference <- str_flatten((vcf_fix$REF != vcf_fix$AA)*1L)
full_ancestor <- str_dup("0", nrow(vcf_fix))

full_cerana <- str_flatten(case_when(
  vcf_fix$Cerana == vcf_fix$AA ~ "0",
  vcf_fix$Cerana == vcf_fix$REF ~ "1",
  vcf_fix$Cerana == vcf_fix$ALT ~ "1",
  TRUE ~ "2"
))

# Count abundance of haplotypes and plot with PCoA
haplotypes_tab <- tibble(full = full_haplotypes_str) |> count(full)
full_dist <- stringdistmatrix(c(full_reference, full_ancestor,
                                full_cerana,
                                haplotypes_tab$full),
                              method = "hamming")
v <- cmdscale(full_dist/nrow(haplotypes_tab), 2)
haplotypes_tab$x <- v[-c(1:3),1]
haplotypes_tab$y <- v[-c(1:3),2]

haplotypes_tab_sorted <- haplotypes_tab %>%
  arrange(desc(n))

haplotypes_tab_sorted <- haplotypes_tab_sorted %>%
  rename(haplotype_count = n)

# Plotting haplotype count (Supplement Figure 1B)
S1B <- ggplot(haplotypes_tab_sorted, aes(x = seq_along(haplotype_count), y = haplotype_count)) +
  geom_bar(stat = "identity", width = 1)
S1B <- S1B + labs(title = "B) Haplotypes Frequency", x = "Haplotypes", y = "Count")
S1B <- S1B + theme_minimal()

# PoCA base plot (Supplement Figure 1C)
S1C <- ggplot(haplotypes_tab, aes(x = x, y = y, size = n)) + 
  geom_point(alpha = 0.5) +
  xlab("Axis 1") + ylab("Axis 2") + 
  ggtitle("C) PCoA of all polymorphisms") +
  scale_size_continuous("Count") +
  coord_fixed()

# Create a data frame for the reference and orthologous species vg sequences
full_reference_xy <- list(x = v[1,1], y = v[1,2])
full_ancestor_xy <- list(x = v[2,1], y = v[2,2])
full_cerana_xy <- list(x = v[3,1], y = v[3,2])

special_points <- data.frame(
  x = c(full_reference_xy$x, full_ancestor_xy$x, full_cerana_xy$x),
  y = c(full_reference_xy$y, full_ancestor_xy$y, full_cerana_xy$y),
  label = c("Reference", "Ancestor", "Outgroup")
)

# Add reference and orthologous points with color aesthetic
S1Cadd <- S1C + 
  geom_point(
    data = special_points, 
    aes(x = x, y = y, color = label),
    size = 1
  ) +
  scale_color_manual(
    name = "",  # Legend title
    values = c("Reference" = "red", "Ancestor" = "blue", "Outgroup" = "green")
  )

# Calculating the explained variance in PCoA
haplo_dist <- stringdistmatrix(haplotypes_tab$full, method = "hamming")
v_full <- cmdscale(haplo_dist, k = 2, eig = TRUE)
eig_vals <- v_full$eig
prop_var <- eig_vals / sum(eig_vals[eig_vals > 0])

n_axes <- 2

var_df <- tibble(
  axis        = factor(1:n_axes,
                       labels = paste0("Axis ", 1:n_axes)),
  prop_var    = prop_var[1:n_axes],
  cum_var     = cumsum(prop_var[1:n_axes])
)

# Plotting explained variance for axis 1 and 2 (Supplement Figure 1D)
S1D <- ggplot(var_df, aes(x = axis, y = prop_var)) +
  geom_col() +
  xlab("PCoA axis") +
  ylab("Proportion of variance explained") +
  ggtitle("D) Variance explained by PCoA axes 1 and 2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Saving Supplement Figure 1
layout_1 <- "
AAAAB
CCCDD
CCCDD
"
FigureS1 <- S1A + S1B + S1C + S1D + plot_layout(design = layout_1)
ggsave(filename = here("..", "figures", "FigureS1.png"), plot = FigureS1, width = 13, height = 6, dpi = 300)


#Plotting distribution of polymorphisms by MAF (Figure 1A)
Fig1A <- ggplot(vcf_fix, aes(x = MAF)) + geom_histogram(bins = 20)
Fig1A <- Fig1A + ggtitle("A) MAF distribution")
Fig1A

# Here we observe a break in the dataset on 0.25, using this as a thershold for "common"
base_o <- vcf_fix$MAF > 0.25

#Filtering out polymorphisms <0.25 (Figure 1B)
Fig1B <- ggplot(vg_positions_df %>% filter(MAF > 0.25), aes(x = POS, y = MAF, color = region)) + geom_point() +
  annotate("rect", xmin = 0, xmax = 31, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 369, xmax = 1087, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 1236, xmax = 2135, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 2213, xmax = 3260, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 3356, xmax = 4981, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 5047, xmax = 5477, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = 5540, xmax = 6094, ymin = 0.25, ymax = 0.55, fill = "purple", alpha = 0.15) +
  annotate("text", x = mean(c(369, 1087)), y = 0.57, label = "Exon 2", size = 3, color = "black") +
  annotate("text", x = mean(c(1236, 2135)), y = 0.57, label = "Exon 3", size = 3, color = "black") +
  annotate("text", x = mean(c(2213, 3260)), y = 0.57, label = "Exon 4", size = 3, color = "black") +
  annotate("text", x = mean(c(3356, 4981)), y = 0.57, label = "Exon 5", size = 3, color = "black") +
  annotate("text", x = mean(c(5047, 5477)), y = 0.57, label = "Exon 6", size = 3, color = "black") +
  annotate("text", x = mean(c(5540, 6094)), y = 0.57, label = "Exon 7", size = 3, color = "black")
Fig1B <- Fig1B + scale_color_manual(values = c("Exon" = "purple", "Intron" = "orange"))
Fig1B <- Fig1B + labs(title = "B) Polymorphisms with MAF > 0.25", x = "Vg Position", y = "MAF", color = "Region")
Fig1B <- Fig1B + theme_minimal() + theme(legend.position = "top")

# create base haplotypes (haplotypes based on polymorphisms >0.25 MAF)
base_haplotypes_mat <- (extract.haps(vcf[base_o], verbose = FALSE) !=
                          vcf_fix[base_o, ]$AA)*1L
colnames(base_haplotypes_mat) <- colnames(base_haplotypes_mat) |>
  str_replace_all(c("_0$" = "_HP1", "_1$" = "_HP2"))
base_haplotypes_str <- apply(base_haplotypes_mat, 2, str_flatten)

# Count base sequences
base_tab <- tibble(base = base_haplotypes_str) |> count(base)
base_tab |> slice_max(n, n = 6) |> knitr::kable()

base_ref <- base_tab |> slice_max(n, n = 1) |> pull(base)

base_dist <- stringdistmatrix(base_tab$base, method = "hamming")
v <- cmdscale(base_dist/nrow(base_tab), 2)
base_tab$x <- v[,1]
base_tab$y <- v[,2]

# Sort the data by haplotype_count in descending order
base_tab <- base_tab %>%
  rename(haplotype_count = n)

base_tab_sorted <- base_tab %>%
  arrange(desc(haplotype_count))

#Plotting haplotype count (Figure 1C)
Fig1C <- ggplot(base_tab_sorted, aes(x = seq_along(haplotype_count), y = haplotype_count)) +
  geom_bar(stat = "identity", width = 1)
Fig1C <- Fig1C + labs(title = "C) Haplotype\nFrequency\n>0.25 MAF", x = "Haplotypes", y = "Count")
Fig1C <- Fig1C + theme_minimal()

#Selecting the two most frequent haplotypes and calculating the hamming distance to all other haplotypes
haplotypes_tab$base <- haplotypes_tab$full |>
  str_sub_all(which(base_o), which(base_o)) |>
  map_chr(str_flatten)

haplotypes_tab$base_dist <- stringdistmatrix(
  haplotypes_tab$base, base_ref, "hamming")[,1]/sum(base_o)

#Identified haplotypes with distance < 0.2 to the two common haplotypes to define the haplogroups calculate the the maximum Euclidean distance
brown_pts <- haplotypes_tab[haplotypes_tab$base_dist > 0.8,]
brown_x_mean <- mean(brown_pts$x)
brown_y_mean <- mean(brown_pts$y)
brown_r <- max(sqrt((brown_pts$x - brown_x_mean)^2 + (brown_pts$y - brown_y_mean)^2))

brown_circle <- data.frame(
  x0 = brown_x_mean,
  y0 = brown_y_mean,
  r  = brown_r
)

green_pts <- haplotypes_tab[haplotypes_tab$base_dist < 0.2, ]
green_x_mean <- mean(green_pts$x)
green_y_mean <- mean(green_pts$y)
green_r <- max(sqrt((green_pts$x - green_x_mean)^2 + (green_pts$y - green_y_mean)^2))

green_circle <- data.frame(
  x0 = green_x_mean,
  y0 = green_y_mean,
  r  = green_r
)

#Coloring the PCoA by "haplogroup distance" and adding the haplogroup circles (Euclidean distance)
Fig1D <- ggplot(haplotypes_tab, aes(x = x, y = y, size = n, color = -1+2*(1-base_dist))) + geom_point()
Fig1D <- Fig1D + xlab("Axis 1") + ylab("Axis 2") + ggtitle("D) PCoA of all polymorphisms")
Fig1D <- Fig1D + scale_size_continuous("Count")
Fig1D <- Fig1D + scale_color_gradient2(
  name = "Haplogroup",
  low = "#BF7F30",  
  mid = "yellow",   
  high = "#358C7C", 
  midpoint = 0
)
Fig1D <- Fig1D + coord_fixed()
Fig1D <- Fig1D + geom_circle(data = brown_circle, aes(x0 = x0, y0 = y0, r = r), inherit.aes = FALSE,
    color = "black",linewidth = 0.5)
Fig1D <- Fig1D + geom_circle(data = green_circle, aes(x0 = x0, y0 = y0, r = r), inherit.aes = FALSE,
                           color = "black",linewidth = 0.5)

#Saving Figure 1
layout_2 <- "
ABBBB
CDDDD
CDDDD
"
Figure1 <- Fig1A + Fig1B + Fig1C + Fig1D + plot_layout(design = layout_2)
ggsave(filename = here("..", "figures", "Figure1.png"), plot = Figure1, width = 9, height = 6, dpi = 300)

#Plotting PCoA with >025 MAF haplotypes (Supplement Figure 2)
S2 <- ggplot(base_tab, aes(x = x, y = y, size = haplotype_count)) + geom_point(alpha = 0.5)
S2 <- S2 + xlab("Axis 1") + ylab("Axis 2")
S2 <- S2 + scale_size_continuous("Count")
S2 <- S2 + coord_fixed()
S2

ggsave(filename = here("..", "figures", "FigureS2.png"), plot = S2, width = 10, height = 5, dpi = 300)

## Recombinant analysis
# We identify recombinant sequences by looking for haplotypes based on polymorphisms >0.25 (base haplotypes) 
# that share fragments of the standard base haplotypes 

#| fig.height=5
base_ref <- base_tab |> slice_max(order_by = haplotype_count, n = 2, with_ties = FALSE)
base_ref_1 <- base_ref$base[1]
base_ref_2 <- base_ref$base[2]
base_anc <- str_dup("0", str_length(base_ref_1))

base_ref_1_v <- str_split_1(base_ref_1, "")
base_ref_2_v <- str_split_1(base_ref_2, "")

base_v <- str_split(base_tab$base, "") |>
  map(function(x) {
    case_when(
      x == base_ref_1_v ~ "1",
      x == base_ref_2_v ~ "2"
    )
  })
base_tab$base2 <- map_chr(base_v, str_flatten)

base_v <- str_split(base_tab$base, "") |>
  map(function(x) {
    case_when(
      x == "0" ~ "0",
      x == base_ref_1_v ~ "1",
      x == base_ref_2_v ~ "2"
    )
  })

base_tab$base3 <- map_chr(base_v, str_flatten)

# Calculate a delta statistics
# https://mol.ax/software/delta
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1894573/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5850291/

n <- str_length(base_ref_1)
df <- expand_grid(i = 1:(n-1))
r1 <- df |> mutate(s = str_c(
  str_sub(base_ref_1, 1, i),
  str_sub(base_ref_2, i+1, n)
))
r2 <- df |> mutate(s = str_c(
  str_sub(base_ref_2, 1, i),
  str_sub(base_ref_1, i+1, n)
))

df <- expand_grid(i = 1:(n-1), j = 1:(n-1)) |> filter(i < j)

r3 <- df |> mutate(s = str_c(
  str_sub(base_ref_1, 1, i),
  str_sub(base_ref_2, i+1, j),
  str_sub(base_ref_1, j+1, n)
))
r4 <- df |> mutate(s = str_c(
  str_sub(base_ref_2, 1, i),
  str_sub(base_ref_1, i+1, j),
  str_sub(base_ref_2, j+1, n)
))

rec_set <- c(base_anc, base_ref_1, base_ref_2, r1$s, r2$s, r3$s, r4$s)

d <- stringdistmatrix(base_tab$base, rec_set, method = "hamming")

base_tab$delta0 <- d[,1] - apply(d[,1:3], 1, min)
base_tab$delta1 <- apply(d[,1:3], 1, min) - apply(d[,1:(2*n+1)], 1, min)
base_tab$delta2 <- apply(d[,1:3], 1, min) - apply(d, 1, min)

haplotypes_tab_orig <- haplotypes_tab 
df <- base_tab |> select(base, base2, base3, delta0, delta1, delta2)
haplotypes_tab <- left_join(haplotypes_tab_orig, df) |>
  mutate(is_recomb = delta2 >= 4)

# Write recombination information to results
full_tab <- enframe(full_haplotypes_str, name = "id", value = "full")
full_tab <- full_tab |> separate_wider_regex(id, c("[^_]+_",
                                                   "apiary" = "..", ".*"), cols_remove = FALSE) |>
  separate_wider_delim(id, "_", names = c("sample", "name", NA),
                       cols_remove = FALSE)
full_tab <- full_tab |> mutate(region = case_match(apiary,
                                                   c("NC", "CA", "MD", "AZ", "MN", "IL") ~ "USA samples",
                                                   c("FN", "SW", "RN", "CS", "DK", "IR", "PL", "TX", "FR") ~ "European Dark honey bee samples",
                                                   c("SL", "IT", "PO", "MK", "MT", "TR") ~ "Other European samples"
))
full_tab <- full_tab |> relocate(id, name, sample, apiary, region)
full_tab <- full_tab |> left_join(haplotypes_tab)

sum(full_tab$is_recomb == TRUE) #The number of recombinant haplotypes

#Plotting PCoA with recombinant haplotypes highlighted (Supplement Figure 3)
S3 <- ggplot(haplotypes_tab, aes(x = x, y = y, size = n, color = base_dist,
                                 shape = is_recomb)) + geom_point()
S3 <- S3 + annotate("point", color = "red", x = full_reference_xy$x, y = full_reference_xy$y)
S3 <- S3 + annotate("point", color = "blue", x = full_ancestor_xy$x, y = full_ancestor_xy$y)
S3 <- S3 + annotate("point", color = "green", x = full_cerana_xy$x, y = full_cerana_xy$y)
S3 <- S3 + xlab("Axis 1") + ylab("Axis 2")
S3 <- S3 + scale_size_continuous("Count")
S3 <- S3 + scale_shape_discrete("Recombinant")
S3 <- S3 + coord_fixed()
S3 <- S3 + scale_color_gradient2(
  name = "Haplogroup",
  low = "#BF7F30",
  mid = "yellow", 
  high = "#358C7C",
  midpoint = 0
)

ggsave(filename = here("..", "figures", "FigureS3.png"), plot = S3, width = 10, height = 5, dpi = 300)


## Sampling Category
full_by_region <- full_tab |> count(full, region, x, y, sort = TRUE)

Fig2A <- ggplot(full_by_region, aes(x = x, y = y, size = n, color = region)) +
  geom_point(alpha = 0.5)
Fig2A <- Fig2A + xlab("Axis 1") + ylab("Axis 2") + ggtitle("A) PCoA per\nsampling category")
Fig2A <- Fig2A + scale_size_continuous("Count")
Fig2A <- Fig2A + scale_color_brewer("Sampling category", palette = "Set1")
Fig2A <- Fig2A + coord_fixed()
Fig2A <- Fig2A + facet_wrap(vars(region), ncol = 1)
Fig2A <- Fig2A + guides(color = "none") + theme(legend.position = "bottom")
Fig2A

full_by_apiary <- full_tab |> count(full, apiary, region, x, y, sort = TRUE)

full_by_apiary_sorted <- full_by_apiary |> 
  mutate(apiary = factor(apiary, levels = unique(apiary[order(region)])))

Fig2B <- ggplot(full_by_apiary_sorted, aes(x = x, y = y, size = n, color = region)) +
  geom_point(alpha = 0.5)
Fig2B <- Fig2B + xlab("Axis 1") + ylab("Axis 2") + ggtitle("B) PCoA per aipary")
Fig2B <- Fig2B + scale_size_continuous("Count")
Fig2B <- Fig2B + scale_color_brewer("Region", palette = "Set1")
Fig2B <- Fig2B + coord_fixed()
Fig2B <- Fig2B + facet_wrap(vars(apiary), ncol = 7)
Fig2B <- Fig2B + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) 
Fig2B

# calculate the distances between each sequence pair
tab <- full_tab |> arrange(apiary)
full_dist <- stringdistmatrix(tab$full, method = "hamming") / 
  nrow(haplotypes_tab)

# construct a table to hold the pairwise distance data
x <- combn(nrow(tab), 2)
a <- x[1, ]
b <- x[2, ]
full_dist_tab <- tibble(
  id_1 = tab$id[a],
  id_2 = tab$id[b],
  apiary_1 = tab$apiary[a],
  apiary_2 = tab$apiary[b],
  dist = full_dist )

# calculate the mean distance between pairs of apiaries
apiary_dist_tab <- summarize(full_dist_tab, d = mean(dist),
                             .by = c(apiary_1, apiary_2))

# convert results into a distance object
labs <- unique(tab$apiary)
apiary_dist <- structure(filter(apiary_dist_tab, apiary_1 != apiary_2)$d,
                         Size = length(labs),
                         Labels = labs,
                         Diag = FALSE,
                         Upper = FALSE,
                         class = "dist")

# plot a hierarchical clustering of the apiaries
plot(hclust(apiary_dist))

#Replot the denogram excluding recombiants
full_tab_norecomb <- full_tab[full_tab$is_recomb == FALSE, ]

tab_norecomb <- full_tab_norecomb |> arrange(apiary)
full_dist_norecomb <- stringdistmatrix(tab_norecomb$full, method = "hamming") / 
  nrow(haplotypes_tab)

# construct a table to hold the pairwise distance data
x2 <- combn(nrow(tab_norecomb), 2)
a2 <- x2[1, ]
b2 <- x2[2, ]
full_dist_tab_norecomb <- tibble(
  id_1_1 = tab_norecomb$id[a2],
  id_2_2 = tab_norecomb$id[b2],
  apiary_1_1 = tab_norecomb$apiary[a2],
  apiary_2_2 = tab_norecomb$apiary[b2],
  dist_norecomb = full_dist_norecomb )

# calculate the mean distance between pairs of apiaries
apiary_dist_tab_norecomb <- summarize(full_dist_tab_norecomb, d = mean(dist_norecomb),
                             .by = c(apiary_1_1, apiary_2_2))

# convert results into a distance object
labs_norecomb <- unique(tab_norecomb$apiary)
apiary_dist_norecomb <- structure(filter(apiary_dist_tab_norecomb, apiary_1_1 != apiary_2_2)$d,
                         Size = length(labs),
                         Labels = labs,
                         Diag = FALSE,
                         Upper = FALSE,
                         class = "dist")

# plot a hierarchical clustering of the apiaries
plot(hclust(apiary_dist_norecomb))
plot(hclust(apiary_dist))

# plot a PCoA analysis of apiaries
v <- cmdscale(apiary_dist, 2)
apiary_xy <- tibble(apiary = rownames(v), x = v[,1], y = v[,2])
apiary_xy <- apiary_xy |> mutate(region = case_match(apiary,
                                                     c("NC", "CA", "MD", "AZ", "MN", "IL") ~ "USA samples",
                                                     c("FN", "SW", "RN", "CS", "DK", "IR", "PL", "TX", "FR") ~ "European Dark honey bee samples",
                                                     c("SL", "IT", "PO", "MK", "MT", "TR") ~ "Other European samples"
))

Fig2C <- ggplot(apiary_xy, aes(x = x, y = y, label = apiary, color = region)) +
  geom_point() + geom_text_repel(show.legend = FALSE)
Fig2C <- Fig2C + coord_fixed()
Fig2C <- Fig2C + ggtitle("C) Mean distance")+ xlab("Axis 1") + ylab("Axis 2")
Fig2C <- Fig2C + scale_color_brewer("Region", palette = "Set1")
Fig2C <- Fig2C + theme(legend.position = "none")

Fig2C

#Fig 2
layout_3 <- "
ABBB
ABBB
ABBB
ACCC
ACCC
"
Figure2 <- Fig2A + Fig2B + Fig2C + plot_layout(design = layout_3)
ggsave(filename = here("..", "figures", "Figure2.png"), plot = Figure2, width = 10, height = 7, dpi = 300)


## Defining HG groups and fetching samples
n95_df <- full_tab[full_tab$base_dist == 0,] #Haplotype n95
n67_df <- full_tab[full_tab$base_dist == 1,] #Haplotype n67

hg1_df <- full_tab[full_tab$base_dist <= 0.2,] #Haplogroup 1
hg2_df <- full_tab[full_tab$base_dist >= 0.8,] #Haplogroup 2

#Numbers for Supplement Figure S5 and S6
hp_pie_df <- data.frame(
  type = c("haplotype n95", "haplotype n67", "haplogroup 1", "haplogroup 2"),
  tot  = c(sum(full_tab$base_dist == 0,   na.rm = TRUE), sum(full_tab$base_dist == 1,   na.rm = TRUE), sum(full_tab$base_dist <= 0.2, na.rm = TRUE), sum(full_tab$base_dist >= 0.8, na.rm = TRUE)),
  amm = c(sum(n95_df$region == "European Dark honey bee samples", na.rm = TRUE), sum(n67_df$region == "European Dark honey bee samples", na.rm = TRUE), sum(hg1_df$region == "European Dark honey bee samples", na.rm = TRUE),sum(hg2_df$region == "European Dark honey bee samples", na.rm = TRUE)),
  amx = c(sum(n95_df$region == "Other European samples", na.rm = TRUE), sum(n67_df$region == "Other European samples", na.rm = TRUE),sum(hg1_df$region == "Other European samples", na.rm = TRUE),sum(hg2_df$region == "Other European samples", na.rm = TRUE)), 
  usa = c(sum(n95_df$region == "USA samples", na.rm = TRUE),sum(n67_df$region == "USA samples", na.rm = TRUE),sum(hg1_df$region == "USA samples", na.rm = TRUE),sum(hg2_df$region == "USA samples", na.rm = TRUE)))

cols <- c("European Dark honey bee samples" = "#E41A1C", "Other European samples" = "#377EB8","USA samples" = "#4DAF4A")

#Plotting Supplement Figure 5, distribution of sampling categories per haplotype and haplogroup
make_pie_from_wide <- function(df_wide, one_type, title_text) {
  df_tmp <- df_wide %>%
    filter(type == one_type) %>%
    select(type, amm, amx, usa) %>%
    pivot_longer(cols = c(amm, amx, usa), names_to = "sampling_category", values_to = "n") %>%
    mutate(
      sampling_category = recode(sampling_category, amm = "European Dark honey bee samples", amx = "Other European samples",usa = "USA samples"),
      prop = n / sum(n),
      pct_label = ifelse(n == 0, "", sprintf("%.1f%%", 100 * prop))
    )
  
  ggplot(df_tmp, aes(x = "", y = n, fill = sampling_category)) +
    geom_col(width = 1, alpha = 0.5) +
    coord_polar(theta = "y") +
    geom_text(aes(label = pct_label), position = position_stack(vjust = 0.5), size = 5) +
    scale_fill_manual(values = cols, name = "") +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom")
}

p1 <- make_pie_from_wide(hp_pie_df, "haplotype n95", "Haplotype n = 95")
p2 <- make_pie_from_wide(hp_pie_df, "haplotype n67", "Haplotype n = 67")
p3 <- make_pie_from_wide(hp_pie_df, "haplogroup 1",  "Haplogroup with n = 95")
p4 <- make_pie_from_wide(hp_pie_df, "haplogroup 2",  "Haplogroup with n = 67")

FigS5 <- (p1 | p2) / (p3 | p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & guides(fill = guide_legend(ncol = 1))

ggsave(filename = here("..", "figures", "FigureS5.png"), plot = FigS5, width = 10, height = 7, dpi = 300)


#Plotting Supplement Figure 6, distribution of haplogroup per sampling category
# haplogroup colors (n95 vs n67)

hg_cols <- c("Haplogroup with n95" = "#358C7C","Haplogroup with n67" = "#BF7F30")

hg_region_df <- hp_pie_df %>%
  filter(type %in% c("haplogroup 1", "haplogroup 2")) %>%
  mutate(haplogroup = recode(type, "haplogroup 1" = "Haplogroup with n95", "haplogroup 2" = "Haplogroup with n67")) %>%
  select(haplogroup, amm, amx, usa) %>%
  pivot_longer(cols = c(amm, amx, usa), names_to = "sampling_category", values_to = "n") %>%
  mutate(
    sampling_category = recode(sampling_category, amm = "European Dark honey bee samples", amx = "Other European samples", usa = "USA samples")) %>%
  group_by(sampling_category) %>%
  mutate(
    prop = n / sum(n),
    pct_label = ifelse(n == 0, "", sprintf("%.1f%%", 100 * prop))
  ) %>%
  ungroup()

make_hg_pie <- function(df_long, one_region, title_text) {
  df_tmp <- df_long %>% filter(sampling_category == one_region)
  
  ggplot(df_tmp, aes(x = "", y = n, fill = haplogroup)) +
    geom_col(width = 1, alpha = 0.9) +
    coord_polar(theta = "y") +
    geom_text(aes(label = pct_label),position = position_stack(vjust = 0.5), size = 5) +
    scale_fill_manual(values = hg_cols, name = "Haplogroup") +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
}

p_amm <- make_hg_pie(hg_region_df, "European Dark honey bee samples", "European Dark honey bee samples")
p_amx <- make_hg_pie(hg_region_df, "Other European samples", "Other European samples")
p_usa <- make_hg_pie(hg_region_df, "USA samples", "USA samples")

FigS6 <- (p_amm | p_amx | p_usa) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & guides(fill = guide_legend(ncol = 2))
ggsave(filename = here("..", "figures", "FigureS6.png"), plot = FigS6, width = 10, height = 5, dpi = 300)

#Chi square statistic:
tab <- matrix(c(hp_pie_df$amm[hp_pie_df$type == "haplogroup 1"],  # European Dark, hg1
    hp_pie_df$amm[hp_pie_df$type == "haplogroup 2"],  # European Dark, hg2
    hp_pie_df$amx[hp_pie_df$type == "haplogroup 1"],  # Other European, hg1
    hp_pie_df$amx[hp_pie_df$type == "haplogroup 2"],  # Other European, hg2
    hp_pie_df$usa[hp_pie_df$type == "haplogroup 1"],  # USA, hg1
    hp_pie_df$usa[hp_pie_df$type == "haplogroup 2"]   # USA, hg2
  ),
  nrow = 3, byrow = TRUE, dimnames = list(c("European Dark honey bee samples", "Other European samples", "USA samples"), c("Haplogroup with n95", "Haplogroup with n67"))
)

test <- chisq.test(tab, correct = FALSE)

cat("Chi-square statistic:", unname(test$statistic), "\n")
cat("p-value:", test$p.value, "\n")
cat("Degrees of freedom:", unname(test$parameter), "\n")

#Number of haplotypes observed per sampling category and aipary for supplement Table
uniquehp_summary <- full_tab %>%
  group_by(region, apiary, base) %>%
  summarise(hg1 = any(base_dist < 0.2, na.rm = TRUE), hg2 = any(base_dist >= 0.8, na.rm = TRUE),.groups = "drop") %>%
  group_by(region, apiary) %>%
  summarise(n_unique_bases = n(), hg1 = sum(hg1),hg2 = sum(hg2), .groups = "drop") %>%
  arrange(region, apiary)

#Soni test 2022

#Import syn and non-syn at aa positions
altaa_df <- read_csv(here("..", "results", "altaa.csv"))
vg_positions_df <- vg_positions_df %>% left_join(altaa_df, by = "POS")

#Get variant information per sample in either hg1 and hg2
gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
sample_ids <- colnames(gt_chr)

alt_on_hap <- function(gt_vec, hap = 1) {
  gt_vec <- as.character(gt_vec)
  out <- rep(NA, length(gt_vec))
  ok <- !is.na(gt_vec) & str_detect(gt_vec, "\\|") & !gt_vec %in% c(".", "./.", ".|.")
  parts <- str_split(gt_vec[ok], "\\|", simplify = TRUE)
  out[ok] <- parts[, hap] == "1"
  out
}

# Assing a samples to group from which region samples are from
apiary_code <- str_extract(sample_ids, "(?<=_)\\w{2}(?=\\d)")
sample_group <- case_when(
  apiary_code %in% c("NC","CA","MD","AZ","MN","IL") ~ "usa",
  apiary_code %in% c("FN","SW","RN","CS","DK","IR","PL","TX","FR") ~ "amm",
  apiary_code %in% c("SL","IT","PO","MK","MT","TR") ~ "amx",
  TRUE ~ NA_character_
)

# Find how many times the snp is observed by group
observed_by_group_simple <- function(df) {
  members <- df %>%
    transmute(
      sample = str_remove(id, "_HP[12]$"),
      hap    = str_extract(id, "HP[12]$"),
      hap_i  = if_else(hap == "HP1", 1L, 2L),
      group  = sample_group[match(sample, sample_ids)]
    ) %>%
    filter(!is.na(group))
  groups <- c("usa", "amm", "amx")
  map_dfc(groups, function(g) {
    mem_g <- members %>% filter(group == g)
    mat_list <- pmap(
      mem_g %>% select(sample, hap_i),
      function(sample, hap_i) {
        j <- match(sample, sample_ids)
        alt_on_hap(gt_chr[, j], hap = hap_i)
      }
    )
    
    mat <- do.call(cbind, mat_list)
    tibble(!!g := apply(mat, 1, \(v) any(v %in% TRUE, na.rm = TRUE)))
  })
}

obs_all <- observed_by_group_simple(full_tab)
vg_positions_df <- bind_cols(vg_positions_df, obs_all)

#Counting private and shared nonsyn and syn SNPs
vg_positions_df_cds <- vg_positions_df[!is.na(vg_positions_df$Exon), ]
vg_positions_df_cds
p_ns_amm = sum(with(vg_positions_df_cds, altaa & (amm) & !(amx | usa)))
p_syn_amm = sum(with(vg_positions_df_cds, !altaa & (amm) & !(amx | usa)))

p_ns_amx = sum(with(vg_positions_df_cds, altaa & (amx) & !(amm | usa)))
p_syn_amx = sum(with(vg_positions_df_cds, !altaa & (amx) & !(amm | usa)))

p_ns_usa = sum(with(vg_positions_df_cds, altaa & (usa) & !(amm | amx)))
p_syn_usa = sum(with(vg_positions_df_cds, !altaa & (usa) & !(amm | amx)))

s_ns_amm_amx = sum(with(vg_positions_df_cds, altaa & (amm) & (amx) & !(usa)))
s_ns_amm_usa = sum(with(vg_positions_df_cds, altaa & (amm) & (usa) & !(amx)))
s_ns_amx_usa = sum(with(vg_positions_df_cds, altaa & (amx) & (usa) & !(amm)))
s_ns_all = sum(with(vg_positions_df_cds, altaa & (amm) & (usa) & (amx)))

s_syn_amm_amx = sum(with(vg_positions_df_cds, !altaa & (amm) & (amx) & !(usa)))
s_syn_amm_usa = sum(with(vg_positions_df_cds, !altaa & (amm) & (usa) & !(amx)))
s_syn_amx_usa = sum(with(vg_positions_df_cds, !altaa & (amx) & (usa) & !(amm)))
s_syn_all = sum(with(vg_positions_df_cds, !altaa & (amm) & (usa) & (amx)))

#Calculate Z (s_ns / s_syn) / (p_ns / p_syn)
z_amm_amx = (s_ns_amm_amx / s_syn_amm_amx) / (p_ns_amm / p_syn_amm)
z_amm_usa = (s_ns_amm_usa / s_syn_amm_usa) / (p_ns_amm / p_syn_amm)
z_amx_amm = (s_ns_amm_amx / s_syn_amm_amx) / (p_ns_amx / p_syn_amx)
z_amx_usa = (s_ns_amx_usa / s_syn_amx_usa) / (p_ns_amx / p_syn_amx)
z_usa_amm = (s_ns_amm_usa / s_syn_amm_usa) / (p_ns_usa / p_syn_usa)
z_usa_amx = (s_ns_amx_usa / s_syn_amx_usa) / (p_ns_usa / p_syn_usa)

#Making venn diagrams for Supplement Figure 8

ns_sets <- list(
  "Dark European\nhoney bee samples" = c(
    paste0("amm_only_", seq_len(p_ns_amm)),
    paste0("amm_amx_",  seq_len(s_ns_amm_amx)),
    paste0("amm_usa_",  seq_len(s_ns_amm_usa)),
    paste0("all_",      seq_len(s_ns_all))
  ),
  "Other European\nsamples" = c(
    paste0("amx_only_", seq_len(p_ns_amx)),
    paste0("amm_amx_",  seq_len(s_ns_amm_amx)),
    paste0("amx_usa_",  seq_len(s_ns_amx_usa)),
    paste0("all_",      seq_len(s_ns_all))
  ),
  "USA samples" = c(
    paste0("usa_only_", seq_len(p_ns_usa)),
    paste0("amm_usa_",  seq_len(s_ns_amm_usa)),
    paste0("amx_usa_",  seq_len(s_ns_amx_usa)),
    paste0("all_",      seq_len(s_ns_all))
  )
)

syn_sets <- list(
  "Dark European\nhoney bee samples" = c(
    paste0("amm_only_", seq_len(p_syn_amm)),
    paste0("amm_amx_",  seq_len(s_syn_amm_amx)),
    paste0("amm_usa_",  seq_len(s_syn_amm_usa)),
    paste0("all_",      seq_len(s_syn_all))
  ),
  "Other European\nsamples" = c(
    paste0("amx_only_", seq_len(p_syn_amx)),
    paste0("amm_amx_",  seq_len(s_syn_amm_amx)),
    paste0("amx_usa_",  seq_len(s_syn_amx_usa)),
    paste0("all_",      seq_len(s_syn_all))
  ),
  "USA samples" = c(
    paste0("usa_only_", seq_len(p_syn_usa)),
    paste0("amm_usa_",  seq_len(s_syn_amm_usa)),
    paste0("amx_usa_",  seq_len(s_syn_amx_usa)),
    paste0("all_",      seq_len(s_syn_all))
  )
)


venn_edge_cols <- c("1" = "#e57373","2" = "#7ea6c9","3" = "#8cc084")
region_fill <- c("1" = "#e57373","2" = "#7ea6c9","3"= "#8cc084", "1/2" = "#b08db3", "1/3" = "#b9a97b","2/3" = "#85b8a6", "1/2/3" = "#b3b3b3")

venn_ns <- process_data(Venn(ns_sets))
venn_syn <- process_data(Venn(syn_sets))

S8A <- ggplot() + geom_polygon(data = venn_regionedge(venn_ns), aes(x = X, y = Y, group = id, fill = id), color = NA, alpha = 0.5) +
  geom_path(data = venn_setedge(venn_ns), aes(x = X, y = Y, group = id, color = id), linewidth = 1) +
  geom_text(data = venn_regionlabel(venn_ns), aes(x = X, y = Y, label = count),size = 5) +
  scale_fill_manual(values = region_fill) +
  scale_color_manual(values = venn_edge_cols) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0, size = 15)) +
  ggtitle("A) Non-synonymous SNPs")
 

S8B <- ggplot() + geom_polygon(data = venn_regionedge(venn_syn), aes(x = X, y = Y, group = id, fill = id), color = NA, alpha = 0.5) +
  geom_path(data = venn_setedge(venn_syn), aes(x = X, y = Y, group = id, color = id), linewidth = 1) +
  geom_text(data = venn_regionlabel(venn_syn), aes(x = X, y = Y, label = count),size = 5) +
  scale_fill_manual(values = region_fill) +
  scale_color_manual(values = venn_edge_cols) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0, size = 15)) +
  ggtitle("B) Synonymous SNPs")

S8B

z_df <- data.frame(
  comparison = c(
    "USA samples vs.\nOther European samples",
    "USA samples vs.\nEuropean Dark honey bee samples",
    "Other European samples vs.\nUSA samples",
    "Other European samples vs.\nEuropean Dark honey bee samples",
    "European Dark honey bee samples vs.\nUSA samples",
    "European Dark honey bee samples vs.\nOther European samples"
  ),
  z = c(z_usa_amx,z_usa_amm,z_amx_usa,z_amx_amm,z_amm_usa,z_amm_amx),
  focal = c("usa", "usa", "amx", "amx", "amm", "amm"))

z_df$focal <- factor(z_df$focal, levels = c("usa", "amx", "amm"))

bar_cols <- c(usa = "#8cc084",amx = "#7ea6c9",amm = "#e57373")

S8C <- ggplot(z_df, aes(x = z, y = comparison, fill = focal)) +
  geom_col() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = bar_cols) +
  labs(title = "C) Z-scores",x = "Z-score",y = NULL) +
  theme_classic() +
  theme(legend.position = "none",plot.title = element_text(hjust = 0, size = 15))

layout4 <- "
AB
AB
CC
"
FigS8 <- S8A + S8B + S8C + plot_layout(design = layout4, heights = c(1, 0.8))
ggsave(filename = here("..", "figures", "FigureS8.png"), plot = FigS8, width = 12, height = 7, dpi = 300)
