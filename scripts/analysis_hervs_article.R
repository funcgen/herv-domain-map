###############################################################################
#                          HERV CONSERVATION ANALYSIS                         #
#                                                                             #
#  Description  : Script to analyze domain conservation in HERV ORFs          #
#                 and generate figures for the Results section of the         #
#                 manuscript                                                  #
#                                                                             #
#  Manuscript   : "A Comprehensive Annotation of Conserved                    # 
#                  Protein Domains in Human Endogenous Retroviruses"          #
#                                                                             #
#                                                                             #
#  Sections     : - Summary statistics of ORFs and domain hits                #
#                 - Domain conservation across gene classes (GAG, POL, ENV)   #
#                 - Chromosomal density of ORFs                               #
#                 - Figure generation for publication                         #
#                                                                             #
#  Author       : Tomàs Montserrat-Ayuso                                      #
#  Institution  : CNAG                                                        #
#  Date         : 16/07/2025                                                  #
#                                                                             #
###############################################################################

# --- Load Libraries ---
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(openxlsx)
library(patchwork)



### 2.1 HIGH PREVALENCE OF CONSERVED DOMAINS IN HERV ORFs ###

# --- File Paths ---
orf_fasta <- "../results/ERV_GyDB_v4_orfs_fixed.fasta"
ann_file  <- "../results/ERV_GyDB_v4_filtered.tsv"
out_plot  <- "figures/distribution_gag_pol_env_orfs.png"

# --- Load ORFs ---
orf_seqs <- readAAStringSet(orf_fasta)
orf_info <- tibble(
  name   = names(orf_seqs),
  length = width(orf_seqs),
  chr    = str_extract(names(orf_seqs), "chr[0-9XYM]+")
)

cat("Total ORFs analyzed:", length(orf_seqs), "\n\n")
print(summary(orf_info$length))

# --- Summary per Chromosome ---
chr_summary <- orf_info %>%
  group_by(chr) %>%
  summarise(
    n_orfs     = n(),
    min_len    = min(length),
    max_len    = max(length),
    median_len = median(length),
    mean_len   = mean(length),
    .groups    = "drop"
  )

# Add chromosome length information
hg38_chr_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
hg38_autosomes <- hg38_chr_lengths[paste0("chr", c(1:22, "X", "Y"))]
chr_df <- data.frame(chr = names(hg38_autosomes), chr_length = as.numeric(hg38_autosomes))

# Merge and compute density
chr_density <- chr_summary %>%
  left_join(chr_df, by = "chr") %>%
  mutate(
    orfs_per_Mb = n_orfs / (chr_length / 1e6),
    bp_per_orf  = chr_length / n_orfs
  )

print(as.data.frame(chr_density))

# --- Load Annotation Results ---
ann_table <- read.table(ann_file, header = TRUE, sep = "\t")

n_annotated_orfs <- length(unique(ann_table$Query))
cat("Number of ORFs with at least one domain hit:", n_annotated_orfs, "\n\n")

# Simplify annotation table
ann_simple <- ann_table %>%
  select(Query, Gene, Domain, Coverage, Score) %>%
  mutate(
    Subfamily   = sub("-int.*", "", Query),
    DomainType  = sub("_.*", "", Domain)
  )

# --- Domain Coverage Threshold ---
coverage_threshold <- 0.4
ann_conserved <- ann_simple %>%
  filter(Coverage >= coverage_threshold)

cat("Number of domain hits with coverage ≥", coverage_threshold, ":", nrow(ann_conserved), "\n\n")

# --- Domain Tally by Gene ---
domain_conservation_counts <- ann_conserved %>%
  count(Gene) %>%
  mutate(
    Gene = factor(Gene, levels = c("gag", "pol", "env", "accessory"))
  )

print(domain_conservation_counts)

# --- Plot: Conserved Domains by Gene Class ---
domain_plot <- ggplot(domain_conservation_counts, aes(x = Gene, y = n)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "#56B4E9") +
  labs(
    title = "Conserved Retroviral Domains in Internal HERV ORFs",
    x = "Gene Class",
    y = "Number of ORFs"
  ) +
  theme_minimal(base_size = 22)
ggsave(
  filename = out_plot,
  plot = domain_plot,
  width = 12, height = 9, bg = "white"
)





### 2.2 SUBSTANTIAL FRACTION OF DOMAINS SHOW HIGH SEQUENCE CONSERVATION ###

# --- Coverage Distribution Plot ---

# Flag domains with full-length coverage
ann_simple <- ann_simple %>%
  mutate(HighCoverage = Coverage > 0.95)

# Plot domain coverage by gene class
coverage_plot <- ggplot(ann_simple, aes(x = Gene, y = Coverage)) +
  geom_violin(fill = "gray90", color = NA, scale = "width", trim = TRUE, alpha = 0.3) +
  geom_jitter(
    aes(color = HighCoverage),
    width = 0.25, height = 0.01,
    size = 1.5, alpha = 0.5
  ) +
  geom_hline(yintercept = 0.40, linetype = "dotdash", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dotdash", color = "grey30", linewidth = 0.6) +
  annotate("text", x = 4.6, y = 0.41, label = "Min. coverage\n(cutoff = 0.4)", 
           hjust = 0, size = 3.5, color = "red") +
  annotate("text", x = 4.6, y = 0.96, label = "Full-length", 
           hjust = 0, size = 4, color = "grey30") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title    = "Domain Coverage in Internal HERV ORFs",
    subtitle = "Thresholds correspond to conservation categories",
    x        = "Gene Class",
    y        = "Domain Coverage",
    color    = "High Coverage"
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#E69F00")) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank()
  )

# Display and save plot
print(coverage_plot)
ggsave(
  filename = "coverage_distribution_gag_pol_env_orfs.png",
  path     = "figures",
  plot     = coverage_plot, 
  width    = 12, 
  height   = 9,
  bg       = "white"
)

# --- Deeper Analysis of Full-Length Domains ---

# Filter full-length domains
ann_full_length <- ann_simple %>%
  filter(Coverage > 0.95)

# Load motif analysis results
motif_hits <- read.table("../results/ERV_GyDB_analysis_conserved_residues.tsv", 
                         header = TRUE, sep = "\t")

# Standardize Protein/ORF identifiers
motif_hits <- motif_hits %>%
  mutate(
    ORFs        = sub("^([^|]+\\|[^|]+)\\|.*$", "\\1", Protein),
    domain_type = sub("^([^|]*\\|){3}([^_|]+).*", "\\2", Protein)
  ) %>%
  filter(ORFs %in% ann_full_length$Query)

# Summary: number of unique full-length domains with motif data
cat("Unique full-length domains with motif data:", length(unique(motif_hits$Protein)), "\n")

# --- Subset by Domain Type ---
motif_RNaseH <- filter(motif_hits, domain_type == "RNaseH")
motif_AP     <- filter(motif_hits, domain_type == "AP")
motif_DUT    <- filter(motif_hits, domain_type == "DUT")
motif_ENV    <- filter(motif_hits, domain_type == "ENV")
motif_RT     <- filter(motif_hits, domain_type == "RT")

cat("RNaseH domains:", length(unique(motif_RNaseH$Protein)), "\n")
cat("AP domains    :", length(unique(motif_AP$Protein)), "\n")
cat("DUT domains   :", length(unique(motif_DUT$Protein)), "\n")
cat("ENV domains   :", length(unique(motif_ENV$Protein)), "\n")

# --- RNaseH: Domains with Both Key Motifs ---
rnaseh_motifs <- motif_RNaseH %>%
  filter(Site.Description %in% c("RNA/DNA hybrid binding site", "active site")) %>%
  distinct(Protein, Site.Description) %>%
  group_by(Protein) %>%
  summarise(n = n()) %>%
  filter(n == 2) %>%
  pull(Protein)

rnaseh_conserved <- filter(motif_RNaseH, Protein %in% rnaseh_motifs)
cat("RNaseH domains with both conserved motifs:", length(unique(rnaseh_conserved$Protein)), "\n")

# --- ENV: Immunosuppressive Regions ---
motif_ENV_immune <- motif_ENV %>%
  filter(Site.Description == "immunosuppressive region")

env_immunosuppressive <- filter(motif_ENV, Protein %in% motif_ENV_immune$Protein)
cat("ENV domains with immunosuppressive region:", length(unique(env_immunosuppressive$Protein)), "\n")




# 
# ann_simple_high_conserved$Query_DomainType = paste0(ann_simple_high_conserved$Query, "_", ann_simple_high_conserved$DomainType)
# ann_simple_high_conserved_env = ann_simple_high_conserved %>% 
#   filter(DomainType == "ENV")
# 
# motif_hits$Query_DomainType = paste0(motif_hits$ORFs, "_", motif_hits$domain_type)
# motif_hits_env = motif_hits %>% 
#   filter(domain_type == "ENV")
# 
# 
# 
# length(which(ann_simple_high_conserved_env$Query_DomainType %in% motif_hits_env$Query_DomainType))
# length(which(motif_hits_env$Query_DomainType %in% ann_simple_high_conserved_env$Query_DomainType))
# 
# 
# 
# View(ann_simple_high_conserved_env[which(!ann_simple_high_conserved_env$Query_DomainType %in% motif_hits_env$Query_DomainType),])










### 2.3 SUBFAMILY-LEVEL PATTERNS OF DOMAIN RETENTION ###

# --- Define Target Subfamilies and Domain Order ---
focus_subfamilies <- c("HERVK", "HERVH", "HERV17", "HERVE")
domain_order <- c("GAG", "RT", "RNaseH", "DUT", "INT", "AP", "ENV", "Accessory")

# --- Prepare Data for Plotting ---
subfamily_plot_df <- ann_simple %>%
  mutate(
    DomainType = recode(DomainType, ORFX = "Accessory"),
    Coverage = as.numeric(Coverage),
    HighCoverage = Coverage > 0.95
  ) %>%
  filter(
    Subfamily %in% focus_subfamilies,
    DomainType %in% domain_order
  ) %>%
  mutate(
    DomainType = factor(DomainType, levels = domain_order),
    Subfamily = factor(Subfamily, levels = focus_subfamilies)
  ) %>%
  distinct(Query, DomainType, .keep_all = TRUE)

# --- Plot: Domain Coverage per Subfamily ---
subfamily_plot <- ggplot(subfamily_plot_df, aes(x = Subfamily, y = Coverage)) +
  geom_violin(fill = "gray90", color = NA, scale = "width", trim = TRUE, alpha = 0.3) +
  geom_jitter(aes(color = HighCoverage), width = 0.25, height = 0.01, size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 0.4, linetype = "dotdash", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dotdash", color = "grey40", linewidth = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "#E69F00"),
    labels = c("≤ 0.95", "> 0.95")
  ) +
  facet_wrap(~ DomainType, ncol = 4, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Domain Coverage in Key HERV Subfamilies",
    subtitle = "Points with coverage > 0.95 are highlighted",
    x = "HERV Subfamily",
    y = "Coverage Score",
    color = "Coverage"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# Display and save plot
print(subfamily_plot)
ggsave(
  filename = "domain_coverage_subfamilies.png",
  path     = "figures",
  plot     = subfamily_plot,
  width    = 12,
  height   = 9,
  bg       = "white"
)

# --- Generate Domain Retention Tables by Subfamily ---

# Table: All annotated domains
domain_table_all <- table(ann_simple$Subfamily, ann_simple$DomainType)

# Table: Conserved domains (coverage ≥ 0.4)
domain_table_conserved <- table(ann_conserved$Subfamily, ann_conserved$DomainType)

# Table: Full-length domains (coverage > 0.95)
domain_table_full_length <- table(ann_full_length$Subfamily, ann_full_length$DomainType)

# Convert to data frames
domain_df_all <- as.data.frame.matrix(domain_table_all) %>%
  tibble::rownames_to_column("Subfamily")

domain_df_conserved <- as.data.frame.matrix(domain_table_conserved) %>%
  tibble::rownames_to_column("Subfamily")

domain_df_full_length <- as.data.frame.matrix(domain_table_full_length) %>%
  tibble::rownames_to_column("Subfamily")

# --- Write Tables to Excel ---
write.xlsx(domain_df_all,          file = "HERV_domain_counts_all.xlsx",              rowNames = FALSE)
write.xlsx(domain_df_conserved,    file = "HERV_domain_counts_all_conserved.xlsx",    rowNames = FALSE)
write.xlsx(domain_df_full_length,  file = "HERV_domain_counts_all_full_length.xlsx",  rowNames = FALSE)






### 2.4 DOMAIN CO-OCCURRENCE AND POTENTIAL FUNCTIONALITY ###

# --- Extract Genomic Locus from ORF Name ---
ann_simple <- ann_simple %>%
  mutate(Locus = sub("_[^_|]+\\|.*$", "", Query))

# --- Function to Identify GAG + POL + ENV Loci ---
identify_top_loci <- function(min_coverage = 0.4) {
  domain_presence <- ann_simple %>%
    filter(Coverage >= min_coverage) %>%
    group_by(Locus) %>%
    summarise(
      Has_GAG   = any(Gene == "gag"),
      Has_POL   = any(Gene == "pol"),
      Has_ENV   = any(Gene == "env"),
      N_domains = n(),
      .groups   = "drop"
    )
  
  # Keep loci with all three major domains
  trios <- domain_presence %>%
    filter(Has_GAG, Has_POL, Has_ENV)
  
  # Return expanded table with all domain info for selected loci
  top_loci <- ann_simple %>%
    filter(Locus %in% trios$Locus) %>%
    arrange(Locus, DomainType, desc(Coverage)) %>%
    select(Locus, Gene, DomainType, Domain, Coverage, Subfamily)
  
  return(top_loci)
}

# --- Moderate Stringency: ≥ 0.4 Coverage ---
top_loci_0.4 <- identify_top_loci(min_coverage = 0.4)
print(top_loci_0.4)

write.xlsx(top_loci_0.4, file = "results/loci_with_trios.xlsx", rowNames = FALSE)

top_loci_0.4_wide <- top_loci_0.4 %>%
  select(Locus, DomainType, Coverage, Subfamily) %>%
  tidyr::pivot_wider(
    names_from  = DomainType,
    values_from = Coverage
  ) %>%
  mutate(
    DomainCount = rowSums(!is.na(select(., -Locus, -Subfamily)))
  )

write.xlsx(top_loci_0.4_wide, file = "results/loci_wide_with_trios.xlsx", rowNames = FALSE)


# --- High Stringency: ≥ 0.8 Coverage ---
top_loci_0.8 <- identify_top_loci(min_coverage = 0.8)
print(top_loci_0.8)

write.xlsx(top_loci_0.8, file = "results/loci_with_trios_top_quality.xlsx", rowNames = FALSE)

# --- Wide Table Format for Manuscript or Supplement ---
top_loci_wide <- top_loci_0.8 %>%
  select(Locus, DomainType, Coverage, Subfamily) %>%
  tidyr::pivot_wider(
    names_from  = DomainType,
    values_from = Coverage
  ) %>%
  mutate(
    DomainCount = rowSums(!is.na(select(., -Locus, -Subfamily)))
  )

write.xlsx(top_loci_wide, file = "results/top_loci_wide_trios.xlsx", rowNames = FALSE)






### FIGURE ###

fig_a <- domain_plot
fig_b <- coverage_plot
fig_c <- subfamily_plot

# Combine plots with patchwork
combined_plot <- (fig_a + fig_b) / fig_c +
  plot_layout(heights = c(1, 1.5)) +  # Make bottom row taller
  plot_annotation(tag_levels = 'a')  # Adds 'a', 'b', 'c' labels automatically

# Save to file
ggsave("figures/Fig2.png", combined_plot,
       width = 28, height = 18, dpi = 300, bg = "white")







