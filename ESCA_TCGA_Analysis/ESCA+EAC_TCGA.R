library(cancereffectsizeR)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

tcga_maf_file <- "TCGA-ESCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "ESCA", filename = tcga_maf_file)
}

# reading the clinical data
tcga_clinical <- fread("clinical.tsv", sep = "\t", header = TRUE)

# Histological codes
# 8070/3 = Squamous cell carcinoma
# 8140/3 = Adenocarcinoma
tcga_clinical_escc <- tcga_clinical[diagnoses.morphology == "8070/3"]
tcga_clinical_eac <- tcga_clinical[diagnoses.morphology == "8140/3"]

# Combine both IDs
combined_ids <- unique(c(tcga_clinical_escc$cases.submitter_id, tcga_clinical_eac$cases.submitter_id))

# preload maf
tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

# Filter the MAF
tcga_maf_combined <- tcga_maf[Unique_Patient_Identifier %in% combined_ids]
nrow(tcga_maf_combined)

# CESAnalysis
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa, maf = tcga_maf_combined, maf_name = "TCGA_ESCA")

# Merge clinical data and remove duplicates
tcga_clinical_combined <- rbind(tcga_clinical_escc, tcga_clinical_eac)
setnames(tcga_clinical_combined, "cases.submitter_id", "Unique_Patient_Identifier")
tcga_clinical_combined <- unique(tcga_clinical_combined, by = "Unique_Patient_Identifier")
cesa <- load_sample_data(cesa, tcga_clinical_combined)

# Use signature exclusions appropriate for a mixed cohort — you can adjust based on your focus
excl_scc <- suggest_cosmic_signature_exclusions("eso-SCC", treatment_naive = TRUE)
excl_eac <- suggest_cosmic_signature_exclusions("eso-AdenoCA", treatment_naive = TRUE)

signature_exclusions <- union(excl_scc, excl_eac)

cesa <- trinuc_mutation_rates(cesa,
  signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)

cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$esophagus)

# view significant genes
dndscv_results <- cesa$dNdScv_results[[1]]
sig_genes <- dndscv_results[qallsubs_cv < .05][order(qallsubs_cv)][1:10]
print(sig_genes)

# top variants + samples
variants_to_check <- cesa$variants[order(-maf_prevalence), variant_id][1:3]
samples_to_check <- cesa$samples$Unique
baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check, samples = samples_to_check)


# cancer effect size analysis
cesa <- ces_variant(cesa = cesa, run_name = "ESCA_combined_selection")

# Plot overall effects
plot_effects(cesa$selection$ESCA_combined_selection)

# Plot effects by gene.....
plot_effects(
  cesa$selection$ESCA_combined_selection,
  group_by = "gene", topn = 20,
  label_individual_variants = FALSE
)

# Count most frequently mutated genes
mut_counts <- cesa$variants %>% count(gene, name = "mutations")
mut_counts %>% arrange(desc(mutations)) %>% head(n = 20)

# Normalize mutation rates by gene length
gr_genes <- readRDS("C:/Users/bodhi/AppData/Local/R/win-library/4.5/ces.refset.hg38/refset/gr_genes.rds")
gene_lengths2 <- as.data.table(gr_genes)[, .(length = sum(width(reduce(GRanges(seqnames, IRanges(start, end)))))), by = gene]

write.csv(gene_lengths2, "gene_lengths2.csv")

normalized <- mut_counts %>%
  inner_join(gene_lengths2, by = "gene") %>%
  mutate(
    length_kb = length / 1e3,
    muts_per_kb = mutations / length_kb
  ) %>%
  arrange(desc(muts_per_kb))

save_cesa(cesa = cesa, file = "esca_combined_analysis.rds")
