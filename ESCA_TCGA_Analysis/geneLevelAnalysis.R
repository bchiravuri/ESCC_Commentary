library(cancereffectsizeR)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)

tcga_maf_file <- "TCGA-ESCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "ESCA", filename = tcga_maf_file)
}

# read clinical data & filter to ESCC (8070/3)
tcga_clinical <- fread("clinical.tsv", sep = "\t", header = TRUE)
tcga_clinical_escc <- tcga_clinical[diagnoses.morphology == "8070/3"]
escc_ids <- unique(tcga_clinical_escc$cases.submitter_id)

# preload MAF
tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
tcga_maf_escc <- tcga_maf[Unique_Patient_Identifier %in% escc_ids]


cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa, maf = tcga_maf_escc, maf_name = "TCGA_ESCC")

# prepare clinical for CESAnalysis
setnames(tcga_clinical_escc, "cases.submitter_id", "Unique_Patient_Identifier")
tcga_clinical_escc <- unique(tcga_clinical_escc, by = "Unique_Patient_Identifier")
cesa <- load_sample_data(cesa, tcga_clinical_escc)

sig_exc <- suggest_cosmic_signature_exclusions(
  cancer_type = "eso-SCC",
  treatment_naive = TRUE
)

# trinucleotide mutation rates
cesa <- trinuc_mutation_rates(
  cesa,
  signature_set        = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = sig_exc
)

# gene mutation rates dNdScv
cesa <- gene_mutation_rates(
  cesa,
  covariates = ces.refset.hg38$covariates$esophagus
)

# define genes of interest & create compound variant table
# intersection of coverage (exome+ & targeted)
cov_all <- c(cesa$coverage_ranges$exome, cesa$coverage_ranges$targeted)
cov_all <- cov_all[names(cov_all) != "exome"]
cov_all <- Reduce(intersect, cov_all)

# list genes to prioritize
genes_of_interest <- c(
  "TP53","NOTCH1","NOTCH2","ERBB4","NFE2L2","PIK3CA", "ARID1A","FAT1",
  "EGFR","ERBB2","FBXW7","FGFR3","RB1","SMAD4","SOX2", "FAM135B", "CSMD3", "EP300"
)

# select variants in those genes & covered regions
variants_all <- select_variants(cesa, genes = genes_of_interest, gr = cov_all)

# filter by impact & prevalence
variants_filtered <- variants_all %>%
  filter(
    # tumor suppressors: allow stop-gain/loss or prevalence >1
    (gene %in% c("TP53","NOTCH1","NOTCH2","ERBB4",
                 "ARID1A","FAT1","FBXW7","RB1","SMAD4", "EP300")
     & (maf_prevalence > 0 |
          (aa_ref != "STOP" & aa_alt == "STOP") |
          (aa_ref == "STOP" & aa_alt != "STOP"))
    )
    |
      # oncogenes: require prevalence >1
      (gene %in% c("NFE2L2","PIK3CA","EGFR","ERBB2","FGFR3","SOX2", "FAM135B", "CSMD3")
       & maf_prevalence > 0)
  )

# exclude genes with low aggregate prevalence (<0) just to define high_prev
high_prev <- variants_filtered %>%
  group_by(gene) %>%
  summarize(total_prev = sum(maf_prevalence)) %>%
  filter(total_prev > 0) %>%
  pull(gene)

variants_for_compound <- variants_filtered %>%
  filter(gene %in% high_prev)

# compound variant table by gene
compound_variants <- define_compound_variants(
  cesa = cesa,
  variant_table = variants_for_compound,
  by = "gene",
  merge_distance = Inf
)

# gene-level CES analysis with compound variants
cesa <- ces_variant(
  cesa = cesa,
  variants = compound_variants,
  run_name = "ESCC_gene_level"
)


# plot of top 20 genes by selection strength
plot_effects(
  cesa$selection$ESCC_gene_level,
  group_by = "variant_name",
  topn = 20,
  label_individual_variants = FALSE
)

# top 10 genes
gene_effects <- cesa$selection$ESCC_gene_level
top_genes <- head(gene_effects[order(-selection_intensity)], 10)
print(top_genes[, .(variant_name, selection_intensity)])

save_cesa(cesa = cesa, file = "escc_gene_level_analysis.rds")


# to clear 
cesa <- clear_effect_output(cesa)

