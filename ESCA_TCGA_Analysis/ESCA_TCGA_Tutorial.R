options(timeout = 600)
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR", dependencies = TRUE, repos = BiocManager::repositories())
BiocManager::install(version = "3.21")
options(timeout = 600)
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@*release")

dir.create("ESCA_TCGA_Analysis")
setwd("ESCA_TCGA_Analysis")

library(cancereffectsizeR)
library(data.table)

tcga_maf_file <- "TCGA-ESCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "ESCA", filename = tcga_maf_file)
}

# reading the clinical data
tcga_clinical <- fread("clinical.tsv")
# setnames(tcga_clinical_escc, "cases.submitter_id", "Unique_Patient_Identifier")
tcga_clinical <- fread("clinical.tsv", sep = "\t", header = TRUE)

# Histological code 8070/3 = Squamous cell carcinoma
tcga_clinical_escc <- tcga_clinical[diagnoses.morphology == "8070/3"]
escc_ids <- unique(tcga_clinical_escc$cases.submitter_id)

# preload maf
tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

# Filter the MAF w/ IDs
tcga_maf_escc <- tcga_maf[Unique_Patient_Identifier %in% escc_ids]
# Check it worked......
nrow(tcga_maf_escc)

#CESAnalysis
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa, maf = tcga_maf_escc, maf_name = "TCGA_ESCC")
setnames(tcga_clinical_escc, "cases.submitter_id", "Unique_Patient_Identifier")

# remove duplicates
tcga_clinical_escc <- unique(tcga_clinical_escc, by = "Unique_Patient_Identifier")
cesa <- load_sample_data(cesa, tcga_clinical_escc)




signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "eso-SCC", treatment_naive = TRUE)

cesa <- trinuc_mutation_rates(cesa,
  signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)


cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$esophagus)


# view significant genes
dndscv_results <- cesa$dNdScv_results[[1]]
sig_genes <- dndscv_results[qallsubs_cv < .05][order(qallsubs_cv)][1:10]
print(sig_genes)


# I think this is how to check top variants and a few sample rates....
variants_to_check <- cesa$variants[order(-maf_prevalence), variant_id][1:3]
samples_to_check <- cesa$samples$Unique_Patient_Identifier[1:3]
baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check, samples = samples_to_check)



# cancer effects PLOTS 
cesa <- ces_variant(cesa = cesa, run_name = "ESCC_selection")

plot_effects(cesa$selection$ESCC_selection)
plot_effects(cesa$selection$ESCC_selection,
             group_by = "gene", topn = 10,
             label_individual_variants = FALSE
)

save_cesa(cesa = cesa, file = "escc_analysis.rds")