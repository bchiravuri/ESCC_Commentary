library(tidyverse)
library(data.table)
library(stringr)      
library(survival)     
library(survminer)  
library(cancereffectsizeR)
setwd("Survival_Analysis")


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


# Get list of ESCC patients with ≥1 TP53 mutation
tp53_patients <- unique(
  substr(
    maf_escc$Tumor_Sample_Barcode[maf_escc$Hugo_Symbol == "TP53"],
    1, 12
  )
)

# Build a data.frame with mutation status
mut_status_df <- escc_clin %>%
  mutate(
    tp53_mut = ifelse(patient_barcode %in% tp53_patients, 1, 0)
  ) %>%
  # compute survival time & event indicator
  mutate(
    time   = ifelse(!is.na(days_to_death),
                    days_to_death,
                    days_to_last_follow_up),
    status = ifelse(vital_status == "Dead", 1, 0)
  )

fit_km <- survfit(Surv(time, status) ~ tp53_mut, data = mut_status_df)

ggsurvplot(
  fit_km,
  data       = mut_status_df,
  pval       = TRUE,
  risk.table = TRUE,
  legend.labs= c("TP53 WT", "TP53 Mut"),
  xlab       = "Time (days)",
  ylab       = "Survival Probability"
)

# If you have age and pathologic_stage in your clinical, include them here:
# mut_status_df <- mut_status_df %>%
#   left_join(clinical_all %>% select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis, pathologic_stage),
#             by = c("patient_barcode" = "bcr_patient_barcode"))

cox_mod <- coxph(
  Surv(time, status) ~ tp53_mut,
  data = mut_status_df
)
summary(cox_mod)
