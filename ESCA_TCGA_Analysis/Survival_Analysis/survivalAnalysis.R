library(tidyverse)
library(data.table)
library(stringr)
library(survival)
library(survminer)
library(cancereffectsizeR)
library(TCGAbiolinks)
setwd("C:/Users/bodhi/ESCC_Commentary/ESCC_Commentary/ESCA_TCGA_Analysis/Survival_Analysis")
dir.create("survival_outputs", showWarnings = FALSE)

# load MAF if not already downloaded
tcga_maf_file <- "TCGA-ESCA-SA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "ESCA", filename = tcga_maf_file)
}

# read clinical data
tcga_clinical <- fread("clinical.tsv", sep = "\t", header = TRUE)

# filter for ESCC histology
tcga_clinical_escc <- tcga_clinical[diagnoses.morphology == "8070/3"]
escc_ids <- unique(tcga_clinical_escc$cases.submitter_id)

# load MAF
tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
tcga_maf_escc <- tcga_maf[Unique_Patient_Identifier %in% escc_ids]

# create CESAnalysis object and load MAF
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa, maf = tcga_maf_escc, maf_name = "TCGA_ESCC")
setnames(tcga_clinical_escc, "cases.submitter_id", "Unique_Patient_Identifier")

# remove duplicate patients
tcga_clinical_escc <- unique(tcga_clinical_escc, by = "Unique_Patient_Identifier")
cesa <- load_sample_data(cesa, tcga_clinical_escc)

# rename columns for survival calculation
setnames(tcga_clinical_escc,
         old = c("demographic.vital_status", 
                 "demographic.days_to_death", 
                 "diagnoses.days_to_last_follow_up"),
         new = c("vital_status", "days_to_death", "days_to_last_follow_up"))

# create OS time and status
tcga_clinical_escc[, time := fifelse(!is.na(days_to_death),
                                     as.numeric(days_to_death),
                                     as.numeric(days_to_last_follow_up))]
tcga_clinical_escc[, status := fifelse(vital_status == "Dead", 1L, 0L)]

# helper functions for survival
days_to_months <- function(x) as.numeric(x) / 30.4375
build_surv_fields <- function(clin) { ... } # unchanged

# function to run survival for one gene
run_expression_survival <- function(expr_dt, clin_dt, gene_symbol, ...) { ... } # unchanged

# load expression matrix
expr_wide <- data.table::fread("TCGA-ESCA.star_fpkm.tsv")
if (!("gene_raw" %in% names(expr_wide))) {
  first_col <- names(expr_wide)[1]; data.table::setnames(expr_wide, first_col, "gene_raw")
}
expr_wide[, ensembl_base := sub("\\..*$", "", gene_raw)]

# reshape to long
expr_long <- melt(expr_wide, id.vars = c("gene_raw","ensembl_base"),
                  variable.name = "barcode", value.name  = "expr",
                  variable.factor = FALSE)
expr_long[, Unique_Patient_Identifier := substr(barcode, 1, 12)]
suppressWarnings(expr_long[, expr := as.numeric(expr)])
expr_long <- expr_long[, .(ensembl_base, Unique_Patient_Identifier, expr)]

# define target genes
ecDNA_genes_headline <- c("COX6C","AZIN1-AS1","MMP12","PVT1")
ecDNA_genes_full     <- c("AZIN1-AS1","MMP12","PVT1","COX6C","SRSF1","MYC","BIRC2","BIRC3","YAP1")
my_genes             <- c("TP53","NOTCH1","PIK3CA")
target_genes         <- c(ecDNA_genes_full, my_genes)

# build symbol→Ensembl mapping
map_all <- data.table::data.table()
if (file.exists("gencode.v36.annotation.gtf.gene.probemap")) {
  pm <- data.table::fread("gencode.v36.annotation.gtf.gene.probemap", header = FALSE)
  m_probemap <- unique(data.table::data.table(
    ensembl_base = sub("\\..*$","", as.character(pm[[1]])),
    gene_symbol  = as.character(pm[[2]])
  ))
  map_all <- unique(rbindlist(list(map_all, m_probemap), use.names = TRUE, fill = TRUE))
}
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressWarnings({
    m2 <- AnnotationDbi::select(org.Hs.eg.db, keys = target_genes,
                                keytype = "SYMBOL", columns = c("ENSEMBL","SYMBOL"))
  })
  if (!is.null(m2) && nrow(m2)) {
    m_org <- unique(data.table::data.table(
      ensembl_base = sub("\\..*$","", m2$ENSEMBL),
      gene_symbol  = m2$SYMBOL
    ))
    map_all <- unique(rbindlist(list(map_all, m_org), use.names = TRUE, fill = TRUE))
  }
}
map_all <- map_all[!is.na(ensembl_base) & !is.na(gene_symbol) & gene_symbol %in% target_genes]
present_ens <- unique(expr_long$ensembl_base)
map_final <- unique(map_all[ensembl_base %in% present_ens])
map_final <- map_final[!duplicated(gene_symbol)]

# check mapping availability
avail <- data.table::data.table(
  gene_symbol = target_genes,
  ensembl_base = map_final$ensembl_base[match(target_genes, map_final$gene_symbol)],
  present = !is.na(map_final$ensembl_base[match(target_genes, map_final$gene_symbol)])
)
print(avail)

# ensure OS fields
if (!("time" %in% names(tcga_clinical_escc))) {
  tcga_clinical_escc[, time := fifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_follow_up))]
}
if (!("status" %in% names(tcga_clinical_escc))) {
  tcga_clinical_escc[, status := fifelse(vital_status == "Dead", 1L, 0L)]
}
tcga_clinical_escc[, time_months := time/30.44]

# load DFS from Xena if available
dfs_files <- c("TCGA-ESCA.survival.tsv", "TCGA-ESCA_survival.tsv", "TCGA-ESCA.clinical_and_survival.tsv")
dfs_file  <- dfs_files[file.exists(dfs_files)][1]
if (!is.na(dfs_file)) {
  xena_surv <- data.table::fread(dfs_file)
  barcode_col <- names(xena_surv)[which.max(grepl("^(sample|id|barcode)$", names(xena_surv), ignore.case = TRUE))]
  if (is.na(barcode_col)) barcode_col <- names(xena_surv)[1]
  xena_surv[, Unique_Patient_Identifier := substr(get(barcode_col), 1, 12)]
  dfi_time_col <- names(xena_surv)[which.max(grepl("^dfi(\\.|_|)time$", names(xena_surv), ignore.case = TRUE))]
  dfi_event_col <- names(xena_surv)[which.max(grepl("^dfi$", names(xena_surv), ignore.case = TRUE))]
  keep_cols <- c("Unique_Patient_Identifier", dfi_time_col, dfi_event_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  xena_keep <- unique(xena_surv[, ..keep_cols])
  setnames(xena_keep, old = c(dfi_time_col, dfi_event_col), new = c("DFI_time", "DFI_event"), skip_absent = TRUE)
  tcga_clinical_escc <- merge(tcga_clinical_escc, xena_keep, by = "Unique_Patient_Identifier", all.x = TRUE)
  tcga_clinical_escc[, time_months := fifelse(!is.na(DFI_time), as.numeric(DFI_time), time_months)]
  tcga_clinical_escc[, status := fifelse(!is.na(DFI_event), as.integer(DFI_event), status)]
}

# process covariates
tcga_clinical_escc[, age_raw := suppressWarnings(as.numeric(`demographic.age_at_index`))]
tcga_clinical_escc[, age := ifelse(mean(age_raw, na.rm = TRUE) > 200, age_raw/365.25, age_raw)]
tcga_clinical_escc[, gender := factor(`demographic.gender`)]
stage_col <- if ("diagnoses.ajcc_pathologic_stage" %in% names(tcga_clinical_escc)) "diagnoses.ajcc_pathologic_stage" else
  if ("diagnoses.uicc_pathologic_stage" %in% names(tcga_clinical_escc)) "diagnoses.uicc_pathologic_stage" else NA_character_
if (!is.na(stage_col)) {
  tcga_clinical_escc[, stage_txt := toupper(gsub("[^IVX]", "", get(stage_col)))]
  tcga_clinical_escc[stage_txt == "", stage_txt := NA_character_]
  tcga_clinical_escc[, stage := factor(stage_txt, levels = c("I","II","III","IV"), ordered = TRUE)]
} else { tcga_clinical_escc[, stage := NA] }

# subset expression to ESCC patients
expr_long <- expr_long[Unique_Patient_Identifier %in% tcga_clinical_escc$Unique_Patient_Identifier]

# run survival for each gene
results_list <- list()
for (sym in target_genes) {
  ens <- map_final$ensembl_base[match(sym, map_final$gene_symbol)]
  if (!is.na(ens) && ens %in% present_ens) {
    tag <- if (sym %in% ecDNA_genes_full) "ecDNA" else "mine"
    results_list[[paste0(tag, "_", sym)]] <- run_expression_survival_ens(expr_long, tcga_clinical_escc, ens, sym, tag)
  }
}

# save results
res_tbl <- Filter(Negate(is.null), results_list)
results_df <- if (length(res_tbl)) data.table::rbindlist(res_tbl, fill = TRUE, use.names = TRUE) else data.table::data.table()
if (nrow(results_df)) {
  out_csv <- file.path("survival_outputs","survival_results.csv")
  data.table::fwrite(results_df, out_csv)
  print(results_df[order(tag, p_logrank), .(tag, gene, ensembl, n, HR_high, p_logrank, HR_high_adj, p_adj)])
}

# list output files
print(list.files("survival_outputs", pattern = "png|csv", full.names = TRUE))
