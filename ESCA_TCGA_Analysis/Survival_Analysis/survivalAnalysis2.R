library(tidyverse)
library(data.table)
library(stringr)
library(survival)
library(survminer)
library(cancereffectsizeR)
library(TCGAbiolinks)
library(UCSCXenaTools)
setwd("C:/Users/bodhi/ESCC_Commentary/ESCC_Commentary/ESCA_TCGA_Analysis/Survival_Analysis")
dir.create("survival_outputs", showWarnings = FALSE)

survivaldata <- read.table("survival_ESCA_survival.txt", header = TRUE, sep = "\t")
colnames(survivaldata)[colnames(survivaldata) == "OS"] <- "vital_status"
colnames(survivaldata)[colnames(survivaldata) == "OS.time"] <- "time"
colnames(survivaldata)[colnames(survivaldata) == "X_PATIENT"] <- "Unique_Patient_Identifier"


# load MAF if not already downloaded
tcga_maf_file <- "TCGA-ESCA-SA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "ESCA", filename = tcga_maf_file)
}

# read clinical data
tcga_clinical <- fread("clinical.tsv", sep = "\t", header = TRUE)

# filter for ESCC histology
tcga_clinical_escc <- tcga_clinical[cases.disease_type == "Squamous Cell Neoplasms"]
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
ecDNA_genes_full     <- c("MMP12","PVT1","COX6C","SRSF1","MYC","BIRC2","BIRC3","YAP1")
my_genes             <- c("TP53","NOTCH1","PIK3CA")
target_genes         <- c(ecDNA_genes_full, my_genes)

# build symbol to Ensembl mapping
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

# subset expression to ESCC patients
expr_long <- expr_long[Unique_Patient_Identifier %in% tcga_clinical_escc$Unique_Patient_Identifier]

expr_long <- expr_long %>% filter(Unique_Patient_Identifier %in% escc_ids)

# merge the two data frames
survivaldata <- as.data.table(survivaldata)
days_to_months <- function(x) as.numeric(x) / 30.4375
survivaldata[, time_months := time/30.44]
small_survivaldata <- survivaldata %>% select(Unique_Patient_Identifier, time_months, vital_status) %>% distinct()
combined <- inner_join(expr_long, small_survivaldata, by = "Unique_Patient_Identifier")
combined_combined <- combined %>% left_join(map_final, by = "ensembl_base")



# loop for genes of interest
results_list <- list()
for (gene in target_genes) {
  message(gene)
  gene_df <- combined_combined %>% filter(gene_symbol == gene)
  if (nrow(gene_df) == 0) {
    warning("No data found for ", gene)
    next
  }
# duplicates
  duplicates <- gene_df %>% count(Unique_Patient_Identifier) %>% arrange(desc(n)) %>% filter(n == 2) %>% pull(Unique_Patient_Identifier)
  
  averages <- gene_df %>% filter(Unique_Patient_Identifier %in% duplicates) %>% select(Unique_Patient_Identifier, expr, ensembl_base) %>%
    pivot_wider(names_from = Unique_Patient_Identifier, values_from = expr, values_fn = mean) %>%
    pivot_longer(cols = !ensembl_base, names_to = "Unique_Patient_Identifier", values_to = "expr")
  
  clean_df <- gene_df %>% mutate(expr = ifelse(Unique_Patient_Identifier %in% duplicates, averages$expr[Unique_Patient_Identifier == Unique_Patient_Identifier], expr)) %>%
    distinct() %>% filter(!is.na(expr))
# split up the data into high and low
  res.cut <- surv_cutpoint(clean_df, time = "time_months", event = "vital_status",
    variables = ("expr"))
  summary(res.cut)
  
  plot(res.cut, "expr", palette = "npg")
  
  res.cat <- surv_categorize(res.cut)
  head(res.cat)
# save plot with number at risk table
  ggsurv <- ggsurvplot(
    fit, data = res.cat,
    risk.table = TRUE,
    conf.int = TRUE, pval = TRUE,
    tables.height = 0.25,             # controls bottom panel height
    tables.theme = theme_cleantable()  # cleaner risk table
  )
  # combine main plot + risk table, then save as ONE image
  combined_plot <- ggarrange(
    ggsurv$plot,
    ggsurv$table + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5)),
    ncol = 1, heights = c(3, 1)
  )
  out_png <- file.path("survival_outputs", paste0("survival_", gene, ".png"))
  ggsave(out_png, combined_plot, width = 7, height = 6, dpi = 300)
  res_tbl <- Filter(Negate(is.null), results_list)
  results_df <- if (length(res_tbl)) data.table::rbindlist(res_tbl, fill = TRUE, use.names = TRUE) else data.table::data.table()
  if (nrow(results_df)) {
    out_csv <- file.path("survival_outputs","survival_results.csv")
    data.table::fwrite(results_df, out_csv)
    print(results_df[order(tag, p_logrank), .(tag, gene, ensembl, n, HR_high, p_logrank, HR_high_adj, p_adj)])
  }
}




 
# COX6C_df <- combined_combined %>% filter(gene_symbol == "COX6C")

# duplicates <- COX6C_df %>% count(Unique_Patient_Identifier) %>% arrange(desc(n)) %>% filter(n == 2) %>% pull(Unique_Patient_Identifier)
# 
# averages <- COX6C_df %>% filter(Unique_Patient_Identifier %in% duplicates) %>% select(Unique_Patient_Identifier, expr, ensembl_base) %>%
#     pivot_wider(names_from = Unique_Patient_Identifier, values_from = expr, values_fn = mean) %>%
#     pivot_longer(cols = !ensembl_base, names_to = "Unique_Patient_Identifier", values_to = "expr")
# 
# clean_df <- COX6C_df %>% mutate(expr = ifelse(Unique_Patient_Identifier %in% duplicates, averages$expr[Unique_Patient_Identifier == Unique_Patient_Identifier], expr)) %>%
#   distinct() %>% filter(!is.na(expr))


# # split up the data into high and low
# res.cut <- surv_cutpoint(clean_df, time = "time_months", event = "vital_status",
#   variables = ("expr"))
# summary(res.cut)
# 
# plot(res.cut, "expr", palette = "npg")
# 
# res.cat <- surv_categorize(res.cut)
# head(res.cat)

# Fit survival curves and visualize
# fit <- survfit(Surv(time_months, vital_status) ~expr, data = res.cat)
# ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval = TRUE)


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
