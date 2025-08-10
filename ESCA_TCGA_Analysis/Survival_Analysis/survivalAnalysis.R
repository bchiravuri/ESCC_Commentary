library(tidyverse)
library(data.table)
library(stringr)      
library(survival)     
library(survminer)  
library(cancereffectsizeR)
library(TCGAbiolinks)
setwd("C:/Users/bodhi/ESCC_Commentary/ESCC_Commentary/ESCA_TCGA_Analysis/Survival_Analysis")
dir.create("survival_outputs", showWarnings = FALSE)


tcga_maf_file <- "TCGA-ESCA-SA.maf.gz"
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

# rename, convenience
setnames(tcga_clinical_escc,
         old = c("demographic.vital_status", 
                 "demographic.days_to_death", 
                 "diagnoses.days_to_last_follow_up"),
         new = c("vital_status", "days_to_death", "days_to_last_follow_up"))


# create time and status things
tcga_clinical_escc[, time := fifelse(
  !is.na(days_to_death),
  as.numeric(days_to_death),
  as.numeric(days_to_last_follow_up)
)]

tcga_clinical_escc[, status := fifelse(
  vital_status == "Dead",
  1L,
  0L
)]

# START
days_to_months <- function(x) as.numeric(x) / 30.4375

# Build survival fields from clinical (you already renamed vital_status / days_to_* earlier)
build_surv_fields <- function(clin) {
  clin <- data.table::copy(clin)
  if (!("time" %in% names(clin))) {
    clin[, time := fifelse(!is.na(days_to_death),
                           as.numeric(days_to_death),
                           as.numeric(days_to_last_follow_up))]
  }
  if (!("status" %in% names(clin))) {
    clin[, status := fifelse(vital_status == "Dead", 1L, 0L)]
  }
  clin[, time_months := days_to_months(time)]
  clin
}

# Core KM/Cox using surv_cutpoint for expression split
run_expression_survival <- function(expr_dt, clin_dt, gene_symbol,
                                    min_prop = 0.1, save_plot = TRUE,
                                    tag = "ecDNA") {
  gdat <- expr_dt[gene == gene_symbol]
  if (nrow(gdat) == 0L) {
    return(data.frame(
      gene = gene_symbol, n = 0, group_high = NA, group_low = NA,
      cutpoint = NA, HR_high = NA, HR_low95 = NA, HR_high95 = NA, p_logrank = NA,
      tag = tag, stringsAsFactors = FALSE
    ))
  }
  
  dat <- merge(
    gdat[, .(Unique_Patient_Identifier, expr)],
    clin_dt[, .(Unique_Patient_Identifier, time_months, status)],
    by = "Unique_Patient_Identifier", all = FALSE
  )
  dat <- dat[!is.na(time_months) & !is.na(status)]
  if (nrow(dat) < 10) {
    return(data.frame(
      gene = gene_symbol, n = nrow(dat), group_high = NA, group_low = NA,
      cutpoint = NA, HR_high = NA, HR_low95 = NA, HR_high95 = NA, p_logrank = NA,
      tag = tag, stringsAsFactors = FALSE
    ))
  }
  
  scp <- surv_cutpoint(
    dat,
    time = "time_months",
    event = "status",
    variables = "expr",
    minprop = min_prop
  )
  cut <- scp$cutpoint$cutpoint[1]
  dat$group <- ifelse(dat$expr > cut, "High", "Low")
  dat$group <- factor(dat$group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(time_months, status) ~ group, data = dat)
  cox <- coxph(Surv(time_months, status) ~ group, data = dat)
  cs  <- summary(cox)
  
  HR   <- unname(cs$coef["groupHigh", "exp(coef)"])
  L95  <- unname(cs$conf.int["groupHigh", "lower .95"])
  U95  <- unname(cs$conf.int["groupHigh", "upper .95"])
  pLR  <- tryCatch({ surv_pvalue(fit)$pval }, error = function(e) NA_real_)
  
  if (save_plot) {
    plt <- ggsurvplot(
      fit,
      data = dat,
      pval = TRUE,
      risk.table = TRUE,
      legend.title = "",
      legend.labs = c(paste0("Low (", sum(dat$group == "Low"), ")"),
                      paste0("High (", sum(dat$group == "High"), ")")),
      xlab = "Months",
      ylab = "Survival probability",
      break.time.by = 12
    )
    hr_lab <- sprintf("HR (High) = %.2f", HR)
    p_lab  <- if (is.na(pLR)) "p = NA" else sprintf("p = %.4g", pLR)
    plt$plot <- plt$plot + annotate("text", x = Inf, y = 0.12, hjust = 1.05, vjust = 0,
                                    label = paste(p_lab, hr_lab, sep = "\n"),
                                    size = 3.8)
    
    out_png <- file.path("survival_outputs",
                         paste0(tag, "_", gene_symbol, "_KM.png"))
    ggsave(out_png, plt$plot, width = 5.2, height = 4, dpi = 300)
    
    out_png_rt <- file.path("survival_outputs",
                            paste0(tag, "_", gene_symbol, "_KM_with_risktable.png"))
    ggsave(out_png_rt, plot = cowplot::ggdraw(plt$plot) +
             cowplot::draw_plot(plt$table + theme_cleantable(), y = -0.32, height = 0.32),
           width = 5.2, height = 6.5, dpi = 300)
  }
  
  data.frame(
    gene = gene_symbol,
    n = nrow(dat),
    group_high = sum(dat$group == "High"),
    group_low  = sum(dat$group == "Low"),
    cutpoint = cut,
    HR_high = HR, HR_low95 = L95, HR_high95 = U95,
    p_logrank = pLR,
    tag = tag,
    stringsAsFactors = FALSE
  )
}

# --- 1) Download TCGA-ESCA FPKM from UCSC Xena (GDC hub) ---------------------

local_fpkm <- "TCGA-ESCA.htseq_fpkm.tsv.gz"

xena_urls <- c(
  # GDC hub (primary)
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-ESCA.htseq_fpkm.tsv.gz",
  # Alt S3 path (sometimes works when the first blocks)
  "https://gdc-hub.s3.amazonaws.com/download/TCGA-ESCA.htseq_fpkm.tsv.gz",
  # Toil hub mirror (may or may not host the exact ESCA FPKM; worth a try)
  "https://toil.xenahubs.net/download/TCGA-ESCA.htseq_fpkm.tsv.gz"
)

read_xena_fpkm <- function() {
  if (file.exists(local_fpkm)) {
    message("Reading local file: ", local_fpkm)
    return(fread(local_fpkm))
  }
  # Try curl -L via data.table pipe to bypass 403 and redirects
  for (u in xena_urls) {
    message("Trying: ", u)
    cmd <- paste("curl -fL --retry 3 --retry-delay 2", shQuote(u))
    dt <- tryCatch(
      fread(cmd = cmd),
      error = function(e) NULL
    )
    if (!is.null(dt)) {
      message("✓ Downloaded via curl")
      return(dt)
    }
  }
  stop("Could not download TCGA-ESCA FPKM from Xena. ",
       "Option A: place '", local_fpkm, "' in your working directory. ",
       "Option B: open the first URL in a browser, download manually, then re-run.")
}

expr_wide <- read_xena_fpkm()


# Expect the first column to be gene identifiers
first_col <- names(expr_wide)[1]
setnames(expr_wide, first_col, "gene_raw")

# Derive gene symbol from "ENSEMBL|SYMBOL" if present
if (any(grepl("\\|", expr_wide$gene_raw))) {
  expr_wide[, gene := sub(".*\\|", "", gene_raw)]
} else {
  expr_wide[, gene := gene_raw]
}
