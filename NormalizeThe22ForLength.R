library(xlsx)
library(cancereffectsizeR)
library(data.table)
library(tidyverse)
library(biomaRt)
library(ggplot2)
options(java.parameters = "-Xmx2g")

# ONCE
names(data)

# set column names and delete extra row 
colnames(data) <- c("Unique_Patient_Identifier")
colnames(data)[2] <- c("Chromosome")
colnames(data)[3] <- c("Start_Position")
colnames(data)[4] <- c("End_Position")
colnames(data)[5] <- c("Reference_Allele")
colnames(data)[6] <- c("Tumor_Seq_Allele1")
colnames(data)[7] <- c("Hugo_Symbol")
colnames(data)[8] <- c("Variant_Classification")
colnames(data)[9] <- c("Variant_Type")
colnames(data)[10] <- c("Func.refGene")
data <- data[-1,]

# make the data into a .maf file
write.table(
  data,
  file = "ESCCdata.maf",
  sep = "\t",
  quote = FALSE,
)

# END ONCE

maf = preload_maf("ESCCdata.maf", refset = "ces.refset.hg38")
maf = maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
cesa <- CESAnalysis(refset = "ces.refset.hg38")
cesa <- load_maf(cesa = cesa, maf = maf, coverage = "genome")

cesa$maf
cesa$variants
cesa$samples
(top_variants <- cesa$variants[order(-maf_prevalence)][1:10, .(variant_name, chr, start, end, maf_prevalence)])

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "Eso-SCC", treatment_naive = TRUE)


cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
                              signature_exclusions = signature_exclusions
)


cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$ESCA)

head(cesa$gene_rates)

dndscv_results <- cesa$dNdScv_results[[1]]
sig_genes <- dndscv_results[qallsubs_cv < .05][order(qallsubs_cv)][1:10]

# take the top 3 variants by MAF prevalence
variants_to_check <- cesa$variants[order(-maf_prevalence), variant_id][1:3]

baseline_mutation_rates(cesa = cesa, variant_ids = variants_to_check)

cesa <- ces_variant(cesa = cesa, run_name = "non_recurrents", variants = cesa$variants)

plot_effects(effects = cesa$selection$non_recurrents)
plot_effects(cesa$selection$non_recurrents,
             group_by = "gene", topn = 10,
             label_individual_variants = FALSE
)

cesa$selection$non_recurrents %>% arrange(desc(included_total))
View(cesa)

#mutations per base pair (normalize for size)
# count mutations per gene
mut_counts <- cesa$variants %>% count(gene, name = "mutations")
# gene_lengths <- data %>% count(Hugo_Symbol) %>% arrange(desc(n))

# compute mutations per kilobase
normalized <- mut_counts %>%
  inner_join(gene_lengths, by = "gene") %>% mutate(length_kb = length / 1e3,
  muts_per_kb = mutations / length_kb) %>%
  arrange(desc(muts_per_kb))


#filter out LINC genes
normalized_filtered <- normalized %>% filter(!str_starts(gene, "LINC"))

# bar graph of mutations per kilobase
normalized_filtered %>%
  arrange(desc(muts_per_kb)) %>%
  slice_head(n = 10) %>%
  ggplot(aes(x = reorder(gene, muts_per_kb), y = muts_per_kb)) +
  geom_col(fill = "seagreen") +
  coord_flip() +
  labs(title = "Top 10 Genes by Mutations per kb",
       x = "Gene", y = "Mutations per kb") +
  theme_minimal()
