library(data.table)
library(cancereffectsizeR)

get_gene_mutations <- function(cesa, genes) {
  # Pull the variants table and subset  genes
  dt <- copy(cesa$variants)[gene %in% genes]
  if (nrow(dt) == 0L) {
    warning("No mutations found for those genes.")
  }
  return(dt[])
}

my_genes <- c("TP53", "NOTCH1", "CSMD3", "EP300", "FAM135B", "TMEM121", "TARP", "DAD1", "TTN", "MUC16")
mutations_in_my_genes <- get_gene_mutations(cesa, my_genes)
print(mutations_in_my_genes)




# check samples
get_variant_samples <- function(cesa, variant_ids) {
  # For each variant_id: samples_with()
  dt_list <- lapply(variant_ids, function(vid) {
    s <- samples_with(cesa, vid)
    if (length(s) == 0L) return(NULL)
    data.table(variant_id = vid, sample_id = s)
  })
  # Combine into one data.table
  out <- rbindlist(dt_list, use.names = TRUE)
  if (nrow(out) == 0L) {
    warning("No samples found for any of those variant_ids.")
  }
  return(out[])
}

mut_dt <- get_gene_mutations(cesa, my_genes) # from prior function

# Pull out the variant_ids
vids <- mut_dt$variant_id

# Build the sample variant table......
variant_samples_dt <- get_variant_samples(cesa, vids)

head(variant_samples_dt)


count_variant_prevalence <- function(variant_samples_dt) {
  variant_samples_dt[
    , .(n_samples = uniqueN(sample_id)), 
    by = variant_id
  ][order(-n_samples)]
}

variant_counts <- count_variant_prevalence(variant_samples_dt)

# view top variants by sample count
print(head(variant_counts, 20))


# COLLAPSING so that it is JUST the sum of all of the times the gene has been mutated
collapse_variant_counts_by_gene <- function(variant_counts, cesa) {
  # each variant_id is now back to its gene
  gene_map <- cesa$variants[, .(variant_id, gene)]
  dt <- merge(variant_counts, gene_map, by = "variant_id", all.x = TRUE)
  # sum the n_samples within each gene
  gene_counts <- dt[, .(combined_n_samples = sum(n_samples)), by = gene][
    order(-combined_n_samples)
  ]
  
  return(gene_counts[])
}

gene_counts <- collapse_variant_counts_by_gene(variant_counts, cesa)
print(gene_counts)
