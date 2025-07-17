library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

cesa@ref_data_dir
list.files("C:/Users/bodhi/AppData/Local/R/win-library/4.5/ces.refset.hg38/refset")

gr_genes <- readRDS("C:/Users/bodhi/AppData/Local/R/win-library/4.5/ces.refset.hg38/refset/gr_genes.rds")

gene_lengths <- as.data.table(gr_genes)[, .(length = sum(width(reduce(GRanges(seqnames, IRanges(start, end)))))), by = gene]

reduce(GRanges(seqnames, IRanges(start, end)))

as.data.frame(gr_genes) %>%
  group_by(gene) %>%
  summarize(gene_length = sum(width(reduce)))



genome_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
seqlevels(gr_genes) <- paste0("chr", seqlevels(gr_genes))
seqlevels(gr_genes) <- intersect(seqlevels(gr_genes), seqlevels(genome_seqinfo))
seqinfo(gr_genes) <- genome_seqinfo[seqlevels(gr_genes)]

seqlevels(gr_genes) <- str_replace(seqlevels(gr_genes), pattern = "chrchr", replacement = "chr")

reduce(gr_genes)
