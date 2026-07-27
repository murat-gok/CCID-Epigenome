#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outfile <- if (length(args) >= 1) args[[1]] else "data/reference/epicv2_probe_annotation_hg38.csv.gz"

required <- c("BiocManager", "minfi", "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in required[-1]) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))

obj_name <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
obj <- get(obj_name, envir = asNamespace(obj_name))
ann <- as.data.frame(minfi::getAnnotation(obj))
ann$cpg_id <- rownames(ann)

pick <- function(candidates, default = NA_character_) {
  hit <- candidates[candidates %in% colnames(ann)]
  if (length(hit) == 0) return(rep(default, nrow(ann)))
  ann[[hit[[1]]]]
}

out <- data.frame(
  cpg_id = ann$cpg_id,
  chromosome = pick(c("chr", "CHR")),
  position = pick(c("pos", "MAPINFO", "Coordinate_36")),
  strand = pick(c("strand", "Strand")),
  gene = pick(c("UCSC_RefGene_Name", "Gene"), ""),
  gene_group = pick(c("UCSC_RefGene_Group", "Gene_Group"), ""),
  island_relation = pick(c("Relation_to_Island", "Relation_to_UCSC_CpG_Island"), ""),
  genome_build = "hg38",
  annotation_source = obj_name,
  stringsAsFactors = FALSE
)
out <- out[!is.na(out$cpg_id) & !is.na(out$chromosome) & !is.na(out$position), ]
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(outfile, open = "wt")
write.csv(out, con, row.names = FALSE, quote = TRUE)
close(con)
cat(sprintf("Wrote %d annotated probes to %s\n", nrow(out), outfile))
