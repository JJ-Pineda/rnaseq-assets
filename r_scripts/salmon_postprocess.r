# Read input arguments (i.e. data directory and design matrix filename)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
design_matrix_fname <- args[2]

# Validate data directory
all_salmon_dirs <- list.files(path = data_dir, pattern = "\\.salmon$")
if (length(all_salmon_dirs) == 0) {
  print("No salmon directories detected using the specified data path...Exiting")
  quit(save = "no")
}

# Read in design matrix CSV
# Expects "sample" (string), "condition" (string), and "reference" (boolean) columns
design_matrix <- read.csv(file.path(data_dir, design_matrix_fname), stringsAsFactors = FALSE)
rownames(design_matrix) <- design_matrix$sample
relevant_salmon_dirs <- design_matrix$sample

# Now load packages
library("AnnotationHub")
library("tximport")
library("tibble")
library("DESeq2")
library("BiocParallel")

# Use 4 cores for any parallelization
register(MulticoreParam(4))

# Load in Ensembl v112 human annotation and build table for transcript/gene conversion
ah <- AnnotationHub(ask = FALSE)  # if "ask = TRUE" (i.e. default), type "yes" when asked to create directory
ahResults <- query(ah, c("Annotation", "EnsDb", "Homo sapiens", "112"))
ensDbHuman112 <- ahResults[[1]]  # Retrieves the resource
k <- keys(ensDbHuman112, keytype = "TXNAME")
tx2gene <- select(ensDbHuman112, k, "GENENAME", "TXNAME")

# Perform transcript count aggregation with tximport and output gene counts
files <- file.path(data_dir, relevant_salmon_dirs, "quant.sf")
names(files) <- relevant_salmon_dirs  # Enables proper sample labeling for the tximport output
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_path <- file.path(data_dir, "salmon_gene_counts.csv")
txi_counts <- as.data.frame(txi$counts)
txi_counts <- rownames_to_column(txi_counts, var = "gene_name")
write.csv(txi_counts, file = txi_path, row.names = FALSE)

# Gear up for DESeq2
ref_condition = design_matrix[design_matrix$reference == TRUE, "condition"][1]
other_conditions = unique(design_matrix[design_matrix$reference == FALSE, "condition"])
design_matrix$condition <- factor(design_matrix$condition)
design_matrix$condition <- relevel(design_matrix$condition, ref = ref_condition)
dds <- DESeqDataSetFromTximport(txi, colData = design_matrix, design = ~ condition)

smallest_group_size <- 1000000
for (cond in other_conditions) {
  group_size <- length(design_matrix[design_matrix$condition == cond])
  if (group_size < smallest_group_size) {
    smallest_group_size <- group_size
  }
}

count_cutoff <- 10
keep <- rowSums(counts(dds) >= count_cutoff) >= smallest_group_size  # Keep rows that have >=10 counts in >=3 samples
dds <- dds[keep,]

# Run DESeq2 (with parallelization)
dds <- DESeq(dds, parallel = TRUE)  # Actually runs DESeq2 (with parallelization)

# Iteratively generate comparison results
for (cond in other_conditions) {
  comparison <- paste("condition", cond, "vs", ref_condition, sep = "_")
  res <- results(dds, name = comparison)
  res_shrunk <- lfcShrink(dds, coef = comparison, type = "apeglm", parallel = TRUE)

  res_fname <- paste0("salmon_deseq2_", cond, "_vs_", ref_condition, ".csv")
  res_shrunk_fname <- paste0("salmon_deseq2_", cond, "_vs_", ref_condition, "_shrunkLFC.csv")

  res <- rownames_to_column(as.data.frame(res), var = "gene_name")
  res_shrunk <- rownames_to_column(as.data.frame(res_shrunk), var = "gene_name")

  write.csv(res, file = file.path(data_dir, res_fname), row.names = FALSE)
  write.csv(res, file = file.path(data_dir, res_shrunk_fname), row.names = FALSE)
}

# # To read in tximport CSV output and convert to DESeqDataSet
# # Note: counts must have integer values
# txi_counts <- read.csv(txi_path, check.names = FALSE)
# rownames(txi_counts) <- txi_counts$gene_name
# txi_counts <- subset(txi_counts, select = -gene_name)
# numeric_cols <- sapply(txi_counts, is.numeric)
# txi_counts[numeric_cols] <- lapply(txi_counts[numeric_cols], as.integer)
#
# dds <- DESeqDataSetFromMatrix(txi_counts, design_matrix, design = ~ condition)
#
# # To read in DESeq CSV output and convert to DESeqResult object
# res_shrunk_df <- read.csv(file.path(data_dir, res_shrunk_fname), check.names = FALSE)
# rownames(res_shrunk_df) <- res_shrunk_df$gene_name
# res_shrunk_df <- subset(res_shrunk_df, select = -gene_name)
# res_shrunk <- DESeqResults(res_shrunk_df)

print("All done!")

