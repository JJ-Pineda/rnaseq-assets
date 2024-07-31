# Load packages
library("AnnotationHub")
library("tximport")
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(4))  # Use 4 cores for parallelization

# The following R code loads in the Ensembl v112 human annotation
ah <- AnnotationHub(ask = FALSE)  # if "ask = TRUE" (i.e. default), type "yes" when asked to create directory
ahResults <- query(ah, c("Annotation", "EnsDb", "Homo sapiens", "112"))
ensDbHuman112 <- ahResults[[1]]  # Retrieves the resource
k <- keys(ensDbHuman112, keytype = "TXNAME")
tx2gene <- select(ensDbHuman112, k, "GENENAME", "TXNAME")

# Ideally the data should reside on a EBS volume
data_dir <- "/home/rstudio/javier/data"
salmon_dirs <- list.files(path = data_dir, pattern = "\\.salmon$")
files <- file.path(data_dir, salmon_dirs, "quant.sf")
names(files) <- salmon_dirs  # Enables proper sample labeling for the tximport output
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Record gene counts as a single dataframe
txi_path <- file.path(data_dir, "salmon_gene_counts.csv")
txi_counts <- as.data.frame(txi$counts)
write.csv(txi_counts, file = txi_path, row.names = TRUE, col.names = TRUE)

# Run DESeq2
sample_df <- data.frame(sample = salmon_dirs)
conditions <- factor(rep(c("Colon", "Control"), each=3))  # Should correspond to "sample_df"
conditions <- relevel(conditions, ref = "Control")  # Sets the reference "level"
sample_df <- cbind(sample_df, condition = conditions)
rownames(sample_df) <- sample_df$sample
dds <- DESeqDataSetFromTximport(txi, colData = sample_df, design = ~ condition)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize  # Keep rows that have >=10 counts in >=3 samples
dds <- dds[keep,]
dds <- DESeq(dds, parallel = TRUE)  # Actually runs DESeq2 (with parallelization)
res <- results(dds, name="condition_Colon_vs_Control")  # Tells DESeq2 what comparison to make
res_shrunk <- lfcShrink(dds, coef="condition_Colon_vs_Control", type="apeglm", parallel = TRUE)  # Shrink LFCs

# Record raw LFC results
output_path <- file.path(data_dir, "salmon_deseq2_results.csv")
write.csv(res, file = output_path, row.names = TRUE, col.names = TRUE)

# Record shrunk LFC results
shrunk_output_path <- file.path(data_dir, "salmon_deseq2_results_shrunk.csv")
write.csv(res_shrunk, file = shrunk_output_path, row.names = TRUE, col.names = TRUE)
