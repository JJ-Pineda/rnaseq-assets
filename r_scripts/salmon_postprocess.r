# The following R code loads in the Ensembl v112 human annotation
library(AnnotationHub)
ah <- AnnotationHub() # type "yes" when asked to create directory
ahResults <- query(ah, c("Annotation", "EnsDb", "Homo sapiens", "112"))
ensDbHuman112 <- ahResults[[1]] # retrieves the resource
k <- keys(ensDbHuman112, keytype = "TXNAME")
tx2gene <- select(ensDbHuman112, k, "GENEID", "TXNAME")

# Ideally the data should reside on a EBS volume
data_dir <- "/home/rstudio/javier/data"
salmon_dirs <- list.files(path = data_dir, pattern = "\\.salmon$")
files <- file.path(getwd(), salmon_dirs, "quant.sf")
names(files) <- salmon_dirs # Enables proper sample labeling for the tximport output
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

