dir.create("data/pbmc3k", recursive = TRUE, showWarnings = FALSE)

url <- "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
destfile <- "data/pbmc3k/pbmc3k.tar.gz"

cat("Downloading PBMC 3k dataset...\n")
download.file(url, destfile, mode = "wb")
cat("Extracting...\n")
untar(destfile, exdir = "data/pbmc3k")
cat("Data downloaded and extracted to data/pbmc3k/\n")
