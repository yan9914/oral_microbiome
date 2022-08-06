library(usethis)
library(devtools)
library(Rcpp)
library(dada2)
library(phyloseq)
library(ggplot2)
library(magrittr)
library(biomformat)
library(ShortRead)
library(stringr)

report_file <- read.table("filereport_read_run_PRJEB14674_tsv.txt", sep = "\t", header = TRUE)
metadata <- read.csv("hillburns_metadata.csv")
metadata <- metadata[metadata$excluded != "Bad metadata" | is.na(metadata$excluded),]
metadata %<>% with(metadata[-unique(c(which(is.na(age)), which(is.na(sex)))),])

path <- "./Hill-Burns_fastq"

setwd(path)
list.files() %>%
  {sapply(with(report_file, run_accession[!sample_accession %in% metadata$Accession]),
          \(x) .[str_detect(., x)])} %>%
  {for(i in seq_along(.)) {unlink(.[i]); cat(paste0("The file is deleted: ",.[i], "\n"))}}
file.rename(
  list.files(".", pattern = ".fastq"),
  list.files(".", pattern = ".fastq") %>% 
    sub(".fastq", "", .) %>%
    sub("out_", "", .) %>% 
    {report_file$sample_accession[match(., report_file$run_accession)]} %>%
    paste0(".fastq")
)
setwd("..")

fnFs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
fnFs <- fnFs[sapply(fnFs, \(x) readLines(x) %>% length %>% `/`(4)) >= 5000]

sample_names <- basename(fnFs) %>% sapply(\(x) sub(".fastq", "", x))
names(sample_names) <- NULL


pdf(file = "qualplot_R1_O.pdf", width = 100, height = 50)
plotQualityProfile(fnFs)
dev.off()

filtFs <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fq.gz"))
names(filtFs) <- sample_names

##check out "Matrix" Version (Should be 1.3-2 or 1.3-3)
packageVersion("Matrix")
remove.packages("Matrix")
devtools::install_version("Matrix", version = "1.3.2", repos = "http://cran.us.r-project.org")

out <- filterAndTrim(fnFs, filtFs, truncLen = 0,
                     maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE) 
filtFs <- filtFs[file.exists(filtFs)]

errF <- learnErrors(filtFs, multithread = TRUE)
png(file = "Errplot_O.png", res = 120)
plotErrors(errF, nominalQ = TRUE)
dev.off()

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

seqtab <- makeSequenceTable(dadaFs)
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus", 
                                    multithread = TRUE, 
                                    verbose = TRUE)
sum(seqtab.nochim)/sum(seqtab)
dim(seqtab)
dim(seqtab.nochim)

taxa <- assignTaxonomy(seqtab.nochim, 
                       "silva_nr_v132_train_set.fa.gz", 
                       multithread = TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

for(i in 1:dim(taxa)[1]) {
  if(is.na(taxa[i,7])) {
    if(!is.na(taxa[i,6])) taxa.print[i,7] <- taxa[i,6]
    else if(!is.na(taxa[i,5])) taxa.print[i,7] <- taxa[i,5]
    else if(!is.na(taxa[i,4])) taxa.print[i,7] <- taxa[i,4]
    else if(!is.na(taxa[i,3])) taxa.print[i,7] <- taxa[i,3]
    else if(!is.na(taxa[i,2])) taxa.print[i,7] <- taxa[i,2]
    else taxa.print[i,7] <- taxa[i,1]
  }else taxa.print[i,7] <- taxa[i,7]
}

sub_df <- `[`(metadata,
              which(`%in%`(metadata$Accession,
                           substr(filtFs,
                                  regexpr("filtered/", filtFs[[1]])[1]+9,
                                  regexpr("_F_filt.fq.gz", filtFs[[1]])[1]-1))),)
rownames(sub_df) <- sub_df$Accession

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(as.data.frame(sub_df)),
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0(taxa.print[,7], seq(ntaxa(ps)))

ps_rar <- rarefy_even_depth(ps, rngseed = 123, sample.size = 5000, replace = FALSE)
save(ps_rar, file = "Hill.rda")
otu = t(as(otu_table(ps_rar), "matrix"))
otu_biom = make_biom(data = otu)
write_biom(otu_biom, "Hill.biom")
writeFasta(getSequences(ps_rar@refseq), file = "Hill.fna")
write.table(ps_rar@sam_data, file = "metadata_Hill.tsv", quote = FALSE, sep = "\t")