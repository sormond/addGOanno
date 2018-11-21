#!/usr/bin/env Rscript

### Shannon Ormond Massey University 2018
### Takes variant (csv) file from VarAFT with refseq ID field and adds GO ID annotation column

# download biomaRt:
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("biomaRt")

# load biomaRt
library("biomaRt")

# load argparse
library("argparse")

parser <- ArgumentParser()
parser$add_argument("-c", "--csv", help="Required argument. VarAFT CSV file to be loaded.")
parser$add_argument("-o", "--out", help="Required argument. Output csv file.")
args <- parser$parse_args()

# load in csv
data <- read.csv(args$csv)

# make list of refseq IDs from loaded csv for annotation
IDs = data["Refseq"]
IDs_vector <- unname(unlist(IDs))
IDs_vector <- gsub("\\..*","",IDs_vector) # removes string after and including '.' so that it's compatible with the biomaRt GO search

# prepare biomaRt for GO annotations
ensembl = useMart("ensembl")
listMarts()
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

# make dataframe with GO terms and for Refseq IDs found in loaded csv
results <- getBM(attributes = c("refseq_mrna", "go_id",
                                "name_1006"), filters = "refseq_mrna",
                 values = IDs_vector, mart = ensembl)

# add GO IDs and GO descriptions to loaded csv
data["go_ids"] <- NA
data["go_descriptions"] <- NA 
for (i in 1:nrow(data)){
  r = data[i, "Refseq"]
  r = gsub("\\..*","",r) # removes string after and including '.' so that it's compatible with the biomaRt GO search
  vector = c()
  go_ids_vector = c()
  go_dspt_vector = c()
  for (j in 1:nrow(results)){
    p = results[j, "refseq_mrna"]
    if (r == p){
      go_id = results[j, "go_id"]
      dspt = results[j, "name_1006"]
      go_ids_vector <- c(go_ids_vector, go_id)
      go_dspt_vector <- c(go_dspt_vector, dspt)
    }
  } 
  go_ids_vector <- paste(go_ids_vector,collapse=" ; ")
  go_dspt_vector <- paste(go_dspt_vector,collapse=" ; ")
  data[i, "go_ids"] = go_ids_vector
  data[i, "go_descriptions"] = go_dspt_vector
}

# write out csv 
write.csv(data, file=args$out)