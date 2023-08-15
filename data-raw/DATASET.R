# required packages

library(readr)
library(dplyr)
library(janitor)
library(here)
library(seqinr)

# Annotation table

# load the tabular information from the samples
sample_annotation <- read_csv(here("data-raw/annotation.csv"))

usethis::use_data(sample_annotation, overwrite = TRUE)

# Normalized protein abundance matrix; FragPipe/TMT-integrator output

prot_abund_mat <- read_tsv(here("data-raw/abundance_protein_MD.tsv")) %>%
                    clean_names() # apply clean_names for consistent column formatting

usethis::use_data(prot_abund_mat, overwrite = TRUE)

# Normalized peptide abundance matrix; FragPipe/TMT-integrator output

pept_abund_mat <- read_tsv(here("data-raw/abundance_peptide_MD.tsv")) %>%
                    clean_names() # apply clean_names for consistent column formatting

usethis::use_data(pept_abund_mat, overwrite = TRUE)

# Tabular results of Protein and Peptide IDs

prot_ident <- read_tsv(here("data-raw/protein.tsv")) %>%
                    clean_names()

usethis::use_data(prot_ident, overwrite = TRUE)

pept_ident <- read_tsv(here("data-raw/peptide.tsv")) %>%
                    clean_names()

usethis::use_data(pept_ident, overwrite = TRUE)

# FASTA file of sequences of identified proteins

fasta <- read.fasta(here("data-raw/protein.fas"),
                    as.string = TRUE,
                    seqtype = "AA")

usethis::use_data(fasta, overwrite = TRUE)

# PSM file of peptide-to-spectrum matches including modified peptide information

psm_tsv <- read.fasta(here("data-raw/psm.tsv"),
                    as.string = TRUE,
                    seqtype = "AA")

usethis::use_data(psm_tsv, overwrite = TRUE)