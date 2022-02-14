#' @title Sample annotation table
#'
#' @description Annotation table of samples and conditions for use case experiment of
#'      polycystic kidney disease in mice.
#'
#' @format A data frame with 30879 rows and 33 variables:
#' \describe{
#'   \item{channel}{TMT channel}
#'   \item{sample_id}{Sample ID}
#'   \item{condition}{Condition}
#' }
#' @source
"sample_annotation"
#' @title Protein sequences of identified proteins
#'
#' @description Protein sequences of identified proteins in the Polycystic Disease mouse model experiment
#'
#' @format A list 5986 elements, each represents a protein sequence identified in the Polycystic Disease mouse model experiment.
#'
#' @source
"fasta"
#' @title Protein abundance data frame
#'
#' @description Normalized values of abundance of identified proteins from the Polycystic Disease mouse model experiment
#'
#' @format A data frame with 5238 rows and 21 variables:
#' \describe{
#'   \item{index}{Protein identifier}
#'   \item{number_psm}{Number of Peptide-to-spectrum matches for the protein}
#'   \item{max_pep_prob}{Max Peptide Probability}
#'   \item{reference_intensity}{Normalized intensity values of virtual normalization channel}
#'   \item{ko/wt}{Normalized intensity values of protein by sample}
#'   \item{mt_*}{Normalized intensity values of empty TMT channels}
#' }
#' @source
"prot_abund_mat"
#' @title Peptide abundance data frame
#'
#' @description Normalized values of abundance of identified proteins from the Polycystic Disease mouse model experiment
#'
#' @format A data frame with 21032 rows and 22 variables:
#' \describe{
#'   \item{index}{Protein and peptide sequence}
#'   \item{gene}{Gene identifier}
#'   \item{protein_id}{Protein identifier}
#'   \item{max_pep_prob}{Max Peptide Probability}
#'   \item{reference_intensity}{Normalized intensity values of virtual normalization channel}
#'   \item{ko/wt}{Normalized intensity values of protein by sample}
#'   \item{mt_*}{Normalized intensity values of empty TMT channels}
#' }
#' @source
"pept_abund_mat"
#' @title Proteins identified and associated modifications
#'
#' @description Protein identification table from the Polycystic Disease mouse model experiment
#'
#' @format A data frame with 5986 rows and 41 variables:
#' \describe{
#'   \item{group}{Protein group}
#'   \item{sub_group}{Protein subgroup}
#'   \item{protein}{Protein header}
#'   \item{protein_id}{Uniprot protein ID}
#'   \item{entry_name}{Uniprot Entry ID}
#'   \item{gene}{Gene name}
#'   \item{length}{Protein length}
#'   \item{percent_coverage}{Coverage of protein identification}
#'   \item{organism}
#'   \item{protein_description}
#'   \item{protein_existence}
#'   \item{protein_description}
#'   \item{protein_probability}
#'   \item{top_peptide_probability}
#'   \item{total_peptides}
#'   \item{unique_peptides}{Number of unique peptides identified}
#'   \item{razor_peptides}
#'   \item{total_spectral_count}{Total spectral counts of identified peptides}
#'   \item{unique_spectral_count}{Spectral counts of unique peptides}
#'   \item{razor_spectral_count}{Spectral counts of razor peptides}
#'   \item{total_intensity}
#'   \item{unique_intensity}
#'   \item{razor_intensity}
#'   \item{razor_assigned_modifications}
#'   \item{razor_observed_modifications}
#'   \item{indistinguishable_proteins}
#'   \item{indistinguishable_proteins}
#'   \item{ko/wt}{Raw intensity values of reporter ions by sample}
#'   \item{mt_*}{Raw intensity values of empty TMT channels}
#' }
#' @source
"prot_ident"
#' @title Peptides identified and associated modifications
#'
#' @description Peptides identification table from the Polycystic Disease mouse model experiment
#'
#' @format A data frame with 30879 rows and 33 variables:
#' \describe{
#'   \item{peptide}{Peptide sequence}
#'   \item{prev_aa}{Amino acid before in the protein sequence}
#'   \item{next_aa}{Amino acid after in the protein sequence}
#'   \item{peptide_length}
#'   \item{probability}{Peptide Probability Score}
#'   \item{spectral_count}{Spectral counts for identified peptide}
#'   \item{intensity}{Raw intensity value of precursor ion}
#'   \item{assigned_modifications}{Chemical modifications assigned to the identified peptide}
#'   \item{observed_modifications}
#'   \item{protein}{Protein header}
#'   \item{protein_id}{Uniprot ID}
#'   \item{entry_name}{Uniprot entry name}
#'   \item{gene}
#'   \item{protein_description}
#'   \item{mapped_genes}
#'   \item{mapped_proteins}{Number of unique peptides identified}
#'   \item{ko/wt}{Raw intensity values of reporter ions by sample}
#'   \item{mt_*}{Raw intensity values of empty TMT channels}
#' }
#' @source
"pept_ident"




