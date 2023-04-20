#' Normalize peptide abundance by protein abundance
#'
#' Takes information of protein intensity and peptide intensity and perform normalization of peptide abundance based on protein abundance.
#'
#' @param peptides Data frame of peptide intensities that will be normalized by protein abundance. Typically gotten from the `purify_n_normalize` function.
#' @param annot Data frame obtained from loading the `annotation.txt` from Fragpipe output. First column should be named `channel` and second column should be named `sample`.
#' @param peptide_annot peptide annotation with specificity information. As it comes from merging the outputs from the functions `annotate_peptides` and `psmtsv_to_modified_peptides`
#' @param summarize_by_specificity
#' @return Data frame of log2 transformed and median-centered peptides abundances corrected by their representative abundance among its correspondent protein abundance based on specific peptides
#' @examples
#'
#' @export
#'
#' @keywords
peptide2protein_normalization <- function(peptides,
                                          annot,
                                          peptide_annot,
                                          summarize_by_specificity = TRUE) {
# pre-process annotation file --------------------------

# we need to modify the sample names in the annotation file to match the
# names of the of the TMT intensities after loading into R.

psm_cols_trim <- colnames(peptides) %>% # get column names of psm file
                 str_remove("\\..*") # eliminate sufix

psm_cols <- colnames(peptides) # get column names of psm file wo eliminate sufix

# create df mapping trimmed sample names and the ones with sufix
trimed2original_sample <- tibble(sample = psm_cols_trim,
                                 sample_trim = psm_cols)

# correct annotation file to match column names for TMT channel intensities
# in the psm file
annot_2 <- annot %>%
                    left_join(., trimed2original_sample) %>% # merge annotation file w untrimmed
                    dplyr::slice(1:nrow(annot)) # keep only sample names present in columns of psm file

# generate protein-quant data frame from peptide abundances ----

## first, we want to be able to differentiate between specific or semi-specific peptides
## ideally, we would want to normalize peptide abundances against protein
## abundances summarized from fully specific peptides.
## This is because we assume that protein abundance based on specific peptides
## better resembles overall abundance of the protein, without the effect of
## proteolytic activity.

if (summarize_by_specificity == TRUE){

specific_peptides <- peptide_annot %>%
                    filter(specificity == "specific") %>%
                    pull(protein_id_modif_pep)
prots_q <- peptides %>%
                    # keep only unique peptides/PSMs
                    filter(`Is Unique` == TRUE,
                           # keep only specific peptides
                           protein_id_modif_pep %in% specific_peptides) %>%
                    dplyr::select(`Protein ID`,
                                  `Modified Peptide`,
                                  annot_2$sample_trim) %>%
                    group_by(`Protein ID`) %>%
                    summarise_if(is.numeric,
                                 sum,
                                 na.rm = TRUE) %>%
                    rename_at(.vars = vars(annot_2$sample_trim),
                              .funs = function(x) paste0(x, "_prot"))

} else if (summarize_by_specificity == FALSE){

prots_q <- peptides %>%
                    filter(`Is Unique` == TRUE) %>% # keep only unique peptides/PSMs
                    dplyr::select(`Protein ID`,
                                  `Modified Peptide`,
                                  annot_2$sample_trim) %>%
                    group_by(`Protein ID`) %>%
                    summarise_if(is.numeric,
                                 sum,
                                 na.rm = TRUE) %>%
                    rename_at(.vars = vars(annot_2$sample_trim),
                              .funs = function(x) paste0(x, "_prot"))

}

# processing files -----

## select quant columns from peptide data frame and add suffix

peptides_q <- peptides %>%
                    dplyr::select(protein_id_modif_pep, `Protein ID`, `Modified Peptide`,
                                  annot_2$sample_trim) %>%
                    rename_at(.vars = vars(annot_2$sample_trim),
                              .funs = function(x) paste0(x, "_peptide")) %>% # add prot suffix to columns
                    mutate(protein_id_modif_pep = paste0(`Protein ID`,"_",`Modified Peptide`))

## merge peptide file with protein file

psm_n_prots <- left_join(peptides_q, prots_q) %>%
                    dplyr::relocate(protein_id_modif_pep)

# function to get the ratio of intensity of Peptide/Protein  -----

## this function is intended to get the fraction of peptide intensity is representative
## of the total protein abundance calculated based on fully specific peptides

pept2prot_ratios <- function(col){

df_q <- psm_n_prots %>%
                    # select id column and columns that match sample names (for map function below)
                    dplyr::select(protein_id_modif_pep, starts_with(col)) %>%
                    # reorganize the column of intensities, so first we have peptide abundances and then proteins
                    dplyr::relocate(protein_id_modif_pep, ends_with("peptide"), ends_with("prot"))
df_q_names <- colnames(df_q) # extract column names
df_q_rat <- df_q %>%
                    # create a column for each sample to get the fraction of the peptide intensity representative of the protein abundance
                    mutate({{col}} := .data[[df_q_names[2]]] / .data[[df_q_names[3]]]) %>%
                    # generate a 'normalized' peptide intensity value, by extracting the peptide abundance associated to it's fraction of the protein abundance
                    mutate("fraction_int_peptide2prot_{col}" := .data[[{{col}}]] * .data[[df_q_names[2]]])
df_q_rat <- df_q_rat %>%
                    # keep only the columns with the fraction of intensities (normalized)
                    dplyr::select(-c(protein_id_modif_pep, ends_with("peptide"), ends_with("prot")))
return(df_q_rat)

}

## PSM intensity/protein ration calculation per sample ----

list_pept2prot_rat <- purrr::map(.x = annot_2$sample,
                                 .f = pept2prot_ratios)

pept2prot_norm1 <- bind_cols(list_pept2prot_rat) %>%
                    mutate(protein_id_modif_pep = psm_n_prots$protein_id_modif_pep) %>%
                    dplyr::relocate(protein_id_modif_pep)

pept2prot_norm_ratio <- pept2prot_norm1 %>%
                    dplyr::select(protein_id_modif_pep, starts_with("fraction_int_"))

## normalizations ---------------------------------------------------

### log2 and median centering of fraction of peptide intensity from peptide / protein ratios -----

mat_ratios2 <- pept2prot_norm_ratio %>%
                    column_to_rownames("protein_id_modif_pep") %>%
                    as.matrix()

mat_ratios2[mat_ratios2 == 0] <- NA

log2_mat_rat2 <- mutate_all(as.data.frame(mat_ratios2),
                            log2)

scaled_mat_rat2 <- scale(log2_mat_rat2,
                         scale = F,
                         center = apply(log2_mat_rat2, 2, median,
                                        na.rm = TRUE) - median(as.matrix(log2_mat_rat2),
                                                               na.rm = TRUE)) %>%
                    abs(.) %>% # log2 values of fractions are negative; take absolute values
                    as.data.frame(.) %>%
                    rownames_to_column("protein_id_modif_pep")

### log2 and median centering of summarized protein abundances -----

mat_ratios3 <- prots_q %>%
                    column_to_rownames("Protein ID") %>%
                    as.matrix()

mat_ratios2[mat_ratios3 == 0] <- NA

log2_mat_rat3 <- mutate_all(as.data.frame(mat_ratios3),
                            log2)

scaled_mat_rat3 <- scale(log2_mat_rat3,
                         scale = F,
                         center = apply(log2_mat_rat3,
                                        2,
                                        median,
                                        na.rm = TRUE) - median(as.matrix(log2_mat_rat3),
                                                               na.rm = TRUE)) %>%
                    abs(.) %>% # log2 values of fractions are negative; take absolute values
                    as.data.frame() %>%
                    rownames_to_column("protein_id")


# create list of results ----

protein_normalized_pepts <- list(protein_normalized_pepts_scaled = scaled_mat_rat2,
                                 protein_normalized_pepts_abundance = pept2prot_norm_ratio,
                                 summarized_protein_abundance = prots_q,
                                 summarized_protein_abundance_scaled = scaled_mat_rat3,
                                 summarize_by_specificity = summarize_by_specificity)

return(protein_normalized_pepts)
}
