#' Proteins from unique PSMs
#'
#' Generate protein-level abundance matrix from unique peptides from `psm.tsv` file after selection of modified peptides by best PSMs.
#'
#' @param psms Data frame. Optimally, the best PSMs obtained after processing the `psm.tsv` with the function `psmtsv_to_modified_peptides`.
#' @param annot Data frame obtained from loading the `annotation.txt` from Fragpipe output. First column should be named `channel` and second column should be named `sample`.
#' @return Data frame of summed protein intensities from unique modified peptides.
#' @examples
#'
#' @export
#'
#' @keywords
#'
proteins_from_psms <- function(psms, annot) {

# pre-process annotation file --------------------------

# we need to modify the sample names in the annotation file to match the
# names of the of the TMT intensities after loading into R.

psm_cols_trim <- colnames(psms) %>% # get column names of psm file
                    str_remove("\\..*") # eliminate sufix

psm_cols <- colnames(psms) # get column names of psm file wo eliminate sufix

# create df mapping trimmed sample names and the ones with sufix
trimed2original_sample <- tibble(sample = psm_cols_trim,
                                 sample_trim = psm_cols)

# correct annotation file to match column names for TMT channel intensities
# in the psm file
annot_2 <- annot %>%
                    left_join(., trimed2original_sample) %>% # merge annotation file w untrimmed
                    dplyr::slice(1:nrow(annot)) # keep only sample names present in columns of psm file

# processing ----------------------

## select quant columns from psm.tsv and add suffix

psms_q <- psms %>%
                    filter(`Is Unique` == TRUE) %>% # keep only unique peptides/PSMs
                    dplyr::select(`Protein ID`,
                                  `Modified Peptide`,
                                  annot_2$sample_trim) %>%
                    group_by(`Protein ID`) %>%
                    summarise_if(is.numeric, sum, na.rm = TRUE)


}
