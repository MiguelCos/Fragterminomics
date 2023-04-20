#' Process `psm.tsv` file into modified peptides
#'
#' Filters the `psm.tsv` based on a 'best PSM' criteria; based on purity and summed TMT reporter intensities. Zero (0) abundance values are substituted by NA (missing values).
#'
#' @param psms Data frame obtained from loading the `psm.tsv` file from Fragpipe output
#' @param annot Data frame obtained from loading the `annotation.txt` from Fragpipe output. First column should be named `channel` and second column should be named `sample`.
#' @param peptide_probability Minimal peptide probability accepted per PSM
#' @param minimal_purity Minimal isobaric purity accepted per PSM
#' @return Data frame of summarized modified peptides based on their 'best PSMs' based on the filtering criteria
#' @examples
#'
#' @export
#'
#' @keywords
#'
psmtsv_to_modified_peptides <- function(psms, annot, peptide_probability, minimal_purity) {

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

# processing and filtering psm.tsv file -----------------------------------

## create reduced dataframe ----
psm_mod1 <- psms %>%
 # select interesting columns for filtering and testing
 dplyr::select(`Peptide`, `Modified Peptide`,Spectrum, Purity,
               Intensity,`PeptideProphet Probability`, `Is Unique`,
               `Quan Usage`, annot_2$sample_trim) %>%
 rowwise() %>%
 # generate column  of summed intensities of TMT channels
 # the selection of the cols representing the TMT channels is based on
 # the annot_2 object
 mutate(summ_int = sum(c_across(cols = annot_2$sample_trim)))

## filter the reduced dataframe ----

psm_mod2 <- psm_mod1 %>%
 # keep only psms with summed TMT intensity > 0
 filter(summ_int > 0) %>%
 # keep only psms with peptideProphet prob > X
 filter(`PeptideProphet Probability` > peptide_probability) %>%
 # keep only psms with purity > minimal_purity
 filter(Purity > minimal_purity) %>%
 # group by modified peptide
 group_by(`Modified Peptide`) %>%
 # if several PSMs to the same modified peptide, select one with best Purity
 slice_max(n = 1,
           order_by = c(Purity)) %>%
 # if a modified peptide was identified with several PSMs and same high purity,
 # then select the PSMs with highest summed TMT intensity
 slice_max(n = 1,
           order_by = summ_int)

## filter the original psm.tsv file ----

best_modified_psms <- psms %>%
 # keep psms after filtering by purity and summed intensity
 filter(Spectrum %in% psm_mod2$Spectrum)

return(best_modified_psms)

}
