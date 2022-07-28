#' Annotate peptides based on their N-terminal (bio)chemical modification.
#'
#' @param peptidestsv A data frame coming from the `peptide.tsv` tabular output from a FragPipe search.
#' @param tmtmass A numeric object. Defining the mass of your TMT labels. 304.2072 for 16plex or 229.1629 for 10/11plex.
#' @param protease_specificity A string. Defining the protease specificity. Examples: "R|K" for Trypsin; "R" for Argc.
#'
#' @return A data frame with annotation of N-terminal modifications.
#' @export
#'
#' @examples
annotate_nterm <- function(peptidestsv, # peptide.tsv table
                           tmtmass = 304.2072, # either 304.2072 for 16plex or 229.1629 for 10/11plex
                           protease_specificity = "R|K") # or R for argc
{

require(dplyr)
require(stringr)
require(magrittr)

if (tmtmass == 304.2072){

                    nterm_tmt <- "N-term\\(304.2072\\)"
                    ktmt <- "K\\(304.2072\\)"

} else if (tmtmass == 229.1629){

                    nterm_tmt <- "N-term\\(229.1629\\)"
                    ktmt <- "K\\(229.1629\\)"

} else if (tmtmass == 304.2071){

                    nterm_tmt <- "N-term\\(304.2071\\)"
                    ktmt <- "K\\(304.20721\\)"

} else {

                    errorCondition("Please check your tmtmass argument input.")

}

nterm_annotated <- peptidestsv %>%
                    mutate(nterm = case_when(str_detect(assigned_modifications, nterm_tmt) ~ "TMT-labelled",
                                             str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "acetylated",
                                             TRUE ~ "free")) %>%
                    mutate(tmt_tag = case_when(str_detect(assigned_modifications, nterm_tmt) ~ "nterm",
                                               str_detect(assigned_modifications, ktmt) ~ "lysine",
                                               str_detect(assigned_modifications, ktmt,
                                                          negate = TRUE) & str_detect(assigned_modifications, "N-term\\(42.0106\\)") ~ "untagged",
                                               str_detect(assigned_modifications, ktmt,
                                                          negate = TRUE) & nterm == "acetylated" ~ "untagged",
                                               str_detect(assigned_modifications, ktmt,
                                                          negate = TRUE) & nterm == "free" ~ "untagged",
                                               TRUE ~ "untagged"))

return(nterm_annotated)

}
