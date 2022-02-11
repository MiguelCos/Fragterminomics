#' Summarize counts of peptide identifications
#'
#' @param annotated_peptide_ids data frame with annotated peptide IDs including specificity and N-terminal modification.
#'
#' @return table in data frame format.
#' @export
#'
#' @examples
summarize_peptide_counts <- function(annotated_peptide_ids){

                    require(dplyr)
                    require(tibble)
                    require(magrittr)

                    to_count_info <- annotated_peptide_ids %>%
                                        dplyr::select(peptide, specificity, nterm, semi_type, tmt_tag)

                    to_count_info_acetyl <- to_count_info %>%
                                        dplyr::filter(nterm == "acetylated")

                    n_semi <- dplyr::count(to_count_info, specificity) %>%
                                        dplyr::rename(feature_type = specificity) %>%
                                        dplyr::mutate(category = "specificity")

                    n_term <- dplyr::count(to_count_info, nterm) %>%
                                        dplyr::rename(feature_type = nterm) %>%
                                        dplyr::mutate(category = "N-term")

                    n_semi_type <- dplyr::count(to_count_info, semi_type) %>%
                                        dplyr::rename(feature_type = semi_type) %>%
                                        dplyr::mutate(category = "Semi type")

                    n_tmt_tag <- dplyr::count(to_count_info, tmt_tag) %>%
                                        dplyr::rename(feature_type = tmt_tag)  %>%
                                        dplyr::mutate(category = "TMT location")

                    n_total <- tibble(feature_type = "Total",
                                      n = nrow(to_count_info),
                                      category = "Total")

                    summary_count <- bind_rows(n_semi,
                                               n_term,
                                               n_semi_type,
                                               n_tmt_tag,
                                               n_total)

                    return(summary_count)

}
#' Title
#'
#' @param annotated_peptide_ids data frame with annotated peptide IDs including specificity and N-terminal modification.
#'
#' @return A list. Each element is a vector of peptide sequences for each interesting feature.
#' @export
#'
#' @examples
list_peptides <- function(annotated_peptide_ids){

                    to_count_info <- annotated_peptide_ids %>%
                                        dplyr::select(peptide,
                                                      specificity,
                                                      nterm,
                                                      semi_type,
                                                      tmt_tag,
                                                      is_terminal)

                    specificity_specific <- filter(to_count_info,
                                                   specificity == "specific") %>%
                                        pull(peptide)

                    specificity_semi <- filter(to_count_info,
                                               specificity == "semi_specific") %>%
                                        pull(peptide)

                    nterm_acetyl <- filter(to_count_info,
                                           nterm == "acetylated") %>%
                                        pull(peptide)

                    nterm_tmtlab <- filter(to_count_info,
                                           nterm == "TMT-labelled") %>%
                                        pull(peptide)

                    nterm_free <- filter(to_count_info,
                                         nterm == "free") %>%
                                        pull(peptide)

                    semitype_Nterm <- filter(to_count_info,
                                             semi_type == "semi_Nterm") %>%
                                        pull(peptide)

                    semitype_Cterm <- filter(to_count_info,
                                             semi_type == "semi_Cterm") %>%
                                        pull(peptide)

                    terminal_yes <- filter(to_count_info,
                                           is_terminal == "terminal") %>%
                                        pull(peptide)

                    terminal_no <- filter(to_count_info,
                                          is_terminal == "not_terminal") %>%
                                        pull(peptide)

                    listed_peptides <- list(specificity_specific = specificity_specific,
                                            specificity_semi = specificity_semi,
                                            nterm_acetyl = nterm_acetyl,
                                            nterm_tmtlab = nterm_tmtlab,
                                            nterm_free = nterm_free,
                                            semitype_Nterm = semitype_Nterm,
                                            semitype_Cterm = semitype_Cterm,
                                            terminal_yes = terminal_yes,
                                            terminal_no = terminal_no)

                    return(listed_peptides)

}

