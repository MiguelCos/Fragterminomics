#' Generete counts of N-terminal TMT-labelling and acetylation per peptide
#'
#' @param nterm_annot data frame. Output of `annotate_nterm` function.
#'
#' @return a list. containing tabular counts and visualizations.
#' @export
#'
#' @examples
count_location_nterm <- function(nterm_annot){

                    require(dplyr)
                    require(tidytext)
                    require(ggplot2)

                    count_info <- nterm_annot %>%
                                        dplyr::select(protein_id, gene,peptide, peptide_length,
                                                      specificity, nterm, semi_type, tmt_tag,
                                                      start_position, end_position, last_aa, aa_after,
                                                      aa_before, following_10_resid, previous_10_resid, protein_length, previous_all_resid) %>%
                                        mutate(start_position_num = as.numeric(start_position)) %>%
                                        mutate(normalized_location = round(start_position_num/protein_length*100))

                    ## keep only acetylated nterm ----

                    count_info_acetylated <- count_info %>%
                                        filter(nterm == "acetylated")

                    ### keep only TMT-labelled semi-tryptic Nterm ----

                    count_info_tmtneo <- count_info %>%
                                        filter(nterm == "TMT-labelled",
                                               semi_type == "semi_Nterm")

                    ### how many times do we have a specific 'AA before' from acetylated N-term ----

                    n_aa_bef_acet <- count_info_acetylated %>%
                                        dplyr::count(aa_before = as.factor(aa_before)) %>%
                                        mutate(`N-term type` = "Acetylated",
                                               `Feature`= "Amino Acid before",
                                               aa_before = fct_reorder(aa_before, n, .desc = FALSE))

                    ### how many times do we have a specific 'normalized location' from acetylated N-term ----


                    n_norm_locat_acet <- count_info_acetylated %>%
                                        dplyr::count(normalized_location = as.factor(normalized_location)) %>%
                                        mutate(`N-term type` = "Acetylated",
                                               `Feature`= "Normalized location",
                                               normalized_location = fct_reorder(normalized_location,
                                                                                 as.numeric(normalized_location),
                                                                                 .desc = FALSE))


                    ### how many times do we have a specific 'AA before' from TMT N-term ----

                    n_aa_bef_tmtneo <- count_info_tmtneo %>%
                                        dplyr::count(aa_before = as.factor(aa_before)) %>%
                                        mutate(`N-term type` = "TMT-labelled",
                                               `Feature`= "Amino Acid before",
                                               aa_before = fct_reorder(aa_before, n, .desc = FALSE))

                    ### how many times do we have a specific 'normalized location' from TMT N-term ----

                    n_norm_locat_tmtneo <- count_info_tmtneo %>%
                                        dplyr::count(normalized_location = as.factor(normalized_location)) %>%
                                        mutate(`N-term type` = "TMT-labelled",
                                               `Feature`= "Normalized location",
                                               normalized_location = fct_reorder(normalized_location, as.numeric(normalized_location), .desc = FALSE))

                    ### Merge count tables of AA before into one for faceted visualization ----

                    counts_aa_before <- bind_rows(n_aa_bef_acet,
                                                  n_aa_bef_tmtneo)

                    ### Merge count tables of normalized location into one ofr faceted visualization ----

                    counts_norm_locat <- bind_rows(n_norm_locat_acet,
                                                   n_norm_locat_tmtneo)


                    ### Plots AA before ----

                    plot_aa_before <- ggplot(counts_aa_before,
                                             aes(x = reorder_within(aa_before, n, `N-term type`), y = n)) +
                                        geom_bar(stat = "identity") +
                                        scale_x_reordered() +
                                        coord_flip() +
                                        labs(y = "Number of IDs", x = "Amino Acid Before")+
                                        facet_wrap(~`N-term type`, scales="free")+
                                        theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                                              axis.text.y = element_text(hjust = 0.5, size = 10),
                                              panel.background = element_blank(),
                                              panel.grid.major = element_line(color = "grey"),
                                              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                                              axis.title=element_text(size=12,face="bold"))

                    ### Plots Normalized location -----

                    plot_norm_locat <- ggplot(counts_norm_locat,
                                              aes(x = forcats::fct_rev(normalized_location), y = n)) +
                                        geom_bar(stat = "identity") +
                                        coord_flip() +
                                        labs(y = "Number of IDs", x = "Normalized start position")+
                                        facet_wrap(~`N-term type`, scales="free")+
                                        theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                                              axis.text.y = element_text(hjust = 0.5, size = 10),
                                              panel.background = element_blank(),
                                              panel.grid.major = element_line(color = "grey"),
                                              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                                              axis.title=element_text(size=12,face="bold"))


                    list_location_nterms <- list(general_peptide_annotation = count_info,
                                                 count_info_acetyl = count_info_acetylated,
                                                 count_info_tmtneo = count_info_tmtneo,
                                                 plot_aa_before = plot_aa_before,
                                                 plot_normalized_location = plot_norm_locat)

                    return(list_location_nterms)


}
