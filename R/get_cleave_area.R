#' Prepare sequences for motif analysis of cleavage areas
#'
#' @param peptide_annotation A data frame. Output from the `annotate_peptides` function.
#'
#' @return A data frame containing information of cleavage areas after proteolytic processing for iceLogo analysis.
#' @export
#'
#' @examples
get_cleave_area <- function(peptide_annotation){

                    require(dplyr)
                    require(stringr)

                    nterm_semi_pepts <- peptide_annotation %>%
                                        filter(semi_type == "semi_Nterm") %>%
                                        dplyr::select(-previous_all_resid) %>%
                                        mutate(previous_10_resid = if_else(is.na(previous_10_resid), true = "", false = previous_10_resid),
                                               following_10_resid = if_else(is.na(following_10_resid), true = "", false = following_10_resid)) %>%
                                        mutate(prev_10_pad = str_pad(previous_10_resid, side = "left", pad = "X", width = 10),
                                               following_10_pad = str_pad(following_10_resid, side = "right", pad = "X", width = 10)) %>%
                                        mutate(cleave_area = paste0(prev_10_pad, Peptide, following_10_pad)) %>%
                                        mutate(cleave_area20 = str_sub(cleave_area,
                                                                       start = 1,
                                                                       end = 20)) %>%
                                        mutate(len_20 = str_count(cleave_area20),
                                               len_clearea = str_count(cleave_area),
                                               len_follow = str_count(following_10_resid),
                                               len_prev = str_count(previous_10_resid),
                                               len_followpad = str_count(prev_10_pad),
                                               len_prevpad = str_count(following_10_pad))

                    cterm_semi_pepts <- peptide_annotation %>%
                                        filter(semi_type == "semi_Cterm") %>%
                                        dplyr::select(-previous_all_resid) %>%
                                        mutate(previous_10_resid = if_else(is.na(previous_10_resid), true = "", false = previous_10_resid),
                                               following_10_resid = if_else(is.na(following_10_resid), true = "", false = following_10_resid)) %>%
                                        mutate(prev_10_pad = str_pad(previous_10_resid, side = "left", pad = "X", width = 10),
                                               following_10_pad = str_pad(following_10_resid, side = "right", pad = "X", width = 10)) %>%
                                        mutate(cleave_area = paste0(prev_10_pad, Peptide, following_10_pad)) %>%
                                        mutate(cleave_area20 = str_sub(cleave_area,
                                                                       -20)) %>%
                                        mutate(len_20 = str_count(cleave_area20),
                                               len_clearea = str_count(cleave_area),
                                               len_follow = str_count(following_10_resid),
                                               len_prev = str_count(previous_10_resid),
                                               len_followpad = str_count(prev_10_pad),
                                               len_prevpad = str_count(following_10_pad))

                     cleave_area20 <- bind_rows(nterm_semi_pepts,
                                               cterm_semi_pepts) %>% 
                     # generate column with descriptive cleavage site
                     mutate(substring1 = str_sub(cleave_area20, start = 1, end = 10),
                            substring2 = str_sub(cleave_area20, start = 11, end = 20)) %>%
                     mutate(cleavage_site = paste(substring1, substring2, sep = " | ")) %>%
                     mutate(short_cleavage_site = substr(cleavage_site, start = 7, stop = nchar(cleavage_site) - 6)) %>%
                     dplyr::select(-substring1, -substring2)

                    list_cleavs <- list(tab_nterm_semi_pepts = nterm_semi_pepts,
                                        tab_cterm_semi_pepts = cterm_semi_pepts,
                                        cleave_area20 = cleave_area20)

}
