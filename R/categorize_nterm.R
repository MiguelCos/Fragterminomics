#' Categorize and count N-terminal peptides based on Uniprot processing annotation
#'
#' @param annotated_peptides data frame of annotated peptides. The output of the `annotate_nterm` function.
#' @param nterm a character vector. Defines the type of modification expected to be found in N-terminal peptides. `c("TMT-labelled", "acetylated")` is default.
#' @param uniprot_features a list of data frames. Each data frame containing Uniprot-annotated features of proteins identified with N-terminal peptides, as resulting from the `feature_to_dataframe` function of the `drawProteins` package.
#'
#'
#' @return a list with 4 elements:
#'         `ntermini_n_processing_features`: data frame of N-terminally identified proteins and their annotated processing features.
#'         `ntermini_category`: data frame of categorized N-terminal peptides according to their matching to annotated processing sites.
#'         `ntermini_category_counts`: data frame of counts of matched features and locations to annotated processing sites.
#'         `ntermini_category_plot`: ggplot2 object. Visualization of counts of matched features and locations.
#'
#'
#' @export
#' @examples
categorize_nterm <- function(annotated_peptides,
                             uniprot_features,
                             distinct = TRUE){

                    require(dplyr)
                    require(ggplot2)

                    # prepare protein_nter data frame

                    protein_semitmt <- annotated_peptides %>%
                                        dplyr::select(protein_id, peptide, nterm, semi_type, specificity,
                                                      is_terminal, last_aa, aa_before, start_position, end_position) %>%
                                        dplyr::filter(str_detect(protein_id,
                                                                 pattern = "Biognosys",
                                                                 negate = TRUE),
                                                      specificity == "semi_specific",
                                                      nterm %in% c("TMT-labelled"))

                    protein_acetyl <- annotated_peptides %>%
                                        dplyr::select(protein_id, peptide, nterm, semi_type, specificity,
                                                      is_terminal, last_aa, aa_before, start_position, end_position) %>%
                                        dplyr::filter(str_detect(protein_id,
                                                                 pattern = "Biognosys",
                                                                 negate = TRUE),
                                                      nterm %in% c("acetylated"))

                    protein_nter <- bind_rows(protein_semitmt,
                                              protein_acetyl)

                    # vector of interesting feature types as annotated by the uniprot API
                    mol_processing_feat <- c("INIT_MET",
                                             "PROPEP",
                                             "SIGNAL",
                                             "TRANSIT")

                    upfeat <- bind_rows(uniprot_features)

                    df_mol_proc_feat <- upfeat %>%
                                        dplyr::filter(type %in% mol_processing_feat, # keep only interesting features
                                               !is.na(length)) %>% # exclude features with missing values
                                        dplyr::rename(protein_id = accession)  # change column name

                    nter_pepts_n_feat <- left_join(protein_nter,
                                                   df_mol_proc_feat,
                                                   by = "protein_id",
                                                   relationship = "many-to-many")

                    # match semi-specific cleavage position of peptide vs end position of annotated processing site

                    init_match <- nter_pepts_n_feat %>%
                                        mutate(matches_start = case_when(semi_type == "semi_Nterm" & abs(as.numeric(start_position) - end) < 4 ~ TRUE,
                                                                         semi_type == "semi_Cterm" & abs(as.numeric(end_position) - end) < 4 ~ TRUE,
                                                                         TRUE ~ FALSE)) %>%
                                        dplyr::filter(!is.na(matches_start)) # eliminate proteins with no processing features annotated

                    # categorize matching locations

                    categ_canon_annot <- init_match %>%
                                        mutate(match_locat = case_when(matches_start == TRUE ~ "start",
                                                                       matches_start == FALSE ~ "none")) %>%
                                        mutate(match_type = case_when(match_locat != "none" & type == "CHAIN" ~ "CHAIN",
                                                                      match_locat != "none" & type == "INIT_MET" ~ "INIT_MET",
                                                                      match_locat != "none" & type == "SIGNAL" ~ "SIGNAL",
                                                                      match_locat != "none" & type == "TRANSIT" ~ "TRANSIT",
                                                                      match_locat != "none" & type == "PROPEP" ~ "PROPEP",
                                                                      match_locat != "none" & type == "PEPTIDE" ~ "PEPTIDE",
                                                                      match_locat == "none" ~ "none"))

                    # select interesting columns

                    categ2_pept_canannot <- categ_canon_annot %>%
                                        dplyr::select(peptide, protein_id,
                                                      match_locat, match_type,
                                                      start_position, end_position,
                                                      begin, end, nterm, specificity, aa_before)

                     # filter to keep one peptide sequence per feature
                     pept_wmatch <- categ2_pept_canannot %>%
                       dplyr::filter(
                         match_type != "none"
                       ) 
                       
                       if(distinct == TRUE){

                         pept_wmatch <- pept_wmatch %>%
                           distinct(
                             peptide,
                             protein_id,
                             match_locat,
                             match_type,
                             .keep_all = TRUE
                           )

                       } else {

                         pept_wmatch <- pept_wmatch
                       
                       }
                                                           
                    pept_womatch <- categ2_pept_canannot %>%
                                        dplyr::filter(!peptide %in% pept_wmatch$peptide) %>%
                                        distinct(peptide, protein_id,
                                                 match_locat, match_type,
                                                 .keep_all = TRUE)

                    categ2_pept_canannot <- bind_rows(pept_wmatch, 
                                                      pept_womatch) %>%
                                        mutate(is_duplicated = duplicated(peptide)) %>%
    # change the definition of the match type
    # to define specific peptides
    mutate(
      match_type = case_when(
        specificity == "specific" ~ "specific_peptide",
        TRUE ~ match_type
      ),
      match_locat = case_when(
        specificity == "specific" ~ "specific_peptide",
        TRUE ~ match_locat
      )
    ) %>%
    mutate(
      match_type = case_when(
        match_type == "INIT_MET" & aa_before == "M" & start_position == 2 ~ "INIT_MET",
        TRUE ~ match_type
      ),
    )

                    # tabular counts of matching locations

                    count_matches <- categ2_pept_canannot %>%
                                        dplyr::count(match_type,
                                                     nterm)

                    # visualize counts of categorized Nterm identifications

                    nterm_categories_plot <- ggplot(count_matches,
                                                    aes(x = reorder(match_type, n), y = n)) +
                                        coord_flip() +
                                        geom_bar(stat = "identity") +
                                        geom_text(aes(label = n), hjust = 1, size = 5) +
                                        facet_grid(nterm~., scales = "free") +
                                        labs(title = "Nr of Nterminal peptides by their annotated category in Uniprot") +
                                        labs(y = "Number of Peptides by matching type",
                                             x = "Type of matching processing information")+
                                        theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                                              axis.text.y = element_text(hjust = 0.5, size = 10),
                                              panel.background = element_blank(),
                                              panel.grid.major = element_line(color = "grey"),
                                              panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),
                                              axis.title=element_text(size=12,face="bold"))

                    # store intermediary outputs into a list object

                    nterm_categorization <- list(ntermini_n_processing_features = nter_pepts_n_feat,
                                                 ntermini_category = categ2_pept_canannot,
                                                 ntermini_category_counts = count_matches,
                                                 ntermini_category_plot = nterm_categories_plot)

                    return(nterm_categorization)

}
