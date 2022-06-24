#' Annotate peptides based on their enzyme specificity and position within protein sequence.
#'
#' @param peptide2protein A data frame with at least two columns: `Peptide` (peptide sequence to annotate) and `Genes` (Uniprot ID of associated protein).
#' @param fasta A list of identified protein sequences, at least matching the ones present in the `peptide2protein` data frame. This list should should be loaded into R with the `read.fasta` function from the `seqinr` package.
#' @param decoy_tag A string with the tag of the decoy sequences.
#' @param specificity A string defining the enzyme specificity.
#'
#' @return a data frame with annotation information of the input peptides.
#' @export
#'
#' @examples
annotate_peptides <- function(peptide2protein,
                              fasta,
                              decoy_tag = "^rev",
                              specificity = "R|K"){

require(dplyr)
require(purrr)
require(stringr)
require(magrittr)

pept2prot <- peptide2protein %>%
                    dplyr::select(Peptide, Genes) %>%
                    mutate(Peptide = as.character(Peptide),
                           Genes = as.character(Genes))

peptide_sequence <- pept2prot$Peptide


if (!is.null(decoy_tag)){
                    fasta <- fasta %>%
                             purrr::discard(str_detect(names(.), decoy_tag))
}

annotation <- data.frame()

for (i in 1:length(peptide_sequence)){

                    # extract the list element that matches the peptide sequence
                    list_elem <- fasta %>%
                                 purrr::keep(str_detect(names(.),
                                             stringr::fixed(pept2prot$Genes[i])))


                    if(is_empty(list_elem)){

                                        protein_id <- NA
                                        protein_description <- NA
                                        protein_length <- NA
                                        Peptide <- peptide_sequence[i]
                                        peptide_position <- NA
                                        start_position <- NA
                                        end_position <- NA
                                        following_10_resid <- NA
                                        previous_10_resid <- NA
                                        previous_all_resid <- NA
                                        aa_after <- NA
                                        aa_before <- NA
                                        last_aa <- NA

                                        temprow <- cbind(protein_id,
                                                         protein_description,
                                                         protein_length,
                                                         Peptide,
                                                         start_position,
                                                         end_position,
                                                         last_aa, aa_after,
                                                         aa_before,
                                                         following_10_resid,
                                                         previous_10_resid,
                                                         previous_all_resid)

                                        annotation <- rbind(annotation,
                                                            temprow)

                                        row.names(annotation) <- NULL

                    } else {

                                        # extract the annotation info stored in the fasta header, check for an Uniprot fasta pattern and extract protein ID and description
                                        if(str_detect(string = attr(list_elem[[1]], "Annot"),
                                                      pattern = "^\\>sp\\|") | str_detect(string = attr(list_elem[[1]],"Annot"),
                                                                                                pattern = "^\\>tr\\|")){

                                                            split1 <- str_split_fixed(attr(list_elem[[1]], "Annot"),
                                                                                pattern = "\\|",
                                                                                n = 3)

                                                            protein_id <- split1[,2] # get protein id for 'uniprot' IDs

                                                            split2 <- str_split_fixed(split1[,3],
                                                                                pattern = "\\s",
                                                                                n = 2)

                                                            protein_description <- str_extract(split2[,2],
                                                                                pattern = ".+?(?=OS=)") %>%
                                                                                str_trim() # get protein description for non-uniprot IDs
                                        } else {

                                                            split1 <- str_split_fixed(attr(list_elem[[1]],
                                                                                                 "Annot"),
                                                                                pattern = "\\s",
                                                                                n = 2)

                                        protein_id  <- split1[,1] # get protein ID

                                        protein_description <- split1[,2] # get protein description for non Uniprot IDs
                    }

                                        # extract protein sequence from the extracted list element

                                        protein_seq <- list_elem[[1]][1] # get the matched protein sequence

                                        protein_length <- str_length(protein_seq) # get the matched protein length

                                        Peptide <- peptide_sequence[i] # get the sequence of the peptide

                                        peptide_position <- str_locate(protein_seq,
                                                                                pattern = Peptide) # get the position of the peptide within the matched protein sequence

                                        start_position <- peptide_position[,1] # start position of the peptide within the protein

                                        end_position <- peptide_position[,2] # end positionj of the peptide within the protein

                                        following_10_resid <- str_sub(protein_seq,
                                                                      start = peptide_position[,2]+1,
                                                                      end = peptide_position[,2]+10) # 10 residues after the end of the peptide sequence

                                        previous_10_resid <- str_sub(protein_seq,
                                                                     # 10 residues before the start of the peptide sequence
                                                                     start = ifelse(peptide_position[,1]-10 > 0,
                                                                                    peptide_position[,1]-10, 1),
                                                                     end = ifelse(peptide_position[,1]-1 > 0,
                                                                                  peptide_position[,1]-1, 0))

                                        previous_all_resid <- str_sub(protein_seq,
                                                                      # all residues before up to annotated protein start
                                                                      start = 2,
                                                                      end = peptide_position[,1]-1)

                                        aa_after <- str_sub(protein_seq,
                                                            # amino acid after
                                                            start = peptide_position[,2]+1,
                                                            end = peptide_position[,2]+1)

                                        aa_before <- str_sub(protein_seq,
                                                             # amino acid before
                                                             start = peptide_position[,1]-1,
                                                             end = peptide_position[,1]-1)


                                        last_aa <- str_sub(Peptide,
                                                           # last amino acid of the identified peptide
                                                           start = str_count(Peptide),
                                                           end = str_count(Peptide))


                                        temprow <- cbind(protein_id,
                                                         protein_description,
                                                         protein_length,
                                                         Peptide,
                                                         start_position,
                                                         end_position,
                                                         last_aa,
                                                         aa_after,
                                                         aa_before,
                                                         following_10_resid,
                                                         previous_10_resid,
                                                         previous_all_resid)

                                        annotation <- rbind(annotation,
                                                            temprow)

                                        row.names(annotation) <- NULL

                                        print(paste(i, peptide_sequence[i], "out of", length(peptide_sequence)))
}
                    }

annotation <-  annotation %>%
                    mutate(semi_type = case_when(str_detect(last_aa, specificity) & str_detect(aa_before, specificity, negate = TRUE) ~ "semi_Nterm",
                                                 str_detect(last_aa, specificity, negate = TRUE) & str_detect(aa_before, specificity) ~ "semi_Cterm",
                                                 str_detect(last_aa, specificity, negate = TRUE) & str_detect(aa_before, specificity, negate = TRUE) ~ "unspecific",
                                                 TRUE ~ "specific"),
                           general_position = case_when(as.numeric(start_position) <= 2 ~ "Nterm_peptide",
                                                        end_position == protein_length ~ "Cterm_peptide",
                                                        TRUE ~ "not_terminal")) %>%
                    mutate(specificity = case_when(semi_type == "semi_Nterm" & general_position == "Nterm_peptide" ~ "specific",
                                                   semi_type == "semi_Cterm" & general_position == "Cterm_peptide" ~ "specific",
                                                   semi_type == "semi_Cterm" & general_position == "not_terminal" ~ "semi_specific",
                                                   semi_type == "semi_Nterm" & general_position == "not_terminal" ~ "semi_specific",
                                                   general_position == "Cterm_peptide" & str_detect(aa_before, specificity) ~ "specific",
                                                   general_position == "Cterm_peptide" & str_detect(aa_before, specificity, negate = TRUE) ~ "semi_specific",
                                                   general_position == "Nterm_peptide" & str_detect(last_aa, specificity) ~ "specific",
                                                   general_position == "Nterm_peptide" & str_detect(last_aa, specificity, negate = TRUE) ~ "semi_specific",
                                                   semi_type == "specific" ~ "specific"),
                           is_terminal = if_else(general_position == "not_terminal",
                                                 true = "not_terminal",
                                                 false = "terminal")) %>%
                    mutate(semi_type = case_when(semi_type == "unspecific" & general_position == "Nterm_peptide" ~ "semi_Nterm",
                                                 semi_type == "unspecific" & general_position == "Cterm_peptide" ~ "semi_Cterm",
                                                 TRUE ~ semi_type))

return(annotation)
}
