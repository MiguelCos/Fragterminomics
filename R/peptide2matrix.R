#' Convert a character vector of peptides to a matrix of amino acids and count the occurrence of each amino acid in each position
#'
#' This function takes a character vector of peptides, where each element represents a peptide of length 20, and converts it into a matrix where each row represents a peptide and each column represents an amino acid position within that peptide. Each element of the matrix is a single character representing the amino acid at that position. Additionally, the function counts the number of times a particular amino acid appears in a particular position and generates a new matrix, in which each column is an amino acid position, and each row represents an amino acid. The intersection represents the number of amino acids that showed up in that position.
#'
#' @param peptides A character vector of peptides
#' @param p_number An integer value. How many elements before and after the cleavage site to show.
#'
#' @return A list containing two elements:
#' - peptide_matrix: a matrix where each row represents a peptide and each column represents an amino acid position within that peptide.
#' Each element of the matrix is a single character representing the amino acid at that position.
#' - amino_acid_count: a matrix where each row represents an amino acid and each column represents an amino acid position within the peptide.
#' The intersection represents the number of times that amino acid showed up in that position.
#'
#' @examples
#' peptides <- c("AKKGSKRHYVWALDMTFEH", "VAVFDELLKIVPNLGSYKR")
#' peptide_matrix_count(peptides)
#'
#' @import purrr
#' @import tidyr
#' @import dplyr
#' @importFrom stringr str_split
#' @export
#'
peptide_matrix_count <- function(peptides,
                                 p_number = 5) {
                    # Split each peptide into individual amino acids using map and str_split
                    split_peptides <- purrr::map(peptides,
                                                 str_split,
                                                 "")

                    # Convert the list of split peptides into a data frame using map_dfc
                    peptide_matrix <- map_dfc(split_peptides,
                                              as.data.frame) %>%
                                        # Transpose the data frame so that each row corresponds to a peptide
                                        t() %>%
                                        suppressMessages()

                    # create character vector that goes from P10 to P1 and then P1' to P10'
                    # assign it to the positions object
                    positions <- c(paste0("P",10:1),
                                   paste0("P",1:10,"'"))

                    # define the order in which the amino acids will be in the matrix
                    # Define the order of amino acids for sorting
                    amino_acid_order <- c("A", "C", "D", "E", "F", "G", "H",
                                          "I", "K", "L", "M", "N", "P", "Q",
                                          "R", "S", "T", "V", "W", "Y")

                    # use the positions object to define the column names of the peptide sequence matrix
                    colnames(peptide_matrix) <- positions

                    rownames(peptide_matrix) <- NULL

                    # Count how many times a particular amino acid appears in a certain position
                    amino_acid_count <- peptide_matrix %>%
                                        as_tibble() %>%
                                        pivot_longer(cols = everything(),
                                                     names_to = "position",
                                                     values_to = "AA") %>%
                                        group_by(position,
                                                 AA) %>%
                                        summarize(n = n(),
                                                  .groups = "drop") %>%
                                        complete(position,
                                                 AA,
                                                 fill = list(n = 0)) %>%
                                        pivot_wider(id_cols = AA,
                                                    names_from = position,
                                                    values_from = n)

                    # Arrange the rows based on the order of the amino acids in the vector above
                    amino_acid_count <- amino_acid_count %>%
                                        as.data.frame() %>%
                                        arrange(factor(rownames(.),
                                                       levels = amino_acid_order)) %>%
                                        column_to_rownames("AA")

                    # Select only columns between P5 and P5' from the count matrix
                    amino_acid_count <- amino_acid_count[,
                                                         c(paste0("P", p_number:1),
                                                           paste0("P", 1:p_number, "'"))] %>%
                                        as.matrix()

                    # Return the final matrices as a list
                    return(list(peptide_matrix = peptide_matrix,
                                amino_acid_count = amino_acid_count))
}
