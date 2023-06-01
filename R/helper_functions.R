#' Pseudo-random sampling of protein IDs from a character vector
#'
#' @param x character vector of protein IDs
#' @param size numeric. Number of proteins to sample
#' @param seed numeric. Random seed number
#'
#' @return character vector of randomly sampled protein IDs.
#' @export
#'
#' @examples
sample_proteins <- function(x, # char vector of protein IDs
                            size, # number of proteins to sample
                            seed = 363 # random seed number
) {

                    require(stringr)
                    require(magrittr)

                    # select/sample N proteins to observe their abundance distribution
                    set.seed(seed)

                    proteins <- x %>%
                                        unique() %>%
                                        str_remove_all(pattern = "Biognosys") # get unique protein IDs

                    which_prots <- sample(proteins,
                                          size = size,
                                          replace = FALSE)

}
