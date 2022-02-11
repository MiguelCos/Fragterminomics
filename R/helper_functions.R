#' Generate a volcano plot from a `limma::topTable` output
#'
#' @param dataarg data frame. Output from `limma::topTable`. It should also contain at least
#'                            a column called `Protein` with protein IDs.
#' @param FC_cutoff
#' @param pval_cutoff
#' @param color_diffex
#' @param color_nondifex
#' @param interesting_proteins
#' @param vert_line_col
#' @param hline_col
#' @param hline_pos
#' @param linetype
#' @param increased_in
#' @param comparison_title
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
our_volcano <- function(dataarg,
                        FC_cutoff = 1,
                        pval_cutoff = 0.05,
                        color_diffex = "red",
                        color_nondifex = "#2a9d8f",
                        interesting_proteins = NULL,
                        vert_line_col = "red",
                        hline_col = "red",
                        hline_pos = 0.05,
                        linetype = "dashed",
                        #vprot_ylim = c(0,2),
                        #vprot_xlim = c(-2,2),
                        increased_in,
                        comparison_title) {

                    require(stringr)
                    require(magrittr)
                    require(dplyr)

                    signi_hits <- filter(dataarg,
                                         adj.P.Val <= pval_cutoff) %>%
                                        pull(Protein)

                    if(is_null(interesting_proteins)){

                                        protlabels <- dataarg %>%
                                                            filter(adj.P.Val <= pval_cutoff)

                    } else {

                                        protlabels <- dataarg %>%
                                                            filter(Gene %in% interesting_proteins)

                    }


                    volcano <- ggplot(data = dataarg,
                                      mapping = aes(x = logFC, y = -log10(adj.P.Val))) +
                                        geom_point()+
                                        geom_point(data = dataarg %>% filter(logFC > FC_cutoff,
                                                                             adj.P.Val < pval_cutoff),
                                                   mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red")+
                                        geom_point(data = dataarg %>% filter(logFC < -FC_cutoff,
                                                                             adj.P.Val < pval_cutoff),
                                                   mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
                                        geom_point(data = dataarg %>% filter(logFC > FC_cutoff,
                                                                             adj.P.Val > pval_cutoff),
                                                   mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f") +
                                        geom_point(data = dataarg %>% filter(logFC < -FC_cutoff,
                                                                             adj.P.Val > pval_cutoff),
                                                   mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "#2a9d8f") +
                                        ggrepel::geom_text_repel(data = protlabels,
                                                                 aes(label = Gene)) +
                                        geom_hline(yintercept = -log10(hline_pos),
                                                   color = "red", linetype = "dashed") +
                                        ggtitle(paste0("Differentially abundant proteins:\n",comparison_title),
                                                paste0("Positive values = Increased in ", increased_in)) +
                                        labs(caption = paste0("Number of significant hits: ",length(signi_hits)))+
                                        theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 14),
                                              panel.background = element_blank(),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                                              axis.title=element_text(size=12,face="bold"))

                    return(volcano)

}

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
