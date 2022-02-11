#' Apply multiple-testing p-value adjustment on specified features
#'
#' @param toptable a data frame in the same format as the output of the `topTable` function from the `limma` package
#' @param interesting_features_table
#' @param method a string to be passed as a method to the `p.adjust` function for multiple-testing correction. "BH" is default.
#'
#' @return a data frame containing a new column `fdr_correction`, indicating if
#'         the the adjusted p-value comes from
#' @export
#'
#' @examples
feature_fdr_correction <- function(toptable,
                                   interesting_features_table,
                                   method = "BH"){

    require(dplyr)

   # merge limma output with feature annotation  ----

   # the right join allows for both filtering to keep only the interesting
   # features from the toptable and merge with the features table

   tab_limma_feature_annot <- right_join(limma_pept_tab,
                                         interesting_features,
                                         by = c("peptide","index")) %>%
                       # apply FDR correction only on the subset
                       mutate(adj.P.Val = p.adjust(p = P.Value,
                                                   method = "BH")) %>%
                       mutate(fdr_correction = 'feature-specific')

   # create a table of non-interesting features containing the globally adjusted p-values ----
   tablimma_subsetout <- filter(limma_pept_tab,
                                !index %in% interesting_features$index) %>%
                       mutate(fdr_correction = 'global')

   output_limma <- bind_rows(tab_limma_feature_annot,
                             tablimma_subsetout)

   return(output_limma)

}
