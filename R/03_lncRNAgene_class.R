#' Long non-coding RNA Gene class
#' 
#' A class to represent long non-coding RNA genes
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about long non-coding RNA 
#' genes, including their lncRNA product.
#' 
#' @details The `gene_product` slot is expected to contain a lncRNA ID and the 
#' corresponding sequence. The validity function ensures that it is correctly 
#' formatted and contain a valid lncRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the long non-coding 
#' RNA genes product and includes:
#'   \itemize{
#'     \item \code{lncrna_id}: a string representing the ID of the lncRNA.
#'     \item \code{lncrna_sequence}: a string representing the sequence of 
#'     the lncRNA.
#'   }
#' 
#' @return This documentation describes the structure of the 
#' \code{LncRNAGene} class.
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @export



setClass(
  "LncRNAGene",
  contains = "Gene",
  prototype = list(
    gene_product = list(
      lncrna_id = NA_character_,
      lncrna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$lncrna_id, 
      object@gene_product$lncrna_sequence, "lncRNA")) {
      return("Invalid lncRNA product.")
    }
    TRUE
  }
)
