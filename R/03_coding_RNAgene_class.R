#' Protein-coding RNA gene class
#' 
#' A class to represent protein-coding genes. 
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about protein-coding 
#' genes, including their protein product. 
#' 
#' @details The `gene_product` slot is expected to contain a protein ID and 
#' the corresponding sequence. The validity function ensures that it is 
#' correctly formatted and contain a valid protein sequence.
#' 
#' @slot gene_product list. This slot is specific for the 
#' protein-coding genes product and includes: 
#'   \itemize{
#'     \item \code{protein_id}: a string representing the ID of the protein.
#'     \item \code{protein_sequence}: a string representing the sequence of 
#'     the protein.
#'   }
#' 
#' @return This documentation describes the structure of the 
#' \code{ProteinCodingGene} class.
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @export


setClass(
  "ProteinCodingGene",
  contains = "Gene",
  prototype = list(
    gene_product = list(
      protein_id = NA_character_,
      protein_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$protein_id, 
      object@gene_product$protein_sequence, "protein")) {
      return("Invalid protein product.")
    }
    TRUE
  }
)