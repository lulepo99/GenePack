#' Housekeeping RNA gene class
#'
#' A virtual class to represent housekeeping RNA genes.
#' This class is a general S4 class inheriting from the \code{Gene} 
#' class and serves as a base class to represent specific types of 
#' housekeeping RNA gene. 
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the virtual 
#' \code{HousekeepingRNAGene} class.
#' 
#' @keywords internal


setClass(
  "HousekeepingRNAGene",
  contains = "Gene",
  representation = "VIRTUAL"
)




#' Transfer RNA gene class
#' 
#' A class to represent transfer RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about tRNA genes, 
#' including their tRNA products.
#' 
#' @details The `gene_product` slot is expected to contain a tRNA ID and the 
#' corresponding sequence. The validity function ensures that it is correctly 
#' formatted and contain a valid tRNA sequence.
#' 
#' @slot gene_product list. This slot is specific for the tRNA genes 
#' product and includes: 
#'   \itemize{
#'     \item \code{trna_id}: a string representing the ID of the tRNA.
#'     \item \code{trna_sequence}: a string representing the sequence of 
#'     the tRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{TRNAGene} class.
#' 
#' @export


setClass(
  "TRNAGene",
  contains = "HousekeepingRNAGene",
  prototype = list(
    gene_product = list(
      trna_id = NA_character_,
      trna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$trna_id, 
      object@gene_product$trna_sequence, "tRNA")) {
      return("Invalid tRNA product.")
    }
    TRUE
  }
)




#' Ribosomal RNA gene class
#' 
#' A class to represent ribosomal RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about rRNA genes, 
#' including their rRNA product.
#' 
#' @details The `gene_product` slot is expected to contain a rRNA ID and the 
#' corresponding sequence. The validity function ensures that it is correctly 
#' formatted and contain a valid rRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the rRNA genes
#' product and includes: 
#'   \itemize{
#'     \item \code{rrna_id}: a string representing the ID of the rRNA.
#'     \item \code{rrna_sequence}: a string representing the sequence of the 
#'     rRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{RRNAGene} class.
#' 
#' @export


setClass(
  "RRNAGene",
  contains = "HousekeepingRNAGene",
  prototype = list(
    gene_product = list(
      rrna_id = NA_character_,
      rrna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$rrna_id, 
      object@gene_product$rrna_sequence, "rRNA")) {
      return("Invalid rRNA product.")
    }
    TRUE
  }
)


