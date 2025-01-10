#' Small non-coding RNA gene class
#'
#' A virtual class to represent small non-coding RNA genes.
#' This class is a general S4 class inheriting from the \code{Gene} class 
#' and serves as a base class to represent specific types of non small 
#' non-coding RNA gene. 
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the virtual 
#' \code{SncRNAGene} class.
#' 
#' @keywords internal


setClass(
  "SncRNAGene",
  contains = "Gene",
  representation = "VIRTUAL"
)




#' MicroRNA gene class
#' 
#' A class to represent microRNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about microRNA genes, 
#' including their microRNA product.
#' 
#' @details The `gene_product` slot is expected to contain a microRNA ID and 
#' the corresponding sequence. The validity function ensures that it is 
#' correctly formatted and contain a valid microRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the microRNA genes 
#' product and includes:
#'   \itemize{
#'     \item \code{microrna_id}: a string representing the ID of the microRNA.
#'     \item \code{microrna_sequence}: a string representing the sequence of 
#'     the microRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{MicroRNAGene} class.
#' 
#' @export


setClass(
  "MicroRNAGene",
  contains = "SncRNAGene",
  prototype = list(
    gene_product = list(
      microrna_id = NA_character_,
      microrna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$microrna_id, 
      object@gene_product$microrna_sequence, "sncRNA")) {
      return("Invalid microRNA product.")
    }
    TRUE
  }
)




#' Short-interfering RNA gene class
#' 
#' A class to represent short-interfering RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about siRNA genes, 
#' including their siRNA product.
#' 
#' @details The `gene_product` slot is expected to contain a siRNA ID and the 
#' corresponding sequence. The validity function ensures that it is correctly 
#' formatted and contain a valid siRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the siRNA genes 
#' product and includes:
#'   \itemize{
#'     \item \code{sirna_id}: a string representing the ID of the siRNA.
#'     \item \code{sirna_sequence}: a string representing the sequence of 
#'     the siRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{SiRNAGene} class.
#' 
#' @export


setClass(
  "SiRNAGene",
  contains = "SncRNAGene",
  prototype = list(
    gene_product = list(
      sirna_id = NA_character_,
      sirna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$sirna_id, 
      object@gene_product$sirna_sequence, "sncRNA")) {
      return("Invalid siRNA product.")
    }
    TRUE
  }
)




#' Piwi-interacting RNA gene class
#' 
#' A class to represent piwi-interacting RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about piRNA genes, including 
#' their piRNA products.
#' 
#' @details The `gene_product` slot is expected to contain a piRNA ID and the 
#' corresponding sequence. The validity function ensures that it is correctly 
#' formatted and contain a valid piRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the piRNA genes 
#' product and includes:
#'   \itemize{
#'     \item \code{pirna_id}: a string representing the ID of the piRNA.
#'     \item \code{pirna_sequence}: a string representing the sequence of 
#'     the piRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{PiRNAGene} class.
#' 
#' @export


setClass(
  "PiRNAGene",
  contains = "SncRNAGene",
  prototype = list(
    gene_product = list(
      pirna_id = NA_character_,
      pirna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$pirna_id, 
      object@gene_product$pirna_sequence, "sncRNA")) {
      return("Invalid piRNA product.")
    }
    TRUE
  }
)




#' Small nuclear RNA gene class
#' 
#' A class to represent small nuclear RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene} 
#' class and is designed to store information about snRNA genes, 
#' including their snRNA products.
#' 
#' @details The `gene_product` slot is expected to contain a snRNA ID and 
#' the corresponding sequence. The validity function ensures that it is 
#' correctly formatted and contain a valid snRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the snRNA genes 
#' product and includes:
#'   \itemize{
#'     \item \code{snrna_id}: a string representing the ID of the snRNA.
#'     \item \code{snrna_sequence}: a string representing the sequence of 
#'     the snRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{SnRNAGene} class.
#' 
#' @export


setClass(
  "SnRNAGene",
  contains = "SncRNAGene",
  prototype = list(
    gene_product = list(
      snrna_id = NA_character_,
      snrna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$snrna_id, 
      object@gene_product$snrna_sequence, "sncRNA")) {
      return("Invalid snRNA product.")
    }
    TRUE
  }
)




#' Small nucleolar RNA gene class
#' 
#' A class to represent small nucleolar RNA genes.
#' This class is a specialized S4 class inheriting from the \code{Gene}
#' class and is designed to store information about snoRNA genes,
#' including their snoRNA products.
#' 
#' @details The `gene_product` slot is expected to contain a snoRNA ID and 
#' the corresponding sequence. The validity function ensures that it is 
#' correctly formatted and contain a valid snoRNA sequence. 
#' 
#' @slot gene_product list. This slot is specific for the snoRNA genes 
#' product and includes:
#'   \itemize{
#'     \item \code{snorna_id}: a string representing the ID of the snoRNA.
#'     \item \code{snorna_sequence}: a string representing the sequence of 
#'     the snoRNA.
#'   }
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @return This documentation describes the structure of the 
#' \code{SnoRNAGene} class.
#' 
#' @export


setClass(
  "SnoRNAGene",
  contains = "SncRNAGene",
  prototype = list(
    gene_product = list(
      snorna_id = NA_character_,
      snorna_sequence = NA_character_
    )
  ),
  validity = function(object) {
    if (!GenePack::validateGeneProduct(object@gene_product$snorna_id, 
      object@gene_product$snorna_sequence, "sncRNA")) {
      return("Invalid snoRNA product.")
    }
    TRUE
  }
)