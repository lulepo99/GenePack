#' Gene class
#'
#' A virtual class to represent genes
#' 
#' @slot id character. Ensembl transcript ID of the gene. In this 
#' implementation, each transcript isoform is considered as an individual gene. 
#' Hence, the Ensembl transcript ID ("ENST") is used instead of the Ensembl 
#' gene ID ("ENSG"). This approach simplifies the management of Gene objects, 
#' avoiding the complexity of handling multiple gene products of different 
#' types within a single Gene object.
#' @slot hugo_symbol character. Hugo symbol of the gene.
#' @slot name character. Optional parameter. Name of the gene.
#' @slot description character. Optional parameter. 
#' Description of the gene.
#' @slot tissue_specificity list. Optional parameter. List of tissues 
#' (character) where the gene is specifically expressed. 
#' @slot gene_structure GenomicRanges::GRanges. Structure of the gene, 
#' including chromosomes ("seqnames", character), start position 
#' ("start", numeric), end position ("end", numeric), 
#' strand ("strand", character). 
#' @slot gene_product list. It represents the product of the gene and 
#' contains its product id (character) and the corresponding 
#' sequence (character).
#' @slot clinical_significance character. Optional parameter. Clinical 
#' relevance of the gene or its product (e.g., association with a disease 
#' or therapeutic targets). 
#' 
#' @details The validity function ensures that
#' \itemize{
#'   \item The \code{id} follows the Ensembl transcript format 
#'   \code{^ENST[0-9]+$}.
#'   \item The \code{hugo_symbol} slot contains only uppercase letters, digits, 
#'   and special characters '-' or '_'.
#'   \item The \code{name}, and \code{description} are non-empty strings, 
#'   if specified.
#'   \item The \code{tissue_specificity} slot is a list of tissues, 
#'   if specified.
#'   \item The \code{gene_structure} slot is a \code{GRanges} object with 
#'   valid chromosome, and strand. The \code{IRanges}
#'   object in \code{GRanges} directly verify that 'start' is less than or 
#'   equal to 'end' when the user creates the GRanges object, so this check is 
#'   not implemented in the package.
#'   \item The \code{gene_product} slot is a list that represents a gene 
#'   product and contains two elements: an ID and a sequence.
#'   \item The \code{clinical_significance} is a non-empty string, 
#'   if specified. 
#' }
#' 
#' @note The \code{gene_structure} and \code{gene_product} arguments are 
#' complex objects and should be created separately by the user before 
#' constructing the Gene object. This approach maximizes flexibility, allowing 
#' the user to better exploit the functionalities provided by the 
#' \code{GenomicRanges} package and to easily reuse the \code{gene_structure} 
#' object in other contexts.  
#' 
#' @return This documentation describes the structure of the virtual 
#' \code{Gene} class.
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @keywords classes
#' @export


setClass(
  "Gene",
  slots = list(
    id = "character",
    hugo_symbol = "character",
    name = "character",
    description = "character",
    tissue_specificity = "list",
    gene_structure = "GRanges",
    gene_product = "list",
    clinical_significance = "character"
  ),
  contains = "VIRTUAL",
  prototype = list(
    id = NA_character_,
    hugo_symbol = NA_character_,
    name = NA_character_,
    description = NA_character_,
    tissue_specificity = list(),
    gene_structure = GenomicRanges::GRanges(),  
    gene_product = list(),
    clinical_significance = NA_character_
  ),
  validity = function(object) {
    if (!is.character(object@id) || !grepl("^ENST[0-9]+$", object@id)) {
      return("Invalid Ensembl transcript ID format. It should start with 'ENST' followed by digits.")
    }
    if (!is.character(object@hugo_symbol) || !grepl("^[A-Z0-9_-]+$", object@hugo_symbol)) {
      return("Invalid HUGO symbol format. It should contain only uppercase letters and digits.")
    }
    if (!is.na(object@name)) {
      if (!is.character(object@name) || nchar(object@name) == 0) {
        return("The 'name' slot must be a non-empty string, if specified.")
      }
    }
    if (!is.na(object@description)) {
      if (!is.character(object@description) || nchar(object@description) == 0) {
        return("The 'description' slot must be a non-empty string, if specified.")
      }
    }
    if (length(object@tissue_specificity) > 0) {
      if (!all(vapply(object@tissue_specificity, is.character, FUN.VALUE = logical(1)))) {
        return("Each tissue in the 'tissue_specificity' slot must be a character string.")
      }
    }
    if (!inherits(object@gene_structure, "GRanges")) {
      return("The 'gene_structure' slot must be a GRanges object.")
    }
    if (!as.character(seqnames(object@gene_structure)) %in% paste0("chr", c(as.character(seq_len(22)), "X", "Y"))) {
      return("Chromosome must start with 'chr' and be between '1' and '22', or 'X' or 'Y'.")
    }
    if (!all(as.character(strand(object@gene_structure)) %in% c("+", "-"))) {
      return("Strand must be '+' or '-'")
    }
    if (!is.list(object@gene_product) || length(object@gene_product) != 2) {
      return("The gene_product argument must be a list with the product ID and the corresponding sequence.")
    }
    if (!is.na(object@clinical_significance)) {
      if (!is.character(object@clinical_significance) || nchar(object@clinical_significance) == 0) {
        return("The clinical_significance argument must be a non-empty string, if specified.")
      }
    }
    TRUE
  }
)
