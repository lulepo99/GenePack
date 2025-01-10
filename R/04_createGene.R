#' Create a Gene object
#'
#' This function creates a \code{Gene} object, serving as a virtual base 
#' class for the creation of more specific gene types. The \code{Gene}  
#' class is not designed to be instantiated directly by users. It forms  
#' the basic structure from which all the specialized gene classes 
#' inherit.
#' 
#' @param id character. Ensembl transcript ID of the gene. In this 
#' implementation, each transcript isoform is considered as an individual gene. 
#' Hence, the Ensembl transcript ID ("ENST") is used instead of the 
#' Ensembl gene ID ("ENSG"). This approach simplifies the management of 
#' Gene objects, avoiding the complexity of handling multiple gene products 
#' of different types within a single Gene object.
#' @param hugo_symbol character. Hugo symbol of the gene.
#' @param name character. Optional parameter. Name of the gene.
#' @param description character. Optional parameter. Description of the gene.
#' @param tissue_specificity list. Optional parameter. List of tissues 
#' (character) where the gene is specifically expressed. 
#' @param gene_structure GenomicRanges::GRanges. Structure of the gene, 
#' including chromosomes ("seqnames", character), start position 
#' ("start", numeric), end position ("end", numeric), 
#' strand ("strand", character). 
#' @param clinical_significance character. Optional parameter. 
#' Clinical relevance of the gene or its product (e.g., association with 
#' a disease or therapeutic targets). 
#' 
#' @return A \code{Gene} object.
#' 
#' @details The validity function within the class definition ensures that
#' \itemize{
#'   \item The \code{id} follows the Ensembl transcript format 
#'   \code{^ENST[0-9]+$}.
#'   \item The \code{hugo_symbol} slot contains only uppercase letters, 
#'   digits, and special characters '-' or '_'.
#'   \item The \code{name}, and \code{description} are non-empty strings, 
#'   if specified.
#'   \item The \code{tissue_specificity} slot is a list of tissues, 
#'   if specified.
#'   \item The \code{gene_structure} slot is a \code{GRanges} object with 
#'   valid chromosome, and strand. The \code{IRanges} object in \code{GRanges} 
#'   directly verify that 'start' is less than or equal to 'end' when the user 
#'   creates the GRanges object, so this check is not implemented 
#'   in the package.
#'   \item The \code{clinical_significance} is a non-empty string, 
#'   if specified. 
#' }
#' 
#' @note The \code{gene_structure} argument is 
#' complex objects and should be created separately by the user before 
#' constructing the Gene object. This approach maximizes flexibility, allowing 
#' the user to better exploit the functionalities provided by the 
#' \code{GenomicRanges} package and to easily 
#' reuse the \code{gene_structure} object in other contexts. 
#' 
#' @keywords internal


createGene <- function(id, hugo_symbol, name = NA_character_, 
                  description = NA_character_, tissue_specificity = list(), 
                  gene_structure, clinical_significance = NA_character_) {
  new("Gene", id = id, hugo_symbol = hugo_symbol, name = name, 
      description = description, 
      tissue_specificity = tissue_specificity, gene_structure = gene_structure, 
      clinical_significance = clinical_significance)
}
