#' Create a Small Nuclear RNA Gene object
#'
#' This function creates a \code{SnRNAGene} object, which represents 
#' a gene that encodes for small nuclear RNAs. The SnRNAGene object 
#' includes base information inherited from the \code{Gene} class 
#' with the specific information regarding its snRNA product.
#' 
#' @param id character. Ensembl transcript ID of the gene. In this 
#' implementation, each transcript isoform
#' is considered as an individual gene. Hence, the Ensembl transcript ID 
#' ("ENST") is used instead of the 
#' Ensembl gene ID ("ENSG"). This approach simplifies the management 
#' of Gene objects, avoiding the 
#' complexity of handling multiple gene products of different types within 
#' a single Gene object.
#' @param hugo_symbol character. Hugo symbol of the gene.
#' @param name character. Name of the gene.
#' @param description character. Description of the gene.
#' @param tissue_specificity list. Optional parameter. List of tissues 
#' (character) where the gene is specifically 
#' expressed. If not specified, its default will be "-".
#' @param gene_structure GenomicRanges::GRanges. Structure of the gene, 
#' including chromosomes ("seqnames", character),
#' start position ("start", numeric), end position ("end", numeric), 
#' strand ("strand", character).
#' @param gene_product list. It represents the product of the gene and 
#' contains its product id (character) and 
#' the corresponding sequence (character).
#' @param clinical_significance character. Optional parameter. 
#' Clinical relevance of the gene or its product 
#' (e.g., association with a disease or therapeutic targets). 
#' If not specified, its default will be "-".
#' 
#' @return A \code{SnRNAGene} object.
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
#'   \item The \code{gene_structure} slot is a \code{GRanges} object 
#'   with valid chromosome, and strand. The \code{IRanges}
#'   object in \code{GRanges} directly verify that 'start' is less 
#'   than or equal to 'end' when the user creates the GRanges
#'   object, so this check is not implemented in the package.
#'   \item The \code{gene_product} slot is a list that represents a gene 
#'   product and contains
#'   two elements: the snrna_id and the snrna_sequence.
#'   \item The \code{clinical_significance} is a non-empty string, 
#'   if specified. 
#' }
#' 
#' @note The \code{gene_structure} and \code{gene_product} arguments are 
#' complex objects and should be created 
#' separately by the user before constructing the Gene object. 
#' This approach maximizes flexibility, allowing 
#' the user to better exploit the functionalities provided by the 
#' \code{GenomicRanges} package and to easily 
#' reuse the \code{gene_structure} in other contexts. 
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' 
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
#'                   ranges = IRanges::IRanges(start = 200, end = 1200),
#'                   strand = "+")
#'                          
#' gene_product <- list(snrna_id = "snRNAID", 
#'                      snrna_sequence = paste0(
#'                          "AUGCGGAUUUACGCCUAACGUUCCG",
#'                          "GAUUACCUCGGAUUCCGA"))
#' 
#' snrna_gene <- createSnRNAGene(id = "ENST000001",
#'               hugo_symbol = "SYMBOL1",
#'               name = "snRNA gene name",
#'               description = "gene description",
#'               tissue_specificity = list("liver", 
#'                                         "small bowel"),
#'               gene_structure = gene_structure,
#'               gene_product = gene_product,
#'               clinical_significance = "association with disease")
#' 
#' @export


createSnRNAGene <- function(id, hugo_symbol, name = NA_character_, 
                            description = NA_character_, 
                            tissue_specificity = list(), 
                            gene_structure, gene_product, 
                            clinical_significance = NA_character_) {
  new("SnRNAGene", id = id, hugo_symbol = hugo_symbol, name = name, 
      description = description, 
      tissue_specificity = tissue_specificity, gene_structure = gene_structure, 
      gene_product = gene_product, 
      clinical_significance = clinical_significance)
}
