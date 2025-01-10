#' Compute the length of a gene product for a gene object
#' 
#' This function computes the length of the product of a gene object. 
#' It provides an integer value that represents the number of nucleotides
#' or amino acids in the sequence, depending on the Gene object class.
#' 
#' @param object Gene object. An object of a specific Gene class 
#' (e.g., "ProteinCodingGene", "LncRNAGene", "SiRNAGene", etc.).
#' @return an integer representing the length of the gene product.
#' @details Based on the class of the Gene object, the function choose the 
#' appropriate method to compute the length of the gene product.
#' 
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
#'                   ranges = IRanges::IRanges(start = 200, end = 1200),
#'                   strand = "+")
#'                          
#' gene_product <- list(lncrna_id = "lncRNAID", 
#'                      lncrna_sequence = paste0(
#'                        "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCG",
#'                        "GCUAAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGAU",
#'                        "CGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCG",
#'                        "AUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUC",
#'                        "GGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGG",
#'                        "CUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC"))
#'                                           
#' lncrna_gene <- createLncRNAGene(id = "ENST000001",
#'                hugo_symbol = "SYMBOL1",
#'                name = "lncRNA gene name",
#'                description = "gene description",
#'                tissue_specificity = list("liver", "small bowel"),
#'                gene_structure = gene_structure,
#'                gene_product = gene_product,
#'                clinical_significance = "association with disease")
#' 
#' lengthProduct(lncrna_gene)
#' 
#' @keywords length, gene, sequence, product
#' 
#' @export


setGeneric("lengthProduct", function(object) {
  standardGeneric("lengthProduct")
})


#' @rdname lengthProduct
#' @export

setMethod("lengthProduct", "Gene", function(object) {
  gene_class <- class(object)
  if (gene_class == "ProteinCodingGene") {
    return(nchar(object@gene_product$protein_sequence))
  } else if (gene_class == "LncRNAGene") {
    return(nchar(object@gene_product$lncrna_sequence))
  } else if (gene_class == "MicroRNAGene") {
    return(nchar(object@gene_product$microrna_sequence))
  } else if (gene_class == "SiRNAGene") {
    return(nchar(object@gene_product$sirna_sequence))
  } else if (gene_class == "PiRNAGene") {
    return(nchar(object@gene_product$pirna_sequence))
  } else if (gene_class == "SnRNAGene") {
    return(nchar(object@gene_product$snrna_sequence))
  } else if (gene_class == "SnoRNAGene") {
    return(nchar(object@gene_product$snorna_sequence))
  } else if (gene_class == "rRNAGene") {
    return(nchar(object@gene_product$rrna_sequence))
  } else if (gene_class == "tRNAGene") {
    return(nchar(object@gene_product$trna_sequence))
  } else {
    stop("Unknown Gene object")
  }
})
