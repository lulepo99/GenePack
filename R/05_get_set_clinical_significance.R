#' Get the clinical significance
#'
#' This function retrieves the clinical significance of the gene.
#' 
#' @param object Gene object
#' @return the clinical significance of the gene, or "-" if it was not 
#' specified at the creation.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class.
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
#' getClinicalSignificance(lncrna_gene)
#' 
#' @export


setGeneric("getClinicalSignificance", function(object) {
  standardGeneric("getClinicalSignificance")
})


#' @rdname getClinicalSignificance
#' @export

setMethod("getClinicalSignificance", "Gene", function(object) {
  object@clinical_significance
})




#' Set the clinical significance
#'
#' This function sets the clinical significance of the gene.
#' 
#' @param object Gene object.
#' @param value the new clinical significance.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. After setting the new 
#' clinical significance, the function checks that the Gene object is still 
#' valid by calling \code{validObject}.
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
#'                clinical_significance = "association with disease1")
#'  
#' setClinicalSignificance(lncrna_gene) <- "association with disease2"
#' 
#' getClinicalSignificance(lncrna_gene)                                 
#' 
#' @export


setGeneric("setClinicalSignificance<-", function(object, value) {
  standardGeneric("setClinicalSignificance<-")
})


#' @rdname setClinicalSignificance-set
#' @export

setMethod("setClinicalSignificance<-", "Gene", function(object, value) {
  object@clinical_significance <- value
  validObject(object)
  object
})
