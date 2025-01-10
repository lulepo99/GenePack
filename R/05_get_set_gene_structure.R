#' Get the gene structure 
#'
#' This function retrieves the gene_structure GRanges object.
#' 
#' @param object Gene object
#' @return the gene_structure object.
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
#' getGeneStructure(lncrna_gene)
#' 
#' @export


setGeneric("getGeneStructure", function(object) {
  standardGeneric("getGeneStructure")
})


#' @rdname getGeneStructure
#' @export

setMethod("getGeneStructure", "Gene", function(object) {
  object@gene_structure
})




#' Set the gene structure
#'
#' This function sets the gene_structure GRanges object.
#' 
#' @param object Gene object.
#' @param value the new gene_structure object.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable to 
#' any object that inherits from the Gene class. After setting the new 
#' gene structure, the function checks that the Gene object is still valid 
#' by calling \code{validObject}. This ensures that all the internal attributes 
#' of the GRanges object are checked according to the validity conditions 
#' specified in the virtual Gene class.
#' 
#' @examples 
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
#'                ranges = IRanges::IRanges(start = 200, end = 1200),
#'                strand = "+")
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
#' new_gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
#'              ranges = IRanges::IRanges(start = 205, end = 1210),
#'              strand = "+")
#'                                  
#' setGeneStructure(lncrna_gene) <- new_gene_structure
#' 
#' getGeneStructure(lncrna_gene)            
#' 
#' @export


setGeneric("setGeneStructure<-", function(object, value) {
  standardGeneric("setGeneStructure<-")
})


#' @rdname setGeneStructure-set
#' @export

setMethod("setGeneStructure<-", "Gene", function(object, value) {
  object@gene_structure <- value
  validObject(object)
  object
})