#' Get the gene description 
#'
#' This function retrieves the description of the gene.
#' 
#' @param object Gene object
#' @return the description of the gene.
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
#' getDescription(lncrna_gene)
#' 
#' @export


setGeneric("getDescription", function(object) {
  standardGeneric("getDescription")
})


#' @rdname getDescription                       
#' @export

setMethod("getDescription", "Gene", function(object) {
  object@description
})




#' Set the gene description
#'
#' This function sets the description of the gene.
#' 
#' @param object Gene object.
#' @param value the new description of the gene.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. After setting the new 
#' description, the function checks that the Gene object is still 
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
#'                clinical_significance = "association with disease")
#'                                  
#' setDescription(lncrna_gene) <- "new gene description"
#' 
#' getDescription(lncrna_gene)                                 
#' 
#' @export


setGeneric("setDescription<-", function(object, value) {
  standardGeneric("setDescription<-")
})


#' @rdname setDescription-set
#' @export

setMethod("setDescription<-", "Gene", function(object, value) {
  object@description <- value
  validObject(object)
  object
})