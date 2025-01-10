#' Get the Hugo symbol
#'
#' This function retrieves the HUGO symbol associated with a Gene object.
#' 
#' @param object Gene object
#' @return the HUGO symbol of the gene.
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
#' getHugoSymbol(lncrna_gene)
#' 
#' @export 


setGeneric("getHugoSymbol", function(object) {
  standardGeneric("getHugoSymbol")
})


#' @rdname getHugoSymbol
#' @export 

setMethod("getHugoSymbol", "Gene", function(object) {
  object@hugo_symbol
})




#' Set the Hugo symbol
#'
#' This function sets the HUGO symbol for a Gene object.
#' 
#' @param object Gene object.
#' @param value the new HUGO symbol.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable to 
#' any object that inherits from the Gene class. After setting the new 
#' HUGO symbol, the function checks that the Gene object is still valid 
#' by calling \code{validObject}.
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
#' setHugoSymbol(lncrna_gene) <- "NEWSYMBOL1"
#' 
#' getHugoSymbol(lncrna_gene)     
#' 
#' @export 


setGeneric("setHugoSymbol<-", function(object, value) {
  standardGeneric("setHugoSymbol<-")
})


#' @rdname setHugoSymbol-set
#' @export 

setMethod("setHugoSymbol<-", "Gene", function(object, value) {
  object@hugo_symbol <- value
  validObject(object)  
  object
})