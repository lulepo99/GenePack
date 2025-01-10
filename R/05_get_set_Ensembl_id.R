#' Get the Ensembl transcript ID
#'
#' This function retrieves the Ensembl transcript ID  associated with 
#' a Gene object.
#' 
#' @param object Gene object
#' @return the Ensembl transcript ID of the gene.
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
#' getID(lncrna_gene)
#' 
#' @export 



setGeneric("getID", function(object) {
  standardGeneric("getID")
})


#' @rdname getID                                  
#' @export 

setMethod("getID", "Gene", function(object) {
  object@id
})




#' Set the Ensembl transcript ID
#'
#' This function sets the Ensembl transcript ID for a Gene object.
#' 
#' @param object Gene object.
#' @param value the new Ensembl transcript ID.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. After setting the new ID, 
#' the function checks that the Gene object is still valid 
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
#' setID(lncrna_gene) <- "ENST000002"
#' 
#' getID(lncrna_gene)                                 
#' 
#' @export 


setGeneric("setID<-", function(object, value) {
  standardGeneric("setID<-")
})


#' @rdname setID-set
#' @export 

setMethod("setID<-", "Gene", function(object, value) {
  validateID(value)
  object@id <- value
  validObject(object)
  object
})