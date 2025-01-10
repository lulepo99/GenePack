#' Get the gene product ID
#'
#' This function retrieves the gene product ID.
#' 
#' @param object Gene object
#' @return the gene product ID.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. The gene product ID is 
#' stored as the second element in the "gene_product" list.
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
#' getProductID(lncrna_gene)
#' 
#' @export


setGeneric("getProductID", function(object) {
  standardGeneric("getProductID")
})


#' @rdname getProductID
#' @export

setMethod("getProductID", "Gene", function(object) {
  object@gene_product[[1]]
})




#' Set the gene product ID
#'
#' This function sets the gene product ID.
#' 
#' @param object Gene object.
#' @param value the new product ID of the gene product.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. After setting the new 
#' product ID, the function checks that the Gene object is still 
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
#' setProductID(lncrna_gene) <- "newlncRNAID"
#' 
#' getProductID(lncrna_gene)    
#' 
#' @export


setGeneric("setProductID<-", function(object, value) {
  standardGeneric("setProductID<-")
})


#' @rdname setProductID-set
#' @export

setMethod("setProductID<-", "Gene", function(object, value) {
  object@gene_product[[1]] <- value
  validObject(object)
  object
})




#' Get the gene product sequence
#'
#' This function retrieves the gene product sequence.
#' 
#' @param object Gene object
#' @return the gene product sequence.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. The gene product sequence 
#' is stored as the second element in the "gene_product" list. 
#' The function is implemented to provide a more readable output, displaying 
#' the sequence in blocks of 80 bases per line.
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
#' getProductSequence(lncrna_gene)
#' 
#' @export


setGeneric("getProductSequence", function(object) {
  standardGeneric("getProductSequence")
})


#' @rdname getProductSequence
#' @export

setMethod("getProductSequence", "Gene", function(object) {
  sequence <- object@gene_product[[2]]  
  
  sequence_length <- nchar(sequence)
  
  for (i in seq(1, sequence_length, by = 80)) {
    cat(substr(sequence, i, min(i + 79, sequence_length)), "\n")
  }
  
  invisible(sequence)
})




#' Set the gene product sequence
#'
#' This function sets the gene product sequence.
#' 
#' @param object Gene object.
#' @param value the new product sequence of the gene product.
#' @return the modified Gene object.
#' @details This function is defined as a generic function to be applicable 
#' to any object that inherits from the Gene class. After setting the new 
#' product sequence, the function checks that the Gene object 
#' is still valid by calling \code{validObject}. 
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
#' setProductSequence(lncrna_gene) <- paste0(
#'                      "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCG",
#'                      "CGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCG",
#'                      "GCUAAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGAU",
#'                      "GGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGG",
#'                      "AUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUC",
#'                      "CUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC")
#' 
#' getProductSequence(lncrna_gene)    
#' 
#' @export


setGeneric("setProductSequence<-", function(object, value) {
  standardGeneric("setProductSequence<-")
})


#' @rdname setProductSequence-set
#' @export

setMethod("setProductSequence<-", "Gene", function(object, value) {
  object@gene_product[[2]] <- value
  validObject(object)
  object
})
