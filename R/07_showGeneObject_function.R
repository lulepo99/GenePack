#' Show Gene object information
#'
#' This function displays the gene object information in a more suitable
#' and detailed way with respect to the standard 'show' method 
#' implemented in R. 
#'
#' @param object Gene object. An object of a specific Gene class 
#' (e.g., "ProteinCodingGene", "LncRNAGene", "SiRNAGene", etc.).
#' @return This function does not return a value. It displays the gene object 
#' information to the console.
#' @details The gene_structure slot is an S4 object of class GRanges, which 
#' is implemented in order to give more flexibility in its creation. 
#' To ensure a more standardized and reliable output format, the R \code{show} 
#' function is used within the \code{showGeneObject} function to display it. 
#' Additionally, the product sequence is formmatted into blocks, with each 
#' line containing 80 bases for better readability.
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
#' showGeneObject(lncrna_gene)
#' 
#' @export


setGeneric("showGeneObject", function(object) {
  standardGeneric("showGeneObject")
})


#' @rdname showGeneObject
#' @export

setMethod("showGeneObject", "Gene", function(object) {
  output <- ""
  
  output <- paste0(output, "Gene type: ", class(object), "\n",
                          "Ensembl ID: ", object@id, "\n",
                          "HUGO Symbol: ", object@hugo_symbol, "\n")
  
  if (!is.na(object@name)) {
    output <- paste0(output, "Name: ", object@name, "\n")
  }
  if (!is.na(object@description)) {
    output <- paste0(output, "Description: ", object@description, "\n")
  }
  if (length(object@tissue_specificity) > 0) {
    output <- paste0(output, "Tissue Specificity: ", 
                     paste(object@tissue_specificity, collapse = ", "), "\n")
  }
  
  output <- paste0(output, "Gene Structure:\n")
  cat(output)
  show(object@gene_structure)
  
  output <- "Gene Product:\n"
  for (elem in names(object@gene_product)) {
    output <- paste0(output, "  ", elem, ": ")
    
    sequence <- object@gene_product[[elem]]
    sequence_length <- nchar(sequence)
    output <- paste0(output, "\n")
    
    for (i in seq(1, sequence_length, by = 80)) {
      output <- paste0(output, "  ", 
                       substr(sequence, i, min(i + 79, sequence_length)), "\n")
    }
  }
  
  if (!is.na(object@clinical_significance)) {
    output <- paste0(output, "Clinical Significance: ", 
                     object@clinical_significance, "\n")
  }
  
  output <- paste0(output, "-----------------\n")
  cat(output)
})

