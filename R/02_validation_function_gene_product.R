#' Validate Gene Product ID
#' 
#' This function validates the gene product ID to ensure it is a 
#' non-empty string.
#' 
#' @param id character. The gene product ID to validate.
#' @return TRUE if the gene product ID is valid, throws an error otherwise.
#' 
#' @keywords internal 


validateID <- function(id) {
  if (!is.character(id) || nchar(id) == 0) {
    stop("Gene product ID must be a non-empty string.")
  }
  TRUE
}




#' Validate Gene Product Sequence
#' 
#' This function validates the gene product sequence, ensuring it is a 
#' string containing only valid nucleotide bases or amino acids, 
#' depending on the Gene object class. Besides, it checks if the gene 
#' product sequence meets specific length requirements based on its type.
#' 
#' @param sequence character. The gene product sequence to validate.
#' @param type character. The gene product type (e.g., "lncRNA", "sncRNA", 
#' "tRNA").
#' @return TRUE if the gene product sequence is valid, throws an error 
#' otherwise.
#' 
#' @details Regarding the length requirements, the function allows some 
#' flexibility due to the variable length ranges defined for some gene types. 
#' A long non-coding RNA gene product must be at least 200 nucleotides
#' in length. All small non coding genes are considered together and their 
#' product must be shorter than 200 nucleotides. A tRNA gene product ranges 
#' from 70 to 90 nucleotides, in general. No checks are defined for rRNA gene 
#' products because of their high variability in nucleotide sequences instead.
#' 
#' @return TRUE if the gene product length is valid, throws an error otherwise.
#' 
#' @keywords internal


validateSequence <- function(sequence, type) {
  if (!is.character(sequence)) {
    stop("Gene product sequence must be a string.")
  }
  
  if (type == "protein") {
    amino_acids <- "ACDEFGHIKLMNPQRSTVWY"
    if (substr(sequence, 1, 1) != "M") {
      stop("Protein sequence must start with Methionine (M).")
    }
    if (!grepl(paste0("^[", amino_acids, "]+$"), sequence)) {
      stop("Protein sequence contains invalid amino acids.")
    }
  } else {
    if (!grepl("^[ACGUNacgun]+$", sequence)) {
      stop("Gene product sequence can only contain the bases A, C, G, U, and", 
           "N (case insensitive).")
    }
    if (type == "lncRNA" && nchar(sequence) < 200) {
      stop("The lncRNA sequence must be at least 200 nucleotides in length.")
    } else if (type == "sncRNA" && nchar(sequence) >= 200) {
      stop("The small non-coding RNA sequence must be shorter than", 
           "200 nucleotides.")
    } else if (type == "tRNA" && (nchar(sequence) < 70 || 
                                  nchar(sequence) > 90)) {
      stop("The tRNA sequence must be between 70 and 90", 
           "nucleotides in length.")
    }
  }
  
  TRUE
}




#' Validate Gene Product
#' 
#' This function validates both the gene product ID and sequence 
#' using internal methods. Users can use it to validate the
#' gene products before creating the gene object. However, the function
#' is also internal in the definition of the gene object classes to 
#' guarantee robustness.
#' 
#' @param id character. The gene product ID to validate.
#' @param sequence character. The gene product sequence to validate.
#' @param type character. The RNA type: "lncRNA", "sncRNA" (including 
#' microRNA, piRNA, snRNA, snoRNA, siRNA), "tRNA", 
#' "rRNA", "protein".
#' @return TRUE if both the ID and sequence are valid, throws an error 
#' otherwise.
#' 
#' @keywords gene product validation function
#' 
#' @examples
#' lncrna_id <- "lncRNAID"
#' 
#' lncrna_sequence <- paste0(
#'                    "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCG",
#'                    "GCUAAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGAU",
#'                    "CGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCG",
#'                    "AUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUC",
#'                    "GGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGG",
#'                    "CUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC") 
#' type <- "lncRNA"
#' 
#' validateGeneProduct(lncrna_id, lncrna_sequence, type)
#' 
#' 
#' @export


validateGeneProduct <- function(id, sequence, type) {
  validateID(id)
  validateSequence(sequence, type)
  TRUE
}
