# This file contains tests for computing the length of the product of a 
# gene object. It shows both an example in which the length is successfully 
# computed and an example where the function throws an error because of an
# unknown gene object.


# Testing the computation of the length of the product of a 
# ProteinCodingGene object.

testthat::test_that(
  "lengthProduct function works for a ProteinCodingGene object", {
  
  codinggene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                          ranges = IRanges::IRanges(start = 200, end = 1200), 
                          strand = "+")

  codinggene_product <- list(protein_id = "proteinID", 
                              protein_sequence = paste0(
            "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVTELTKVEPPEVPVNGP",
            "AAANLLGGQMKKLGFLHSGTAKSVTCTYSPALNKLTENSSEKGSITRRRCFEGIKTM"))

  proteincoding_gene <- GenePack::createProteinCodingGene(
                          id = "ENST00005859745",
                          hugo_symbol = "SYMBOL1",
                          name = "protein-coding gene name",
                          description = "gene description",
                          tissue_specificity = list("liver", "small bowel"),
                          gene_structure = codinggene_structure,
                          gene_product = codinggene_product,
                          clinical_significance = "association with disease1")

  computed_length <- GenePack::lengthProduct(proteincoding_gene)
  expected_length <- nchar(codinggene_product$protein_sequence)

  testthat::expect_equal(computed_length, expected_length)
})


# Testing how the function handles the presence of an unknown object.

testthat::test_that(
  "lengthProduct function throws an error with an unknown object", {
  
  unknown_object <-  list(unknown_id = "unknownID", 
                           unknown_sequence = 
                          "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVTELTKVEPPEVPVNGPAA")
  
  testthat::expect_error(GenePack::lengthProduct(unknown_object))
})

