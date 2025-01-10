# This file contains a test for displaying information of a 
# ProteinCodingGene object. The test serves as an example and can 
# be adapted for all the other Gene object types.

testthat::test_that(
  "showGeneObject function works for a ProteinCodingGene object", {
  
  codinggene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                          ranges = IRanges::IRanges(start = 200, end = 1200), 
                          strand = "+")
  
  codinggene_product <- list(protein_id = "proteinID", 
                             protein_sequence = paste0(
                  "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVTELTKVEPPEVPVNGPAAANLLGGQ",
                  "MKKLGFLHSGTAKSVTCTYSPALNKLTENSSEKGSITRRRCFEGIKTM"))
  
  proteincoding_gene <- GenePack::createProteinCodingGene(
                      id = "ENST00005859745",
                      hugo_symbol = "SYMBOL1",
                      name = "protein-coding gene name",
                      description = "gene description",
                      tissue_specificity = list("liver", "small bowel"),
                      gene_structure = codinggene_structure,
                      gene_product = codinggene_product,
                      clinical_significance = "association with bowel cancer")
  
  testthat::expect_output(GenePack::showGeneObject(proteincoding_gene), 
                          "Gene type: ProteinCodingGene")
  testthat::expect_output(GenePack::showGeneObject(proteincoding_gene), 
                          "HUGO Symbol: SYMBOL1")
  testthat::expect_output(GenePack::showGeneObject(proteincoding_gene), 
                          "Tissue Specificity: liver, small bowel")
  testthat::expect_output(GenePack::showGeneObject(proteincoding_gene), 
                          "Gene Product:\n  protein_id: \n  proteinID\n")
})

