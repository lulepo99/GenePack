# Testing the creation of a ProteinCodingGene object with its related slots. 

testthat::test_that(
  "createProteinCodingGene creates a valid ProteinCodingGene object", {
  
  codinggene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                          ranges = IRanges::IRanges(start = 200, end = 1200), 
                          strand = "+")
  
  codinggene_product <- list(protein_id = "proteinID", 
            protein_sequence = paste0(
              "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVTELTKVEPPEVPVNGPAAANLL",
              "GGQMKKLGFLHSGTAKSVTCTYSPALNKLTENSSEKGSITRRRCFEGIKTM"))
  
  proteincoding_gene <- GenePack::createProteinCodingGene(
                        id = "ENST00005859745",
                        hugo_symbol = "SYMBOL1",
                        name = "protein-coding gene name",
                        description = "gene description",
                        tissue_specificity = list("liver", "small bowel"),
                        gene_structure = codinggene_structure,
                        gene_product = codinggene_product,
                        clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(proteincoding_gene, "ProteinCodingGene")
  
  testthat::expect_type(proteincoding_gene@id, "character")
  testthat::expect_type(proteincoding_gene@hugo_symbol, "character")
  testthat::expect_type(proteincoding_gene@name, "character")
  testthat::expect_type(proteincoding_gene@description, "character")
  testthat::expect_type(proteincoding_gene@tissue_specificity, "list")
  testthat::expect_s4_class(proteincoding_gene@gene_structure, "GRanges")
  testthat::expect_type(proteincoding_gene@gene_product$protein_id, "character")
  testthat::expect_type(proteincoding_gene@gene_product$protein_sequence, 
                        "character")
  testthat::expect_type(proteincoding_gene@clinical_significance, "character")
  
})

