# This file contains tests for the creation of Gene objects belonging to 
# specialized classes derived from the SncRNAGene virtual class. 
# The tested classes include RRNAGene, and TRNAGene.

# Testing the creation of a RRNAGene object with its related slots. 

testthat::test_that("createRRNAGene creates a valid RRNAGene object", {
  
  rrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                        ranges = IRanges::IRanges(start = 200, end = 1200), 
                        strand = "+")
  
  rrnagene_product <- list(rrna_id = "rrnaID", 
              rrna_sequence = paste0("GAGCCCGGCGUCCCGUCUCAGACCCGGCCGACAGGAG",
                                     "GGGUGGAGUGGGUGUGUGCGUG",
                                     "UGGAGGAGGGUGUGGUGGUGCG",
                                     "GUGUGGUGGGCGUGUGAGUGUG",
                                     "UGGCGCGUGGAGGAG"))
  
  rrna_gene <- GenePack::createRRNAGene(id = "ENST00005859745",
                      hugo_symbol = "SYMBOL1",
                      name = "rRNA gene name",
                      description = "gene description",
                      tissue_specificity = list("liver", "small bowel"),
                      gene_structure = rrnagene_structure,
                      gene_product = rrnagene_product,
                      clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(rrna_gene, "RRNAGene")
  
  testthat::expect_type(rrna_gene@id, "character")
  testthat::expect_type(rrna_gene@hugo_symbol, "character")
  testthat::expect_type(rrna_gene@name, "character")
  testthat::expect_type(rrna_gene@description, "character")
  testthat::expect_type(rrna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(rrna_gene@gene_structure, "GRanges")
  testthat::expect_type(rrna_gene@gene_product$rrna_id, "character")
  testthat::expect_type(rrna_gene@gene_product$rrna_sequence, "character")
  testthat::expect_type(rrna_gene@clinical_significance, "character")
  
})


# Testing the creation of a TRNAGene object with its related slots. 

testthat::test_that("createTRNAGene creates a valid TRNAGene object", {
  
  trnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                        ranges = IRanges::IRanges(start = 200, end = 1200), 
                        strand = "+")
  
  trnagene_product <- list(trna_id = "trnaID", 
                           trna_sequence = paste0(
                             "GGCUACCCUGUUCGAACCUCGAUUAACACAGGCUUCUUUU",
                             "CGACUGGUCGAACCCUGACAACCUCGAUUAACACAG"))
  
  trna_gene <- GenePack::createTRNAGene(id = "ENST00005859746",
                  hugo_symbol = "SYMBOL1",
                  name = "tRNA gene name",
                  description = "gene description",
                  tissue_specificity = list("liver", "small bowel"),
                  gene_structure = trnagene_structure,
                  gene_product = trnagene_product,
                  clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(trna_gene, "TRNAGene")
  
  testthat::expect_type(trna_gene@id, "character")
  testthat::expect_type(trna_gene@hugo_symbol, "character")
  testthat::expect_type(trna_gene@name, "character")
  testthat::expect_type(trna_gene@description, "character")
  testthat::expect_type(trna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(trna_gene@gene_structure, "GRanges")
  testthat::expect_type(trna_gene@gene_product$trna_id, "character")
  testthat::expect_type(trna_gene@gene_product$trna_sequence, "character")
  testthat::expect_type(trna_gene@clinical_significance, "character")
  
})

