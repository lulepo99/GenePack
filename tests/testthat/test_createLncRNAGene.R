# Test the creation of a LncRNAGene object with its related slots. 

testthat::test_that("createLncRNAGene creates a valid LncRNAGene object", {
  
  lncrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                          ranges = IRanges::IRanges(start = 200, end = 1200), 
                          strand = "+")
  
  lncrnagene_product <- list(lncrna_id = "lncRNAID", 
                lncrna_sequence = paste0(
                "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUA",
                "CGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUAC",
                "GAUCGAUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUCGGCUUAGCGUACG",
                "GUAGGCUUAACUGCGGACGAUCGAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGA",
                "UCGAUCGC"))
  
  lncrna_gene <- GenePack::createLncRNAGene(id = "ENST00005859745",
              hugo_symbol = "SYMBOL1",
              name = "lncRNA gene name",
              description = "gene description",
              tissue_specificity = list("liver", "small bowel"),
              gene_structure = lncrnagene_structure,
              gene_product = lncrnagene_product,
              clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(lncrna_gene, "LncRNAGene")
  
  testthat::expect_type(lncrna_gene@id, "character")
  testthat::expect_type(lncrna_gene@hugo_symbol, "character")
  testthat::expect_type(lncrna_gene@name, "character")
  testthat::expect_type(lncrna_gene@description, "character")
  testthat::expect_type(lncrna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(lncrna_gene@gene_structure, "GRanges")
  testthat::expect_type(lncrna_gene@gene_product$lncrna_id, "character")
  testthat::expect_type(lncrna_gene@gene_product$lncrna_sequence, "character")
  testthat::expect_type(lncrna_gene@clinical_significance, "character")
  
})

