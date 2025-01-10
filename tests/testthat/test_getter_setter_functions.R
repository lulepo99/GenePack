# This file contains tests for the getter and setter functions associated 
# with the arguments of Gene objects created 
# using the GenePack package. The arguments include id, hugo_symbol, 
# name, species, tissue_specificity, 
# gene_structure, gene_product, and clinical_significance. The file also 
# shows how the functions handle errors.

testthat::test_that("Getter and Setter functions work correctly", {
  
  gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                    ranges = IRanges::IRanges(start = 200, end = 1200),
                    strand = "+")
  
  gene_product <- list(lncrna_id = "lncRNAID", 
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
                 gene_structure = gene_structure,
                 gene_product = gene_product,
                 clinical_significance = "association with small bowel cancer")
  
  # Testing for getID and setID.
  
  expect_equal(GenePack::getID(lncrna_gene), "ENST00005859745")
  
  GenePack::setID(lncrna_gene) <- "ENST00005859750"
  expect_equal(GenePack::getID(lncrna_gene), "ENST00005859750")
  
  testthat::expect_error(GenePack::setID(lncrna_gene) <- "EN00005859751a",      
    paste("Invalid Ensembl transcript ID format.", 
          "It should start with 'ENST' followed by digits."))
  
  
  # Testing for getHugoSymbol and setHugoSymbol.
  
  expect_equal(GenePack::getHugoSymbol(lncrna_gene), "SYMBOL1")
  
  GenePack::setHugoSymbol(lncrna_gene) <- "NEWSYMBOL1"
  expect_equal(GenePack::getHugoSymbol(lncrna_gene), "NEWSYMBOL1")
  
  testthat::expect_error(GenePack::setHugoSymbol(
    lncrna_gene) <- "Newsymbol1",      
    paste("Invalid HUGO symbol format.", 
          "It should contain only uppercase letters and digits."))
  
  
  # Testing for getName and setName.
  
  expect_equal(GenePack::getName(lncrna_gene), "lncRNA gene name")
  
  GenePack::setName(lncrna_gene) <- "newgenename"
  expect_equal(GenePack::getName(lncrna_gene), "newgenename")
  
  testthat::expect_error(GenePack::setName(lncrna_gene) <- "",      
    "The 'name' slot must be a non-empty string, if specified.")
  
  
  # Testing for getDescription and setDescription.
  
  expect_equal(GenePack::getDescription(lncrna_gene), "gene description")
  
  GenePack::setDescription(lncrna_gene) <- "new gene description"
  expect_equal(GenePack::getDescription(lncrna_gene), "new gene description")
  
  testthat::expect_error(GenePack::setDescription(lncrna_gene) <- "",      
    "The 'description' slot must be a non-empty string, if specified.")
  
  
  # Testing for getTissues and setTissues.
  
  expect_equal(GenePack::getTissues(lncrna_gene), list("liver", "small bowel"))
  
  GenePack::setTissues(lncrna_gene) <- list("liver", "large bowel")
  expect_equal(GenePack::getTissues(lncrna_gene), list("liver", "large bowel"))
  
  testthat::expect_error(GenePack::setTissues(lncrna_gene) <- list(1, 2),      
    "Each tissue in the 'tissue_specificity' slot must be a character string.")
  
  
  # Testing for getGeneStructure and setGeneStructure.
  
  expect_equal(GenePack::getGeneStructure(lncrna_gene), gene_structure)
  
  new_gene_structure1 <- GenomicRanges::GRanges(seqnames = "chr1",
                         ranges = IRanges(start = 205, end = 1210),
                         strand = "+")
  
  GenePack::setGeneStructure(lncrna_gene) <- new_gene_structure1
  expect_equal(GenePack::getGeneStructure(lncrna_gene), new_gene_structure1)
  
  new_gene_structure2 <- GenomicRanges::GRanges(seqnames = "1",
                         ranges = IRanges(start = 205, end = 1210),
                         strand = "+")
  
  testthat::expect_error(GenePack::setGeneStructure(lncrna_gene
                            ) <- new_gene_structure2,      
        paste("Chromosome must start with 'chr' and be between '1'", 
              "and '22', or 'X' or 'Y'."))
  
  
  # Testing for gene product: getProductID and setProductID.
  
  expect_equal(GenePack::getProductID(lncrna_gene), gene_product[[1]])
  
  GenePack::setProductID(lncrna_gene) <- "newlncRNAID"
  expect_equal(GenePack::getProductID(lncrna_gene), "newlncRNAID")
  
  testthat::expect_error(GenePack::setProductID(lncrna_gene) <- "",      
                         "Gene product ID must be a non-empty string.")
  
  
  # Testing for gene product: getProductSequence and setProductSequence.
  # capture.output is used to capture the output produced by the cat function, 
  # which is normally displayed in the console.
  
  sequence <- capture.output({
    GenePack::getProductSequence(lncrna_gene)
  })
  sequence <- gsub("\\s+", "", paste(sequence, collapse = ""))
  expect_equal(sequence, gsub("\\s+", "", gene_product[[2]]))
  
  seq <- paste0(
    "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUACGGUAGG",
    "CUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUUAGCGUAC",
    "GUAGGCUCAACUGCGUACGAUCGAUCGGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGGCUUAG",
    "CGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC")
  
  GenePack::setProductSequence(lncrna_gene) <- seq
  new_sequence <- capture.output({ 
    GenePack::getProductSequence(lncrna_gene)
  })
  new_sequence <- gsub("\\s+", "", paste(new_sequence, collapse = ""))
  expect_equal(new_sequence, gsub("\\s+", "", seq))
  
  testthat::expect_error(GenePack::setProductSequence(lncrna_gene
                                            ) <- "AUGCUUAGCGUACGGUAGGC",
    "The lncRNA sequence must be at least 200 nucleotides in length."
  )

  
  # Testing for getClinicalSignificance and setClinicalSignificance
  
  expect_equal(GenePack::getClinicalSignificance(lncrna_gene), 
               "association with small bowel cancer")
  
  GenePack::setClinicalSignificance(lncrna_gene
                        ) <- "association with large bowel cancer"
  expect_equal(GenePack::getClinicalSignificance(lncrna_gene), 
               "association with large bowel cancer")
  
  testthat::expect_error(GenePack::setClinicalSignificance(
    lncrna_gene) <- "",      
    paste("The clinical_significance argument must be", 
          "a non-empty string, if specified."))
  
})
  
  
  
  