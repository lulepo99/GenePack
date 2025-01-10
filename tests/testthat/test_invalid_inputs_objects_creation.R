# This file contains tests that show how the create functions handle most 
# of invalid inputs for each attribute 
# during the gene objects creation. A LncRNAGene object is chosen 
# as an example for all types of gene objects.

testthat::test_that("createLncRNAGene handle correctly input errors", {
  
  
  # Testing for invalid transcript ID
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "INVALID_ID", 
          hugo_symbol = "SYMBOL1", 
          name = "lncRNA gene name",
          description = "gene description", 
          tissue_specificity = list("liver", "small bowel"),
          gene_structure = lncrnagene_structure, 
          gene_product = lncrnagene_product,
          clinical_significance = "association with bowel cancer")},
    paste("Invalid Ensembl transcript ID format.", 
          "It should start with 'ENST' followed by digits."))
  
  
  # Testing for invalid Hugo Symbol
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
          hugo_symbol = "Invalidsymbol", 
          name = "lncRNA gene name",
          description = "gene description", 
          tissue_specificity = list("liver", "small bowel"),
          gene_structure = lncrnagene_structure, 
          gene_product = lncrnagene_product,
          clinical_significance = "association with bowel cancer")},
    paste("Invalid HUGO symbol format.", 
          "It should contain only uppercase letters and digits."))
  
  
  # Testing for invalid gene name
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
                hugo_symbol = "SYMBOL1", 
                name = "",
                description = "gene description", 
                tissue_specificity = list("liver", "small bowel"),
                gene_structure = lncrnagene_structure, 
                gene_product = lncrnagene_product,
                clinical_significance = "association with bowel cancer")},
    "The 'name' slot must be a non-empty string, if specified.")
  
  
  # Testing for invalid gene description
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
            hugo_symbol = "SYMBOL1", 
            name = "lncRNA gene name",
            description = "", 
            tissue_specificity = list("liver", "small bowel"),
            gene_structure = lncrnagene_structure, 
            gene_product = lncrnagene_product,
            clinical_significance = "association with bowel cancer")},
    "The 'description' slot must be a non-empty string, if specified.")
  
  
  # Testing for invalid tissue_specificity argument
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
                hugo_symbol = "SYMBOL1", 
                name = "lncRNA gene name",
                description = "gene description", 
                tissue_specificity = list(1, 2),
                gene_structure = lncrnagene_structure, 
                gene_product = lncrnagene_product,
                clinical_significance = "association with bowel cancer")},
    "Each tissue in the 'tissue_specificity' slot must be a character string.")
  
  
  # Testing for invalid gene_structure argument (invalid 'seqnames')
  
  testthat::expect_error({
    
    lncrnagene_structure <- GenomicRanges::GRanges(seqnames = "1",
                            ranges = IRanges::IRanges(start = 200, end = 1200), 
                            strand = "+")
    
    lncrnagene_product <- list(lncrna_id = "lncRNAID", 
                               lncrna_sequence = paste0(
            "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUA",
            "CGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUAC",
            "GAUCGAUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUCGGCUUAGCGUACG",
            "GUAGGCUUAACUGCGGACGAUCGAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGA",
            "UCGAUCGC"))
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
          hugo_symbol = "SYMBOL1", 
          name = "lncRNA gene name",
          description = "gene description", 
          tissue_specificity = list("liver", "small bowel"),
          gene_structure = lncrnagene_structure, 
          gene_product = lncrnagene_product,
          clinical_significance = "association with bowel cancer")},
    paste("Chromosome must start with 'chr' and be between '1'", 
          "and '22', or 'X' or 'Y'."))
  
  
  # Testing for invalid gene_product argument, example 1 (invalid product ID)
  
  testthat::expect_error({
    
    lncrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                            ranges = IRanges::IRanges(start = 200, end = 1200), 
                            strand = "+")
    
    lncrnagene_product <- list(lncrna_id = "", 
                               lncrna_sequence = paste0(
          "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUA",
          "CGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUAC",
          "GAUCGAUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUCGGCUUAGCGUACG",
          "GUAGGCUUAACUGCGGACGAUCGAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGA",
          "UCGAUCGC"))
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
              hugo_symbol = "SYMBOL1", 
              name = "lncRNA gene name",
              description = "gene description", 
              tissue_specificity = list("liver", "small bowel"),
              gene_structure = lncrnagene_structure, 
              gene_product = lncrnagene_product,
              clinical_significance = "association with bowel cancer")},
    "Gene product ID must be a non-empty string.")
  
  
  # Testing for invalid gene_product argument, example 2 
  # (invalid product sequence)
  
  testthat::expect_error({
    
    lncrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                            ranges = IRanges::IRanges(start = 200, end = 1200), 
                            strand = "+")
    
    lncrnagene_product <- list(lncrna_id = "lncRNAID", 
                               lncrna_sequence = paste0(
              "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCG",
              "CGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUG"))
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
              hugo_symbol = "SYMBOL1", 
              name = "lncRNA gene name",
              description = "gene description", 
              tissue_specificity = list("liver", "small bowel"),
              gene_structure = lncrnagene_structure, 
              gene_product = lncrnagene_product,
              clinical_significance = "association with bowel cancer")},
    "The lncRNA sequence must be at least 200 nucleotides in length.")
    
  
  # Testing for invalid clinical_significance argument
  
  testthat::expect_error({
    
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
    
    GenePack::createLncRNAGene(id = "ENST00005859745", 
            hugo_symbol = "SYMBOL1", 
            name = "lncRNA gene name",
            description = "gene description", 
            tissue_specificity = list("liver", "small bowel"),
            gene_structure = lncrnagene_structure, 
            gene_product = lncrnagene_product,
            clinical_significance = "")},
    paste("The clinical_significance argument must be", 
          "a non-empty string, if specified."))

})  
  
  
  
  
  
  
  