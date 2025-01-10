# This file contains tests for the creation of Gene objects belonging to 
#specialized classes derived from the 
# SncRNAGene virtual class. The tested classes include 
# MicroRNAGene, SiRNAGene, PiRNAGene, SnRNAGene, 
# and SnoRNAGene.  


# Testing the creation of a MicroRNAGene object with its related slots.

testthat::test_that("createMicroRNAGene creates a valid MicroRNAGene object", {
  
  micrornagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                            ranges = IRanges::IRanges(start = 200, end = 1200), 
                            strand = "+")
  
  micrornagene_product <- list(microrna_id = "micrornaID", 
                               microrna_sequence = "UGAGGUAGUAGGUUGUAUAGUU")
  
  microrna_gene <- GenePack::createMicroRNAGene(id = "ENST00005859745",
                    hugo_symbol = "SYMBOL1",
                    name = "microRNA gene name",
                    description = "gene description",
                    tissue_specificity = list("liver", "small bowel"),
                    gene_structure = micrornagene_structure,
                    gene_product = micrornagene_product,
                    clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(microrna_gene, "MicroRNAGene")
  
  testthat::expect_type(microrna_gene@id, "character")
  testthat::expect_type(microrna_gene@hugo_symbol, "character")
  testthat::expect_type(microrna_gene@name, "character")
  testthat::expect_type(microrna_gene@description, "character")
  testthat::expect_type(microrna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(microrna_gene@gene_structure, "GRanges")
  testthat::expect_type(microrna_gene@gene_product$microrna_id, "character")
  testthat::expect_type(microrna_gene@gene_product$microrna_sequence, 
                        "character")
  testthat::expect_type(microrna_gene@clinical_significance, "character")
  
})


# Testing the creation of a SiRNAGene object with its related slots.

testthat::test_that("createSiRNAGene creates a valid SiRNAGene object", {
  
  sirnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                         ranges = IRanges::IRanges(start = 200, end = 1200), 
                         strand = "+")
  
  sirnagene_product <- list(sirna_id = "sirnaID", 
                            sirna_sequence = "GGAUACUGGGAUCUAGAGAUU")
  
  sirna_gene <- GenePack::createSiRNAGene(id = "ENST00005859746",
                hugo_symbol = "SYMBOL1",
                name = "siRNA gene name",
                description = "gene description",
                tissue_specificity = list("liver", "small bowel"),
                gene_structure = sirnagene_structure,
                gene_product = sirnagene_product,
                clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(sirna_gene, "SiRNAGene")
  
  testthat::expect_type(sirna_gene@id, "character")
  testthat::expect_type(sirna_gene@hugo_symbol, "character")
  testthat::expect_type(sirna_gene@name, "character")
  testthat::expect_type(sirna_gene@description, "character")
  testthat::expect_type(sirna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(sirna_gene@gene_structure, "GRanges")
  testthat::expect_type(sirna_gene@gene_product$sirna_id, "character")
  testthat::expect_type(sirna_gene@gene_product$sirna_sequence, "character")
  testthat::expect_type(sirna_gene@clinical_significance, "character")
  
})


# Testing the creation of a PiRNAGene object with its related slots.

testthat::test_that("createPiRNAGene creates a valid PiRNAGene object", {
  
  pirnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                         ranges = IRanges::IRanges(start = 200, end = 1200), 
                         strand = "+")
  
  pirnagene_product <- list(pirna_id = "pirnaID", 
                            pirna_sequence = "AAGCAGAGCUUUUUGCCAGGUGCUAAAGUC")
  
  pirna_gene <- GenePack::createPiRNAGene(id = "ENST00005859747",
                  hugo_symbol = "SYMBOL1",
                  name = "piRNA gene name",
                  description = "gene description",
                  tissue_specificity = list("liver", "small bowel"),
                  gene_structure = pirnagene_structure,
                  gene_product = pirnagene_product,
                  clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(pirna_gene, "PiRNAGene")
  
  testthat::expect_type(pirna_gene@id, "character")
  testthat::expect_type(pirna_gene@hugo_symbol, "character")
  testthat::expect_type(pirna_gene@name, "character")
  testthat::expect_type(pirna_gene@description, "character")
  testthat::expect_type(pirna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(pirna_gene@gene_structure, "GRanges")
  testthat::expect_type(pirna_gene@gene_product$pirna_id, "character")
  testthat::expect_type(pirna_gene@gene_product$pirna_sequence, "character")
  testthat::expect_type(pirna_gene@clinical_significance, "character")
  
})


# Testing the creation of a SnRNAGene object with its related slots.

testthat::test_that("createSnRNAGene creates a valid SnRNAGene object", {
  
  snrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                         ranges = IRanges::IRanges(start = 200, end = 1200), 
                         strand = "+")
  
  snrnagene_product <- list(snrna_id = "snrnaID", 
                snrna_sequence = paste0(
    "GGGAGACUGACUGUCAUUUGAAGACGAACACUGAGGGACCUUCGUGCAGGACGCCU",
    "GCUGACCCUGCAGCGGCUCCUGCCCAAGGAGGUUUGCAGGUGGCAGGAUUGGAGGU",
    "GCUAGGUGCCCGGGAGGAGACAGGCAGGGGAGUUGGGGGCGGGAGCGCAGGGUGGG"))
  
  snrna_gene <- GenePack::createSnRNAGene(id = "ENST00005859748",
                hugo_symbol = "SYMBOL1",
                name = "snRNA gene name",
                description = "gene description",
                tissue_specificity = list("liver", "small bowel"),
                gene_structure = snrnagene_structure,
                gene_product = snrnagene_product,
                clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(snrna_gene, "SnRNAGene")
  
  testthat::expect_type(snrna_gene@id, "character")
  testthat::expect_type(snrna_gene@hugo_symbol, "character")
  testthat::expect_type(snrna_gene@name, "character")
  testthat::expect_type(snrna_gene@description, "character")
  testthat::expect_type(snrna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(snrna_gene@gene_structure, "GRanges")
  testthat::expect_type(snrna_gene@gene_product$snrna_id, "character")
  testthat::expect_type(snrna_gene@gene_product$snrna_sequence, "character")
  testthat::expect_type(snrna_gene@clinical_significance, "character")
  
})


# Testing the creation of a SnoRNAGene object with its related slots.

testthat::test_that("createSnoRNAGene creates a valid SnoRNAGene object", {
  
  snornagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                          ranges = IRanges::IRanges(start = 200, end = 1200), 
                          strand = "+")
  
  snornagene_product <- list(snorna_id = "snornaID", 
          snorna_sequence = paste0(
          "CUAUGUGGCGUGUGUGCUAGACGAACGAUGGCUAUAGAUAGGAAUCGCGUCU",
          "CAGAUAGAGGAUUGGAACUGAGUCCGAUAGGGUCACGGAUAGAGCCCGGAUG",
          "GGCGGGAACGCUCAGAGAGGCUCGGUGGUUCUA"))
  
  snorna_gene <- GenePack::createSnoRNAGene(id = "ENST00005859749",
                 hugo_symbol = "SYMBOL1",
                 name = "snoRNA gene name",
                 description = "gene description",
                 tissue_specificity = list("liver", "small bowel"),
                 gene_structure = snornagene_structure,
                 gene_product = snornagene_product,
                 clinical_significance = "association with disease1")
  
  testthat::expect_s4_class(snorna_gene, "SnoRNAGene")
  
  testthat::expect_type(snorna_gene@id, "character")
  testthat::expect_type(snorna_gene@hugo_symbol, "character")
  testthat::expect_type(snorna_gene@name, "character")
  testthat::expect_type(snorna_gene@description, "character")
  testthat::expect_type(snorna_gene@tissue_specificity, "list")
  testthat::expect_s4_class(snorna_gene@gene_structure, "GRanges")
  testthat::expect_type(snorna_gene@gene_product$snorna_id, "character")
  testthat::expect_type(snorna_gene@gene_product$snorna_sequence, "character")
  testthat::expect_type(snorna_gene@clinical_significance, "character")
  
})
