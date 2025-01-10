## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(GenePack)

## -----------------------------------------------------------------------------
lncrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                                                 ranges = IRanges::IRanges(start = 200, end = 1200), 
                                                 strand = "+")

## -----------------------------------------------------------------------------
lncrnagene_product <- list(lncrna_id = "lncRNAID", 
                             lncrna_sequence = 
                                paste0("AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUA",
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
                                          clinical_significance = "association with bowel cancer")

## -----------------------------------------------------------------------------
codinggene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                                                 ranges = IRanges::IRanges(start = 200, end = 1200), 
                                                 strand = "+")
  
codinggene_product <- list(protein_id = "proteinID", 
                           protein_sequence = 
                             paste0("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVTELTKVEPPEVPVNGPAAANLL",
                                    "GGQMKKLGFLHSGTAKSVTCTYSPALNKLTENSSEKGSITRRRCFEGIKTM"))
  
proteincoding_gene <- GenePack::createProteinCodingGene(id = "ENST00005859745",
                                            hugo_symbol = "SYMBOL1",
                                            name = "protein-coding gene name",
                                            description = "gene description",
                                            tissue_specificity = list("liver", "small bowel"),
                                            gene_structure = codinggene_structure,
                                            gene_product = codinggene_product,
                                            clinical_significance = "association with bowel cancer")

## -----------------------------------------------------------------------------
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
                                                clinical_significance = "association with bowel cancer")

## -----------------------------------------------------------------------------
rrnagene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                                               ranges = IRanges::IRanges(start = 200, end = 1200), 
                                               strand = "+")
  
rrnagene_product <- list(rrna_id = "rrnaID", 
                         rrna_sequence = 
                           paste0("GAGCCCGGCGUCCCGUCUCAGACCCGGCCGACAGGAGGGGUGGAGUGGGUGUGUGCGUG",
                                  "UGGAGGAGGGUGUGGUGGUGCGGUGUGGUGGGCGUGUGAGUGUGUGGCGCGUGGAGGAG"))
  
rrna_gene <- GenePack::createRRNAGene(id = "ENST00005859745",
                                                    hugo_symbol = "SYMBOL1",
                                                    name = "rRNA gene name",
                                                    description = "gene description",
                                                    tissue_specificity = list("liver", "small bowel"),
                                                    gene_structure = rrnagene_structure,
                                                    gene_product = rrnagene_product,
                                                    clinical_significance = "association bowel cancer")

## -----------------------------------------------------------------------------
GenePack::getID(lncrna_gene)

GenePack::setID(lncrna_gene) <- "ENST00005859750"
GenePack::getID(lncrna_gene)


## -----------------------------------------------------------------------------
GenePack::getHugoSymbol(lncrna_gene)

GenePack::setHugoSymbol(lncrna_gene) <- "NEWSYMBOL1"
GenePack::getHugoSymbol(lncrna_gene)

## -----------------------------------------------------------------------------
GenePack::getGeneStructure(lncrna_gene)


new_gene_structure2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                               ranges = IRanges::IRanges(start = 205, end = 1210),
                                               strand = "+")

GenePack::setGeneStructure(lncrna_gene) <- new_gene_structure2
GenePack::getGeneStructure(lncrna_gene)

## -----------------------------------------------------------------------------
GenePack::getProductID(lncrna_gene)

GenePack::setProductID(lncrna_gene) <- "newlncRNAID"
GenePack::getProductID(lncrna_gene)

## -----------------------------------------------------------------------------
GenePack::getProductSequence(lncrna_gene)


seq <- paste0("AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUAAUCGGCUUAGCGUACGGUAGG",
                "CUUAACUGCGUACGAUCGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGGCUUAGCGUAC",
                "GUAGGCUCAACUGCGUACGAUCGAUCGGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGGCUUAG",
                "CGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC")

GenePack::setProductSequence(lncrna_gene) <- seq
GenePack::getProductSequence(lncrna_gene)

## -----------------------------------------------------------------------------
computed_length <- GenePack::lengthProduct(lncrna_gene)

## -----------------------------------------------------------------------------
GenePack::showGeneObject(lncrna_gene)

## -----------------------------------------------------------------------------
sessionInfo()

