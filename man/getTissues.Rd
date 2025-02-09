% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_get_set_tissue_specificity.R
\name{getTissues}
\alias{getTissues}
\alias{getTissues,Gene-method}
\title{Get the tissue specificity list}
\usage{
getTissues(object)

\S4method{getTissues}{Gene}(object)
}
\arguments{
\item{object}{Gene object}
}
\value{
the tissue specificity of the gene, or "-" if it was not specified
at the creation.
}
\description{
This function retrieves the list of tissues where the gene is
specifically expressed.
}
\details{
This function is defined as a generic function to be applicable to
any object that inherits from the Gene class.
}
\examples{
gene_structure <- GenomicRanges::GRanges(seqnames = "chr1",
                  ranges = IRanges::IRanges(start = 200, end = 1200),
                  strand = "+")
                         
gene_product <- list(lncrna_id = "lncRNAID", 
                     lncrna_sequence = paste0(
                       "AUGCUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCG",
                       "GCUAAUCGGCUUAGCGUACGGUAGGCUUAACUGCGUACGAU",
                       "CGAUCGGCUAAGCGUACGGUAGGCUUAACUGCGUACGAUCG",
                       "AUCGGCUUAGCGUACGUAGGCUCAACUGCGUACGAUCGAUC",
                       "GGCUUAGCGUACGGUAGGCUUAACUGCGGACGAUCGAUCGG",
                       "CUUAGCGUACGGUAGGCUUAACUGCGUACGAUCGAUCGC"))
                                          
lncrna_gene <- createLncRNAGene(id = "ENST000001",
               hugo_symbol = "SYMBOL1",
               name = "lncRNA gene name",
               description = "gene description",
               tissue_specificity = list("liver", "small bowel"),
               gene_structure = gene_structure,
               gene_product = gene_product,
               clinical_significance = "association with disease")

}
