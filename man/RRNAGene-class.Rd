% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_housekeepingRNAgene_classes.R
\docType{class}
\name{RRNAGene-class}
\alias{RRNAGene-class}
\title{Ribosomal RNA gene class}
\value{
This documentation describes the structure of the
\code{RRNAGene} class.
}
\description{
A class to represent ribosomal RNA genes.
This class is a specialized S4 class inheriting from the \code{Gene}
class and is designed to store information about rRNA genes,
including their rRNA product.
}
\details{
The \code{gene_product} slot is expected to contain a rRNA ID and the
corresponding sequence. The validity function ensures that it is correctly
formatted and contain a valid rRNA sequence.
}
\section{Slots}{

\describe{
\item{\code{gene_product}}{list. This slot is specific for the rRNA genes
product and includes:
\itemize{
\item \code{rrna_id}: a string representing the ID of the rRNA.
\item \code{rrna_sequence}: a string representing the sequence of the
rRNA.
}}
}}

