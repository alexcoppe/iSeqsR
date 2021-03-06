% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variants.table.R
\name{variants.recurrence.by.gene.name}
\alias{variants.recurrence.by.gene.name}
\title{Summarises genes by the number of patients in which are found mutated}
\usage{
variants.recurrence.by.gene.name(variants.table)
}
\arguments{
\item{variants.table}{A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and gene_name columns.}
}
\value{
A tbl_df with 2 columns: gene name, number of patients in which the
}
\description{
Given variants data obtained with read.variants.tsv, returns a data frame with the number of patients in which a gene is found mutated.
}
