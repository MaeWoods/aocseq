# aocseq 0.1.0: An R package for mixed cell type annotation following the activation of cells & sequencing.

**aocseq** is a suite of statistical tools that are publically available for the analysis of multimodal immunosequencing data. For example, the package can be used to trace recently stimulated T cells in the blood with their gene expression transcripts. 

The aim of the tool is to adapt probability and control theory for the classification of cells to provide 1) A sample breakdown of mixed cell types and 2) To detect, score and rank cells with optimal or desirable characteristics. Characteristics are defined by the gene expression of cells that are considered desirable by their response to perturbation. Version 0.1.0 can be used to quantify cells that are activated through the binding of their T cell receptors (TCRs) to the cognate antigens. Thus providing a method to identify full length TCRs with specificity directed toward a target antigen, that can then be tested for use as transgenes. In a later update, functionality for time series analysis will be included to analyze changes in clone frequency in the TCR repertoire in the format of nucleic acid sequence productive frequency, the corresponding amino acid sequences and the total number of T cells. Currently, the data can be from multiple conditions or treatments, hashtagged and include multiple time points. 
# Biological motivation
The TCR sequence is important to study because T cell specificity depends on the **CDR3** region of the TCR and this varys in sequence between different people and single cells. T cells can be grouped into clonotypes that share a common CDR3 beta chain and this way, used to estimate the frequency of target specific T cells in the blood. Harnessing this heterogeneity in sequence between T cells for a quantitative analysis of adaptive immunology has broad applicability in immunology and immunotherapy because the clearing of infection and cancer depends on availability of immune cells (including T cells) with capacity to mount a response. 
Immunosequencing is a PCR-based based method that exploits the capacity of high-throughput sequencing technology to characterize tens of thousands of TCR CDR3 chains simultaneously and **aocseq** has been developed to analyse and annotate this data.

Specifically, this is a software package of statistical tools that can be used to trace, analyse, annotate and query clonotypes subject to amplification or reduction in frequency following antigen stimulation or between experimental conditions. **aocseq** has initially been applied to Virus specific T cells (VSTs) because these clinical blood products contain non specific bystander T cells in addition to potent virus specific clonotypes. However, the tool can be adapted to model other barcoded time series frequency data and the accompanying vignette demonstrates how the tool could be used to track clonotypes *in vivo*, using Adaptive's ImmuneAccess database. 

# Installation and running the software: 

install_github("MaeWoods/aocseq");

library("aocseq")

Steps to run the software are illustrated below in the flow chart and functions are defined in the documentation. The fuctionality circled with the dashed line is available in version 0.1.0 of the software. See the **vignettes** for further help loading and running the software.


<img src="aocseq.png" width="100%"></img>

Flow chart showing aocseq usage  “Created in Lucidchart, www.lucidchart.com”.

# Getting started: 
Initial preprocessing of gene expression arrays with cell type annotation is provided in the function CombineData. aocseq can be used to analyze data in the form of VDJ enrichment, immunoprofiling, Multiplex data, Hash tagged data and Cell surface labelling. Examples are provided below with accompanying vignettes

To read in gene expression data and analyze clonotypes:
[Getting started with aocseq](./html/GettingStartedWithaocseq.svg)



If importing immunoseq data without single cell RNA or hashtagged sequencing, outputs of aocseq are compatible with the environments single cell experiment, Seurat and scanpy but are stored in S4 objects that can be added to an S4 object of class Seurat. If importing an Adaptive TCR immunoseq assay, total numbers $(n)$ of T cells for each condition and time point must be included as a vector. The numbers should be listed in sequential order of time points and the ordering of the conditions should not change between time points, for example for $j$ conditions $n_{1}-n_{j}$ over $k$ time points $n_{1}(1)-n_{j}(k)$, the input vector should be in the form $$N=(n_{1}(1),n_{2}(1),...,n_{j}(1),...,n_{1}(k),n_{2}(k),...,n_{j}(k)).$$ The input files are the track rearrangements files for both the nucleic acid and amino acid sequences and the rearrangements file. For alternative data, matrices must be included in a specific format. Two matrices are required for input, one with the CDR3 sequence in amino acids and another with the CDR3 sequence in nucleic acids. Rows of the matrix must correspond to unique TCRBeta CDR3s and columns of the matrix should correspond to the TCR repertoire for each sample so that elements of the matrix are the productive frequency of each rearrangement for each sample.

Columns should be organised in the same way as the cell number input, i.e. by treatment or donor first for the initial time point, so that for $q$ unique CDR3 sequences, $j$ conditions and $k$ time points, the input matrix $M$ has $q$ rows and $j\mbox{k}$ columns

$$M=\begin{bmatrix}a_{11}(1) & a_{12}(1) & ... & a_{1j}(1) & ... & a_{11}(k) & a_{12}(k) & ... & a_{1j}(k) \\\ a_{21}(1) & a_{22}(1) & ... & a_{2j}(1) & ... & a_{21}(k) & a_{22}(k) & ... & a_{2j}(k) \\\ \vdots & \vdots & & \vdots & \ddots & \vdots & \vdots& & \vdots \\\a_{q1}(1) & a_{q2}(1) & ... & a_{qj}(1) & ... & a_{q1}(k) & a_{q2}(k) & ... & a_{qj}(k) \end{bmatrix},$$

where $a_{qj}(k)$ is the frequency of the $q\mbox{th}$ TCR in the $j\mbox{th}$ condition at time point $k$.


# Documentation: 
The aocseq package contains documentation and a set of vignettes are being developed to demonstrate the processing of hashtagged data, immunosequencing data and single cell RNA sequencing data.
