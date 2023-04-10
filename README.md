# Trace sequencing

**traceseq** is software for analysis of immunosequencing data and can be used to predict and trace recently stimulated T cells in the blood with their gene expression transcripts. The package was developed by the department of Cell and Gene Therapy at Baylor College of Medicine. If you use **traceseq** please consider citing one of the following: ...

The aim of the tool is to provide a method to rapidly quantify the number of activated cells that can bind a target, either in an immunotherapy product or after its infusion when administered as therapy. The tool is designed to be compatible with 10X genomics Chromium Single Cell Immune Profiling data after FASTQ files have been processed and if cells are barcoded by their T cell receptor (TCR), unique clones with target specificity can be identified by their gene expression. Thus the tool in addition serves as a method to identify noval neoantigen specific TCRs. The tool is compatible with multiple sequencing data formats including Adaptive biotechnologies' T cell receptor (TCR) immunoseq assay (https://www.immunoseq.com), or any time series CDR3 TCRbeta dataset of a TCR repertoire in the format of nucleic acid sequence productive frequency, the corresponding amino acid sequences and the total number of T cells at each time point. The data can be from multiple conditions or treatments and multiple time points. 
# Biological motivation
The TCR sequence is important to study because T cell specificity depends on the **CDR3** region of the TCR and this varys in sequence between different people and single cells. T cells can be grouped into clonotypes that share a common CDR3 beta chain and this way, used to estimate the frequency of target specific T cells in the blood. Harnessing this heterogeneity in sequence between T cells for a quantitative analysis of adaptive immunology has broad applicability in immunology and immunotherapy because the clearing of infection and cancer depends on availability of immune cells (including T cells) with capacity to mount a response. 
Immunosequencing is a PCR-based based method that exploits the capacity of high-throughput sequencing technology to characterize tens of thousands of TCR CDR3 chains simultaneously and **traceseq** has been developed to efficiently analyse and annotate this data.

Specifically, this is a software package of statistical tools that can be used to trace, analyse, annotate and query clonotypes subject to amplification or reduction in frequency following antigen stimulation or between experimental conditions. **traceseq** has initially been applied to Virus specific T cells (VSTs) because these clinical blood products contain non specific bystander T cells in addition to potent virus specific clonotypes. However, the tool can be adapted to model other barcoded time series frequency data and the accompanying vignette demonstrates how the tool could be used to track clonotypes *in vivo*, using Adaptive's ImmuneAccess database. 

# Installation and running the software: 

install_github("MaeWoods/traceseq");

library("traceseq")

Steps to run the software are illustrated below in the flow chart and functions are documented in detail in the manual. Sequence and time series data are imported along with total T cell numbers at each time point.  **traceseq** is different to alternative TCR frequencing tracking methods because instead of setting a difference in frequency to label cells *a priori*, or modelling the probability of cell capture as a binomial distribution, **VRTrack** is Bayesian, in the sense that the frequency of all clonotypes at all time points are included in a statistical model to fit the average expansion of all clonotypes within a product over time. From this model, individual clonotypes can be classified as It is common for clonotypes to fall below the limit of detection in time series immunosequencing experiments model and the tool is designed to model this by including drop out events.

Results of VRTrack provide additional to extrapolate the CDR3 sequences of T cell clonotypes with the greatest expansion, annotates these clonotypes with sequence metadata so that the frequency of clonotypes that share the same amino acid sequence (homoplastic frequency) can be jointly merged with expansion, queried in online databases and included as additional metadata in a single cell RNA sequencing (scRNAseq) experiment.

<img src="VRTrack.png" width="100%"></img>

Flow chart showing VRTrack usage.

# Data inputs: 
Total numbers $(n)$ of T cells for each condition and time point must be included as a vector. The numbers should be listed in sequential order of time points and the ordering of the conditions should not change between time points, for example for $j$ conditions $n_{1}-n_{j}$ over $k$ time points $n_{1}(1)-n_{j}(k)$, the input vector should be in the form $$N=(n_{1}(1),n_{2}(1),...,n_{j}(1),...,n_{1}(k),n_{2}(k),...,n_{j}(k)).$$ If importing an Adaptive TCR immunoseq assay, the input files are the track rearrangements files for both the nucleic acid and amino acid sequences and the rearrangements file. For alternative data, matrices must be included in a specific format. Two matrices are required for input, one with the CDR3 sequence in amino acids and another with the CDR3 sequence in nucleic acids. Rows of the matrix must correspond to unique TCRBeta CDR3s and columns of the matrix should correspond to the TCR repertoire for each sample so that elements of the matrix are the productive frequency of each rearrangement for each sample.

Columns should be organised in the same way as the cell number input, i.e. by treatment or donor first for the initial time point, so that for $q$ unique CDR3 sequences, $j$ conditions and $k$ time points, the input matrix $M$ has $q$ rows and $j\mbox{k}$ columns

$$M=\begin{bmatrix}a_{11}(1) & a_{12}(1) & ... & a_{1j}(1) & ... & a_{11}(k) & a_{12}(k) & ... & a_{1j}(k) \\\ a_{21}(1) & a_{22}(1) & ... & a_{2j}(1) & ... & a_{21}(k) & a_{22}(k) & ... & a_{2j}(k) \\\ \vdots & \vdots & & \vdots & \ddots & \vdots & \vdots& & \vdots \\\a_{q1}(1) & a_{q2}(1) & ... & a_{qj}(1) & ... & a_{q1}(k) & a_{q2}(k) & ... & a_{qj}(k) \end{bmatrix},$$

where $a_{qj}(k)$ is the frequency of the $q$th TCR in the $j$th condition at time point $k$.


# Documentation: 
traceseq can be used following the flow chart above, to read in the data and create a clonal object with CDR3 annotations and time series productive frequency, run the function CreateClonalObject().

cp paste the function definitions

To model clonotype expansion, a suite of statistical models are provided to 
