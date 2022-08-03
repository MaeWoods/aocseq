# VRTrack
VRTrack is a software package written in R and Rcpp for the prediction of antigen specific T cells from the masses.

The tool is compatible with Adaptive biotechnologies' T cell receptors (TCR) immunoseq assay (https://www.immunoseq.com), or any time series CDR3 TCRbeta dataset of a TCR repertoire in the format of nucleic acid sequence productive frequency, the corresponding amino acid sequences and the total number of T cells at each time point. The data can be from multiple conditions or treatments and multiple time points.

# Biological motivation
The TCR sequence is important to study because the specificity of a T cell depends on its variable CDR3 region of the TCR. The CDR3 region can be used to group cells into clonotypes that share a common beta sequence and this way, clonotypes serve to approximate the potential frequency of target specific T cells in a blood product. This frequency therefore has implications in immunotherapy and in vivo data analysis because the immune response to a pathogen depends on the availability of immune cells with capacity to mitigate infection and the number of antigen specific T cells contributes to this. Variable region Track (VRTrack) is a software package that demonstrates statistical methods that can be used to analyse antigen specific clonotypes that are present at low frequency and require amplification by antigen stimulation. VRTrack is applied to Virus specific T cells (VSTs) that contain a variable number of cells that are specific for the target and this is caused by a mix of bystander and potent VST clonotypes. Antigen specific clonotypes are present at low frequency and require amplification by antigen stimulation. However, activation signals can promote expansion of bystander clonotypes. Therefore, given time series expansion, in addition to using frequency doubling metrics or the probability that two observations of a clonotype frequency fall outside the 95% confidence interval of a binomial distribution, we can use the change in frequency of all clonotypes in a stistical model to predict the average expansion of all clonotypes within a product This tool uses a statistical model that includes drop out events to extrapolate the CDR3 sequences of T cell clonotypes with the greatest expansion, annotates these clonotypes with sequence metadata so that the frequency of clonotypes that share the same amino acid sequence (homoplastic frequency) can be jointly merged with expansion, queried in online databases and included as additional metadata in a single cell RNA sequencing (scRNAseq) experiment.

<img src="VRTrack.png" width="90%"></img>

Flow chart showing VRTrack usage.

# Data inputs: 
Total numbers of T cells for each condition and time point must beTCR immunoseq assay (https://www.immunoseq.com), or matrices with a specific format. If using matrices, two matrices are required for input, one with the CDR3 sequence in amino acids and another with the CDR3 sequence in nucleic acids. Rows of the matrix must correspond to unique TCRBeta CDR3s and columns of the matrix should correspond to the TCR repertoire for each sample so that elements of the matrix are the productive frequency of each rearrangement for each sample.

Columns should be organised by treatment or donor first for the initial time point, the order of these columns must be identical for each time point, so that the matrix has N rows and (T x S) columns, where N is the total number of unique CDR3s, T is the number of time points and S are the number of treatments or conditions.

# Installation: 

install_github("MaeWoods/VRTrack");

library("VRTrack")

# Documentation: 
VRTrack can be used following the flow chart above, to read in the data and create a clonal object with CDR3 annotations and time series productive frequency, run the function CreateClonalObject().

cp paste the function definitions

To model clonotype expansion, a suite of statistical models are provided to 
