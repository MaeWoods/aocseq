# VRTrack
VRTrack is a software package written in R and Rcpp for the prediction of antigen specific T cells from the masses.

VRTrack is compatible with Adaptive biotechnologies' TCR immunoseq assay (https://www.immunoseq.com), or any time series TCR repertoire CDR3 TCRbeta data of nucleic acid sequence productive frequency, their corresponding amino acid sequences and the total number of T cells.

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
