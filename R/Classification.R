#'This is a copy of DEsingle that was lifted from the orginal package, see https://bioconductor.org/packages/release/bioc/html/DEsingle.html
#'and included in the aocseq package with modifications for package compatibility
#' DEsingle: Detecting differentially expressed genes from scRNA-seq data
#' This function is used to detect differentially expressed genes between two specified groups of cells in a raw read counts matrix of single-cell RNA-seq (scRNA-seq) data. It takes a non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object which contains the read counts matrix. The rows of the matrix are genes and columns are samples/cells.
#' @param group A vector of factor which specifies the two groups to be compared, corresponding to the columns in the counts matrix.
#' @param parallel If FALSE (default), no parallel computation is used; if TRUE, parallel computation using \code{BiocParallel}, with argument \code{BPPARAM}.
#' @param BPPARAM An optional parameter object passed internally to \code{\link{bplapply}} when \code{parallel=TRUE}. If not specified, \code{\link{bpparam}()} (default) will be used.
#' @return
#' A data frame containing the differential expression (DE) analysis results, rows are genes and columns contain the following items:
#' \itemize{
#'   \item theta_1, theta_2, mu_1, mu_2, size_1, size_2, prob_1, prob_2: MLE of the zero-inflated negative binomial distribution's parameters of group 1 and group 2.
#'   \item total_mean_1, total_mean_2: Mean of read counts of group 1 and group 2.
#'   \item foldChange: total_mean_1/total_mean_2.
#'   \item norm_total_mean_1, norm_total_mean_2: Mean of normalized read counts of group 1 and group 2.
#'   \item norm_foldChange: norm_total_mean_1/norm_total_mean_2.
#'   \item chi2LR1: Chi-square statistic for hypothesis testing of H0.
#'   \item pvalue_LR2: P value of hypothesis testing of H20 (Used to determine the type of a DE gene).
#'   \item pvalue_LR3: P value of hypothesis testing of H30 (Used to determine the type of a DE gene).
#'   \item FDR_LR2: Adjusted P value of pvalue_LR2 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item FDR_LR3: Adjusted P value of pvalue_LR3 using Benjamini & Hochberg's method (Used to determine the type of a DE gene).
#'   \item pvalue: P value of hypothesis testing of H0 (Used to determine whether a gene is a DE gene).
#'   \item pvalue.adj.FDR: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method (Used to determine whether a gene is a DE gene).
#'   \item Remark: Record of abnormal program information.
#' }
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{DEtype}}, for the classification of differentially expressed genes found by \code{\link{DEsingle}}.
#'
#' \code{\link{TestData}}, a test dataset for DEsingle.
#'
#' @examples
#' # Load test data for DEsingle
#' data(TestData)
#'
#' # Specifying the two groups to be compared
#' # The sample number in group 1 and group 2 is 50 and 100 respectively
#' group <- factor(c(rep(1,50), rep(2,100)))
#'
#' # Detecting the differentially expressed genes
#' results <- DEsingle(counts = counts, group = group)
#'
#' # Dividing the differentially expressed genes into 3 categories
#' results.classified <- DEtype(results = results, threshold = 0.05)
#'
#' @import stats
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom maxLik maxLik
#' @importFrom pscl zeroinfl
#' @importMethodsFrom Matrix colSums
#' @concept Statistical inference
#'
#' @export

DEsingle <- function(counts, group, goi, parallel = FALSE, BPPARAM = bpparam()){

  # Handle SingleCellExperiment
  if(class(counts)[1] == "SingleCellExperiment"){
    if(!require(SingleCellExperiment))
      stop("To use SingleCellExperiment as input, you should install the package firstly")
    counts <- counts(counts)
  }

  # Invalid input control
  if(!is.matrix(counts) & !is.data.frame(counts) & class(counts)[1] != "dgCMatrix")
    stop("Wrong data type of 'counts'")
  if(sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'");gc();
  if(sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'");gc();
  if(all(counts == 0))
    stop("All elements of 'counts' are zero");gc();
  if(any(colSums(counts) == 0))
    warning("Library size of zero detected in 'counts'");gc();

  if(!is.factor(group))
    stop("Data type of 'group' is not factor")
  if(length(levels(group)) != 2)
    stop("Levels number of 'group' is not two")
  if(table(group)[1] < 2 | table(group)[2] < 2)
    stop("Too few samples (< 2) in a group")
  if(ncol(counts) != length(group))
    stop("Length of 'group' must equal to column number of 'counts'")

  if(!is.logical(parallel))
    stop("Data type of 'parallel' is not logical")
  if(length(parallel) != 1)
    stop("Length of 'parallel' is not one")

  # Preprocessing
  counts <- round(as.matrix(counts))
  #counts <- qq
  storage.mode(counts) <- "integer"
  if(any(rowSums(counts) == 0))
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
  if(length(goi)==1){
    counts <- t(as.matrix(counts[rowSums(counts) != 0,]))
  }else{
  counts <- (as.matrix(counts[rowSums(counts) != 0,]))
  }
  geneNum <- nrow(counts)
  sampleNum <- ncol(counts)
  gc()

  # Normalization
  message("Normalizing the data")
  GEOmean <- rep(NA,geneNum)
  for (i in 1:geneNum)
  {
    gene_NZ <- counts[i,counts[i,] > 0]
    GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE) / length(gene_NZ))
  }

  counts_norm <- matrix(data=0,nrow=geneNum,ncol=sampleNum)
  for (j in 1:sampleNum)
  {
    counts_norm[,j] = unname(counts[,j])/GEOmean[i]
  }
  counts_norm <- ceiling(counts_norm)
  remove(GEOmean, gene_NZ, i, j)
  gc()

  # Cache totalMean and foldChange for each gene
  totalMean_1 <- mean(counts[1,1:length(subset(group,group=="1"))])
  totalMean_2 <- mean(counts[1,(length(subset(group,group=="1"))+1):length(group)])
  foldChange <- totalMean_1/totalMean_2
  All_Mean_FC <- cbind(totalMean_1, totalMean_2, foldChange)

  # Memory management
  remove(counts, totalMean_1, totalMean_2, foldChange)
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()


  # Function of testing homogeneity of two ZINB populations
  CallDE <- function(i){

    # Memory management
    if(i %% 100 == 0)
      gc()

    # Function input and output
    counts_1 <- counts_norm[i, 1:length(subset(group,group=="1"))]
    counts_2 <- counts_norm[i, (length(subset(group,group=="1"))+1):length(group)]
    results_gene <- data.frame(row.names = row.names(counts_norm)[i], theta_1 = NA, theta_2 = NA, mu_1 = NA, mu_2 = NA, size_1 = NA, size_2 = NA, prob_1 = NA, prob_2 = NA, total_mean_1 = NA, total_mean_2 = NA, foldChange = NA, norm_total_mean_1 = NA, norm_total_mean_2 = NA, norm_foldChange = NA, chi2LR1 = NA, pvalue_LR2 = NA, pvalue_LR3 = NA, FDR_LR2 = NA, FDR_LR3 = NA, pvalue = NA, pvalue.adj.FDR = NA, Remark = NA)

    # Log likelihood functions
    logL <- function(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2){
      logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1, pstr0 = theta_1, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2, pstr0 = theta_2, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL2 <- function(param){
      theta_resL2 <- param[1]
      size_1_resL2 <- param[2]
      prob_1_resL2 <- param[3]
      size_2_resL2 <- param[4]
      prob_2_resL2 <- param[5]
      logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL2NZ <- function(param){
      theta_resL2 <- 0
      size_1_resL2 <- param[1]
      prob_1_resL2 <- param[2]
      size_2_resL2 <- param[3]
      prob_2_resL2 <- param[4]
      logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3 <- function(param){
      theta_1_resL3 <- param[1]
      size_resL3 <- param[2]
      prob_resL3 <- param[3]
      theta_2_resL3 <- param[4]
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3NZ1 <- function(param){
      theta_1_resL3 <- 0
      size_resL3 <- param[1]
      prob_resL3 <- param[2]
      theta_2_resL3 <- param[3]
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3NZ2 <- function(param){
      theta_1_resL3 <- param[1]
      size_resL3 <- param[2]
      prob_resL3 <- param[3]
      theta_2_resL3 <- 0
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3AZ1 <- function(param){
      theta_1_resL3 <- 1
      size_resL3 <- param[1]
      prob_resL3 <- param[2]
      theta_2_resL3 <- param[3]
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3AZ2 <- function(param){
      theta_1_resL3 <- param[1]
      size_resL3 <- param[2]
      prob_resL3 <- param[3]
      theta_2_resL3 <- 1
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3NZ1AZ2 <- function(param){
      theta_1_resL3 <- 0
      size_resL3 <- param[1]
      prob_resL3 <- param[2]
      theta_2_resL3 <- 1
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    logL3NZ2AZ1 <- function(param){
      theta_1_resL3 <- 1
      size_resL3 <- param[1]
      prob_resL3 <- param[2]
      theta_2_resL3 <- 0
      logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    judgeParam <- function(param){
      if((param >= 0) & (param <= 1))
        res <- TRUE
      else
        res <- FALSE
      res
    }

    counts_1[is.na(counts_1)] = 0
    counts_2[is.na(counts_2)] = 0
    # MLE of parameters of ZINB counts_1
    if(sum(counts_1 == 0) > 0){
      if(sum(counts_1 == 0) == length(counts_1)){
        theta_1 <- 1
        mu_1 <- 0
        size_1 <- 1
        prob_1 <- size_1/(size_1 + mu_1)
      }else{
        options(show.error.messages = FALSE)
        zinb_try <- try(zeroinfl(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        options(show.error.messages = TRUE)
        if('try-error' %in% class(zinb_try)){
          zinb_try_twice <- try(gamlssML(counts_1, family="ZINBI"), silent=TRUE)
          if('try-error' %in% class(zinb_try_twice)){
            print("MLE of ZINB failed!");
            results_gene[1,"Remark"] <- "ZINB failed!"
            return(results_gene)
          }else{
            zinb_1 <- zinb_try_twice
            theta_1 <- zinb_1$nu;names(theta_1) <- NULL
            mu_1 <- zinb_1$mu;names(mu_1) <- NULL
            size_1 <- 1/zinb_1$sigma;names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL

          }

        }else{
          zinb_1 <- zinb_try
          theta_1 <- plogis(zinb_1$coefficients$zero);names(theta_1) <- NULL
          mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
          size_1 <- zinb_1$theta;names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(counts_1, "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(counts_1~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_1), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_gene[1,"Remark"] <- "NB failed!"
              return(results_gene)
            }else{
              nb_1 <- nb_try_fourth
              theta_1 <- 0
              mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
              size_1 <- nb_1$theta;names(size_1) <- NULL
              prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
            }
          }else{
            nb_1 <- nb_try_again
            theta_1 <- 0
            mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
            size_1 <- 1/nb_1@coef["invk"];names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
          }
        }else{
          nb_1 <- nb_try_twice
          theta_1 <- 0
          mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
          size_1 <- nb_1$estimate["size"];names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }else{
        nb_1 <- nb_try
        theta_1 <- 0
        mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
        size_1 <- nb_1$theta;names(size_1) <- NULL
        prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
      }
    }



    # MLE of parameters of ZINB counts_2
    if(sum(counts_2 == 0) > 0){
      if(sum(counts_2 == 0) == length(counts_2)){
        theta_2 <- 1
        mu_2 <- 0
        size_2 <- 1
        prob_2 <- size_2/(size_2 + mu_2)
      }else{
        options(show.error.messages = FALSE)
        zinb_try <- try(zeroinfl(formula = counts_2 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        options(show.error.messages = TRUE)
        if('try-error' %in% class(zinb_try)){
          zinb_try_twice <- try(gamlssML(counts_2, family="ZINBI"), silent=TRUE)
          if('try-error' %in% class(zinb_try_twice)){
            print("MLE of ZINB failed!");
            results_gene[1,"Remark"] <- "ZINB failed!"
            return(results_gene)
          }else{
            zinb_2 <- zinb_try_twice
            theta_2 <- zinb_2$nu;names(theta_2) <- NULL
            mu_2 <- zinb_2$mu;names(mu_2) <- NULL
            size_2 <- 1/zinb_2$sigma;names(size_2) <- NULL
            prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
          }
        }else{
          zinb_2 <- zinb_try
          theta_2 <- plogis(zinb_2$coefficients$zero);names(theta_2) <- NULL
          mu_2 <- exp(zinb_2$coefficients$count);names(mu_2) <- NULL
          size_2 <- zinb_2$theta;names(size_2) <- NULL
          prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
        }
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(counts_2, "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(counts_2~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_2), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_gene[1,"Remark"] <- "NB failed!"
              return(results_gene)
            }else{
              nb_2 <- nb_try_fourth
              theta_2 <- 0
              mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
              size_2 <- nb_2$theta;names(size_2) <- NULL
              prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
            }
          }else{
            nb_2 <- nb_try_again
            theta_2 <- 0
            mu_2 <- exp(nb_2@coef["logmu"]);names(mu_2) <- NULL
            size_2 <- 1/nb_2@coef["invk"];names(size_2) <- NULL
            prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
          }
        }else{
          nb_2 <- nb_try_twice
          theta_2 <- 0
          mu_2 <- nb_2$estimate["mu"];names(mu_2) <- NULL
          size_2 <- nb_2$estimate["size"];names(size_2) <- NULL
          prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
        }
      }else{
        nb_2 <- nb_try
        theta_2 <- 0
        mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
        size_2 <- nb_2$theta;names(size_2) <- NULL
        prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
      }
    }

    # Restricted MLE under H0 (MLE of c(counts_1, counts_2))
    if(sum(c(counts_1, counts_2) == 0) > 0){
      options(show.error.messages = FALSE)
      zinb_try <- try(gamlssML(c(counts_1, counts_2), family="ZINBI"), silent=TRUE)
      options(show.error.messages = TRUE)
      if('try-error' %in% class(zinb_try)){
        zinb_try_twice <- try(zeroinfl(formula = c(counts_1, counts_2) ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(zinb_try_twice)){
          print("MLE of ZINB failed!");
          results_gene[1,"Remark"] <- "ZINB failed!"
          return(results_gene)
        }else{
          zinb_res <- zinb_try_twice
          theta_res <- plogis(zinb_res$coefficients$zero);names(theta_res) <- NULL
          mu_res <- exp(zinb_res$coefficients$count);names(mu_res) <- NULL
          size_res <- zinb_res$theta;names(size_res) <- NULL
          prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
        }
      }else{
        zinb_res <- zinb_try
        theta_res <- zinb_res$nu;names(theta_res) <- NULL
        mu_res <- zinb_res$mu;names(mu_res) <- NULL
        size_res <- 1/zinb_res$sigma;names(size_res) <- NULL
        prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(c(counts_1, counts_2), "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(c(counts_1, counts_2)~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(c(counts_1, counts_2)), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_gene[1,"Remark"] <- "NB failed!"
              return(results_gene)
            }else{
              nb_res <- nb_try_fourth
              theta_res <- 0
              mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
              size_res <- nb_res$theta;names(size_res) <- NULL
              prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
            }
          }else{
            nb_res <- nb_try_again
            theta_res <- 0
            mu_res <- exp(nb_res@coef["logmu"]);names(mu_res) <- NULL
            size_res <- 1/nb_res@coef["invk"];names(size_res) <- NULL
            prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
          }
        }else{
          nb_res <- nb_try_twice
          theta_res <- 0
          mu_res <- nb_res$estimate["mu"];names(mu_res) <- NULL
          size_res <- nb_res$estimate["size"];names(size_res) <- NULL
          prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
        }
      }else{
        nb_res <- nb_try
        theta_res <- 0
        mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
        size_res <- nb_res$theta;names(size_res) <- NULL
        prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
      }
    }

    # # LRT test of H0
    chi2LR1 <- 2 *(logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_res, size_res, prob_res, counts_2, theta_res, size_res, prob_res))
    pvalue <- 1 - pchisq(chi2LR1, df = 3)

    # Format output
    results_gene[1,"theta_1"] <- theta_1
    results_gene[1,"theta_2"] <- theta_2
    results_gene[1,"mu_1"] <- mu_1
    results_gene[1,"mu_2"] <- mu_2
    results_gene[1,"size_1"] <- size_1
    results_gene[1,"size_2"] <- size_2
    results_gene[1,"prob_1"] <- prob_1
    results_gene[1,"prob_2"] <- prob_2
    results_gene[1,"norm_total_mean_1"] <- mean(counts_1)
    results_gene[1,"norm_total_mean_2"] <- mean(counts_2)
    results_gene[1,"norm_foldChange"] <- results_gene[1,"norm_total_mean_1"] / results_gene[1,"norm_total_mean_2"]
    results_gene[1,"chi2LR1"] <- chi2LR1
    results_gene[1,"pvalue"] <- pvalue

    # Restricted MLE of logL2 and logL3 under H20 and H30 when pvalue <= 0.05
    if(pvalue <= 0.05){
      if(sum(c(counts_1, counts_2) == 0) > 0){
        options(warn=-1)
        # Restricted MLE of logL2
        A <- matrix(rbind(c(1, 0, 0, 0, 0), c(-1, 0, 0, 0, 0), c(0, 0, 1, 0 ,0), c(0, 0, -1, 0 ,0), c(0, 0, 0, 0 ,1), c(0, 0, 0, 0 ,-1)), 6, 5)
        B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10, 1e-10, 1+1e-10)
        mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 0.5, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
        if('try-error' %in% class(mleL2)){
          mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 0, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
        }
        if('try-error' %in% class(mleL2)){
          mleL2 <- try(maxLik(logLik = logL2, start = c(theta_resL2 = 1, size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
        }
        if('try-error' %in% class(mleL2)){
          A <- matrix(rbind(c(0, 1, 0, 0), c(0, -1, 0, 0), c(0, 0, 0 ,1), c(0, 0, 0 ,-1)), 4, 4)
          B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10)
          mleL2 <- maxLik(logLik = logL2NZ, start = c(size_1_resL2 = 1, prob_1_resL2 = 0.5, size_2_resL2 = 1, prob_2_resL2 = 0.5), constraints=list(ineqA=A, ineqB=B))
          theta_resL2 <- 0
          size_1_resL2 <- mleL2$estimate["size_1_resL2"];names(size_1_resL2) <- NULL
          prob_1_resL2 <- mleL2$estimate["prob_1_resL2"];names(prob_1_resL2) <- NULL
          size_2_resL2 <- mleL2$estimate["size_2_resL2"];names(size_2_resL2) <- NULL
          prob_2_resL2 <- mleL2$estimate["prob_2_resL2"];names(prob_2_resL2) <- NULL
        }else{
          theta_resL2 <- mleL2$estimate["theta_resL2"];names(theta_resL2) <- NULL
          size_1_resL2 <- mleL2$estimate["size_1_resL2"];names(size_1_resL2) <- NULL
          prob_1_resL2 <- mleL2$estimate["prob_1_resL2"];names(prob_1_resL2) <- NULL
          size_2_resL2 <- mleL2$estimate["size_2_resL2"];names(size_2_resL2) <- NULL
          prob_2_resL2 <- mleL2$estimate["prob_2_resL2"];names(prob_2_resL2) <- NULL
        }

        # Restricted MLE of logL3
        if((sum(counts_1 == 0) > 0) & (sum(counts_2 == 0) > 0)){
          # logL3
          if(sum(counts_1 == 0) == length(counts_1)){
            A <- matrix(rbind(c(0, 1, 0), c(0, -1, 0), c(0, 0 ,1), c(0, 0 ,-1)), 4, 3)
            B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3AZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- 1
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
          }else if(sum(counts_2 == 0) == length(counts_2)){
            A <- matrix(rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 0 ,1), c(0, 0 ,-1)), 4, 3)
            B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3AZ2, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- 1
          }else{
            A <- matrix(rbind(c(1, 0, 0, 0), c(-1, 0, 0, 0), c(0, 0, 1, 0), c(0, 0, -1, 0), c(0, 0, 0 ,1), c(0, 0, 0 ,-1)), 6, 4)
            B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10, 1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
          }
        }else if(sum(counts_1 == 0) == 0){
          # logL3
          if(sum(counts_2 == 0) == length(counts_2)){
            A <- matrix(rbind(c(0, 1), c(0, -1)), 2, 2)
            B <- c(1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3NZ1AZ2, start = c(size_resL3 = 1, prob_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- 0
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- 1
          }else{
            A <- matrix(rbind(c(0, 1, 0), c(0, -1, 0), c(0, 0 ,1), c(0, 0 ,-1)), 4, 3)
            B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3NZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5, theta_2_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- 0
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- mleL3$estimate["theta_2_resL3"];names(theta_2_resL3) <- NULL
          }
        }else if(sum(counts_2 == 0) == 0){
          # logL3
          if(sum(counts_1 == 0) == length(counts_1)){
            A <- matrix(rbind(c(0, 1), c(0, -1)), 2, 2)
            B <- c(1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3NZ2AZ1, start = c(size_resL3 = 1, prob_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- 1
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- 0
          }else{
            A <- matrix(rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 0 ,1), c(0, 0 ,-1)), 4, 3)
            B <- c(1e-10, 1+1e-10, 1e-10, 1+1e-10)
            mleL3 <- maxLik(logLik = logL3NZ2, start = c(theta_1_resL3 = 0.5, size_resL3 = 1, prob_resL3 = 0.5), constraints=list(ineqA=A, ineqB=B))
            theta_1_resL3 <- mleL3$estimate["theta_1_resL3"];names(theta_1_resL3) <- NULL
            size_resL3 <- mleL3$estimate["size_resL3"];names(size_resL3) <- NULL
            prob_resL3 <- mleL3$estimate["prob_resL3"];names(prob_resL3) <- NULL
            theta_2_resL3 <- 0
          }
        }
        options(warn=0)
      }else{
        # Restricted MLE of logL2
        theta_resL2 <- 0
        size_1_resL2 <- size_1
        prob_1_resL2 <- prob_1
        size_2_resL2 <- size_2
        prob_2_resL2 <- prob_2

        # Restricted MLE of logL3
        theta_1_resL3 <- 0
        size_resL3 <- size_res
        prob_resL3 <- prob_res
        theta_2_resL3 <- 0
      }

      # Judge parameters
      if(!(judgeParam(theta_resL2) & judgeParam(prob_1_resL2) & judgeParam(prob_2_resL2)))
        results_gene[1,"Remark"] <- "logL2 failed!"
      if(!(judgeParam(theta_1_resL3) & judgeParam(theta_2_resL3) & judgeParam(prob_resL3)))
        results_gene[1,"Remark"] <- "logL3 failed!"

      # LRT test of H20 and H30
      chi2LR2 <- 2 *(logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_resL2, size_1_resL2, prob_1_resL2, counts_2, theta_resL2, size_2_resL2, prob_2_resL2))
      pvalue_LR2 <- 1 - pchisq(chi2LR2, df = 1)
      chi2LR3 <- 2 *(logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_1_resL3, size_resL3, prob_resL3, counts_2, theta_2_resL3, size_resL3, prob_resL3))
      pvalue_LR3 <- 1 - pchisq(chi2LR3, df = 2)

      # Format output
      results_gene[1,"pvalue_LR2"] <- pvalue_LR2
      results_gene[1,"pvalue_LR3"] <- pvalue_LR3
    }else{
      results_gene[is.na(results_gene)] = 1
    }

    # Return results_gene
    return(results_gene)
  }


  # Call DEG gene by gene
  if(!parallel){
    results <- matrix(data=NA, nrow = geneNum, ncol = 22, dimnames = list(row.names(counts_norm), c("theta_1", "theta_2", "mu_1", "mu_2", "size_1", "size_2", "prob_1", "prob_2", "total_mean_1", "total_mean_2", "foldChange", "norm_total_mean_1", "norm_total_mean_2", "norm_foldChange", "chi2LR1", "pvalue_LR2", "pvalue_LR3", "FDR_LR2", "FDR_LR3", "pvalue", "pvalue.adj.FDR", "Remark")))
    results <- as.data.frame(results)
    for(i in 1:geneNum){
      cat("\r",paste0("DEsingle is analyzing ", i," of ",geneNum," expressed genes"))
      results[i,] <- CallDE(i)
    }
  }else{
    message("DEsingle is analyzing ", geneNum, " expressed genes in parallel")
    results <- do.call(rbind, bplapply(1:geneNum, CallDE, BPPARAM = BPPARAM))
  }

  # Format output results
  results[, c("total_mean_1", "total_mean_2", "foldChange")] <- All_Mean_FC
  results[,"FDR_LR2"] <- p.adjust(results[,"pvalue_LR2"], method="fdr")
  results[,"FDR_LR3"] <- p.adjust(results[,"pvalue_LR3"], method="fdr")
  results[,"pvalue.adj.FDR"] <- p.adjust(results[,"pvalue"], method="fdr")
  results <- results[order(results[,"chi2LR1"], decreasing = TRUE),]

  # Abnormity control
  if(exists("lastFuncGrad") & exists("lastFuncParam"))
    remove(lastFuncGrad, lastFuncParam, envir=.GlobalEnv)
  if(sum(!is.na(results[,"Remark"])) != 0)
    cat(paste0("\n\n ",sum(!is.na(results[,"Remark"])), " gene failed.\n\n"))

  return(results)


}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Isolation forest solvers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# construct_tree
#' Build a binary tree from a numerical data frame
#' 
#' This function will build a binary tree from a numerical data frame and return
#' a data frame containing each unique data point from the input data and the
#' height at which each data point becomes "isolated" (reaches an external node)
#' or the max_height.
#' 
#' @param data a numerical data frame
#' @param current_height a parameter which tracks the current height of a given construct_tree instance, allows for recursive calling of the function
#' @param max_height the maximum height the tree will grow to before stopping
#' @param kurtosis splits the data based on the kurtosis over the maximum gene distribution across the cells in the input data.
#' 
#' @return a data frame containing the data for each unique point in the input data and the height at which that point was isolated
#' @concept Routine functions
#'
#' @export

construct_tree<-function(data, current_height, max_height, kurtosis=TRUE){
  
  if(current_height==0){
    output_df<<-unique(data)
    output_df$height<<-rep(-1, nrow(output_df))
    if(kurtosis){
      output_df$ID<<-do.call(paste0, output_df[,1:ncol(data),drop=FALSE])
    }
  }
  if(current_height == max_height|nrow(unique(data))==1){
    for (i in 1:nrow(unique(data))) {
      if(kurtosis){
        data$ID<-do.call(paste0, data)
        output_df$height[output_df$ID%in%data$ID]<<-current_height
      }
      else{
        output_df$height[do.call(paste0, output_df[,1:ncol(data),drop=FALSE]) %in% do.call(paste0, unique(data)[i,1:ncol(data),drop=FALSE])]<<-current_height
      }
    }
  }
  else{
    
    if(kurtosis){
      data_pared<-data[unlist(lapply(data, function(x){length(unique(x))!=1}))]
      sample_attribute<-data_pared[kurtosis(data_pared)==max(kurtosis(data_pared), na.rm = TRUE)]
      split_value<-runif(1, min=min(sample_attribute), max=max(sample_attribute))
      split_left<-data[sample_attribute[,1]<=split_value, ,drop=FALSE]
      split_right<-data[sample_attribute[,1]>split_value, , drop=FALSE]
      
      if(nrow(split_left)!=0){
        left<-construct_tree(split_left, current_height+1, max_height)
      }
      if(nrow(split_right)!=0){
        right<-construct_tree(split_right, current_height+1, max_height)
      }
    }
    else{
      sample_attribute<-data[sample(1:ncol(data), 1)]
      split_value<-runif(1, min=min(sample_attribute), max=max(sample_attribute))
      split_left<-data[sample_attribute[,1]<=split_value, ,drop=FALSE]
      split_right<-data[sample_attribute[,1]>split_value, , drop=FALSE]
      if(nrow(split_left)!=0){
        left<-construct_tree(split_left, current_height+1, max_height,kurtosis=FALSE)
      }
      if(nrow(split_right)!=0){
        right<-construct_tree(split_right, current_height+1, max_height,kurtosis=FALSE)
      }
    }
    
  }
  if(all(output_df$height!=-1)){
    return(output_df)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Iso_forest
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build many trees on samples from data and average the heights
#' 
#' @param data a numerical data frame
#' @param num_trees the number of trees to build and average over
#' @param max_height the maximum height each tree will grow to before stopping
#' @param subsample_count the size of the random sample that a tree will be built on. Cannot be larger than the number of rows in the data
#' @param kurtosis_param splits the data based on the kurtosis over the maximum gene distribution across the cells in the input data.
#' 
#' @return a data frame containing the data for each unique point in the input data and the average height at which that point was isolated
#' @concept Rutine functions
#'
#' @export

Iso_forest<-function(data, num_trees, max_height, subsample_count=nrow(data),kurtosis_param=TRUE){
  results_df<-unique(data)
  results_df$avg_height<-rep(-1, nrow(results_df))
  tree_counter<-0
  while(tree_counter<num_trees) {
    random_subsample<-data[sample(1:nrow(data), min(subsample_count, nrow(data))), ,drop=FALSE]
    tree<-construct_tree(random_subsample, 0, max_height,kurtosis=kurtosis_param)
    for(i in 1:nrow(tree)){
      matched_logical<-do.call(paste0, results_df[,1:ncol(data), drop=FALSE]) %in% do.call(paste0, tree[i,1:ncol(data), drop=FALSE])
      if(results_df$avg_height[matched_logical]==-1){
        results_df$avg_height[matched_logical]<-tree[i, 'height']
      }
      else{
        results_df$avg_height[matched_logical]<-mean(results_df$avg_height[matched_logical],tree[i, 'height'])
      }
    }
    tree_counter<-tree_counter+1
  }
  if(any(results_df$avg_height==-1)){
    print('Some data points were never sampled. Increase num_trees or subsample_count.')
  }
  return(results_df)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Aocseq functions - uses inbuilt version of DEsingle
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Combine single cell data and annotate with cell labels based on functionality
#'
#' This function will read in Seurat objects that have been processed with aocseq,
#' and accompanying annotation tables to classify clonotypes using distance scores.
#' Takes as input gene expression from scRNAseq and VDJ enrichment. Outputs a classification table
#' that lists samples and each clonotype classified.
#'
#' @param Clonal_Obs A seurat object. Pre-processed with aocseq.
#' @param Clonotypes Table. An annotation table from the aocseq package.
#' @param method Vector of "ZINB", "Mahalanobis", "z-score" and "Taxi cab". "ZINB" can be used without a reference dataset.
#' @param goi Gene. Gene of interest (goi). For use with "ZINB" classification.
#' @param percentile Double. Percentile cutoff for distance based scoreing metrics.
#' @param path Directory for output tables.
#' @param verbose Print progress bars and output.
#'
#' @return A Seurat object list containing metadata and VDJ annotations.
#' @concept annotation
#' @concept Single cell analysis
#'
#' @export
ClassifyCellTypes <- function(
    Clonal_Obs,
    Clonotypes,
    within.sample=FALSE,
    method="ZINB",
    goi=c("IFNG"),
    percentile=0.01,
    path="",
    path.glist="",
    reference=matrix(nrow=1,ncol=1),
    verbose=TRUE
){

  if(method=="ZINB"){
  path=paste(paste(path,"/",sep=""),paste(paste(goi[1],"ClassificationTable",sep=""),".csv",sep=""),sep="")
  print(path)
  CloneList=subset(Clonotypes,avg>3)$`Clone (nucleic)`
  ClassArray=data.frame(matrix(data="-",ncol=((4*(length(Clonal_Obs)-1))+2),nrow=length(CloneList)))
  names.sample=0
  if(within.sample==FALSE){
  for(j in 1:(length(Clonal_Obs)-1)){
  names.sample=append(names.sample,c(
    paste(levels(factor(Clonal_Obs[[1]]@meta.data$orig.ident)),".status",sep=""),
    paste(levels(factor(Clonal_Obs[[1+j]]@meta.data$orig.ident)),".status",sep=""),
    paste(levels(factor(Clonal_Obs[[1]]@meta.data$orig.ident)),".abundance",sep=""),
    paste(levels(factor(Clonal_Obs[[1+j]]@meta.data$orig.ident)),".abundance",sep="")
  ))

  }
  }else{
    ClassArray=data.frame(matrix(data="-",ncol=4+2,nrow=length(CloneList)))

      names.sample=append(names.sample,c(
        paste("sample",".status",sep=""),
        paste(levels(factor(Clonal_Obs[[1]]@meta.data$orig.ident)),".status",sep=""),
        paste("sample",".abundance",sep=""),
        paste(levels(factor(Clonal_Obs[[1]]@meta.data$orig.ident)),".abundance",sep="")
      ))


  }
  names.sample=names.sample[-1]
  names.sample=append(names.sample,c("cdr3_na","phenotype"))
  names(ClassArray) <- names.sample

  df_speed=Clonal_Obs[[1]]@meta.data


  if(within.sample==FALSE){
  for(g in 1:length(CloneList)){

    LE=dim(ClassArray)[[2]]
    ClassArray[g,(LE-1)]=CloneList[g]
    cd4prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD4cells=="1")$CD4cells)
    cd8prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD8cells=="1")$CD8cells)
    lenCln=length(subset(df_speed,cdr3_na == CloneList[g])$CD8cells)
    if((lenCln>5)&(cd4prop>cd8prop)){
      ClassArray[g,(LE)]="CD4"
    }else if((lenCln>5)&(cd8prop>cd4prop)){
      ClassArray[g,(LE)]="CD8"
    }else{
      ClassArray[g,(LE)]="-"
    }

    for(j in 1:(length(Clonal_Obs)-1)){
    TCR1=CloneList[g]
    GOIUS=match(goi,row.names(Clonal_Obs[[1]]@assays$RNA@counts))
    GOISTIM=match(goi,row.names(Clonal_Obs[[1+j]]@assays$RNA@counts))
    if((length(intersect(Clonal_Obs[[1]]@meta.data$cdr3_na,TCR1))>0)&&
       (length(intersect(Clonal_Obs[[1+j]]@meta.data$cdr3_na,TCR1))>0)){

      if((dim(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1))[2]>2)&&
         (dim(subset(Clonal_Obs[[1+j]],cdr3_na %in% TCR1))[2]>2)){

        if(length(goi)==1){
          mats=cbind((t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOIUS,]))),
                     (t(as.matrix(subset(Clonal_Obs[[1+j]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM,]))))
          grouper=c(rep(1,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOIUS[1],])))[2]),rep(2,dim(t(as.matrix(subset(Clonal_Obs[[1+j]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM[1],])))[2]))

          if(sum(mats!=0)>5){
            output1=DEsingle(mats,factor(grouper),goi)
            ThetaIDX=which.min(output1$pvalue_LR2)
            MuIDX=which.min(output1$pvalue_LR3)
            if((min(output1$pvalue_LR2)<0.05)&&(min(output1$pvalue_LR2)>10^-(5))){
              if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                ClassArray[g,4*(j-1)+2]=="-"
              }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                ClassArray[g,4*(j-1)+1]="-"
                ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
              }
            }else if(min(output1$pvalue_LR2)<10^-(5)){
              if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                ClassArray[g,4*(j-1)+2]="-"
              }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                ClassArray[g,4*(j-1)+1]="-"
                ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
              }
            }
            if((min(output1$pvalue_LR3)<0.05)&&(min(output1$pvalue_LR3)>10^-(5))){
              if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                ClassArray[g,4*(j-1)+4]="-"
              }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                ClassArray[g,4*(j-1)+3]="-"
                ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
              }
            }else if(min(output1$pvalue_LR3)<10^-(5)){
              if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                ClassArray[g,4*(j-1)+4]="-"
              }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                ClassArray[g,4*(j-1)+3]="-"
                ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
              }
            }
          }
          }
        else{
        mats=cbind((as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOIUS,])),
                   (as.matrix(subset(Clonal_Obs[[1+j]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM,])))
        grouper=c(rep(1,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOIUS[1],])))[2]),rep(2,dim(t(as.matrix(subset(Clonal_Obs[[1+j]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM[1],])))[2]))
        if(sum(mats!=0)>5){
          output1=DEsingle(mats,factor(grouper),goi)
          ThetaIDX=which.min(output1$pvalue_LR2)
          MuIDX=which.min(output1$pvalue_LR3)
          if((min(output1$pvalue_LR2)<0.05)&&(min(output1$pvalue_LR2)>10^-(5))){
            if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
              ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
              ClassArray[g,4*(j-1)+2]=="-"
            }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
              ClassArray[g,4*(j-1)+1]="-"
              ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
            }
          }else if(min(output1$pvalue_LR2)<10^-(5)){
            if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
              ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
              ClassArray[g,4*(j-1)+2]="-"
            }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
              ClassArray[g,4*(j-1)+1]="-"
              ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
            }
          }
          if((min(output1$pvalue_LR3)<0.05)&&(min(output1$pvalue_LR3)>10^-(5))){
            if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
              ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
              ClassArray[g,4*(j-1)+4]="-"
            }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
              ClassArray[g,4*(j-1)+3]="-"
              ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
            }
          }else if(min(output1$pvalue_LR3)<10^-(5)){
            if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
              ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
              ClassArray[g,4*(j-1)+4]="-"
            }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
              ClassArray[g,4*(j-1)+3]="-"
              ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
            }
          }
        }
}
      }
    }
    rm(output1)
    }
    print(g)
  }
  }else{
    for(g in 1:length(CloneList)){

      LE=dim(ClassArray)[[2]]
      ClassArray[g,(LE-1)]=CloneList[g]
      cd4prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD4cells=="1")$CD4cells)
      cd8prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD8cells=="1")$CD8cells)
      lenCln=length(subset(df_speed,cdr3_na == CloneList[g])$CD8cells)
      if((lenCln>5)&(cd4prop>cd8prop)){
        ClassArray[g,(LE)]="CD4"
      }else if((lenCln>5)&(cd8prop>cd4prop)){
        ClassArray[g,(LE)]="CD8"
      }else{
        ClassArray[g,(LE)]="-"
      }

      for(j in 1:1){
        TCR1=CloneList[g]
        GOIUS=match(goi,row.names(Clonal_Obs[[1]]@assays$RNA@counts))
        GOISTIM=GOIUS
        TCRS=setdiff(levels(factor(Clonal_Obs[[1]]@meta.data$cdr3_na)),TCR1)
        if((length(intersect(Clonal_Obs[[1]]@meta.data$cdr3_na,TCR1))>0)&&
           (length(intersect(Clonal_Obs[[1]]@meta.data$cdr3_na,TCRS))>0)){

          if((dim(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1))[2]>2)&&
             (dim(subset(Clonal_Obs[[1]],cdr3_na %in% TCRS))[2]>2)){

            if(length(goi)==1){
              mats=cbind((t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCRS)@assays$RNA@counts[GOIUS,]))),
                         (t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM,]))))
              grouper=c(rep(1,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCRS)@assays$RNA@counts[GOIUS,])))[2]),rep(2,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM,])))[2]))

              if(sum(mats!=0)>5){
                output1=DEsingle(mats,factor(grouper),goi)
                ThetaIDX=which.min(output1$pvalue_LR2)
                MuIDX=which.min(output1$pvalue_LR3)
                if((min(output1$pvalue_LR2)<0.05)&&(min(output1$pvalue_LR2)>10^-(5))){
                  if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                    ClassArray[g,4*(j-1)+2]=="-"
                  }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]="-"
                    ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
                  }
                }else if(min(output1$pvalue_LR2)<10^-(5)){
                  if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                    ClassArray[g,4*(j-1)+2]="-"
                  }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]="-"
                    ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
                  }
                }
                if((min(output1$pvalue_LR3)<0.05)&&(min(output1$pvalue_LR3)>10^-(5))){
                  if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                    ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                    ClassArray[g,4*(j-1)+4]="-"
                  }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                    ClassArray[g,4*(j-1)+3]="-"
                    ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
                  }
                }else if(min(output1$pvalue_LR3)<10^-(5)){
                  if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                    ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                    ClassArray[g,4*(j-1)+4]="-"
                  }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                    ClassArray[g,4*(j-1)+3]="-"
                    ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
                  }
                }
              }
            }else{
              mats=cbind((as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCRS)@assays$RNA@counts[GOIUS,])),
                         (as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM,])))
              grouper=c(rep(1,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCRS)@assays$RNA@counts[GOIUS[1],])))[2]),rep(2,dim(t(as.matrix(subset(Clonal_Obs[[1]],cdr3_na %in% TCR1)@assays$RNA@counts[GOISTIM[1],])))[2]))
              if(sum(mats!=0)>5){
                output1=DEsingle(mats,factor(grouper),goi)
                ThetaIDX=which.min(output1$pvalue_LR2)
                MuIDX=which.min(output1$pvalue_LR3)
                if((min(output1$pvalue_LR2)<0.05)&&(min(output1$pvalue_LR2)>10^-(5))){
                  if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                    ClassArray[g,4*(j-1)+2]=="-"
                  }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]="-"
                    ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
                  }
                }else if(min(output1$pvalue_LR2)<10^-(5)){
                  if(output1$theta_1[ThetaIDX]<output1$theta_2[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]=output1$pvalue_LR2[ThetaIDX]
                    ClassArray[g,4*(j-1)+2]="-"
                  }else if(output1$theta_2[ThetaIDX]<output1$theta_1[ThetaIDX]){
                    ClassArray[g,4*(j-1)+1]="-"
                    ClassArray[g,4*(j-1)+2]=output1$pvalue_LR2[ThetaIDX]
                  }
                }
                if((min(output1$pvalue_LR3)<0.05)&&(min(output1$pvalue_LR3)>10^-(5))){
                  if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                    ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                    ClassArray[g,4*(j-1)+4]="-"
                  }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                    ClassArray[g,4*(j-1)+3]="-"
                    ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
                  }
                }else if(min(output1$pvalue_LR3)<10^-(5)){
                  if(output1$mu_1[MuIDX]>output1$mu_2[MuIDX]){
                    ClassArray[g,4*(j-1)+3]=output1$pvalue_LR3[MuIDX]
                    ClassArray[g,4*(j-1)+4]="-"
                  }else if(output1$mu_2[MuIDX]>output1$mu_1[MuIDX]){
                    ClassArray[g,4*(j-1)+3]="-"
                    ClassArray[g,4*(j-1)+4]=output1$pvalue_LR3[MuIDX]
                  }
                }
              }
            }
          }
        }
        rm(output1)
      }
      print(g)
    }
  }
  write.csv(ClassArray,path)
  }
  else{
    if(method=="Mahalanobis"){
      inputdistance=0
    }
    else if(method=="z-score"){
      inputdistance=2
    }
    else if(method=="Taxi cab"){
      inputdistance=1
    }
    else if(method=="isoForest"){
      inputdistance=3
      ClassifyCells(cell.data,cell.arrayPFlog1pPF,signature.ref,Glist,distance=inputdistance)
    }
    path=paste(path,paste(paste(goi,"ClassificationTable",sep=""),".csv",sep=""),sep="")
    CloneList=subset(Clonotypes,avg>3)$Clone..nucleic.
    ClassArray=data.frame(matrix(ncol=(((length(Clonal_Obs)))+2),nrow=length(CloneList)))
    names.sample=0
    for(j in 1:(length(Clonal_Obs)-1)){
      names.sample=append(names.sample,c(
        paste(levels(factor(Clonal_Obs[[j]]@meta.data$orig.ident)),".Rscore",sep="")))

    }
    names.sample=names.sample[-1]
    names.sample=append(names.sample,c("cdr3_na","phenotype"))
    names(ClassArray) <- names.sample

    df_speed=Clonal_Obs[[1]]@meta.data

    for(g in 1:length(CloneList)){
      LE=dim(ClassArray)[[2]]
      ClassArray[g,(LE-1)]=CloneList[g]
      cd4prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD4cells=="1")$CD4cells)
      cd8prop=length(subset(df_speed,cdr3_na == CloneList[g] & CD8cells=="1")$CD8cells)
      lenCln=length(subset(df_speed,cdr3_na == CloneList[g])$CD8cells)
      if((lenCln>5)&(cd4prop>cd8prop)){
        ClassArray[g,(LE)]="CD4"
      }
      else if((lenCln>5)&(cd8prop>cd4prop)){
        ClassArray[g,(LE)]="CD8"
      }
      else{
        ClassArray[g,(LE)]="-"
      }
      print(g)

      for(j in 1:(length(Clonal_Obs)-1)){
    array.name=1
    Clonal_Obs[[j]]=AddDistances(Clonal_Obs[[j]],array.name,reference,path.glist,distance = inputdistance)
    ClassArray[g,j]=mean(subset(Clonal_Obs[[j]],`cdr3_na` %in% CloneList[g])@meta.data$Mdist)
      }
    }

    ###Find the percentile cutoff for the distance
    Thresh=quantile(ClassArray[,c.index],percentile)
    for(g in 1:length(CloneList)){

      for(j in 1:length(Clonal_Obs)){
        if(ClassArray[g,j]<Thresh){
      ClassArray[g,j]="Specific"
        }
        else{
          ClassArray[g,j]="Bystander"
        }
      }

    }

    write.csv(ClassArray,path)
  }
}

