#' iNterpolation and EXTrapolation with Beta Diversity for TD, PD, and FD
#'
#' \code{iNEXTbeta3D()} computes standardized 3D gamma and alpha diversity estimates at common sample sizes, or gamma, alpha, and beta diversity, as well as four dissimilarity indices, at common sample coverage levels.
#' Estimates can be obtained at default or user-specified sample sizes or coverage levels. See Chao et al. (2023b) for the theory.
#'
#' @param data For \code{datatype = "abundance"}, data for a single dataset can be provided as a species-by-assemblage abundance matrix/data.frame, with species in rows and assemblages in columns. Data for multiple datasets should be provided as a list of such matrices/data.frames, with each element representing one dataset.\cr
#' For \code{datatype = "incidence_raw"}, data for a single dataset should be provided as a list of species-by-sampling-unit incidence matrices/data.frames, with each element representing one assemblage. Data for multiple datasets should be provided as a list of such lists, with each outer-list element representing one dataset.
#'
#' @param diversity Selection of the diversity type: \code{diversity = "TD"} for taxonomic diversity, \code{diversity = "PD"} for phylogenetic diversity, and \code{diversity = "FD"} for functional diversity.
#'
#' @param q A numerical vector specifying the diversity orders. Default is \code{c(0, 1, 2)}.
#'
#' @param datatype Type of input data: individual-based abundance data (\code{datatype = "abundance"}) or sampling-unit-based incidence data (\code{datatype = "incidence_raw"}). For incidence raw data, each species-by-sampling-unit matrix contains only 0 (non-detection) or 1 (detection).
#'
#' @param base Standardization base: coverage-based rarefaction and extrapolation for gamma, alpha, and beta diversity, as well as four classes of dissimilarity indices (\code{base = "coverage"}), or size-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}). Default is \code{base = "coverage"}.
#'
#' @param level A numeric vector specifying the target sample coverage levels (when \code{base = "coverage"}) or sample sizes (when \code{base = "size"}) used to compute standardized diversity and dissimilarity estimates. For coverage-based standardization, the specified coverage levels must be between 0 and 1. Setting \code{level = 1} requests estimates at complete coverage (i.e., asymptotic estimates), when supported by the estimator.\cr
#' By default, when \code{base = "coverage"}, the function initially uses coverage levels ranging from 0.5 to 1 in increments of 0.025, together with the estimated gamma and alpha coverage levels at the reference sample size and twice the reference sample size. For \code{q = 0}, estimates corresponding to sample sizes greater than twice the reference sample size are excluded. Thus, the extrapolation limits correspond to \code{SC(2n, gamma)} for gamma diversity and \code{SC(2n, alpha)} for alpha and beta diversity, as well as dissimilarity measures. For \code{q = 1} and \code{q = 2}, estimates can be extended to complete coverage when supported by the data. The reference coverage for beta diversity and dissimilarity is defined as that for alpha diversity.\cr
#' When \code{base = "size"}, the function computes standardized 3D gamma and alpha diversity estimates at sample-size knots spanning rarefaction, the reference sample size, and extrapolation up to twice the reference sample size.
#'
#' @param nboot A non-negative integer specifying the number of bootstrap replications used to assess sampling uncertainty and construct confidence intervals. Bootstrap replications can be time-consuming. Set \code{nboot = 0} to skip the bootstrap procedures. Default is \code{nboot = 10}. For more accurate uncertainty estimates, a larger number of replications (e.g., 100 or 200) is recommended.
#'
#' @param conf A number between 0 and 1 specifying the confidence level. Default is \code{conf = 0.95}.
#'
#' @param PDtree Required when \code{diversity = "PD"}. A phylogenetic tree (an R \code{phylo} object) spanning all observed species across the input datasets. The tree can be imported from a file in Newick format.
#'
#' @param PDreftime Argument only for \code{diversity = "PD"} and \code{PDtype = "meanPD"}. A numeric value specifying a fixed reference time for phylogenetic diversity. If \code{PDreftime = NULL} (default), the root age of the phylogenetic tree spanning the pooled observed species is used.
#'
#' @param PDtype Argument only for \code{diversity = "PD"}. Specifies the phylogenetic diversity measure: \code{PDtype = "meanPD"} for the effective number of equally divergent lineages at a specified reference time, or \code{PDtype = "AUC"} for an overall measure obtained by integrating and averaging mean phylogenetic diversity over all reference times from zero to the root age of the phylogenetic tree. Default is \code{PDtype = "AUC"}.
#'
#' @param PDcut_number Argument only for \code{diversity = "PD"} and \code{PDtype = "AUC"}. An integer specifying the number of equally spaced integration points used to compute the AUC over the interval from zero to the root age of the phylogenetic tree. Default is \code{PDcut_number = 30}, and the minimum is 3.
#'
#' @param FDdistM Required when \code{diversity = "FD"}. A species pairwise distance matrix covering all species observed across the input datasets.
#'
#' @param FDtype Argument only for \code{diversity = "FD"}. Specifies the functional diversity measure: \code{FDtype = "tau_value"} for FD at a specified threshold of functional distinctiveness, or \code{FDtype = "AUC"} for an overall FD measure obtained by integrating the tau-profile over threshold values from 0 to 1. Default is \code{FDtype = "AUC"}.
#'
#' @param FDtau Argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}. A numeric value between 0 and 1 specifying the threshold of functional distinctiveness used to compute FD. If \code{FDtau = NULL} (default), the threshold is set to the mean pairwise functional distance between two individuals randomly selected from the pooled data across all datasets (i.e., Rao's quadratic entropy).
#'
#' @param FDcut_number Argument only for \code{diversity = "FD"} and \code{FDtype = "AUC"}. An integer specifying the number of equally spaced integration points used to compute the AUC over the interval from 0 to 1. Default is \code{FDcut_number = 30}, and the minimum is 3.
#'
#' @param by_pair A logical variable specifying whether diversity decomposition is performed for all assemblages as a whole or for every pair of assemblages. If \code{by_pair = FALSE}, gamma, alpha, and beta diversity, as well as dissimilarity, are computed for all assemblages within each dataset considered as a whole. If \code{by_pair = TRUE}, the corresponding estimates are computed for every pair of assemblages. Default is \code{by_pair = FALSE}.
#'
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import iNEXT.3D
#' @import future.apply
#' @import tibble
#' @import dplyr
#' @import tidytree
#' @importFrom tidyr gather
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom stats rmultinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' @importFrom grDevices hcl
#' @importFrom utils combn
#'
#' @return
#' For \code{base = "coverage"}, \code{iNEXTbeta3D()} returns an \code{iNEXTbeta3D} object containing seven data frames for each dataset: three for gamma, alpha, and beta diversity, and four for dissimilarity measures (\code{1-C}, \code{1-U}, \code{1-V}, and \code{1-S}).
#' For \code{base = "size"}, the object contains two data frames for each dataset, corresponding to gamma and alpha diversity.
#' Size-based beta diversity and dissimilarity measures are not provided.
#'
#' For \code{base = "coverage"}, the output columns are described below:
#'
#' \describe{
#'   \item{Dataset}{The name of the dataset.}
#'   \item{Pair}{The pair of assemblages being compared; this column is included only when \code{by_pair = TRUE}.}
#'   \item{Order.q}{The diversity order q.}
#'   \item{SC}{The target standardized sample coverage value.}
#'   \item{Size / mT}{The corresponding sample size for gamma and alpha diversity: \code{Size} for abundance data and \code{mT} for incidence data. This column is not included for beta diversity or dissimilarity measures.}
#'   \item{Gamma / Alpha / Beta / Dissimilarity}{The estimated diversity or dissimilarity at the target standardized coverage value.
#'     For \code{PDtype = "meanPD"}, the gamma output includes \code{Gamma_PD} and \code{Gamma_meanPD}, whereas the alpha output includes \code{Alpha_PD} and \code{Alpha_meanPD}.
#'     These quantities satisfy \code{Gamma_PD = Gamma_meanPD * Reftime} and \code{Alpha_PD = Alpha_meanPD * Reftime}. The PD-based and meanPD-based beta diversity and dissimilarity estimates are identical after multiplicative decomposition; therefore, only one set of estimates is reported.}
#'   \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample. For beta diversity and dissimilarity measures, the reference coverage is defined as that for alpha diversity.}
#'   \item{s.e.}{The standard error of the estimated diversity or dissimilarity.}
#'   \item{LCL, UCL}{The bootstrap lower and upper confidence limits for the estimated diversity or dissimilarity at the specified confidence level (default: \code{conf = 0.95}).}
#'   \item{Diversity}{The diversity type:
#'     \code{"TD"} = taxonomic diversity;
#'     \code{"meanPD"} = mean phylogenetic diversity;
#'     \code{"PD_AUC"} = phylogenetic diversity (AUC);
#'     \code{"FD_tau"} = functional diversity (specified tau);
#'     \code{"FD_AUC"} = functional diversity (AUC).}
#'   \item{Reftime}{The reference time used for phylogenetic diversity when \code{PDtype = "meanPD"}.}
#'   \item{Tau}{The threshold of functional distinctiveness used for functional diversity when \code{FDtype = "tau_value"}.}
#' }
#'
#' Similar output columns are provided for \code{base = "size"}, where applicable.
#'  
#'  
#' @examples
#' 
#' ## (R/E Analysis) Taxonomic diversity for abundance data
#' # Coverage-based standardized TD estimates and related statistics
#' # See Example 1 in the iNEXT.beta3D vignette for details and graphical output.
#' 
#' data(Brazil_rainforests)
#' output_TDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'TD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10)
#' output_TDc_abun
#'
#' 
#'
#'
#' # Size-based standardized TD estimates and related statistics
#' 
#' data(Brazil_rainforests)
#' output_TDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'TD', 
#'                               datatype = 'abundance', base = "size", nboot = 10)
#' output_TDs_abun
#'
#'
#' 
#' 
#' \donttest{
#' ## (R/E Analysis) Taxonomic diversity for incidence data
#' # Coverage-based standardized TD estimates and related statistics (all assemblages as a whole)
#' # See Example 2 in the iNEXT.beta3D vignette for details and graphical output.
#' 
#' data(Second_growth_forests)
#' output_TDc_inci = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                               datatype = 'incidence_raw', base = "coverage", nboot = 10)
#' output_TDc_inci
#' 
#' 
#' 
#' 
#' # Size-based standardized TD estimates and related statistics (all assemblages as a whole)
#' 
#' data(Second_growth_forests)
#' output_TDs_inci = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                               datatype = 'incidence_raw', base = "size", nboot = 10)
#' output_TDs_inci
#' 
#' 
#' 
#' 
#' # Coverage-based standardized TD estimates and related statistics for every pair of assemblages
#' # See Example 3 in the iNEXT.beta3D vignette for details and graphical output.
#' 
#' data(Second_growth_forests)
#' output_TDc_inci_pair = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                                      datatype = 'incidence_raw', base = "coverage", 
#'                                      nboot = 10, by_pair = TRUE)
#' output_TDc_inci_pair
#' 
#'
#' # Size-based standardized TD estimates and related statistics for every pair of assemblages
#' 
#' data(Second_growth_forests)
#' output_TDs_inci_pair = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                                      datatype = 'incidence_raw', base = "size", 
#'                                      nboot = 10, by_pair = TRUE)
#' output_TDs_inci_pair
#' 
#' 
#' }
#' 
#' \dontrun{
#' ## (R/E Analysis) Phylogenetic diversity for abundance data
#' # Coverage-based standardized PD estimates and related statistics 
#' # See Example 4 in the iNEXT.beta3D vignette for details and graphical output.
#' 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output_PDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10, 
#'                               PDtree = Brazil_tree, PDtype = 'AUC', PDcut_number = 30) 
#' output_PDc_abun
#' 
#' 
#' 
#' 
#' # Size-based standardized PD estimates and related statistics 
#' 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output_PDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                               datatype = 'abundance', base = "size", nboot = 10, 
#'                               PDtree = Brazil_tree, PDtype = 'AUC', PDcut_number = 30) 
#' output_PDs_abun
#' 
#' 
#' 
#' 
#' ## (R/E Analysis) Functional diversity for abundance data when all thresholds from 0 to 1 
#' ## are considered
#' # Coverage-based standardized FD estimates and related statistics 
#' # See Example 5 in the iNEXT.beta3D vignette for details and graphical output.
#' 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output_FDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10, 
#'                               FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' output_FDc_abun
#' 
#' 
#' 
#' 
#' # Size-based standardized FD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output_FDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                               datatype = 'abundance', base = "size", nboot = 10, 
#'                               FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' output_FDs_abun
#' 
#' }
#' 
#' 
#' @references
#' Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zelený, D., Colwell, R. K., and Magurran, A. E. (2023b). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs, e1588. https://doi.org/10.1002/ecm.1588
#' @export
iNEXTbeta3D = function(data, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', 
                       base = 'coverage', level = NULL, nboot = 10, conf = 0.95, 
                       PDtree = NULL, PDreftime = NULL, PDtype = 'AUC', PDcut_number = 30,
                       FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 30,
                       by_pair = FALSE) {
  
  ## Check parameter setting
  if (is.na(pmatch(diversity, c("TD", "PD", "FD")))) stop("invalid diversity")
  
  if (!inherits(q, "numeric"))
    stop("invlid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  if (min(q) < 0){
    warning("ambigous of order q, we only compute postive q.", call. = FALSE)
    q <- q[q >= 0]
  }
  
  if (datatype == "incidence" | datatype == "incidence_freq") stop('Please use datatype = "incidence_raw".')  
  if(is.na(pmatch(datatype, c("abundance", "incidence_raw"))))
    stop("invalid datatype")
  
  if (is.na(pmatch(base, c("size", "coverage")))) stop("invalid datatype")
  
  if (! (is.null(level) | inherits(level, c("numeric", "integer", "double"))))
    stop("invlid class of level, level should be a postive value/vector of numeric object", call. = FALSE)
  
  
  if ((nboot < 0) | (is.numeric(nboot) == F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  if (! (is.null(PDreftime) | inherits(PDreftime, "numeric")))
    stop("invalid class of reference time, PDreftime should be a postive value of numeric object.", call. = FALSE)
  if (length(PDreftime) > 1)
    stop("PDreftime can only accept a value instead of a vector.", call. = FALSE)
  
  if(is.na(pmatch(PDtype, c("meanPD", "AUC"))))
    stop("Incorrect phylogenetic diversity type. Please use 'meanPD', or 'AUC'.", call. = FALSE)
  
  if (PDtype == "AUC") {if (!is.numeric(PDcut_number) || length(PDcut_number) != 1 || is.na(PDcut_number) || PDcut_number < 3 || PDcut_number %% 1 != 0) {stop("Invalid PDcut_number. PDcut_number must be an integer of at least 3.",call. = FALSE)}}
  
  if (FDtype == "tau_values") stop('Please try FDtype = "tau_value".')  
  if(is.na(pmatch(FDtype, c("AUC", "tau_value"))))
    stop("Incorrect type of functional diversity type, please use either 'AUC' or 'tau_value'", call. = FALSE)
  
  if (! (is.null(FDtau) | inherits(FDtau, "numeric")))
    stop("invalid class of tau value, FDtau should be a postive value between zero and one.", call. = FALSE)
  if (length(FDtau) > 1)
    stop("FDtau only accept a value instead of a vector.", call. = FALSE)
  
  if (FDtype == "AUC") {
    if (!is.numeric(FDcut_number) ||
        length(FDcut_number) != 1 ||
        is.na(FDcut_number) ||
        !is.finite(FDcut_number) ||
        FDcut_number < 3 ||
        FDcut_number %% 1 != 0) {
      stop("FDcut_number must be an integer of at least 3.", call. = FALSE)
    }
  }
  
  ##
  
  
  if (datatype == 'abundance') {
    
    if ( inherits(data, "data.frame") | inherits(data, "matrix") ) data = list(Dataset_1 = data)
    
    if ( inherits(data, "list") ){
      
      if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
      if (by_pair == FALSE) Ns = sapply(data, ncol) else Ns = sapply(data, function(x) 2)
      data_list = data
      
    }
    
  }
  
  if (datatype == 'incidence_raw') {
    
    if (!inherits(data, "list"))
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( inherits(data, "list") & (inherits(data[[1]], "data.frame") | inherits(data[[1]], "matrix")) ) data = list(Dataset_1 = data)
    
    if (! ( inherits(data[[1]][[1]], "data.frame") | inherits(data[[1]][[1]], "matrix") ) )
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( sum( sapply(1:length(data), function(i) ( length(unique(sapply(data[[i]], nrow))) != 1 | 
                                                   length(unique(sapply(data[[i]], ncol))) != 1 ) ) ) > 0 )
      stop("Number of species (row) or sampling units (column) should be the same within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
    
    data = lapply(data, function(x) {
      
      if (is.null(rownames(x[[1]]))) tmp = x else {
        
        nT = ncol(x[[1]])
        if (nT <= 3) stop("Number of sampling units of some datasets is too less. Please add more sampling units data.", call. = FALSE)
        
        tmp = lapply(x, function(i) data.frame(i) %>% rownames_to_column(var = "Species"))
        
        tmp1 = tmp[[1]]
        for (i in 2:length(tmp)) {
          tmp1 = full_join(tmp1, tmp[[i]], by = "Species")
        }
        tmp1 = tmp1 %>% column_to_rownames(var = "Species")
        
        if (nrow(tmp1) != nrow(x[[1]]))
          stop("Species names (rownames) should be matched within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
        
        if (sum(!as.matrix(tmp1) %in% c(0,1)) > 0)
          stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
        
        tmp = lapply(1:length(tmp), function(i) tmp1[, ((i-1)*nT+1):(i*nT)])
        names(tmp) = if (is.null(data)) paste0("Assemblage_", 1:length(data)) else names(x)
      }
      
      return(tmp)
    })
    
    if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
    if (by_pair == FALSE) Ns = sapply(data, length) else Ns = sapply(data, function(x) 2)
    data_list = data
  }
  
  
  if (datatype == 'abundance') {
    
    pool.name <- lapply(data_list, function(x) rownames(x)) %>% unlist %>% unique
    
  } else if (datatype == 'incidence_raw') {
    
    pool.name <- lapply(data_list, function(x) lapply(x, function(y) rownames(y))) %>% unlist %>% unique
    
  }
  
  if (diversity == "PD") {
    
    if (sum(c(duplicated(PDtree$tip.label), duplicated(PDtree$node.label[PDtree$node.label!=""])))>0)
      stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
    
    if ( is.null(pool.name) )
      stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% PDtree$tip.label) != length(pool.name))
      stop("Data and tree tip label contain unmatched species", call. = FALSE)
  }
  
  if (diversity == "FD") {
    
    if (is.null(rownames(FDdistM)))
      stop('The species names are not provided in distance matrix.', call. = FALSE)
    
    if( is.null(pool.name) )
      stop("Row names of data must be the species names that match row names in distance matrix and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% rownames(FDdistM)) != length(pool.name))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  ##
  
  
  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  
  trunc = ifelse(is.null(level), T, F)
  
  if (base == 'coverage' ) {
    
    if ( is.null(level) ) {
      
      level <- lapply(1:length(data_list), function(i) {
        
        if (by_pair == FALSE) {
          
          level = seq(0.5, 1, 0.025)
          
          if(datatype == "abundance") {
            
            n = sum(data_list[[i]])
            data_gamma = rowSums(data_list[[i]])
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_list[[i]]) %>% as.vector
            
            ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
            ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
            ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
            ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
            
          }else if(datatype == "incidence_raw"){
            
            n = unique(sapply(data_list[[i]], ncol))
            gamma = Reduce('+', data_list[[i]])
            gamma[gamma>1] = 1
            data_gamma_raw = gamma
            data_gamma_freq = c(n, rowSums(gamma))
            
            data_alpha_freq = sapply(data_list[[i]], rowSums) %>% c(n, .)
            
            ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
            ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
            ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
            ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
            
          }
          
          c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
          
        } else {
          
          level = seq(0.5, 1, 0.025)
          
          if(datatype == "abundance") {
            
            two_idx <- combn(1:ncol(data_list[[i]]), 2)
            
            lapply(1:ncol(two_idx), function(j) {
              
              n = sum(data_list[[i]][,two_idx[,j]])
              data_gamma = rowSums(data_list[[i]][,two_idx[,j]])
              data_gamma = data_gamma[data_gamma>0]
              data_alpha = as.matrix(data_list[[i]][,two_idx[,j]]) %>% as.vector
              
              ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
              ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
              ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
              ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
              
              c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
              
            })
            
            
          }else if(datatype == "incidence_raw"){
            
            two_idx <- combn(1:length(data_list[[i]]), 2)
            
            lapply(1:ncol(two_idx), function(j) {
              
              n = unique(sapply(data_list[[i]][two_idx[,j]], ncol))
              gamma = Reduce('+', data_list[[i]][two_idx[,j]])
              gamma[gamma>1] = 1
              data_gamma_raw = gamma
              data_gamma_freq = c(n, rowSums(gamma))
              
              data_alpha_freq = sapply(data_list[[i]][two_idx[,j]], rowSums) %>% c(n, .)
              
              ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
              ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
              ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
              ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
              
              c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
              
            })
            
          }
          
        }
        
      })
      
    } else {
      
      if (inherits(level, "numeric") | inherits(level, "integer") | inherits(level, "double")) {
        level <- list(level = level)
      }
      
      if ( (length(level) != length(data_list)) & (by_pair == FALSE) ) level <- lapply(1:length(data_list), function(x) level[[1]])
      
      
      
      if ( (length(level) != length(data_list)) & (by_pair == TRUE) ) level <- lapply(1:length(data_list), function(x) {
        
        if (datatype == "abundance") two_idx <- combn(1:ncol(data_list[[x]]), 2)
        if (datatype == "incidence_raw") two_idx <- combn(1:length(data_list[[x]]), 2)
        
        lapply(1:ncol(two_idx), function(y) level[[1]])
      })
      
      }
    
    } else if ( base == 'size' ) {
    
      if ( is.null(level) ) {
        
        if (by_pair == FALSE) {
          
          if (datatype == "abundance") {
            
            endpoint <- sapply(data_list, function(x) 2*sum(x))
            
          } else if (datatype == "incidence_raw") {
            
            endpoint <- sapply(data_list, function(x) 2*ncol(x[[1]]))
            
          }
          
          level <- lapply(1:length(data_list), function(i) {
            
            if(datatype == "abundance") {
              
              ni <- sum(data_list[[i]])
              
            }else if(datatype == "incidence_raw"){
              
              ni <- ncol(data_list[[i]][[1]])
            }
            
            mi <- floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint[i], length.out = 20)))
          })
          
        } else {
          
          level <- lapply(1:length(data_list), function(i) {
            
            if(datatype == "abundance") {
              
              two_idx <- combn(1:ncol(data_list[[i]]), 2)
              
              mi <- lapply(1:ncol(two_idx), function(j) {
                
                endpoint <- 2*sum(data_list[[i]][,two_idx[,j]])
                
                ni <- sum(data_list[[i]][,two_idx[,j]])
                
                floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint, length.out = 20)))
                
              })
              
            }else if(datatype == "incidence_raw"){
              
              two_idx <- combn(1:length(data_list[[i]]), 2)
              
              mi <- lapply(1:ncol(two_idx), function(j) {
                
                endpoint <- 2*ncol(data_list[[i]][[two_idx[1,j]]])
                
                ni <- ncol(data_list[[i]][[two_idx[1,j]]])
                
                floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint, length.out = 20)))
                
              })
            }
            
          })
          
        }
        
      } else {
        
        if (inherits(level, "numeric") | inherits(level, "integer") | inherits(level, "double")) {
          level <- list(level = level)
        }
        
        if ( (length(level) != length(data_list)) & (by_pair == FALSE) ) level <- lapply(1:length(data_list), function(x) level[[1]])
        
        if ( (length(level) != length(data_list)) & (by_pair == TRUE) ) level <- lapply(1:length(data_list), function(x) {
          
          if (datatype == "abundance") two_idx <- combn(1:ncol(data_list[[x]]), 2)
          if (datatype == "incidence_raw") two_idx <- combn(1:length(data_list[[x]]), 2)
          
          lapply(1:ncol(two_idx), function(y) level[[1]])
        })
          
      }
  }
  
  if (diversity == 'FD' & FDtype == 'tau_value' & is.null(FDtau) == T) {
    
    if (datatype == 'abundance') {
      
      tt <- lapply(data_list, rowSums)
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    } else if (datatype == 'incidence_raw') {
      
      tt <- lapply(data_list, function(x) {tt = Reduce('+', x); tt[tt > 1] = 1; rowSums(tt) })
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    }
    
    FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
  }
  
  if (diversity == 'PD') {
    
    if (datatype == "abundance") 
      
      if (length(data_list) > 1) {
        
        pool.data = data_list[[1]] %>% data.frame %>% rownames_to_column()
        for (i in 2:length(data_list)) 
          pool.data = full_join(pool.data, data_list[[i]] %>% data.frame %>% rownames_to_column(), 'rowname')
        pool.data[is.na(pool.data)] = 0
        pool.data = pool.data %>% column_to_rownames() %>% rowSums
        
      } else pool.data = do.call(cbind, data_list) %>% rowSums
    
    if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums
    
    pool.name = names(pool.data[pool.data>0])
    tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
    mytree = ape::drop.tip(PDtree, tip)
    
    # H_max = get.rooted.tree.height(mytree)
    H_max = max(ape::node.depth.edgelength(mytree))
    
    
    if (PDtype == "meanPD") {
      if (is.null(PDreftime)) {
        
        reft = H_max
        
      } else if (PDreftime <= 0) {
        
        stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
        
      } else {
        
        reft = PDreftime
  
      }
    }
    
    if (PDtype == "AUC") {
      
      PD_AUC_max = H_max

      PDcut = seq(from = 1e-08,to = PD_AUC_max,length.out = PDcut_number)
      
      PDwidth = diff(PDcut)
      
      reft = PD_AUC_max
    }
    
  }
    

  for_each_dataset = function(data, dataset_name, N, level, by_pair = FALSE) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
      
      # level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level<1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all datasets must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
      
      # level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level < 1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq, i, datatype='incidence_freq'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq, i, datatype='incidence_freq'))
      
      # if (sum(m_gamma < 1 | m_alpha < 1) > 0) {
      #   index = m_gamma < 1 | m_alpha < 1
      #   level = level[!index]
      #   m_gamma = m_gamma[!index]
      #   m_alpha = m_alpha[!index]
      # }
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        alpha = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        alpha = lapply(1:length(level), function(i){
          estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        
        
      }
      
      gamma = (cbind(SC = rep(level, each=length(q)), gamma[,-c(1,3,7,8,9)]) %>% 
                 mutate(Method = ifelse(SC>=ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,4,1,3)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      
      alpha = (cbind(SC = rep(level, each = length(q)), alpha[,-c(1,3,7,8,9)]) %>% 
                 mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,4,1,3)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Estimate = alpha$Estimate / N
      
      
      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate), 
      #                                Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q==1, log(Estimate)/log(N), (Estimate^(1-Order.q) - 1)/(N^(1-Order.q)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q==1, log(Estimate)/log(N), (Estimate^(Order.q-1) - 1)/(N^(Order.q-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
            
            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
            
            gamma = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
            alpha = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
          }
          
          if (datatype == 'incidence_raw') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            
            raw = lapply(1:ncol(bootstrap_population), function(j){
              
              lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            gamma = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
            alpha = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
          }
          
          gamma = gamma$qTD
          
          alpha = alpha$qTD
          alpha = alpha / N
          
          beta = gamma/alpha
          
          Order.q = rep(q, length(level))
          
          beta = data.frame(Estimate=beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1,log(Estimate)/log(N),(Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1,log(Estimate)/log(N),(Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
      if (nboot>0 & datatype == 'abundance') {
        gamma.se = estimate3D(data_gamma, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        
        alpha.se = estimate3D(data_alpha, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N
        
      } else if (nboot>0 & datatype == 'incidence_raw') {
        
        gamma.se = estimate3D(data_gamma_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        
        alpha.se = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N
        
      } else {
        gamma.se = alpha.se = NA
      }
      
      se[1:( length(level) * length(q) ), 'gamma'] = gamma.se
      
      se[1:( length(level) * length(q) ), 'alpha'] = alpha.se
      
    }
    
    if (diversity == 'PD') {
      
      PD_by_time = function(data, tree, reft_i, datatype, m_gamma, m_alpha) {
        
        if (datatype == "abundance") {
          
          data_gamma_i = rowSums(data)
          data_gamma_i = data_gamma_i[data_gamma_i > 0]
          
          # Gamma
          aL_gamma = iNEXT.3D:::phyBranchAL_Abu(phylo = tree, data = data_gamma_i, rootExtend = TRUE, refT = reft_i)
          
          aL_gamma$treeNabu$branch.length = aL_gamma$BLbyT[, 1]
          
          aL_table_gamma = aL_gamma$treeNabu %>%
            dplyr::select(branch.abun, branch.length, tgroup)
          
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length),
                                       m = m_gamma, q = q, nt = n, reft = reft_i, cal = "PD") %>% t() %>% as.vector()
          
          
          # Alpha
          aL_table_alpha = NULL
          
          for (k in seq_len(N)) {
            
            x = data[data[, k] > 0, k]
            
            names(x) = rownames(data)[data[, k] > 0]
            
            aL_alpha = iNEXT.3D:::phyBranchAL_Abu(phylo = tree, data = x, rootExtend = TRUE, refT = reft_i)
            
            aL_alpha$treeNabu$branch.length = aL_alpha$BLbyT[, 1]
            
            aL_table_k = aL_alpha$treeNabu %>%
              dplyr::select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table_k)
          }
          
          alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length),
                                        m = m_alpha, q = q, nt = n, reft = reft_i, cal = "PD") / N
          ) %>% t() %>% as.vector()
        }
        
        
        if (datatype == "incidence_raw") {
          
          gamma_raw = Reduce("+", data)
          gamma_raw[gamma_raw > 1] = 1
          
          # Gamma
          aL_gamma = iNEXT.3D:::phyBranchAL_Inc(phylo = tree, data = as.matrix(gamma_raw), datatype = "incidence_raw",
                                                rootExtend = TRUE, refT = reft_i)
          
          aL_gamma$treeNabu$branch.length = aL_gamma$BLbyT[, 1]
          
          aL_table_gamma = aL_gamma$treeNabu %>%
            dplyr::select(branch.abun, branch.length, tgroup)
          
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun,Lis = as.matrix(aL_table_gamma$branch.length),
                                       m = m_gamma, q = q, nt = n, reft = reft_i, cal = "PD") %>% t() %>% as.vector()
          
          
          # Alpha
          aL_table_alpha = NULL
          
          for (k in seq_len(N)) {
            
            x = data[[k]]
            
            aL_alpha = iNEXT.3D:::phyBranchAL_Inc(phylo = tree, data = x, datatype = "incidence_raw", rootExtend = TRUE, refT = reft_i)
            
            aL_alpha$treeNabu$branch.length = aL_alpha$BLbyT[, 1]
            
            aL_table_k = aL_alpha$treeNabu %>% dplyr::select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table_k)
          }
          
          alpha = (
            iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length),
                                 m = m_alpha, q = q, nt = n, reft = reft_i, cal = "PD") / N) %>% t() %>% as.vector()
        }
        
        gamma = gamma / reft_i
        alpha = alpha / reft_i
        
        return(data.frame(gamma = gamma, alpha = alpha))
      }
      
      PD_AUC_by_data = function(data, tree, datatype, m_gamma, m_alpha) {
        
        gamma_alpha_over_time = lapply(
          PDcut,
          function(reft_i) {
            PD_by_time(
              data = data,
              tree = tree,
              reft_i = reft_i,
              datatype = datatype,
              m_gamma = m_gamma,
              m_alpha = m_alpha
            )
          }
        )
        
        gamma_over_time = sapply(
          gamma_alpha_over_time,
          function(x) x$gamma
        )
        
        alpha_over_time = sapply(
          gamma_alpha_over_time,
          function(x) x$alpha
        )
        
        if (is.vector(gamma_over_time)) {
          gamma_over_time = matrix(gamma_over_time, nrow = 1)
        }
        
        if (is.vector(alpha_over_time)) {
          alpha_over_time = matrix(alpha_over_time, nrow = 1)
        }
        
        # Beta: first gamma / alpha at each reference time
        beta_over_time = gamma_over_time / alpha_over_time
        
        # Left and right endpoints for trapezoidal integration
        gamma_left  = gamma_over_time[, -ncol(gamma_over_time), drop = FALSE]
        gamma_right = gamma_over_time[, -1, drop = FALSE]
        
        alpha_left  = alpha_over_time[, -ncol(alpha_over_time), drop = FALSE]
        alpha_right = alpha_over_time[, -1, drop = FALSE]
        
        beta_left   = beta_over_time[, -ncol(beta_over_time), drop = FALSE]
        beta_right  = beta_over_time[, -1, drop = FALSE]
        
        width_matrix = matrix(
          PDwidth,
          nrow = nrow(gamma_over_time),
          ncol = length(PDwidth),
          byrow = TRUE
        )
        
        gamma = rowSums(
          ((gamma_left + gamma_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        alpha = rowSums(
          ((alpha_left + alpha_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        beta = rowSums(
          ((beta_left + beta_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        return(
          list(
            gamma = gamma,
            alpha = alpha,
            beta = beta
          )
        )
      }
      

      
      
      if(PDtype == "AUC"){
        
        PD_auc <- PD_AUC_by_data(data = data, tree = PDtree, datatype = datatype, m_gamma = m_gamma, m_alpha = m_alpha)
        
        gamma = data.frame(Order.q = rep(q, each = length(level)), Estimate = PD_auc$gamma,
                           SC = rep(level, times = length(q)), Size = rep(m_gamma, times = length(q)))
        
        alpha = data.frame(Order.q = rep(q, each = length(level)), Estimate = PD_auc$alpha,
                           SC = rep(level, times = length(q)), Size = rep(m_alpha, times = length(q)))
        

      }else{
        
              if (datatype == 'abundance') {
        
        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_gamma, length(q)))
        
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        
        qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, reft = reft, cal = "PD")
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))
        
      }
      
              if (datatype == 'incidence_raw') {
        
        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_gamma, length(q)))
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[[i]]
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, reft = reft, cal = "PD")/N) %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))
        
        
      }
        
      }
      
   
      
      gamma = (gamma %>% 
                 mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      gamma$Order.q = as.numeric(gamma$Order.q)
      
      
      alpha = (alpha %>% 
                 mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Order.q = as.numeric(alpha$Order.q)
      
      if (PDtype == 'meanPD') {
        gamma$PD = gamma$Estimate
        alpha$PD = alpha$Estimate
        
        gamma$Estimate = gamma$PD / reft
        alpha$Estimate = alpha$PD / reft
      }
      
      beta = alpha
      
      if (PDtype == "AUC") {
        
        beta$Estimate = PD_auc$beta
        
      } else {
        
        beta$Estimate = gamma$Estimate / alpha$Estimate
        
      }
      
      
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate), 
      #                                Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          
          if (datatype == "abundance") {
            tree_bt = PDtree
            bootstrap_population = bootstrap_population_multiple_assemblage(data, 
                                                                                           data_gamma, "abundance")
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
            if (nrow(p_bt) > nrow(data) & sum(unseen_p) > 
                0) {
              unseen = unseen_p[which(rowSums(unseen_p) > 
                                        0), ]
              unseen = matrix(unseen, ncol = ncol(unseen_p))
              p_bt = rbind(p_bt[(1:nrow(data)), ], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_", 
                                                                      i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              bootstrap_sample = sapply(1:ncol(data), 
                                        function(k) rmultinom(n = 1, size = sum(data[, 
                                                                                     k]), prob = p_bt[, k]))
              x_bt = bootstrap_sample
              rownames(x_bt) = rownames(p_bt)
              if (sum(x_bt[-(1:nrow(data)), ]) > 0) {
                g0_hat = apply(data, 2, function(x) {
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, 
                                                  data = x, rootExtend = T, refT = reft)
                  aL$treeNabu$branch.length = aL$BLbyT[, 
                                                       1]
                  aL = aL$treeNabu %>% select(branch.abun, 
                                              branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 
                                          1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 
                                          2] %>% sum
                  g0_hat = ifelse(g2 > ((g1 * f2)/(2 * 
                                                     f1)), ((n - 1)/n) * (g1^2/(2 * g2)), 
                                  ((n - 1)/n) * (g1 * (f1 - 1)/(2 * 
                                                                  (f2 + 1))))
                  if (is.na(g0_hat)) {
                    g0_hat <- 0
                  }
                  g0_hat
                })
                te = (x_bt[1:nrow(data), ] * (data == 
                                                0)) > 0
                used_length = sapply(1:ncol(data), function(i) {
                  if (sum(te[, i]) == 0) 
                    return(0)
                  else {
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, 
                                               data = x_bt[1:nrow(data), i], 
                                               rootExtend = T, refT = reft)$treeNabu %>% 
                      subset(label %in% names(which(te[, 
                                                       i] == TRUE))) %>% select(branch.length) %>% 
                      sum
                  }
                })
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                unseen_sample = x_bt[-(1:nrow(data)), 
                ]
                if (is.vector(unseen_sample)) 
                  unseen_sample = matrix(unseen_sample, 
                                         ncol = ncol(x_bt))
                L0_hat = sapply(1:length(g0_hat), function(i) if (sum(unseen_sample[, 
                                                                                    i] > 0) > 0) 
                  (g0_hat[i]/nrow(unseen))
                  else 0)
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), 
                                         ncol(unseen_sample), byrow = T) * 
                                    unseen_sample))/rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 
                               0)] = 0
                for (i in 1:length(L0_hat)) {
                  tip = list(edge = matrix(c(2, 1), 
                                           1, 2), tip.label = unseen_name[i], 
                             edge.length = L0_hat[i], Nnode = 1)
                  class(tip) = "phylo"
                  tree_bt = tree_bt + tip
                }
              }
              else {
                x_bt = x_bt[1:nrow(data), ]
                p_bt = p_bt[1:nrow(data), ]
              }
            }
            else {
              p_bt = p_bt[1:nrow(data), ]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, 
                                                                size = sum(data[, k]), prob = p_bt[, 
                                                                                                   k]))
              rownames(x_bt) = rownames(data)
            }
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 
                                                          0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% 
              as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 
                                                          0]
            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, 
                                                                                i, datatype = "abundance"))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, 
                                                                                i, datatype = "abundance"))
            
            
            if (PDtype == "AUC") {
              
              PD_auc_bt = PD_AUC_by_data(data = x_bt, tree = tree_bt, datatype = "abundance",
                                         m_gamma = m_gamma, m_alpha = m_alpha)
              
              gamma = PD_auc_bt$gamma
              alpha = PD_auc_bt$alpha
              beta  = PD_auc_bt$beta
              
            }else{
              
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, 
                                              data = bootstrap_data_gamma, rootExtend = T, 
                                              refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[, 1]
              aL_table_gamma = aL$treeNabu %>% select(branch.abun, 
                                                      branch.length, tgroup)
              gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, 
                                                     Lis = as.matrix(aL_table_gamma$branch.length), 
                                                     m = m_gamma, q = q, nt = n, reft = reft, 
                                                     cal = "PD") %>% t)
              aL_table_alpha = c()
              
              for (i in 1:N) {
                x = x_bt[, i]
                names(x) = rownames(p_bt)
                x = x[x_bt[, i] > 0]
                aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, 
                                                data = x, rootExtend = T, refT = reft)
                aL$treeNabu$branch.length = aL$BLbyT[, 
                                                     1]
                aL_table = aL$treeNabu %>% select(branch.abun, 
                                                  branch.length, tgroup)
                aL_table_alpha = rbind(aL_table_alpha, 
                                       aL_table)
              }
              
              alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, 
                                                      Lis = as.matrix(aL_table_alpha$branch.length), 
                                                      m = m_alpha, q = q, nt = n, reft = reft, 
                                                      cal = "PD")/N) %>% t)
              
            }
            
            
          }
          if (datatype == "incidence_raw") {
            tree_bt = PDtree
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, 
                                                                                           data_gamma_freq, "incidence")
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])), ] %>% 
              matrix(ncol = N)
            if (nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p) > 
                0) {
              unseen = unseen_p[which(rowSums(unseen_p) > 
                                        0), ]
              unseen = matrix(unseen, ncol = ncol(unseen_p))
              p_bt = rbind(p_bt[(1:nrow(data[[1]])), 
              ], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_", 
                                                                      i))
              rownames(p_bt) = c(rownames(data[[1]]), 
                                 unseen_name)
              raw = lapply(1:ncol(p_bt), function(j) {
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, 
                                                        size = 1, prob = p_bt[i, j])) %>% 
                  do.call(rbind, .)
              })
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              if (lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])), 
              ]) %>% do.call(sum, .) > 0) {
                R0_hat = sapply(data, function(x) {
                  nT = ncol(x)
                  Q1 = sum(rowSums(x) == 1)
                  Q2 = sum(rowSums(x) == 2)
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, 
                                                  data = x, datatype = "incidence_raw", 
                                                  rootExtend = T, refT = reft)
                  aL$treeNabu$branch.length = aL$BLbyT[, 
                                                       1]
                  aL = aL$treeNabu %>% select(branch.abun, 
                                              branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 
                                          1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 
                                          2] %>% sum
                  R0_hat = ifelse(R2 > ((R1 * Q2)/(2 * 
                                                     Q1)), ((nT - 1)/nT) * (R1^2/(2 * 
                                                                                    R2)), ((nT - 1)/nT) * (R1 * (Q1 - 
                                                                                                                   1)/(2 * (Q2 + 1))))
                  if (is.na(R0_hat)) {
                    R0_hat <- 0
                  }
                  R0_hat
                })
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]), 
                ] * (sapply(data, rowSums) == 0)) > 
                  0
                used_length = sapply(1:N, function(i) {
                  if (sum(te[, i]) == 0) 
                    return(0)
                  else {
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, 
                                               data = raw[[i]][1:nrow(data[[1]]), 
                                               ], datatype = "incidence_raw", 
                                               rootExtend = T, refT = reft)$treeNabu %>% 
                      subset(label %in% names(which(te[, 
                                                       i] == TRUE))) %>% select(branch.length) %>% 
                      sum
                  }
                })
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])), 
                ]
                if (is.vector(unseen_sample)) 
                  unseen_sample = matrix(unseen_sample, 
                                         ncol = N)
                L0_hat = sapply(1:length(R0_hat), function(i) if (sum(unseen_sample[, 
                                                                                    i] > 0) > 0) 
                  (R0_hat[i]/nrow(unseen))
                  else 0)
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), 
                                         ncol(unseen_sample), byrow = T) * 
                                    unseen_sample))/rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 
                               0)] = 0
                for (i in 1:length(L0_hat)) {
                  tip = list(edge = matrix(c(2, 1), 
                                           1, 2), tip.label = unseen_name[i], 
                             edge.length = L0_hat[i], Nnode = 1)
                  class(tip) = "phylo"
                  tree_bt = tree_bt + tip
                }
              }
              else raw = lapply(raw, function(i) i[1:nrow(data[[1]]), 
              ])
            }
            else {
              p_bt = p_bt[1:nrow(data[[1]]), ]
              raw = lapply(1:ncol(p_bt), function(j) {
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, 
                                                        size = 1, prob = p_bt[i, j])) %>% 
                  do.call(rbind, .)
              })
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
            }
            gamma = Reduce("+", raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            bootstrap_data_alpha_freq = sapply(raw, 
                                               rowSums) %>% c(n, .)
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 
                                                                    0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 
                                                                    0]
            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma_freq, 
                                                                                i, datatype = "incidence_freq"))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha_freq, 
                                                                                i, datatype = "incidence_freq"))
            
            
            if (PDtype == "AUC") {
              
              PD_auc_bt = PD_AUC_by_data(data = raw, tree = tree_bt, datatype = "incidence_raw",
                                         m_gamma = m_gamma, m_alpha = m_alpha)
              
              gamma = PD_auc_bt$gamma
              alpha = PD_auc_bt$alpha
              beta  = PD_auc_bt$beta
              
            }else{
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, 
                                              data = bootstrap_data_gamma_raw, datatype = "incidence_raw", 
                                              rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[, 1]
              aL_table_gamma = aL$treeNabu %>% select(branch.abun, 
                                                      branch.length, tgroup)
              gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, 
                                                     Lis = as.matrix(aL_table_gamma$branch.length), 
                                                     m = m_gamma, q = q, nt = n, reft = reft, 
                                                     cal = "PD") %>% t)
              aL_table_alpha = c()
              for (i in 1:N) {
                x = raw[[i]]
                aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, 
                                                data = x, datatype = "incidence_raw", 
                                                rootExtend = T, refT = reft)
                aL$treeNabu$branch.length = aL$BLbyT[, 
                                                     1]
                aL_table = aL$treeNabu %>% select(branch.abun, 
                                                  branch.length, tgroup)
                aL_table_alpha = rbind(aL_table_alpha, 
                                       aL_table)
              }
              alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, 
                                                      Lis = as.matrix(aL_table_alpha$branch.length), 
                                                      m = m_alpha, q = q, nt = n, reft = reft, 
                                                      cal = "PD")/N) %>% t)
            }
            
          }
          
          
          
          if (PDtype == 'meanPD') {
            gamma = gamma / reft
            alpha = alpha / reft
          } 
          
          if (PDtype != "AUC") {
            beta = gamma / alpha
          }
          
          # gamma = c(gamma, rep(0, length(q)))
          # alpha = c(alpha, rep(0, length(q)))
          # beta = c(beta, beta_obs)
          # 
          # Order.q = rep(q, each = length(level) + 1)[under_max_alpha]
          Order.q = rep(q, each = length(level))
          
          beta = data.frame(Estimate = beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl) 
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, level, datatype, m_gamma, m_alpha) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(zik)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            aik = (1 - dij/tau) %*% as.matrix(zik)
          }
          positive_id = rowSums(aik)>0
          
          
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = gamma_x/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          # gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          ## Add MLE
          gamma_a_MLE = rowSums(aik)[positive_id]
          gamma_v_MLE = gamma_x/gamma_a_MLE
          
          ai_vi_gamma_MLE = list(ai = data.frame(gamma_a_MLE), vi = data.frame(gamma_v_MLE))
          
          gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n, ai_vi_gamma_MLE) %>% as.vector
          ##
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          # alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
          
          ## Add MLE
          alpha_a_MLE = as.vector(aik)
          
          alpha_v_MLE = rep(gamma_v_MLE, N)
          alpha_v_MLE = alpha_v_MLE[alpha_a_MLE>0]
          alpha_a_MLE = alpha_a_MLE[alpha_a_MLE>0]
          
          ai_vi_alpha_MLE = list(ai = data.frame(alpha_a_MLE), vi = data.frame(alpha_v_MLE))
          
          alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n, ai_vi_alpha_MLE)/N) %>% as.vector
          ##
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(gamma_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          }
          
          gamma_a[gamma_a > n] = n
          gamma_a_MLE = gamma_a
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = gamma_Y/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          # gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          ## Add MLE
          gamma_v_MLE = gamma_Y/gamma_a_MLE
          
          ai_vi_gamma_MLE = list(ai = data.frame(gamma_a_MLE), vi = data.frame(gamma_v_MLE))
          
          gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n, ai_vi_gamma_MLE) %>% as.vector
          ##
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            alpha_a = (1 - dij/1) %*% as.matrix(alpha_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          }
          
          alpha_a[alpha_a > n] = n
          alpha_a_MLE = alpha_a
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          # alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
          
          ## Add MLE
          alpha_a_MLE = as.vector(alpha_a_MLE)
          
          alpha_v_MLE = rep(gamma_v_MLE, N)
          alpha_v_MLE = alpha_v_MLE[alpha_a_MLE > 0]
          alpha_a_MLE = alpha_a_MLE[alpha_a_MLE > 0]
          
          ai_vi_alpha_MLE = list(ai = data.frame(alpha_a_MLE), vi = data.frame(alpha_v_MLE))
          
          alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n, ai_vi_alpha_MLE)/N) %>% as.vector
          ##
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', m_gamma = m_gamma, m_alpha = m_alpha)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC>=ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each=length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC>=ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each=length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
          ## Observed Beta ##
          # output_obs = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
          # gamma_obs = output_obs$gamma
          # alpha_obs = output_obs$alpha
          # obs_beta = gamma_obs/alpha_obs
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype='incidence_raw', m_gamma=m_gamma, m_alpha=m_alpha)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC>=ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
          ## Observed Beta ##
          # output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_gamma_freq[1])
          # gamma_obs = output_obs$gamma
          # alpha_obs = output_obs$alpha
          # obs_beta = gamma_obs/alpha_obs
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        cut = seq(0.00000001, 1, length.out = FDcut_number)
        # cut = seq(0, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', m_gamma = m_gamma, m_alpha = m_alpha)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = data.frame(SC = level, beta) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          ## Observed Beta ##
          # obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
          #   
          #   FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
          #   
          # })
          # 
          # obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
          # 
          # if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)
          # 
          # obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', m_gamma = m_gamma, m_alpha = m_alpha)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = data.frame(SC = level, beta) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          ## Observed Beta ##
          # obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
          #   
          #   FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_alpha_freq[1])
          #   
          # })
          # 
          # obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
          # 
          # if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)
          # 
          # obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
          
        }
        
      }
      
      gamma = gamma[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha = alpha[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      beta = beta[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% 
      #   rbind(., data.frame(Estimate = obs_beta, Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess, workers=7)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype='abundance', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha
              
              beta=gamma/alpha
              
              ## Observed Beta ##
              # output_obs = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
              # gamma_obs = output_obs$gamma
              # alpha_obs = output_obs$alpha
              # beta_obs = gamma_obs/alpha_obs
              
            }
            
            if (FDtype == 'AUC'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', m_gamma = m_gamma, m_alpha = m_alpha)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
              ## Observed Beta ##
              # gamma_alpha_over_tau = lapply(cut, function(tau) {
              #   
              #   FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
              #   
              # })
              # 
              # gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              # 
              # alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              # 
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta_obs = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            # p_bt = p_bt[rowSums(p_bt)>0,]
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha
              
              beta = gamma/alpha
              
              ## Observed Beta ##
              # output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
              # gamma_obs = output_obs$gamma
              # alpha_obs = output_obs$alpha
              # beta_obs = gamma_obs/alpha_obs
              
            }
            
            if (FDtype == 'AUC'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', m_gamma = m_gamma, m_alpha = m_alpha)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
              ## Observed Beta ##
              # gamma_alpha_over_tau = lapply(cut, function(tau) {
              #   
              #   FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
              #   
              # })
              # 
              # gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              # 
              # alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              # 
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta_obs = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (FDtype != "AUC") {
            beta = gamma / alpha
          }
          
          # gamma = c(gamma, rep(0, length(q)))
          # alpha = c(alpha, rep(0, length(q)))
          # beta = c(beta, beta_obs)
          # 
          # Order.q = rep(q, each=length(level) + 1)[under_max_alpha]
          Order.q = rep(q, each=length(level))
          
          beta = data.frame(Estimate = beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "PD" & PDtype == "AUC") index = "PD_AUC"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, any_of("PD"), Estimate, Method, s.e., LCL, UCL, Diversity)
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, any_of("PD"), Estimate, Method, s.e., LCL, UCL, Diversity)
    
    if (diversity == "PD" & PDtype == "meanPD") {
      gamma = gamma %>%
        rename(Gamma_meanPD = Estimate, Gamma_PD = PD)
      
      alpha = alpha %>%
        rename(Alpha_meanPD = Estimate, Alpha_PD = PD)
      
    } else {
      gamma = gamma %>% rename(Gamma = Estimate)
      
      alpha = alpha %>% rename(Alpha = Estimate)
    }
    
    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, Estimate,
                                      Method, s.e., LCL, UCL, Diversity) %>% rename("Beta" = "Estimate")
    
    C = C %>% mutate(        s.e. = se$C,
                             LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, Estimate,
                                      Method, s.e., LCL, UCL, Diversity) %>% rename("Dissimilarity" = "Estimate")
    
    
    U = U %>% mutate(        s.e. = se$U,
                             LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, Estimate,
                                      Method, s.e., LCL, UCL, Diversity) %>% rename("Dissimilarity" = "Estimate")
    
    V = V %>% mutate(        s.e. = se$V,
                             LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, Estimate,
                                      Method, s.e., LCL, UCL, Diversity) %>% rename("Dissimilarity" = "Estimate")
    
    S = S %>% mutate(        s.e. = se$S,
                             LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% select(Dataset, Order.q, SC, Size, Estimate,
                                      Method, s.e., LCL, UCL, Diversity) %>% rename("Dissimilarity" = "Estimate")
    
    
    if (trunc) {
      
      gamma = gamma %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      alpha = alpha %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      beta  = beta  %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      C    =  C    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      U    =  U    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      V    =  V    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      S    =  S    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
    }
    
    if (diversity == "PD" & PDtype == "meanPD") {
      
      gamma = gamma %>% mutate(Reftime = reft)
      
      alpha = alpha %>% mutate(Reftime = reft)
      
      beta  = beta  %>% mutate(Reftime = reft)
      
      C     =  C    %>% mutate(Reftime = reft)
      
      U     =  U    %>% mutate(Reftime = reft)
      
      V     =  V    %>% mutate(Reftime = reft)
      
      S     =  S    %>% mutate(Reftime = reft)
    }
    
    if (diversity == "FD" & FDtype == "tau_value") {
      
      gamma = gamma %>% mutate(Tau = FDtau)
      
      alpha = alpha %>% mutate(Tau = FDtau)
      
      beta  = beta  %>% mutate(Tau = FDtau)
      
      C     =  C    %>% mutate(Tau = FDtau)
      
      U     =  U    %>% mutate(Tau = FDtau)
      
      V     =  V    %>% mutate(Tau = FDtau)
      
      S     =  S    %>% mutate(Tau = FDtau)
    }
    
    if (datatype == "abundance") {
      
      alpha$Method[alpha$SC == ref_gamma_max] = 
        beta$Method[beta$SC == ref_gamma_max] = 
        C$Method[C$SC == ref_gamma_max] = 
        U$Method[U$SC == ref_gamma_max] = 
        V$Method[V$SC == ref_gamma_max] = 
        S$Method[S$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha_max] = 
        alpha$Method[alpha$SC == ref_alpha_max] = 
        beta$Method[beta$SC == ref_alpha_max] = 
        C$Method[C$SC == ref_alpha_max] = 
        U$Method[U$SC == ref_alpha_max] = 
        V$Method[V$SC == ref_alpha_max] = 
        S$Method[S$SC == ref_alpha_max] = "Extrap_SC(2n, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"
      
      
      alpha$Method[alpha$SC == ref_gamma] = 
        beta$Method[beta$SC == ref_gamma] = 
        C$Method[C$SC == ref_gamma] = 
        U$Method[U$SC == ref_gamma] = 
        V$Method[V$SC == ref_gamma] = 
        S$Method[S$SC == ref_gamma] = "Observed_SC(n, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha] = 
        alpha$Method[alpha$SC == ref_alpha] = 
        beta$Method[beta$SC == ref_alpha] = 
        C$Method[C$SC == ref_alpha] = 
        U$Method[U$SC == ref_alpha] = 
        V$Method[V$SC == ref_alpha] = 
        S$Method[S$SC == ref_alpha] = "Observed_SC(n, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma] = "Observed_SC(n, gamma)"
      
      
    }
    
    if (datatype == "incidence_raw") {
      
      alpha$Method[alpha$SC == ref_gamma_max] = 
        beta$Method[beta$SC == ref_gamma_max] = 
        C$Method[C$SC == ref_gamma_max] = 
        U$Method[U$SC == ref_gamma_max] = 
        V$Method[V$SC == ref_gamma_max] = 
        S$Method[S$SC == ref_gamma_max] = "Extrap_SC(2T, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha_max] = 
        alpha$Method[alpha$SC == ref_alpha_max] = 
        beta$Method[beta$SC == ref_alpha_max] = 
        C$Method[C$SC == ref_alpha_max] = 
        U$Method[U$SC == ref_alpha_max] = 
        V$Method[V$SC == ref_alpha_max] = 
        S$Method[S$SC == ref_alpha_max] = "Extrap_SC(2T, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma_max] = "Extrap_SC(2T, gamma)"
      
      
      alpha$Method[alpha$SC == ref_gamma] = 
        beta$Method[beta$SC == ref_gamma] = 
        C$Method[C$SC == ref_gamma] = 
        U$Method[U$SC == ref_gamma] = 
        V$Method[V$SC == ref_gamma] = 
        S$Method[S$SC == ref_gamma] = "Observed_SC(T, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha] = 
        alpha$Method[alpha$SC == ref_alpha] = 
        beta$Method[beta$SC == ref_alpha] = 
        C$Method[C$SC == ref_alpha] = 
        U$Method[U$SC == ref_alpha] = 
        V$Method[V$SC == ref_alpha] = 
        S$Method[S$SC == ref_alpha] = "Observed_SC(T, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma] = "Observed_SC(T, gamma)"
      
      
      gamma = gamma %>% rename("mT" = "Size")
      alpha = alpha %>% rename("mT" = "Size")
      beta  = beta  %>% rename("mT" = "Size")
      C     = C     %>% rename("mT" = "Size")
      U     = U     %>% rename("mT" = "Size")
      V     = V     %>% rename("mT" = "Size")
      S     = S     %>% rename("mT" = "Size")
    }
    
    if (by_pair) {
      
      if (datatype == "abundance") ab_name = colnames(data) else ab_name = names(data)
      
      gamma = gamma %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      alpha = alpha %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      beta  = beta  %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      C    =  C    %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      U    =  U    %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      V    =  V    %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      S    =  S    %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
    }
    
    beta = beta %>% select(-any_of(c("Size", "mT", "PD")))
    C    = C    %>% select(-any_of(c("Size", "mT", "PD")))
    U    = U    %>% select(-any_of(c("Size", "mT", "PD")))
    V    = V    %>% select(-any_of(c("Size", "mT", "PD")))
    S    = S    %>% select(-any_of(c("Size", "mT", "PD")))
    
    
    list(gamma = gamma, alpha = alpha, beta = beta, `1-C` = C, `1-U` = U, `1-V` = V, `1-S` = S)
    
  }
  
  for_each_dataset.size = function(data, dataset_name, N, level, by_pair = FALSE) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all datasets must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
        alpha = estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
        alpha = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
      }
      
      se = cbind(gamma$s.e., alpha$s.e. / N)
      colnames(se) = c("gamma", "alpha")
      se = as.data.frame(se)
      # se[is.na(se)] = 0
      
      gamma = (cbind(Size = rep(level, each=length(q)), gamma[,-c(1,3,8,9)]) %>% 
                 mutate(Method = ifelse(Size>=ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,3,4,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      alpha = (cbind(Size = rep(level, each = length(q)), alpha[,-c(1,3,8,9)]) %>% 
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,3,4,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Estimate = alpha$Estimate / N
      
    }
    
    if (diversity == 'PD') {
      
      PD_by_time = function(data, tree, reft_i, datatype, m_gamma, m_alpha) {
        
        if (datatype == "abundance") {
          
          data_gamma_i = rowSums(data)
          data_gamma_i = data_gamma_i[data_gamma_i > 0]
          
          # Gamma
          aL_gamma = iNEXT.3D:::phyBranchAL_Abu(phylo = tree, data = data_gamma_i, rootExtend = TRUE, refT = reft_i)
          
          aL_gamma$treeNabu$branch.length = aL_gamma$BLbyT[, 1]
          
          aL_table_gamma = aL_gamma$treeNabu %>%
            dplyr::select(branch.abun, branch.length, tgroup)
          
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length),
                                       m = m_gamma, q = q, nt = n, reft = reft_i, cal = "PD") %>% t() %>% as.vector()
          
          
          # Alpha
          aL_table_alpha = NULL
          
          for (k in seq_len(N)) {
            
            x = data[data[, k] > 0, k]
            
            names(x) = rownames(data)[data[, k] > 0]
            
            aL_alpha = iNEXT.3D:::phyBranchAL_Abu(phylo = tree, data = x, rootExtend = TRUE, refT = reft_i)
            
            aL_alpha$treeNabu$branch.length = aL_alpha$BLbyT[, 1]
            
            aL_table_k = aL_alpha$treeNabu %>%
              dplyr::select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table_k)
          }
          
          alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length),
                                        m = m_alpha, q = q, nt = n, reft = reft_i, cal = "PD") / N
          ) %>% t() %>% as.vector()
        }
        
        
        if (datatype == "incidence_raw") {
          
          gamma_raw = Reduce("+", data)
          gamma_raw[gamma_raw > 1] = 1
          
          # Gamma
          aL_gamma = iNEXT.3D:::phyBranchAL_Inc(phylo = tree, data = as.matrix(gamma_raw), datatype = "incidence_raw",
                                                rootExtend = TRUE, refT = reft_i)
          
          aL_gamma$treeNabu$branch.length = aL_gamma$BLbyT[, 1]
          
          aL_table_gamma = aL_gamma$treeNabu %>%
            dplyr::select(branch.abun, branch.length, tgroup)
          
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun,Lis = as.matrix(aL_table_gamma$branch.length),
                                       m = m_gamma, q = q, nt = n, reft = reft_i, cal = "PD") %>% t() %>% as.vector()
          
          
          # Alpha
          aL_table_alpha = NULL
          
          for (k in seq_len(N)) {
            
            x = data[[k]]
            
            aL_alpha = iNEXT.3D:::phyBranchAL_Inc(phylo = tree, data = x, datatype = "incidence_raw", rootExtend = TRUE, refT = reft_i)
            
            aL_alpha$treeNabu$branch.length = aL_alpha$BLbyT[, 1]
            
            aL_table_k = aL_alpha$treeNabu %>% dplyr::select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table_k)
          }
          
          alpha = (
            iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length),
                                 m = m_alpha, q = q, nt = n, reft = reft_i, cal = "PD") / N) %>% t() %>% as.vector()
        }
        
        gamma = gamma / reft_i
        alpha = alpha / reft_i
        
        return(data.frame(gamma = gamma, alpha = alpha))
      }
      
      PD_AUC_by_data = function(data, tree, datatype, m_gamma, m_alpha) {
        
        gamma_alpha_over_time = lapply(
          PDcut,
          function(reft_i) {
            PD_by_time(
              data = data,
              tree = tree,
              reft_i = reft_i,
              datatype = datatype,
              m_gamma = m_gamma,
              m_alpha = m_alpha
            )
          }
        )
        
        gamma_over_time = sapply(
          gamma_alpha_over_time,
          function(x) x$gamma
        )
        
        alpha_over_time = sapply(
          gamma_alpha_over_time,
          function(x) x$alpha
        )
        
        if (is.vector(gamma_over_time)) {
          gamma_over_time = matrix(gamma_over_time, nrow = 1)
        }
        
        if (is.vector(alpha_over_time)) {
          alpha_over_time = matrix(alpha_over_time, nrow = 1)
        }
        
        # Beta: first gamma / alpha at each reference time
        beta_over_time = gamma_over_time / alpha_over_time
        
        # Left and right endpoints for trapezoidal integration
        gamma_left  = gamma_over_time[, -ncol(gamma_over_time), drop = FALSE]
        gamma_right = gamma_over_time[, -1, drop = FALSE]
        
        alpha_left  = alpha_over_time[, -ncol(alpha_over_time), drop = FALSE]
        alpha_right = alpha_over_time[, -1, drop = FALSE]
        
        beta_left   = beta_over_time[, -ncol(beta_over_time), drop = FALSE]
        beta_right  = beta_over_time[, -1, drop = FALSE]
        
        width_matrix = matrix(
          PDwidth,
          nrow = nrow(gamma_over_time),
          ncol = length(PDwidth),
          byrow = TRUE
        )
        
        gamma = rowSums(
          ((gamma_left + gamma_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        alpha = rowSums(
          ((alpha_left + alpha_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        beta = rowSums(
          ((beta_left + beta_right) / 2) * width_matrix
        ) / sum(PDwidth)
        
        return(
          list(
            gamma = gamma,
            alpha = alpha,
            beta = beta
          )
        )
      }
      
      
      if (PDtype == "AUC") {
        
        PD_auc <- PD_AUC_by_data(data = data, tree = PDtree, datatype = datatype, m_gamma = level, m_alpha = level)
        
        
        if (datatype == "abundance") {
          

          gamma = data.frame(Order.q = rep(q, each = length(level)),
                             Estimate = PD_auc$gamma,
                             Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance",level),length(q)),
                             Size = rep(level, times = length(q)))
          
          alpha = data.frame(Order.q = rep(q, each = length(level)),
                             Estimate = PD_auc$alpha,
                             Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha,"abundance",level),length(q)),
                             Size = rep(level,times = length(q)))
          
          
        }
        
        
        if (datatype == "incidence_raw") {
          
          gamma = data.frame(Order.q = rep(q, each = length(level)),
            Estimate = PD_auc$gamma,
            Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq",level),length(q)),
            Size = rep(level, times = length(q)))
          
          alpha = data.frame(Order.q = rep(q, each = length(level)),
            Estimate = PD_auc$alpha,
            Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq,"incidence_freq",level),length(q)),
            Size = rep(level,times = length(q)))
        }
        
      } else {
        

        if (datatype == 'abundance') {
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>% 
            set_colnames(q) %>% gather(Order.q, Estimate) %>% 
            mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)), Size = rep(level, length(q)))
          
          
          aL_table_alpha = c()
          
          for (i in 1:N){
            
            x = data[data[,i]>0,i]
            names(x) = rownames(data)[data[,i]>0]
            
            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table)
            
          }
          
          
          qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD")
          qPDm = qPDm/N
          alpha = qPDm %>% t %>% as.data.frame %>% 
            set_colnames(q) %>% gather(Order.q, Estimate) %>% 
            mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)), Size = rep(level, length(q)))
          
        }
        
        if (datatype == 'incidence_raw') {
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD") %>% t %>% as.data.frame %>% 
            set_colnames(q) %>% gather(Order.q, Estimate) %>% 
            mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)), Size = rep(level, length(q)))
          
          aL_table_alpha = c()
          
          for (i in 1:N){
            
            x = data[[i]]
            
            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table)
            
          }
          
          alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, reft = reft, cal = "PD")/N) %>% t %>% as.data.frame %>% 
            set_colnames(q) %>% gather(Order.q, Estimate) %>% 
            mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)), Size = rep(level, length(q)))
          
          
        }
        
        
      }
      
      

      
      gamma = (gamma %>% 
                 mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      gamma$Order.q = as.numeric(gamma$Order.q)
      
      
      alpha = (alpha %>% 
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Order.q = as.numeric(alpha$Order.q)
      
      if (PDtype == 'meanPD') {
        gamma$PD = gamma$Estimate
        alpha$PD = alpha$Estimate
        
        gamma$Estimate = gamma$PD / reft
        alpha$Estimate = alpha$PD / reft
      }
      
      # beta = alpha
      # beta$Estimate = gamma$Estimate/alpha$Estimate
      # 
      # C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == "abundance") {
            tree_bt = PDtree
            bootstrap_population = bootstrap_population_multiple_assemblage(data, 
                                                                                           data_gamma, "abundance")
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
            if (nrow(p_bt) > nrow(data) & sum(unseen_p) > 
                0) {
              unseen = unseen_p[which(rowSums(unseen_p) > 
                                        0), ]
              unseen = matrix(unseen, ncol = ncol(unseen_p))
              p_bt = rbind(p_bt[(1:nrow(data)), ], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_", 
                                                                      i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              bootstrap_sample = sapply(1:ncol(data), 
                                        function(k) rmultinom(n = 1, size = sum(data[, 
                                                                                     k]), prob = p_bt[, k]))
              x_bt = bootstrap_sample
              rownames(x_bt) = rownames(p_bt)
              if (sum(x_bt[-(1:nrow(data)), ]) > 0) {
                g0_hat = apply(data, 2, function(x) {
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, 
                                                  data = x, rootExtend = T, refT = reft)
                  aL$treeNabu$branch.length = aL$BLbyT[, 
                                                       1]
                  aL = aL$treeNabu %>% select(branch.abun, 
                                              branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 
                                          1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 
                                          2] %>% sum
                  g0_hat = ifelse(g2 > ((g1 * f2)/(2 * 
                                                     f1)), ((n - 1)/n) * (g1^2/(2 * g2)), 
                                  ((n - 1)/n) * (g1 * (f1 - 1)/(2 * 
                                                                  (f2 + 1))))
                  g0_hat
                })
                te = (x_bt[1:nrow(data), ] * (data == 
                                                0)) > 0
                used_length = sapply(1:ncol(data), function(i) {
                  if (sum(te[, i]) == 0) 
                    return(0)
                  else {
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, 
                                               data = x_bt[1:nrow(data), i], 
                                               rootExtend = T, refT = reft)$treeNabu %>% 
                      subset(label %in% names(which(te[, 
                                                       i] == TRUE))) %>% select(branch.length) %>% 
                      sum
                  }
                })
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                unseen_sample = x_bt[-(1:nrow(data)), 
                ]
                if (is.vector(unseen_sample)) 
                  unseen_sample = matrix(unseen_sample, 
                                         ncol = ncol(x_bt))
                L0_hat = sapply(1:length(g0_hat), function(i) if (sum(unseen_sample[, 
                                                                                    i] > 0) > 0) 
                  (g0_hat[i]/nrow(unseen))
                  else 0)
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), 
                                         ncol(unseen_sample), byrow = T) * 
                                    unseen_sample))/rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 
                               0)] = 0
                for (i in 1:length(L0_hat)) {
                  tip = list(edge = matrix(c(2, 1), 
                                           1, 2), tip.label = unseen_name[i], 
                             edge.length = L0_hat[i], Nnode = 1)
                  class(tip) = "phylo"
                  tree_bt = tree_bt + tip
                }
              }
              else {
                x_bt = x_bt[1:nrow(data), ]
                p_bt = p_bt[1:nrow(data), ]
              }
            }
            else {
              p_bt = p_bt[1:nrow(data), ]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, 
                                                                size = sum(data[, k]), prob = p_bt[, 
                                                                                                   k]))
              rownames(x_bt) = rownames(data)
            }
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 
                                                          0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% 
              as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 
                                                          0]
            
            
            if (PDtype == "AUC") {
              
              PD_auc_bt = PD_AUC_by_data(data = x_bt, tree = tree_bt, datatype = "abundance",
                                         m_gamma = level, m_alpha = level)
              
              gamma = PD_auc_bt$gamma
              alpha = PD_auc_bt$alpha
              
            }else{
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, 
                                              data = bootstrap_data_gamma, rootExtend = T, 
                                              refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[, 1]
              aL_table_gamma = aL$treeNabu %>% select(branch.abun, 
                                                      branch.length, tgroup)
              gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, 
                                                     Lis = as.matrix(aL_table_gamma$branch.length), 
                                                     m = level, q = q, nt = n, reft = reft, 
                                                     cal = "PD") %>% t)
              aL_table_alpha = c()
              for (i in 1:N) {
                x = x_bt[, i]
                names(x) = rownames(p_bt)
                x = x[x_bt[, i] > 0]
                aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, 
                                                data = x, rootExtend = T, refT = reft)
                aL$treeNabu$branch.length = aL$BLbyT[, 
                                                     1]
                aL_table = aL$treeNabu %>% select(branch.abun, 
                                                  branch.length, tgroup)
                aL_table_alpha = rbind(aL_table_alpha, 
                                       aL_table)
              }
              alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, 
                                                      Lis = as.matrix(aL_table_alpha$branch.length), 
                                                      m = level, q = q, nt = n, reft = reft, 
                                                      cal = "PD")/N) %>% t)
              
            }
            
            
          }
          if (datatype == "incidence_raw") {
            tree_bt = PDtree
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, 
                                                                                           data_gamma_freq, "incidence")
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])), ] %>% 
              matrix(ncol = N)
            if (nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p) > 
                0) {
              unseen = unseen_p[which(rowSums(unseen_p) > 
                                        0), ]
              unseen = matrix(unseen, ncol = ncol(unseen_p))
              p_bt = rbind(p_bt[(1:nrow(data[[1]])), 
              ], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_", 
                                                                      i))
              rownames(p_bt) = c(rownames(data[[1]]), 
                                 unseen_name)
              raw = lapply(1:ncol(p_bt), function(j) {
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, 
                                                        size = 1, prob = p_bt[i, j])) %>% 
                  do.call(rbind, .)
              })
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              if (lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])), 
              ]) %>% do.call(sum, .) > 0) {
                R0_hat = sapply(data, function(x) {
                  nT = ncol(x)
                  Q1 = sum(rowSums(x) == 1)
                  Q2 = sum(rowSums(x) == 2)
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, 
                                                  data = x, datatype = "incidence_raw", 
                                                  rootExtend = T, refT = reft)
                  aL$treeNabu$branch.length = aL$BLbyT[, 
                                                       1]
                  aL = aL$treeNabu %>% select(branch.abun, 
                                              branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 
                                          1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 
                                          2] %>% sum
                  R0_hat = ifelse(R2 > ((R1 * Q2)/(2 * 
                                                     Q1)), ((nT - 1)/nT) * (R1^2/(2 * 
                                                                                    R2)), ((nT - 1)/nT) * (R1 * (Q1 - 
                                                                                                                   1)/(2 * (Q2 + 1))))
                  R0_hat
                })
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]), 
                ] * (sapply(data, rowSums) == 0)) > 
                  0
                used_length = sapply(1:N, function(i) {
                  if (sum(te[, i]) == 0) 
                    return(0)
                  else {
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, 
                                               data = raw[[i]][1:nrow(data[[1]]), 
                                               ], datatype = "incidence_raw", 
                                               rootExtend = T, refT = reft)$treeNabu %>% 
                      subset(label %in% names(which(te[, 
                                                       i] == TRUE))) %>% select(branch.length) %>% 
                      sum
                  }
                })
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])), 
                ]
                if (is.vector(unseen_sample)) 
                  unseen_sample = matrix(unseen_sample, 
                                         ncol = N)
                L0_hat = sapply(1:length(R0_hat), function(i) if (sum(unseen_sample[, 
                                                                                    i] > 0) > 0) 
                  (R0_hat[i]/nrow(unseen))
                  else 0)
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), 
                                         ncol(unseen_sample), byrow = T) * 
                                    unseen_sample))/rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 
                               0)] = 0
                for (i in 1:length(L0_hat)) {
                  tip = list(edge = matrix(c(2, 1), 
                                           1, 2), tip.label = unseen_name[i], 
                             edge.length = L0_hat[i], Nnode = 1)
                  class(tip) = "phylo"
                  tree_bt = tree_bt + tip
                }
              }
              else raw = lapply(raw, function(i) i[1:nrow(data[[1]]), 
              ])
            }
            else {
              p_bt = p_bt[1:nrow(data[[1]]), ]
              raw = lapply(1:ncol(p_bt), function(j) {
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, 
                                                        size = 1, prob = p_bt[i, j])) %>% 
                  do.call(rbind, .)
              })
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
            }
            gamma = Reduce("+", raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            bootstrap_data_alpha_freq = sapply(raw, 
                                               rowSums) %>% c(n, .)
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 
                                                                    0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 
                                                                    0]
            
            
            if (PDtype == "AUC") {
              
              PD_auc_bt = PD_AUC_by_data(data = raw, tree = tree_bt, datatype = "incidence_raw",
                                         m_gamma = level, m_alpha = level)
              
              gamma = PD_auc_bt$gamma
              alpha = PD_auc_bt$alpha
              
            }else{
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, 
                                              data = bootstrap_data_gamma_raw, datatype = "incidence_raw", 
                                              rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[, 1]
              aL_table_gamma = aL$treeNabu %>% select(branch.abun, 
                                                      branch.length, tgroup)
              gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, 
                                                     Lis = as.matrix(aL_table_gamma$branch.length), 
                                                     m = level, q = q, nt = n, reft = reft, 
                                                     cal = "PD") %>% t)
              aL_table_alpha = c()
              for (i in 1:N) {
                x = raw[[i]]
                aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, 
                                                data = x, datatype = "incidence_raw", 
                                                rootExtend = T, refT = reft)
                aL$treeNabu$branch.length = aL$BLbyT[, 
                                                     1]
                aL_table = aL$treeNabu %>% select(branch.abun, 
                                                  branch.length, tgroup)
                aL_table_alpha = rbind(aL_table_alpha, 
                                       aL_table)
              }
              alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, 
                                                      Lis = as.matrix(aL_table_alpha$branch.length), 
                                                      m = level, q = q, nt = n, reft = reft, 
                                                      cal = "PD")/N) %>% t)
              
            }
            
          }
          
          if (PDtype == 'meanPD') {
            gamma = gamma / reft
            alpha = alpha / reft
          } 
          
          # beta = gamma/alpha
          # 
          # Order.q = rep(q, each = length(level))
          # 
          # beta = data.frame(Estimate = beta, Order.q)
          # 
          # C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          # 
          # beta = beta$Estimate
          # 
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl) 
        # plan(sequential)
        
      } else {
        
        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)
        
        se = matrix(NA, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, datatype, m_gamma, m_alpha) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(zik)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            aik = (1 - dij/tau) %*% as.matrix(zik)
          }
          positive_id = rowSums(aik)>0
          
          
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = gamma_x/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          # gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          ## Add MLE
          gamma_a_MLE = rowSums(aik)[positive_id]
          gamma_v_MLE = gamma_x/gamma_a_MLE
          
          ai_vi_gamma_MLE = list(ai = data.frame(gamma_a_MLE), vi = data.frame(gamma_v_MLE))
          
          gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n, ai_vi_gamma_MLE) %>% as.vector
          ##
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          # alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
          
          ## Add MLE
          alpha_a_MLE = as.vector(aik)
          
          alpha_v_MLE = rep(gamma_v_MLE, N)
          alpha_v_MLE = alpha_v_MLE[alpha_a_MLE>0]
          alpha_a_MLE = alpha_a_MLE[alpha_a_MLE>0]
          
          ai_vi_alpha_MLE = list(ai = data.frame(alpha_a_MLE), vi = data.frame(alpha_v_MLE))
          
          alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n, ai_vi_alpha_MLE)/N) %>% as.vector
          ##
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(gamma_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          }
          
          gamma_a[gamma_a > n] = n
          gamma_a_MLE = gamma_a
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = gamma_Y/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          # gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          ## Add MLE
          gamma_v_MLE = gamma_Y/gamma_a_MLE
          
          ai_vi_gamma_MLE = list(ai = data.frame(gamma_a_MLE), vi = data.frame(gamma_v_MLE))
          
          gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n, ai_vi_gamma_MLE) %>% as.vector
          ##
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            alpha_a = (1 - dij/1) %*% as.matrix(alpha_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          }
          
          alpha_a[alpha_a > n] = n
          alpha_a_MLE = alpha_a
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          # alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
          
          ## Add MLE
          alpha_a_MLE = as.vector(alpha_a_MLE)
          
          alpha_v_MLE = rep(gamma_v_MLE, N)
          alpha_v_MLE = alpha_v_MLE[alpha_a_MLE > 0]
          alpha_a_MLE = alpha_a_MLE[alpha_a_MLE > 0]
          
          ai_vi_alpha_MLE = list(ai = data.frame(alpha_a_MLE), vi = data.frame(alpha_v_MLE))
          
          alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n, ai_vi_alpha_MLE)/N) %>% as.vector
          ##
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(data, FDdistM, FDtau, datatype = 'abundance', m_gamma = level, m_alpha = level)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, datatype='incidence_raw', m_gamma=level, m_alpha=level)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        cut = seq(0.00000001, 1, length.out = FDcut_number)
        # cut = seq(0, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
          # beta = data.frame(Size = level, beta) %>% 
          #   mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
          #          Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
          # beta = data.frame(Size = level, beta) %>% 
          #   mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
          #          Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
        }
        
      }
      
      gamma = gamma[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      alpha = alpha[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      # beta = beta[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      # 
      # C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess, workers=7)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, datatype='abundance', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha
              
              # beta=gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            # p_bt = p_bt[rowSums(p_bt)>0,]
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha
              
              # beta = gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              if (is.vector(gamma_over_tau)) gamma_over_tau = matrix(gamma_over_tau, nrow = 1)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              if (is.vector(alpha_over_tau)) alpha_over_tau = matrix(alpha_over_tau, nrow = 1)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          # beta = gamma/alpha
          # 
          # Order.q = rep(q, each=length(level))
          # 
          # beta = data.frame(Estimate = beta, Order.q)
          # 
          # C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          # 
          # beta = beta$Estimate
          # 
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)
        
        se = matrix(NA, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "PD" & PDtype == "AUC") index = "PD_AUC"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, Size) %>% select(Dataset, Order.q, SC, Size, any_of("PD"), Estimate, Method, s.e., LCL, UCL, Diversity)
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, Size) %>% select(Dataset, Order.q, SC, Size, any_of("PD"), Estimate, Method, s.e., LCL, UCL, Diversity)
    
    if (diversity == "PD" & PDtype == "meanPD") {
      gamma = gamma %>%
        rename(Gamma_meanPD = Estimate, Gamma_PD = PD)
      
      alpha = alpha %>%
        rename(Alpha_meanPD = Estimate, Alpha_PD = PD)
      
    } else {
      gamma = gamma %>% rename(Gamma = Estimate)
      
      alpha = alpha %>% rename(Alpha = Estimate)
    }
    
    if (datatype == 'incidence_raw') {
      
      gamma = gamma %>% rename("mT" = "Size")
      alpha = alpha %>% rename("mT" = "Size")
    }
    
    if (diversity == "PD" & PDtype == "meanPD") {
      
      gamma = gamma %>% mutate(Reftime = reft)
      
      alpha = alpha %>% mutate(Reftime = reft)
    }
    
    if (diversity == "FD" & FDtype == "tau_value") {
      
      gamma = gamma %>% mutate(Tau = FDtau)
      
      alpha = alpha %>% mutate(Tau = FDtau)
    }
    
    if (by_pair) {
      
      if (datatype == "abundance") ab_name = colnames(data) else ab_name = names(data)
      
      gamma = gamma %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
      alpha = alpha %>% mutate('Pair' = paste(ab_name[1], ab_name[2], sep = " vs. "), .after = "Dataset")
      
    }
    
    
    list(gamma = gamma, alpha = alpha)
    
  }
  
  
  if (by_pair == FALSE) {
    
    if (base == 'coverage') output = lapply(1:length(data_list), function(i) for_each_dataset(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))
    
    if (base == 'size') output = lapply(1:length(data_list), function(i) for_each_dataset.size(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))
    
  }
  
  
  ## Added by KaiHsiang 2024/11/11
  if (by_pair == TRUE) {
    
    if (base == 'coverage') output = lapply(1:length(data_list), function(i) {
      
      if (datatype == "abundance") {
        
        two_idx <- combn(1:ncol(data_list[[i]]), 2)
        
        pair_out = lapply(1:ncol(two_idx), function(j) for_each_dataset(data = data_list[[i]][,two_idx[,j]], dataset_name = dataset_names[i], N = Ns[i], 
                                                                        level = level[[i]][[j]], by_pair = TRUE))
      }
      
      if (datatype == "incidence_raw") {
        
        two_idx <- combn(1:length(data_list[[i]]), 2)
        
        pair_out = lapply(1:ncol(two_idx), function(j) for_each_dataset(data = data_list[[i]][two_idx[,j]], dataset_name = dataset_names[i], N = Ns[i], 
                                                                        level = level[[i]][[j]], by_pair = TRUE))
      }
      
      
      lapply(1:7, function(k) lapply(pair_out, function(x) x[[k]]) %>% do.call(rbind,.)) %>% set_names(names(pair_out[[1]]))
    })
    
    if (base == 'size') output = lapply(1:length(data_list), function(i) {
      
      if (datatype == "abundance") {
        
        two_idx <- combn(1:ncol(data_list[[i]]), 2)
        
        pair_out = lapply(1:ncol(two_idx), function(j) for_each_dataset.size(data = data_list[[i]][,two_idx[,j]], dataset_name = dataset_names[i], N = Ns[i], 
                                                                             level = level[[i]][[j]], by_pair = TRUE))
      }
      
      if (datatype == "incidence_raw") {
        
        two_idx <- combn(1:length(data_list[[i]]), 2)
        
        pair_out = lapply(1:ncol(two_idx), function(j) for_each_dataset.size(data = data_list[[i]][two_idx[,j]], dataset_name = dataset_names[i], N = Ns[i], 
                                                                             level = level[[i]][[j]], by_pair = TRUE))
      }
      
      
      lapply(1:2, function(k) lapply(pair_out, function(x) x[[k]]) %>% do.call(rbind,.)) %>% set_names(names(pair_out[[1]]))
    })
    
  }
  
  names(output) = dataset_names
  
  class(output) <- c("iNEXTbeta3D")
  return(output)
  
  }






#' ggplot2 Extension for iNEXTbeta3D Objects
#'
#' \code{ggiNEXTbeta3D()} is a plotting function based on \code{ggplot2} that generates sample-size-based and coverage-based rarefaction and extrapolation (R/E) sampling curves from \code{iNEXTbeta3D()} output. It supports results generated with both \code{by_pair = FALSE} and \code{by_pair = TRUE}.
#'
#' @param output An object returned by \code{iNEXTbeta3D()}.
#' @param type Argument used only when \code{base = "coverage"}. \code{type = "B"} generates R/E sampling curves for gamma, alpha, and beta diversity, whereas \code{type = "D"} generates R/E sampling
#' curves for four dissimilarity indices. The \code{type} argument is not required when plotting size-based R/E sampling curves for gamma and alpha diversity.
#'
#' @return When the input was generated with \code{by_pair = FALSE}, the function returns a single \code{ggplot} object that can be customized using \code{ggplot2}. 
#' When the input was generated with \code{by_pair = TRUE}, the function returns a named list of \code{ggplot} objects, with one plot for each dataset and different pairs of assemblages distinguished by color.
#'
#' @examples
#' \donttest{
#' ## (Graphic Display) Taxonomic diversity for abundance data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' # See Example 1 in the iNEXT.beta3D vignette for details.
#' 
#' data(Brazil_rainforests)
#' output_TDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'TD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10)
#' 
#' ggiNEXTbeta3D(output_TDc_abun, type = 'B')
#' ggiNEXTbeta3D(output_TDc_abun, type = 'D')
#' }
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' 
#' data(Brazil_rainforests)
#' output_TDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'TD', 
#'                               datatype = 'abundance', base = "size", nboot = 10)
#' 
#' ggiNEXTbeta3D(output_TDs_abun)
#' 
#' \dontrun{
#' ## (Graphic Display) Taxonomic diversity for incidence data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' # See Example 2 in the iNEXT.beta3D vignette for details.
#' 
#' data(Second_growth_forests)
#' output_TDc_inci = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                               datatype = 'incidence_raw', base = "coverage", nboot = 10)
#' 
#' ggiNEXTbeta3D(output_TDc_inci, type = 'B')
#' ggiNEXTbeta3D(output_TDc_inci, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' 
#' data(Second_growth_forests)
#' output_TDs_inci = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', 
#'                               datatype = 'incidence_raw', base = "size", nboot = 10)
#' 
#' ggiNEXTbeta3D(output_TDs_inci)
#'
#'
#'
#' # Coverage-based rarefaction and extrapolation sampling curves for every pair of assemblages
#' # See Example 3 in the iNEXT.beta3D vignette for details.
#'
#' output_TDc_inci_pair = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD',
#'                                    datatype = 'incidence_raw', base = "coverage",
#'                                    nboot = 10, by_pair = TRUE)
#'
#' plot_TDc_inci_pair = ggiNEXTbeta3D(output_TDc_inci_pair, type = 'B')
#' names(plot_TDc_inci_pair)
#' plot_TDc_inci_pair[["CR_2005_2011_2017"]]
#' plot_TDc_inci_pair[["JE_2005_2011_2017"]]
#'
#'
#' 
#' ## (Graphic Display) Phylogenetic diversity for abundance data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' # See Example 4 in the iNEXT.beta3D vignette for details.
#' 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output_PDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10, 
#'                               PDtree = Brazil_tree, PDtype = 'AUC', PDcut_number = 30) 
#' 
#' ggiNEXTbeta3D(output_PDc_abun, type = 'B')
#' ggiNEXTbeta3D(output_PDc_abun, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output_PDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                               datatype = 'abundance', base = "size", nboot = 10, 
#'                               PDtree = Brazil_tree, PDtype = 'AUC', PDcut_number = 30) 
#' 
#' ggiNEXTbeta3D(output_PDs_abun)
#' 
#' 
#' ## (Graphic Display) Functional diversity for abundance data when all threshold levels 
#' ## from 0 to 1 are considered
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' # See Example 5 in the iNEXT.beta3D vignette for details.
#' 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output_FDc_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                               datatype = 'abundance', base = "coverage", nboot = 10, 
#'                               FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggiNEXTbeta3D(output_FDc_abun, type = 'B')
#' ggiNEXTbeta3D(output_FDc_abun, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output_FDs_abun = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                               datatype = 'abundance', base = "size", nboot = 10, 
#'                               FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggiNEXTbeta3D(output_FDs_abun)
#' }
#' 
#' 
#' @export
ggiNEXTbeta3D = function(output, type = 'B'){
  
  transp = 0.4
  
  add_alpha_size <- function(y, tmp) {
    
    alpha_size <- y[["alpha"]]
    
    if ("mT" %in% colnames(alpha_size)) {
      alpha_size <- alpha_size %>% rename(Size = mT)
    }
    
    join_vars <- c("Dataset", "Order.q", "SC")
    
    if ("Pair" %in% colnames(tmp) && "Pair" %in% colnames(alpha_size)) {
      join_vars <- c("Dataset", "Pair", "Order.q", "SC")
    }
    
    alpha_size <- alpha_size %>%
      select(all_of(join_vars), Size) %>%
      distinct()
    
    tmp %>%
      left_join(alpha_size, by = join_vars)
  }
  
  get_double_extrapolation <- function(df){
    group_vars <- c("Dataset","Order.q","div_type")
    if("Pair" %in% colnames(df)) group_vars <- c("Dataset","Pair","Order.q","div_type")
    obs <- df %>% filter(Method=="Observed") %>% select(all_of(group_vars), obs_size=Size) %>% distinct()
    df %>% filter(Method=="Extrapolation") %>% left_join(obs, by=group_vars) %>% filter(!is.na(obs_size), round(Size)==round(2*obs_size))
  }
  
  make_single_plot <- function(df, xvar, xlab, ylab){
    lty <- c(Rarefaction="solid", Extrapolation="dashed")
    df$Method <- factor(df$Method, levels=c("Rarefaction","Extrapolation","Observed"))

    if("Pair" %in% colnames(df)) df$PlotGroup <- df$Pair else df$PlotGroup <- df$Dataset
    n_group <- length(unique(df$PlotGroup))
    cbPalette <- c("#D55E00", "#0072B2", "#330066", "#C4A484", "#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73")
    if (n_group > length(cbPalette)) {
      cbPalette <- c(
        cbPalette,
        ggplotColors(n_group - length(cbPalette))
      )
    } else {
      cbPalette <- cbPalette[seq_len(n_group)]
    }
    
    double_extrapolation <- get_double_extrapolation(df)
    
    if(xvar=="SC"){
      fig <- ggplot(df, aes(x=SC, y=Estimate, col=PlotGroup, group=interaction(PlotGroup,div_type,Order.q,Method)))
      } else {
        fig <- ggplot(df, aes(x=Size, y=Estimate, col=PlotGroup, group=interaction(PlotGroup,div_type,Order.q,Method)))
      }
    fig <- fig +
      geom_line(data=subset(df, Method!="Observed"), aes(linetype=Method), size=1.1) +
      scale_linetype_manual(values=lty) +
      geom_point(data=subset(df, Method=="Observed" & div_type=="Gamma"), shape=19, size=3) +
      geom_point(data=subset(df, Method=="Observed" & div_type!="Gamma"), shape=1, size=3, stroke=1.5) +
      geom_point(data=subset(double_extrapolation, div_type=="Gamma"), shape=17, size=3) +
      geom_point(data=subset(double_extrapolation, div_type!="Gamma"), shape=2, size=3, stroke=1.5) +
      scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette) +
      facet_grid(div_type ~ Order.q, scales="free") + theme_bw() +
      theme(legend.position="bottom", legend.title=element_blank(), strip.text=element_text(size=15,face="bold"),
            axis.title=element_text(hjust=.5,size=15,face="bold"), axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12), legend.box="vertical", legend.text=element_text(size=13)) +
      labs(x=xlab, y=ylab) + guides(linetype=guide_legend(keywidth=2.5))
    
    if(sum(!is.na(df$LCL))>0) fig <- fig + geom_ribbon(aes(ymin=LCL,ymax=UCL, fill=PlotGroup,
                                                           group=interaction(PlotGroup,div_type,Order.q),col=NULL), alpha=transp)
    
    return(fig)
  }

  
  if ((length(output[[1]]) == 7 & type == "B") | (length(output[[1]]) == 2)) {
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD_AUC') { ylab = "Phylogenetic diversity (AUC)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional diversity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
  }
  
  if (length(output[[1]]) == 7 & type == "D") {
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD_AUC') { ylab = "Phylogenetic dissimilarity (AUC)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional dissimilarity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }
  }
  
  if (length(output[[1]]) == 7) {
    if (type == 'B'){
      
      gamma = lapply(output, function(y) {
        
        tmp = y[["gamma"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = tmp[tmp == 'Observed_SC(T, gamma)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = tmp[tmp == 'Extrap_SC(2T, gamma)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>%  mutate(Estimate = if ("Gamma_meanPD" %in% colnames(.)) Gamma_meanPD else Gamma) %>% select(-any_of(c("Gamma", "Gamma_meanPD", "Gamma_PD"))) %>% mutate(div_type = "Gamma") %>% as_tibble()
      
      alpha = lapply(output, function(y) {
        
        tmp = y[["alpha"]]
        if ('mT' %in% colnames(tmp)) tmp = tmp %>% rename("Size" = "mT")
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% mutate(Estimate = if ("Alpha_meanPD" %in% colnames(.)) Alpha_meanPD else Alpha) %>% select(-any_of(c("Alpha", "Alpha_meanPD", "Alpha_PD"))) %>% mutate(div_type = "Alpha") %>% as_tibble()
      
      beta = lapply(output, function(y) {
        
        tmp = y[["beta"]]
      
        tmp = add_alpha_size(y, tmp)
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      })  %>% do.call(rbind,.) %>% rename("Estimate" = "Beta") %>% mutate(div_type = "Beta")  %>% as_tibble()
      
      
      # # Dropping out the points extrapolated over double reference size
      # gamma1 = data.frame() ; alpha1 = data.frame() ; beta1 = data.frame()
      # 
      # for(i in 1:length(unique(gamma$Dataset))){
      #   
      #   Gamma <- gamma %>% filter(Dataset==unique(gamma$Dataset)[i]) ; ref_size = unique(Gamma[Gamma$Method=="Observed",]$Size)
      #   Gamma = Gamma %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   Alpha <- alpha %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Alpha = Alpha %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   Beta <- beta %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Beta = Beta %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   gamma1 = rbind(gamma1,Gamma) ; alpha1 = rbind(alpha1,Alpha) ; beta1 = rbind(beta1,Beta)
      #   
      # }
      # 
      # gamma = gamma1 ; alpha = alpha1 ; beta= beta1
      
      df = rbind(gamma, alpha, beta)
      for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
      
      id_obs = which(df$Method == 'Observed')
      
      if (length(id_obs) > 0) {
        
        for (i in 1:length(id_obs)) {
          
          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'
          
          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'
          
          df = rbind(df, new, newe)
          
        }
      }
      
      
    }
    
    if (type == 'D'){
      
      C = lapply(output, function(y) {
        
        tmp = y[["1-C"]]
        tmp = add_alpha_size(y, tmp)
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-CqN") %>% as_tibble()
      
      U = lapply(output, function(y) {
        
        tmp = y[["1-U"]]
        tmp = add_alpha_size(y, tmp)
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-UqN") %>% as_tibble()
      
      V = lapply(output, function(y) {
        
        tmp = y[["1-V"]]
        tmp = add_alpha_size(y, tmp)
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-VqN") %>% as_tibble()
      
      S = lapply(output, function(y) {
        
        tmp = y[["1-S"]]
        tmp = add_alpha_size(y, tmp)
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-SqN") %>% as_tibble()
      
      # C = C %>% filter(Method != 'Observed')
      # U = U %>% filter(Method != 'Observed')
      # V = V %>% filter(Method != 'Observed')
      # S = S %>% filter(Method != 'Observed')
     
      
      # # Dropping out the points extrapolated over double reference size
      # c1 = data.frame() ; u1 = data.frame() ; v1 = data.frame() ; s1 = data.frame()
      # 
      # for(i in 1:length(unique(C$Dataset))){
      #   
      #   CC <- C %>% filter(Dataset==unique(C$Dataset)[i]) ; ref_size = unique(CC[CC$Method=="Observed",]$Size)
      #   CC = CC %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   UU <- U %>% filter(Dataset==unique(C$Dataset)[i]) ; UU = UU %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   VV <- V %>% filter(Dataset==unique(C$Dataset)[i]) ; VV = VV %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   SS <- S %>% filter(Dataset==unique(C$Dataset)[i]) ; SS = SS %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   c1 = rbind(c1,CC) ; u1 = rbind(u1,UU) ; v1 = rbind(v1,VV) ; s1 = rbind(s1,SS)
      #   
      # }
      # 
      # C = c1 ; U = u1 ; V = v1 ; S = s1
      
      
      df = rbind(C, U, V, S)
      for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
      
      id_obs = which(df$Method == 'Observed')
      
      if (length(id_obs) > 0) {
        
        for (i in 1:length(id_obs)) {
          
          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'
          
          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'
          
          df = rbind(df, new, newe)
          
        }
      }

      
    }
    
    df$Method <- factor(df$Method, levels=c("Rarefaction","Extrapolation","Observed"))
    
    if("Pair" %in% colnames(df)){
      dataset_names <- unique(df$Dataset)
      fig <- lapply(dataset_names, function(dd){
        df_dataset <- df %>% filter(Dataset==dd)
        make_single_plot(df_dataset, "SC", "Sample coverage", ylab) + ggtitle(dd)
      })
      names(fig) <- dataset_names
    } else {
      fig <- make_single_plot(df, "SC", "Sample coverage", ylab)
    }
    
  } else if (length(output[[1]]) == 2) {
    
    gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(Estimate = if ("Gamma_meanPD" %in% colnames(.)) Gamma_meanPD else Gamma) %>%
      select(-any_of(c("Gamma", "Gamma_meanPD", "Gamma_PD"))) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(Estimate = if ("Alpha_meanPD" %in% colnames(.)) Alpha_meanPD else Alpha) %>%
      select(-any_of(c("Alpha", "Alpha_meanPD", "Alpha_PD"))) %>% mutate(div_type = "Alpha") %>% as_tibble()
    
    if ('mT' %in% colnames(gamma)) {
      
      xlab = 'Number of sampling units'
      colnames(gamma)[colnames(gamma) == 'mT'] = 'Size'
      colnames(alpha)[colnames(alpha) == 'mT'] = 'Size'
      
    } else xlab = 'Number of individuals'
    
    df = rbind(gamma, alpha)
    for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha"))
    
    id_obs = which(df$Method == 'Observed')
    
    if (sum(id_obs) > 0) {
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$Size = new$Size - 0.0001
        new$Method = 'Rarefaction'
        
        newe = df[id_obs[i],]
        newe$Size = newe$Size + 0.0001
        newe$Method = 'Extrapolation'
        
        df = rbind(df, new, newe)
        
      }
    }
    
    df$Method <- factor(df$Method, levels=c("Rarefaction","Extrapolation","Observed"))
    
    if("Pair" %in% colnames(df)){
      dataset_names <- unique(df$Dataset)
      fig <- lapply(dataset_names, function(dd){
        df_dataset <- df %>% filter(Dataset==dd)
        make_single_plot(df_dataset, "Size", xlab, ylab) + ggtitle(dd)
      })
      names(fig) <- dataset_names
    } else {
      fig <- make_single_plot(df, "Size", xlab, ylab)
    }
    
  }
  
  return(fig)
}



coverage_to_size = function(x, C, datatype = 'abundance'){
  
  if (datatype == 'abundance'){
    
    n <- sum(x)
    refC <- iNEXT.3D:::Coverage(x, 'abundance', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'abundance', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  } else {
    
    m <- NULL
    n <- max(x)
    refC <- iNEXT.3D:::Coverage(x, 'incidence_freq', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'incidence_freq', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  }
  
  return(mm)
}

bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
  
  if (datatype == 'abundance'){
    
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
    
    output = apply(data, 2, function(x){
      
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      
      if(length(p_i_hat) != length(x)){
        
        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)
        
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
        
        p_i_hat
        
      } else {
        
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
        
      }
    })
    
  }
  
  if (datatype == 'incidence'){
    
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
    
    output = apply(data, 2, function(x){
      
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      
      if(length(pi_i_hat) != (length(x) - 1)){
        
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = unique(pi_i_hat_unobs)
        
        pi_i_hat
        
      } else {
        
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
        
      }
    })
    
  }
  
  return(output)
  
}

Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){
  
  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }
  
  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  
  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat == 0) {
    d = dij
  } else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(d)
  
}

FD.m.est_0 = function (ai_vi, m, q, nT) {
  EFD = function(m, qs, obs, asy, beta, av) {
    m = m - nT
    out <- sapply(1:length(qs), function(i) {
      if (qs[i] != 2) {
        obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
      }
      else if (qs[i] == 2 & beta[i] != 0 ) {
        V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[, 1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] * (av[, 1] - 1)/(nT * (nT - 1)))))
      }
      else if( qs[i] == 2 & beta[i] == 0 ){
        asy[i]
      } 
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
  asy <- iNEXT.3D:::FD_est(ai_vi, q, nT, ai_vi)
  obs <- iNEXT.3D:::FD_mle(ai_vi, q)
  out <- sapply(1:ncol(ai_vi$ai), function(i) {
    ai <- ai_vi$ai[, i]
    # ai[ai < 1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
    RFD_m = iNEXT.3D:::RFD(av, nT, nT - 1, q, V_bar)
    beta <- rep(0, length(q))
    asy_i <- asy[, i]
    obs_i <- obs[, i]
    
    
    # obs_i <- sapply(1:length(q), function(j) {
    #   max(RFD_m[j], obs_i[j])
    # })
    RFD_m <- sapply(1:length(q), function(j) {
      min(RFD_m[j], obs_i[j])
    })
    
    asy_i <- sapply(1:length(q), function(j) {
      max(asy_i[j], obs_i[j], RFD_m[j])
    })
    
    
    beta0plus <- which(asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])
    
    if (sum(m < nT) != 0) {
      int.m = sort(unique(c(floor(m[m < nT]), ceiling(m[m < nT]))))
      mRFD = rbind(int.m, sapply(int.m, function(k) iNEXT.3D:::RFD(av, nT, k, q, V_bar)))
      
      if (0 %in% int.m) mRFD[,mRFD[1,] == 0] = 0
    }
    
    sapply(m, function(mm) {
      if (mm < nT) {
        if (mm == round(mm)) {
          mRFD[-1, mRFD[1, ] == mm]
        } else {
          (ceiling(mm) - mm) * mRFD[-1, mRFD[1, ] == floor(mm)] + (mm - floor(mm)) * mRFD[-1, mRFD[1, ] == ceiling(mm)]
        }
      }
      else if (mm == nT) {
        obs_i
      }
      else if (mm == Inf) {
        asy_i
      }
      else {
        EFD(m = mm, qs = q, obs = obs_i, asy = asy_i, beta = beta, av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out, ncol = ncol(ai_vi$ai))
}


# Generate Color Palette for ggplot2
#
# This function creates a color palette suitable for ggplot2 visualizations by evenly spacing colors in the HCL color space. The function ensures that the colors are well-distributed and visually distinct, making it ideal for categorical data where each category needs to be represented by a different color.
#
# @param g An integer indicating the number of distinct colors to generate. This value should be a positive integer, with higher values resulting in a broader range of colors.
# @return A vector of color codes in hexadecimal format, suitable for use in ggplot2 charts and plots. The length of the vector will match the input parameter `g`.
# @examples
# # Generate a palette of 5 distinct colors
# ggplotColors(5)
#
# # Use the generated colors in a ggplot2 chart
# library(ggplot2)
# df <- data.frame(x = 1:5, y = rnorm(5), group = factor(1:5))
# ggplot(df, aes(x, y, color = group)) +
#   geom_point() +
#   scale_color_manual(values = ggplotColors(5))
#
ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}



#' Data Information for Reference Samples
#'
#' \code{DataInfobeta3D()} provides basic information about (1) the reference sample from each individual assemblage, (2) the gamma reference sample from the pooled assemblage, and (3) the alpha reference sample from the joint assemblage, for taxonomic (TD), phylogenetic (PD), and functional (FD) diversity.
#' @param data (a) For \code{datatype = "abundance"}, data for a single dataset can be provided as a species-by-assemblage abundance matrix/data.frame, with species in rows and assemblages in columns. Data for multiple datasets should be provided as a list of such matrices/data.frames, with each element representing one dataset.\cr
#' (b) For \code{datatype = "incidence_raw"}, data for a single dataset should be provided as a list of species-by-sampling-unit incidence matrices/data.frames, with each element representing one assemblage. Data for multiple datasets should be provided as a list of such lists, with each outer-list element representing one dataset.
#' @param diversity Selection of the diversity type: \code{diversity = "TD"} for taxonomic diversity, \code{diversity = "PD"} for phylogenetic diversity, and \code{diversity = "FD"} for functional diversity.
#' @param datatype Type of input data: individual-based abundance data (\code{datatype = "abundance"}) or sampling-unit-based incidence data (\code{datatype = "incidence_raw"}). For incidence raw data, each species-by-sampling-unit matrix contains only 0 (non-detection) or 1 (detection).
#' @param PDtree Required when \code{diversity = "PD"}. A phylogenetic tree (an R \code{phylo} object) spanning all observed species across the input datasets. The tree can be imported from a file in Newick format.
#' @param PDreftime Argument only for \code{diversity = "PD"} and \code{PDtype = "meanPD"}. A numeric value specifying a fixed reference time for phylogenetic diversity. If \code{PDreftime = NULL} (default), the root age of the phylogenetic tree spanning the pooled observed species is used.
#' @param PDtype Argument only for \code{diversity = "PD"}. Specifies the phylogenetic diversity measure: \code{PDtype = "meanPD"} for the effective number of equally divergent lineages at a specified reference time, or \code{PDtype = "AUC"} for an overall measure obtained by integrating and averaging mean phylogenetic diversity across reference times from zero to the root age of the phylogenetic tree. Default is \code{PDtype = "meanPD"}.
#' @param FDdistM Required when \code{diversity = "FD"}. A species pairwise distance matrix covering all species observed across the input datasets.
#' @param FDtype Argument only for \code{diversity = "FD"}. Specifies the functional diversity measure: \code{FDtype = "tau_value"} for FD at a specified threshold of functional distinctiveness, or \code{FDtype = "AUC"} for an overall FD measure obtained by integrating the tau-profile over threshold values from 0 to 1. Default is \code{FDtype = "AUC"}.
#' @param FDtau Argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}. A numeric value between 0 and 1 specifying the threshold of functional distinctiveness used to compute FD. If \code{FDtau = NULL} (default), the threshold is set to the mean pairwise functional distance between two individuals randomly selected from the pooled data across all datasets (i.e., Rao's quadratic entropy).
#' @param by_pair A logical variable specifying whether data information is computed for all assemblages as a whole or for every pair of assemblages. If \code{by_pair = FALSE}, data information is provided for all assemblages within each dataset considered as a whole; if \code{by_pair = TRUE}, the corresponding information is provided for every pair of assemblages. Default is \code{by_pair = FALSE}.
#'
#' @return A data frame containing basic information about reference samples.\cr\cr For abundance data, basic information shared by TD, PD, and FD includes the dataset name (\code{Dataset}), the pair of assemblages being compared (\code{Pair}, included only when \code{by_pair = TRUE}), individual/pooled/joint assemblages (\code{Assemblage}), sample size (\code{n}), observed species richness (\code{S.obs}), estimated sample coverage of the reference sample (\code{SC(n)}), and estimated sample coverage of twice the reference sample size (\code{SC(2n)}). Additional information for each diversity dimension is described below.\cr\cr
#'
#' (1) TD: The first five species abundance frequency counts in the reference sample (\code{f1}--\code{f5}).\cr\cr
#' (2) PD (\code{PDtype = "meanPD"}): The observed phylogenetic diversity (\code{PD.obs}), representing the observed total branch length for the assemblage under consideration; the observed mean phylogenetic diversity (\code{meanPD.obs}), calculated as \code{PD.obs / Reftime}; the numbers of singletons (\code{f1*}) and doubletons (\code{f2*}) in the node/branch abundance set; the corresponding total branch lengths (\code{g1} and \code{g2}); and the reference time (\code{Reftime}). If \code{PDreftime = NULL}, the root age of the phylogenetic tree spanning the pooled observed species is used.\cr\cr
#' (3) PD (\code{PDtype = "AUC"}): The integration end time (\code{Tmax}) for AUC-based phylogenetic diversity, corresponding to the root age of the phylogenetic tree spanning the pooled observed species. The integration starts at zero.\cr\cr
#' (4) FD (\code{FDtype = "AUC"}): The minimum distance (\code{dmin}) among all non-diagonal elements in the species pairwise distance matrix; the mean functional distance (\code{dmean}) between two individuals randomly selected from the assemblage under consideration; and the maximum pairwise distance (\code{dmax}) among species in the distance matrix.\cr\cr
#' (5) FD (\code{FDtype = "tau_value"}): The numbers of singletons (\code{a1*}) and doubletons (\code{a2*}) in the functional abundances at the specified threshold of functional distinctiveness (\code{Tau}), as well as the total contributions of functional singletons (\code{h1}) and doubletons (\code{h2}) at that threshold. If \code{FDtau = NULL}, the threshold is set to the mean pairwise functional distance between two individuals randomly selected from the pooled data across all datasets.\cr\cr
#'
#' For incidence raw data, basic information shared by TD, PD, and FD includes the dataset name (\code{Dataset}), the pair of assemblages being compared (\code{Pair}, included only when \code{by_pair = TRUE}), individual/pooled/joint assemblages (\code{Assemblage}), the number of sampling units (\code{T}), the total number of incidences (\code{U}), observed species richness (\code{S.obs}), estimated sample coverage of the reference sample (\code{SC(T)}), and estimated sample coverage of twice the reference sample size (\code{SC(2T)}).\cr\cr
#'
#' For TD, the output also includes the first five species incidence frequency counts in the reference sample (\code{Q1}--\code{Q5}). For PD and FD, additional output information is analogous to that described above for abundance data.



#' 
#' @examples
#' ## (Data Information) Taxonomic diversity for abundance data
#' data(Brazil_rainforests)
#' info_TD_abun = DataInfobeta3D(data = Brazil_rainforests[1 :2], diversity = 'TD', 
#'                               datatype = 'abundance')
#' info_TD_abun
#' 
#' 
#' 
#' ## (Data Information) Taxonomic diversity for incidence data (all assemblages as a whole)
#' data(Second_growth_forests)
#' info_TD_inci = DataInfobeta3D(data = Second_growth_forests, diversity = 'TD',
#'                               datatype = 'incidence_raw')
#' info_TD_inci
#' 
#' 
#' ## (Data Information) Taxonomic diversity for incidence data for every pair of assemblages
#' data(Second_growth_forests)
#' info_TD_inci_pair = DataInfobeta3D(data = Second_growth_forests, diversity = 'TD',
#'                                    datatype = 'incidence_raw', by_pair = TRUE)
#' info_TD_inci_pair
#' 
#' \donttest{
#'
#' ## (Data Information) Mean phylogenetic diversity for abundance data 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' info_meanPD_abun = DataInfobeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                                   datatype = 'abundance', PDtype = 'meanPD', PDtree = Brazil_tree,
#'                                   PDreftime = NULL)
#' info_meanPD_abun
#' 
#' ## (Data Information) Phylogenetic diversity (AUC) for abundance data 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' info_PDAUC_abun = DataInfobeta3D(data = Brazil_rainforests[1:2], diversity = 'PD', 
#'                                   datatype = 'abundance', PDtype = 'AUC',
#'                                   PDtree = Brazil_tree)
#' info_PDAUC_abun 
#' 
#' ## (Data Information) Functional diversity for abundance data under a specified threshold level
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' info_FDtau_abun = DataInfobeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                                  datatype = 'abundance', FDdistM = Brazil_distM, 
#'                                  FDtype = 'tau_value', FDtau = NULL)
#' info_FDtau_abun
#' 
#' ## (Data Information) Functional diversity for abundance data when all threshold levels
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' info_FDAUC_abun = DataInfobeta3D(data = Brazil_rainforests[1:2], diversity = 'FD', 
#'                                  datatype = 'abundance', FDdistM = Brazil_distM, FDtype = 'AUC')
#' info_FDAUC_abun
#' }
#' 
#' @export
DataInfobeta3D = function(data, diversity = 'TD', datatype = 'abundance', PDtype = 'meanPD', PDtree = NULL, PDreftime = NULL,
                          FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, by_pair = FALSE) {
  
  if (is.na(pmatch(diversity, c("TD", "PD", "FD")))) stop("invalid diversity")
  
  if (datatype == "incidence" | datatype == "incidence_freq") stop('Please try datatype = "incidence_raw".')  
  if(is.na(pmatch(datatype, c("abundance", "incidence_raw"))))
    stop("invalid datatype")
  
  if (! (is.null(PDreftime) | inherits(PDreftime, "numeric")))
    stop("invalid class of reference time, PDreftime should be a postive value of numeric object.", call. = FALSE)
  if (length(PDreftime) > 1)
    stop("PDreftime can only accept a value instead of a vector.", call. = FALSE)
  
  if(is.na(pmatch(PDtype, c("meanPD", "AUC"))))
    stop("Incorrect phylogenetic diversity type. Please use 'meanPD' or 'AUC'.",
         call. = FALSE)
  
  if (FDtype == "tau_values") stop('Please try FDtype = "tau_value".')  
  if(is.na(pmatch(FDtype, c("AUC", "tau_value"))))
    stop("Incorrect type of functional diversity type, please use either 'AUC' or 'tau_value'", call. = FALSE)
  
  if (! (is.null(FDtau) | inherits(FDtau, "numeric")))
    stop("invalid class of tau value, FDtau should be a postive value between zero and one.", call. = FALSE)
  if (length(FDtau) > 1)
    stop("FDtau only accept a value instead of a vector.", call. = FALSE)
  
  
  if (datatype == 'abundance') {
    
    if ( inherits(data, "data.frame") | inherits(data, "matrix") ) data = list(Dataset_1 = data)
    
    if ( inherits(data, "list")) {
      
      data = lapply(data, as.matrix)
      if (is.null(names(data))) names(data) = paste0("Dataset_", 1:length(data))
    }
  }
  
  if (datatype == 'incidence_raw') {
    
    if (!inherits(data, "list"))
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( inherits(data, "list") & (inherits(data[[1]], "data.frame") | inherits(data[[1]], "matrix")) ) data = list(Dataset_1 = data)
    
    if (! ( inherits(data[[1]][[1]], "data.frame") | inherits(data[[1]][[1]], "matrix") ) )
      stop("Invalid data format for incidence raw data. Please refer to example of iNEXTbeta3D.", call. = FALSE)
    
    if ( sum( sapply(1:length(data), function(i) ( length(unique(sapply(data[[i]], nrow))) != 1 | 
                                                   length(unique(sapply(data[[i]], ncol))) != 1 ) ) ) > 0 )
      stop("Number of species (row) or sampling units (column) should be the same within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
    
    data = lapply(data, function(x) {
      
      if (is.null(rownames(x[[1]]))) tmp = x else {
        
        nT = ncol(x[[1]])
        if (nT <= 3) stop("Number of sampling units of some datasets is too less. Please add more sampling units data.", call. = FALSE)
        
        tmp = lapply(x, function(i) data.frame(i) %>% rownames_to_column(var = "Species"))
        
        tmp1 = tmp[[1]]
        for (i in 2:length(tmp)) {
          tmp1 = full_join(tmp1, tmp[[i]], by = "Species")
        }
        tmp1 = tmp1 %>% column_to_rownames(var = "Species")
        
        if (nrow(tmp1) != nrow(x[[1]]))
          stop("Species names (rownames) should be matched within each dataset. Please check you data or refer to example of iNEXTbeta3D.", call. = FALSE)
        
        if (sum(!as.matrix(tmp1) %in% c(0,1)) > 0)
          stop("The data for datatype = 'incidence_raw' can only contain values zero (undetected) or one (detected). Please transform values to zero or one.", call. = FALSE)
        
        tmp = lapply(1:length(tmp), function(i) tmp1[, ((i-1)*nT+1):(i*nT)])
        names(tmp) = if (is.null(data)) paste0("Assemblage_", 1:length(data)) else names(x)
      }
      
      return(tmp)
    })
    
    if (is.null(names(data))) names(data) = paste0("Dataset_", 1:length(data))
  }
  
  
  if (datatype == 'abundance') {
    
    pool.name <- lapply(data, function(x) rownames(x)) %>% unlist %>% unique
    
  } else if (datatype == 'incidence_raw') {
    
    pool.name <- lapply(data, function(x) lapply(x, function(y) rownames(y))) %>% unlist %>% unique
    
  }
  
  if (diversity == "PD") {
    
    if (sum(c(duplicated(PDtree$tip.label), duplicated(PDtree$node.label[PDtree$node.label!=""])))>0)
      stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
    
    if ( is.null(pool.name) )
      stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% PDtree$tip.label) != length(pool.name))
      stop("Data and tree tip label contain unmatched species", call. = FALSE)
    
    
    if (datatype == "abundance") {
      
      pool.data = lapply(data, rowSums)
      
      pool.data = lapply(pool.data, function(x) {
        data.frame(value = x) %>%
          rownames_to_column(var = "Species")
      })
      
      pdata = pool.data[[1]]
      
      if (length(pool.data) > 1) {
        for (i in 2:length(pool.data)) {
          pdata = full_join(pdata, pool.data[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      
      pdata = pdata %>%
        column_to_rownames(var = "Species") %>%
        rowSums()
    }
    
    
    if (datatype == "incidence_raw") {
      
      pool.data = lapply(data, function(x) {
        tmp = Reduce("+", x)
        tmp[tmp > 1] = 1
        rowSums(tmp)
      })
      
      pool.data = lapply(pool.data, function(x) {
        data.frame(value = x) %>%
          rownames_to_column(var = "Species")
      })
      
      pdata = pool.data[[1]]
      
      if (length(pool.data) > 1) {
        for (i in 2:length(pool.data)) {
          pdata = full_join(pdata, pool.data[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      
      pdata = pdata %>%
        column_to_rownames(var = "Species") %>%
        rowSums()
    }
    
    
    pool.name = names(pdata[pdata > 0])
    
    tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
    
    mytree = ape::drop.tip(PDtree, tip)
    
    H_max = max(ape::node.depth.edgelength(mytree))
    
  }
  
  if (diversity == "FD") {
    
    if (is.null(rownames(FDdistM)))
      stop('The species names are not provided in distance matrix.', call. = FALSE)
    
    if( is.null(pool.name) )
      stop("Row names of data must be the species names that match row names in distance matrix and thus can not be empty.", call. = FALSE)
    
    if (sum(pool.name %in% rownames(FDdistM)) != length(pool.name))
      stop("Data and distance matrix contain unmatched species", call. = FALSE)
  }
  
  if (diversity == "TD") {
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        if (is.null(colnames(data[[i]]))) colnames(data[[i]]) = paste('Assemblage_', 1:ncol(data[[i]]), sep = "")
        
        single = DataInfo3D(data[[i]], datatype = "abundance")
        
        if (by_pair == FALSE) {
          
          multiple = DataInfo3D(list("Pooled assemblage" = rowSums(data[[i]]),
                                     "Joint assemblage" = as.vector(data[[i]])), datatype = "abundance")
          
        } else {
          
          two_idx <- combn(1:ncol(data[[i]]), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            ab_name = colnames(data[[i]])[two_idx[,j]]
            
            DataInfo3D(list("Pooled assemblage" = rowSums(data[[i]][,two_idx[,j]]),
                            "Joint assemblage" = as.vector(data[[i]][,two_idx[,j]])), datatype = "abundance") %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.) %>% .[,1:11])
        
        }) %>% do.call(rbind,.)
    }
    
    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        if (is.null(names(data[[i]]))) names(data[[i]]) = paste('Assemblage_', 1:length(data[[i]]), sep = "")
        
        single = DataInfo3D(data[[i]], datatype = "incidence_raw")
        
        if (by_pair == FALSE) {
          
          x = lapply(data[[i]], as.matrix)
          data_gamma = Reduce('+', x)
          data_gamma[data_gamma > 1] = 1
          
          data_alpha = do.call(rbind, x)
          
          multiple = DataInfo3D(list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                                     "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))), 
                                datatype = "incidence_freq")
          
        } else {
          
          two_idx <- combn(1:length(data[[i]]), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            ab_name = names(data[[i]])[two_idx[,j]]
            
            x = lapply(data[[i]][two_idx[,j]], as.matrix)
            data_gamma = Reduce('+', x)
            data_gamma[data_gamma > 1] = 1
            
            data_alpha = do.call(rbind, x)
            
            DataInfo3D(list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                            "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))), 
                       datatype = "incidence_freq") %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.) %>% .[,1:12])
        
        }) %>% do.call(rbind,.)
    }
    
  }
  
  if (diversity == "PD" & PDtype == "meanPD") {
    
    if(is.null(PDreftime)) PDreftime = H_max else if (PDreftime <= 0) { 
      stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
      }
    
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
        
        single = DataInfo3D(x, diversity = "PD", datatype = "abundance", PDtree = PDtree, PDreftime = PDreftime)
        
        single = single %>%
          mutate(meanPD.obs = PD.obs / PDreftime, .after = PD.obs)
        
        if (by_pair == FALSE) {
          
          data_gamma = rowSums(x)
          data_gamma = data_gamma[data_gamma > 0]
          
          aL_table_gamma = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = PDreftime)
          aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
          aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          N = ncol(x)
          aL_table_alpha = c()
          
          for (j in 1:N) {
            
            y = x[x[,j] > 0, j]
            names(y) = rownames(x)[x[,j]>0]
            
            aL_table = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = y, rootExtend = T, refT = PDreftime)
            aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
            aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table)
            
          }
          
          output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
            
            ai = y$branch.abun
            Li = y$branch.length
            I1 <- which(ai == 1 & Li > 0)
            I2 <- which(ai == 2 & Li > 0)
            S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
            PD_obs <- sum(Li)
            meanPD_obs <- PD_obs / PDreftime
            f1 <- length(I1)
            f2 <- length(I2)
            g1 <- sum(Li[I1])
            g2 <- sum(Li[I2])
            c(S.obs, PD_obs, meanPD_obs, f1, f2, g1, g2)
            
          }) %>% t()
          
          Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
          mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
          
          multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                            'n' = sum(x), 'S.obs' = output[,1], 'SC(n)' = Chat, 'SC(2n)' = mul_C2n, 'PD.obs' = output[,2],
                            'meanPD.obs' = output[,3], 'f1*' = output[,4], 'f2*' = output[,5], 'g1' = output[,6], 'g2' = output[,7],
                            'Reftime' = PDreftime)
          
        } else {
          
          two_idx <- combn(1:ncol(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            x = x[,two_idx[,j]]
            ab_name = colnames(x)
            
            data_gamma = rowSums(x)
            data_gamma = data_gamma[data_gamma > 0]
            
            aL_table_gamma = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = PDreftime)
            aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
            aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            N = 2
            aL_table_alpha = c()
            
            for (j in 1:N) {
              
              y = x[x[,j] > 0, j]
              names(y) = rownames(x)[x[,j]>0]
              
              aL_table = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = y, rootExtend = T, refT = PDreftime)
              aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
              aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
              
              ai = y$branch.abun
              Li = y$branch.length
              I1 <- which(ai == 1 & Li > 0)
              I2 <- which(ai == 2 & Li > 0)
              S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
              PD_obs <- sum(Li)
              meanPD_obs <- PD_obs / PDreftime
              f1 <- length(I1)
              f2 <- length(I2)
              g1 <- sum(Li[I1])
              g2 <- sum(Li[I2])
              c(S.obs, PD_obs, meanPD_obs, f1, f2, g1, g2)
              
            }) %>% t()
            
            Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
            mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
            
            tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                   'n' = sum(x), 'S.obs' = output[,1], 'SC(n)' = Chat, 'SC(2n)' = mul_C2n, 'PD.obs' = output[,2],
                   'meanPD.obs' = output[,3], 'f1*' = output[,4], 'f2*' = output[,5], 'g1' = output[,6], 'g2' = output[,7],
                   'Reftime' = PDreftime) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      }
    
    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        x = lapply(data[[i]], as.matrix)
        if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
        
        single = DataInfo3D(x, diversity = "PD", datatype = "incidence_raw", PDtree = PDtree, PDreftime = PDreftime)
        
        single = single %>%
          mutate(meanPD.obs = PD.obs / PDreftime, .after = PD.obs)
        
        if (by_pair == FALSE) {
          
          data_gamma = Reduce('+', x)
          data_gamma[data_gamma > 1] = 1
          
          aL_table_gamma = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma), datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
          aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
          aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          N = length(x)
          aL_table_alpha = c()
          
          for (j in 1:N) {
            
            aL_table = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x[[j]], datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
            aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
            aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            aL_table_alpha = rbind(aL_table_alpha, aL_table)
            
          }
          
          output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
            
            ai = y$branch.abun
            Li = y$branch.length
            I1 <- which(ai == 1 & Li > 0)
            I2 <- which(ai == 2 & Li > 0)
            S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
            PD_obs <- sum(Li)
            meanPD_obs <- PD_obs / PDreftime
            f1 <- length(I1)
            f2 <- length(I2)
            g1 <- sum(Li[I1])
            g2 <- sum(Li[I2])
            c(S.obs, PD_obs, meanPD_obs, f1, f2, g1, g2)
            
          }) %>% t()
          
          Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
          mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
          multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                            'T' = ncol(x[[1]]), 'U' = c(sum(data_gamma), sum(sapply(x, rowSums))), 
                            'S.obs' = output[,1], 'SC(T)' = Chat, 'SC(2T)' = mul_C2T, 'PD.obs' = output[,2],
                            'meanPD.obs' = output[,3], 'Q1*' = output[,4], 'Q2*' = output[,5], 'R1' = output[,6], 'R2' = output[,7],
                            'Reftime' = PDreftime)
          
        } else {
          
          two_idx <- combn(1:length(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            x = x[two_idx[,j]]
            ab_name = names(x)
            
            data_gamma = Reduce('+', x)
            data_gamma[data_gamma > 1] = 1
            
            aL_table_gamma = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma), datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
            aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
            aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            N = length(x)
            aL_table_alpha = c()
            
            for (j in 1:N) {
              
              aL_table = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x[[j]], datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
              aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
              aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
              
              ai = y$branch.abun
              Li = y$branch.length
              I1 <- which(ai == 1 & Li > 0)
              I2 <- which(ai == 2 & Li > 0)
              S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
              PD_obs <- sum(Li)
              meanPD_obs <- PD_obs / PDreftime
              f1 <- length(I1)
              f2 <- length(I2)
              g1 <- sum(Li[I1])
              g2 <- sum(Li[I2])
              c(S.obs, PD_obs, meanPD_obs, f1, f2, g1, g2)
              
            }) %>% t()
            
            Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
            mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
            
            tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                   'T' = ncol(x[[1]]), 'U' = c(sum(data_gamma), sum(sapply(x, rowSums))), 
                   'S.obs' = output[,1], 'SC(T)' = Chat, 'SC(2T)' = mul_C2T, 'PD.obs' = output[,2],
                   'meanPD.obs' = output[,3], 'Q1*' = output[,4], 'Q2*' = output[,5], 'R1' = output[,6], 'R2' = output[,7],
                   'Reftime' = PDreftime) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      
    }
    
    rownames(output) = NULL
  }
  
  if (diversity == "PD" & PDtype == "AUC") {
    

    PDreftime = H_max

    
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        
        if (is.null(colnames(x)))colnames(x) = paste("Assemblage_", 1:ncol(x), sep = "")
        
        single = DataInfo3D(x, datatype = "abundance")
        single = single[, c("Assemblage","n","S.obs","SC(n)","SC(2n)")]
        
        single$Tmax = PDreftime
        
        if (by_pair == FALSE) {
          
          multiple = DataInfo3D(list("Pooled assemblage" = rowSums(x), "Joint assemblage" = as.vector(x)), datatype = "abundance")
          
          multiple = multiple[, c("Assemblage","n","S.obs","SC(n)","SC(2n)")]
          
          multiple$Tmax = PDreftime
          
        } else {
          
          two_idx = combn(1:ncol(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            xx = x[, two_idx[, j], drop = FALSE]
            
            ab_name = colnames(xx)
            
            tmp = DataInfo3D(list("Pooled assemblage" = rowSums(xx), "Joint assemblage" = as.vector(xx)), datatype = "abundance")
            
            tmp = tmp[, c("Assemblage","n","S.obs","SC(n)","SC(2n)")]
            
            tmp$Tmax = PDreftime
            
            tmp = cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "), tmp)
            
            return(tmp)
            
          }) %>% do.call(rbind, .)
          
          single = cbind(Pair = "", single)
        }
        
        
        out = rbind(single, multiple) %>%
          cbind(Dataset = names(data)[i], .)
        
        return(out)
        
      }) %>% do.call(rbind, .)
    }
    
    

    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        x = lapply(data[[i]], as.matrix)
        
        if (is.null(names(x))) {
          names(x) = paste("Assemblage_", 1:length(x), sep = "")
        }
        
        single = DataInfo3D(x, datatype = "incidence_raw")
        
        single = single[, c("Assemblage", "T", "U", "S.obs", "SC(T)", "SC(2T)")]
        
        single$Tmax = PDreftime
        
        
        ## Multiple assemblages
        if (by_pair == FALSE) {
          
          data_gamma = Reduce("+", x)
          data_gamma[data_gamma > 1] = 1
          
          data_alpha = do.call(rbind, x)
          
          multiple = DataInfo3D(
            list("Pooled assemblage" =c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                 "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))),
                  datatype = "incidence_freq")
          
          multiple = multiple[, c("Assemblage", "T", "U", "S.obs", "SC(T)", "SC(2T)")]
          
          multiple$Tmax = PDreftime
          
        } else {
          
          two_idx = combn(1:length(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            xx = x[two_idx[, j]]
            
            ab_name = names(xx)
            
            data_gamma = Reduce("+", xx)
            data_gamma[data_gamma > 1] = 1
            
            data_alpha = do.call(rbind, xx)
            
            tmp = DataInfo3D(
              list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                   "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))),
                    datatype = "incidence_freq")
            
            tmp = tmp[, c("Assemblage", "T", "U", "S.obs", "SC(T)", "SC(2T)")]
            
            tmp$Tmax = PDreftime
            
            tmp = cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "), tmp)
            
            return(tmp)
            
          }) %>% do.call(rbind, .)
          
          single = cbind(Pair = "", single)
        }
        
        out = rbind(single, multiple) %>%
          cbind(Dataset = names(data)[i], .)
        
        return(out)
        
      }) %>% do.call(rbind, .)
    }
    
    
    rownames(output) = NULL
  }
    
  if (diversity == "FD" & FDtype == "tau_value") {
    
    FDdistM = as.matrix(FDdistM)
    
    if (is.null(FDtau) == T) {
      
      if (datatype == 'abundance') {
        
        tmp <- lapply(data, rowSums)
        tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
        pdata = tmp[[1]]
        
        if (length(tmp) > 1) {
          for(i in 2:length(tmp)){
            pdata = full_join(pdata, tmp[[i]], by = "Species")
          }
        }
        
        pdata[is.na(pdata)] = 0
        pdata = pdata %>% column_to_rownames("Species")
        pdata = rowSums(pdata)
        
        order_sp <- match(names(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- matrix(pdata/sum(pdata), ncol = 1)
        
      } else if (datatype == 'incidence_raw') {
        
        tmp <- lapply(data, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) })
        tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
        pdata = tmp[[1]]
        
        if (length(tmp) > 1) {
          for(i in 2:length(tmp)){
            pdata = full_join(pdata, tmp[[i]], by = "Species")
          }
        }
        
        pdata[is.na(pdata)] = 0
        pdata = pdata %>% column_to_rownames("Species")
        pdata = rowSums(pdata)
        
        order_sp <- match(names(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- matrix(pdata/sum(pdata), ncol = 1)
        
      }
      
      FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
    }
    
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
        
        single = DataInfo3D(x, diversity = "FD", datatype = "abundance", FDtype = "tau_values", FDtau = FDtau, FDdistM = FDdistM)
        
        
        if (by_pair == FALSE) {
          
          dij = FDdistM
          dij = dij[rownames(dij) %in% rownames(x), colnames(dij) %in% rownames(x)]
          order_sp <- match(rownames(x), rownames(dij))
          dij <- dij[order_sp, order_sp]
          
          dij = dij[rowSums(x)>0, rowSums(x)>0]
          x = x[rowSums(x)>0,]
          
          
          if (FDtau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(x)
          } else {
            dij[which(dij > FDtau, arr.ind = T)] = FDtau
            aik = (1 - dij/FDtau) %*% as.matrix(x)
          }
          
          positive_id = rowSums(aik) > 0
          
          
          gamma_x = rowSums(x)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = gamma_x/gamma_a
          
          alpha_a = as.vector(aik)
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          alpha_v = rep(gamma_v, ncol(x))
          
          Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
          mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
          
          multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                            'n' = sum(x), 
                            'S.obs' = c(sum(rowSums(x) > 0), sum(as.vector(x) > 0)), 
                            'SC(n)' = Chat, 'SC(2n)' = mul_C2n,
                            'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                            'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                            'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                            'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                            'Tau' = FDtau)
          
        } else {
          
          two_idx <- combn(1:ncol(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            x = x[,two_idx[,j]]
            ab_name = colnames(x)
            
            dij = FDdistM
            dij = dij[rownames(dij) %in% rownames(x), colnames(dij) %in% rownames(x)]
            order_sp <- match(rownames(x), rownames(dij))
            dij <- dij[order_sp, order_sp]
            
            dij = dij[rowSums(x)>0, rowSums(x)>0]
            x = x[rowSums(x)>0,]
            
            
            if (FDtau == 0) {
              dij[dij > 0] <- 1
              aik = (1 - dij/1) %*% as.matrix(x)
            } else {
              dij[which(dij > FDtau, arr.ind = T)] = FDtau
              aik = (1 - dij/FDtau) %*% as.matrix(x)
            }
            
            positive_id = rowSums(aik) > 0
            
            
            gamma_x = rowSums(x)[positive_id]
            gamma_a = rowSums(aik)[positive_id]
            gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
            gamma_v = gamma_x/gamma_a
            
            alpha_a = as.vector(aik)
            alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
            alpha_v = rep(gamma_v, ncol(x))
            
            Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
            mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
            
            tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                   'n' = sum(x), 
                   'S.obs' = c(sum(rowSums(x) > 0), sum(as.vector(x) > 0)), 
                   'SC(n)' = Chat, 'SC(2n)' = mul_C2n,
                   'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                   'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                   'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                   'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                   'Tau' = FDtau) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
    }
    
    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        x = lapply(data[[i]], as.matrix)
        if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
        
        single = DataInfo3D(x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_values", FDtau = FDtau, FDdistM = FDdistM)
        
        nT = ncol(x[[1]])
        
        if (by_pair == FALSE) {
          
          data_gamma = Reduce('+', x)
          data_gamma[data_gamma > 1] = 1
          data_gamma_freq = rowSums(data_gamma)
          
          data_2D = sapply(x, rowSums)
          
          ##
          dij = FDdistM
          dij = dij[rownames(dij) %in% names(data_gamma_freq), colnames(dij) %in% names(data_gamma_freq)]
          order_sp <- match(names(data_gamma_freq), rownames(dij))
          dij <- dij[order_sp, order_sp]
          
          dij = dij[data_gamma_freq > 0, data_gamma_freq > 0]
          data_gamma_freq = data_gamma_freq[data_gamma_freq > 0]
          
          
          if (FDtau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(data_gamma_freq)
          } else {
            dij[which(dij > FDtau, arr.ind = T)] = FDtau
            gamma_a = (1 - dij/FDtau) %*% as.matrix(data_gamma_freq)
          }
          
          gamma_a[gamma_a > nT] = nT
          gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
          gamma_v = data_gamma_freq/gamma_a
          
          
          ##
          dij = FDdistM
          dij = dij[rownames(dij) %in% rownames(data_2D), colnames(dij) %in% rownames(data_2D)]
          order_sp <- match(rownames(data_2D), rownames(dij))
          dij <- dij[order_sp, order_sp]
          
          dij = dij[rowSums(data_2D) > 0, rowSums(data_2D) > 0]
          data_2D = data_2D[rowSums(data_2D) > 0,]
          
          
          if (FDtau == 0) {
            dij[dij > 0] <- 1
            alpha_a = (1 - dij/1) %*% as.matrix(data_2D)
          } else {
            dij[which(dij > FDtau, arr.ind = T)] = FDtau
            alpha_a = (1 - dij/FDtau) %*% as.matrix(data_2D)
          }
          
          alpha_a[alpha_a > nT] = nT
          alpha_a = as.vector(alpha_a)
          alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
          alpha_v = rep(gamma_v, length(x))
          
          Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
          mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
          
          multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                            'T' = rep(ncol(data_gamma), 2), 
                            'U' = c(sum(data_gamma_freq), sum(data_2D)),
                            'S.obs' = c(sum(rowSums(data_gamma) > 0), sum(data_2D > 0)), 
                            'SC(T)' = Chat, 'SC(2T)' = mul_C2T,
                            'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                            'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                            'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                            'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                            'Tau' = FDtau)
          
        } else {
          
          two_idx <- combn(1:length(x), 2)
          
          multiple = lapply(1:ncol(two_idx), function(j) {
            
            x = x[two_idx[,j]]
            ab_name = names(x)
            
            data_gamma = Reduce('+', x)
            data_gamma[data_gamma > 1] = 1
            data_gamma_freq = rowSums(data_gamma)
            
            data_2D = sapply(x, rowSums)
            
            ##
            dij = FDdistM
            dij = dij[rownames(dij) %in% names(data_gamma_freq), colnames(dij) %in% names(data_gamma_freq)]
            order_sp <- match(names(data_gamma_freq), rownames(dij))
            dij <- dij[order_sp, order_sp]
            
            dij = dij[data_gamma_freq > 0, data_gamma_freq > 0]
            data_gamma_freq = data_gamma_freq[data_gamma_freq > 0]
            
            
            if (FDtau == 0) {
              dij[dij > 0] <- 1
              gamma_a = (1 - dij/1) %*% as.matrix(data_gamma_freq)
            } else {
              dij[which(dij > FDtau, arr.ind = T)] = FDtau
              gamma_a = (1 - dij/FDtau) %*% as.matrix(data_gamma_freq)
            }
            
            gamma_a[gamma_a > nT] = nT
            gamma_a = ifelse(gamma_a < 1 & gamma_a > 0, 1, round(gamma_a))
            gamma_v = data_gamma_freq/gamma_a
            
            
            ##
            dij = FDdistM
            dij = dij[rownames(dij) %in% rownames(data_2D), colnames(dij) %in% rownames(data_2D)]
            order_sp <- match(rownames(data_2D), rownames(dij))
            dij <- dij[order_sp, order_sp]
            
            dij = dij[rowSums(data_2D) > 0, rowSums(data_2D) > 0]
            data_2D = data_2D[rowSums(data_2D) > 0,]
            
            
            if (FDtau == 0) {
              dij[dij > 0] <- 1
              alpha_a = (1 - dij/1) %*% as.matrix(data_2D)
            } else {
              dij[which(dij > FDtau, arr.ind = T)] = FDtau
              alpha_a = (1 - dij/FDtau) %*% as.matrix(data_2D)
            }
            
            alpha_a[alpha_a > nT] = nT
            alpha_a = as.vector(alpha_a)
            alpha_a = ifelse(alpha_a < 1 & alpha_a > 0, 1, round(alpha_a))
            alpha_v = rep(gamma_v, length(x))
            
            
            Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
            mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
            
            tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                   'T' = rep(ncol(data_gamma), 2), 
                   'U' = c(sum(data_gamma_freq), sum(data_2D)),
                   'S.obs' = c(sum(rowSums(data_gamma) > 0), sum(data_2D > 0)), 
                   'SC(T)' = Chat, 'SC(2T)' = mul_C2T,
                   'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                   'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                   'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                   'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                   'Tau' = FDtau) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
          single = single %>% cbind(Pair = "",.)
        }
        
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      
    }
    
    rownames(output) = NULL
  }
  
  if (diversity == "FD" & FDtype == "AUC") {
    
    FDdistM = as.matrix(FDdistM)
    
    # if (datatype == 'abundance') {
    #   
    #   tmp <- lapply(data, rowSums)
    #   tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
    #   pdata = tmp[[1]]
    #   for(i in 2:length(tmp)){
    #     pdata = full_join(pdata, tmp[[i]], by = "Species")
    #   }
    #   pdata[is.na(pdata)] = 0
    #   pdata = pdata %>% column_to_rownames("Species")
    #   pdata = rowSums(pdata)
    #   
    #   order_sp <- match(names(pdata),rownames(FDdistM))
    #   FDdistM <- FDdistM[order_sp,order_sp]
    #   pdata <- matrix(pdata/sum(pdata), ncol = 1)
    #   
    # } else if (datatype == 'incidence_raw') {
    #   
    #   tmp <- lapply(data, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) })
    #   tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
    #   pdata = tmp[[1]]
    #   for(i in 2:length(tmp)){
    #     pdata = full_join(pdata, tmp[[i]], by = "Species")
    #   }
    #   pdata[is.na(pdata)] = 0
    #   pdata = pdata %>% column_to_rownames("Species")
    #   pdata = rowSums(pdata)
    #   
    #   order_sp <- match(names(pdata),rownames(FDdistM))
    #   FDdistM <- FDdistM[order_sp,order_sp]
    #   pdata <- matrix(pdata/sum(pdata), ncol = 1)
    #   
    # }
    # 
    # dmin <- min(FDdistM[FDdistM > 0])
    # dmean <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
    # dmax <- max(FDdistM[FDdistM > 0])
    # Tau = c(dmin, dmean, dmax)
    
    output = lapply(1:length(data), function(i) {
      
      x = data[[i]]
      
      if (by_pair == FALSE) {
        
        if (datatype == "abundance") {
          
          if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
          
          single = DataInfo3D(x, diversity = "FD", datatype = "abundance", FDtype = "AUC", FDdistM = FDdistM)
          
          pdata = data.frame(rowSums(x))
          
          order_sp <- match(rownames(pdata),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          pdata <- pdata/sum(pdata)
          
          dmean <- sum ( (as.matrix(pdata) %*% t(as.matrix(pdata)) ) * FDdistM) # dmean
          
          dmin <- min(FDdistM[lower.tri(FDdistM)])
          dmax <- max(FDdistM[FDdistM > 0])
          
          mul_C2n = sapply(list("Pooled assemblage" = rowSums(x),
                                "Joint assemblage" = as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
          
          multiple <- cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = rowSums(x),
                                                   "Joint assemblage" = as.vector(x)), "abundance")[,1:5], 
                            "dmin" = dmin,
                            "dmean" = dmean,
                            "dmax" = dmax)
        }
        
        if (datatype == "incidence_raw") {
          
          if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
          
          single = DataInfo3D(x, diversity = "FD", datatype = "incidence_raw", FDtype = "AUC", FDdistM = FDdistM)
          
          x = lapply(x, as.matrix)
          
          tmp = Reduce('+', x)
          tmp[tmp > 1] = 1
          pdata = data.frame(rowSums(tmp))
          
          order_sp <- match(rownames(pdata),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          pdata <- pdata/sum(pdata)
          
          dmean <- sum ( (as.matrix(pdata) %*% t(as.matrix(pdata)) ) * FDdistM) # dmean
          
          dmin <- min(FDdistM[lower.tri(FDdistM)])
          dmax <- max(FDdistM[FDdistM > 0])
          
          data_gamma = Reduce('+', x)
          data_gamma[data_gamma > 1] = 1
          
          data_alpha = do.call(rbind, x)
          
          mul_C2T = sapply(list(data_gamma, data_alpha), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
          
          multiple <- cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                                                   "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))), "incidence_freq")[,1:6],
                            "dmin" = dmin,
                            "dmean" = dmean,
                            "dmax" = dmax)
        }
        
      } else {
        
        if (datatype == "abundance") {
          
          if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
          
          single = DataInfo3D(x, diversity = "FD", datatype = "abundance", FDtype = "AUC", FDdistM = FDdistM)
          
          two_idx <- combn(1:ncol(x), 2)
          
          multiple <- lapply(1:ncol(two_idx), function(j) {
            
            x = x[,two_idx[,j]]
            ab_name = colnames(x)
            
            pdata = data.frame(rowSums(x))
            
            order_sp <- match(rownames(pdata),rownames(FDdistM))
            FDdistM <- FDdistM[order_sp,order_sp]
            pdata <- pdata/sum(pdata)
            
            dmean <- sum ( (as.matrix(pdata) %*% t(as.matrix(pdata)) ) * FDdistM) # dmean
            
            dmin <- min(FDdistM[lower.tri(FDdistM)])
            dmax <- max(FDdistM[FDdistM > 0])
            
            mul_C2n = sapply(list("Pooled assemblage" = rowSums(x),
                                  "Joint assemblage" = as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
            
            cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = rowSums(x),
                                         "Joint assemblage" = as.vector(x)), "abundance")[,1:5], 
                  "dmin" = dmin,
                  "dmean" = dmean,
                  "dmax" = dmax) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
        }
        
        if (datatype == "incidence_raw") {
          
          if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
          
          single = DataInfo3D(x, diversity = "FD", datatype = "incidence_raw", FDtype = "AUC", FDdistM = FDdistM)
          
          two_idx <- combn(1:length(x), 2)
          
          multiple <- lapply(1:ncol(two_idx), function(j) {
            
            x = lapply(x, as.matrix)
            
            x = x[two_idx[,j]]
            ab_name = names(x)
            
            tmp = Reduce('+', x)
            tmp[tmp > 1] = 1
            pdata = data.frame(rowSums(tmp))
            
            order_sp <- match(rownames(pdata),rownames(FDdistM))
            FDdistM <- FDdistM[order_sp,order_sp]
            pdata <- pdata/sum(pdata)
            
            dmean <- sum ( (as.matrix(pdata) %*% t(as.matrix(pdata)) ) * FDdistM) # dmean
            
            dmin <- min(FDdistM[lower.tri(FDdistM)])
            dmax <- max(FDdistM[FDdistM > 0])
            
            data_gamma = Reduce('+', x)
            data_gamma[data_gamma > 1] = 1
            
            data_alpha = do.call(rbind, x)
            
            mul_C2T = sapply(list(data_gamma, data_alpha), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
            
            cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                                         "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))), "incidence_freq")[,1:6],
                  "dmin" = dmin,
                  "dmean" = dmean,
                  "dmax" = dmax) %>% 
              cbind(Pair = paste(ab_name[1], ab_name[2], sep = " vs. "),.)
            
          }) %>% do.call(rbind,.)
          
        }
        
        single = single %>% cbind(Pair = "",.)
        
      }
      
      out = rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.)
      
      return(out)
      
      }) %>% do.call(rbind,.)
    
    rownames(output) = NULL
    }
  
  return(output)
  
  }




#' Printing \code{iNEXTbeta3D} object
#' 
#' \code{print.iNEXTbeta3D}: Print method for objects inheriting from class \code{"iNEXTbeta3D"}
#' @param x an \code{iNEXTbeta3D} object computed by \code{iNEXTbeta3D}.
#' @param ... additional arguments.
#' @return Prints the diversity and dissimilarity estimation tables contained in the input \code{iNEXTbeta3D} object to the console and invisibly returns \code{NULL}.
#' @export
print.iNEXTbeta3D <- function(x, ...){
  
  res <- lapply(x, function(y) {
    
    tmp = lapply(1:length(y), function(i) {
      
      # index1 = which(y[[i]]$SC == min(y[[i]]$SC))
      # index2 = which( (y[[i]]$Order.q != 0 & y[[i]]$SC == max(y[[i]]$SC)) | (y[[i]]$Order.q == 0 & y[[i]]$SC == max(y[[i]][y[[i]]$Order.q == 0, "SC"])) )
      # index3 = which(y[[i]]$Method %in% c("Observed", "Observed_SC(n, alpha)"))
      # 
      # index = sort(unique(c(index1, index1 + 1, index2, index2 - 1, index3 - 1, index3, index3 + 1)))
      # if (length(index) > 0) out = y[[i]][index,] else out = y[[i]]
      
      out = y[[i]]
      
      num_cols <- sapply(out, is.numeric)
      out[, num_cols] <- lapply(out[, num_cols, drop = FALSE], round, 4)

      ## Added by KaiHsiang 2024/11/11
      if ("Pair" %in% colnames(out)) {
        
        
        lapply(unique(out$Order.q), function(q) rbind(matrix(c("", "", paste("Order q = ", q, sep = ""), rep("", ncol(out)-3)), nrow = 1) %>% 
                                                        set_colnames(colnames(out)), 
                                                      out %>% filter(Order.q == q))) %>% do.call(rbind,.)
        
      } else {
        
        
        lapply(unique(out$Order.q), function(q) rbind(matrix(c("", paste("Order q = ", q, sep = ""), rep("", ncol(out)-2)), nrow = 1) %>% 
                                                        set_colnames(colnames(out)), 
                                                      out %>% filter(Order.q == q))) %>% do.call(rbind,.)
      }
      
    })
    
    names(tmp) = names(y)
    
    return(tmp)
    }
  )
  
  print(res)
  # cat("\n")
  # cat("NOTE: The above output only shows seven rows for each table; call any table in iNEXT.beta3D.object to view complete output.\n", sep = "")
  return(invisible())
}


## ========== no visible global function definition for R CMD check ========== ##

## ========== no visible global variable binding for R CMD check ========== ##
utils::globalVariables(c(".", "branch.abun", "branch.length", "tgroup", "Method",
  "Size", "SC", "Estimate", "Dataset", "div_type", "LCL", "UCL",
  "Order.q", "s.e.", "label", "Alpha", "Alpha_meanPD", "Diversity", "Gamma", "Gamma_meanPD",
  "PD", "PD.obs", "PlotGroup", "mT", "obs_size"
))



