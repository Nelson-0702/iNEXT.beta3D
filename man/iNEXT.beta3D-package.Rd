\name{iNEXT.beta3D-package}
\alias{iNEXT.beta3D-package}
\alias{iNEXT.beta3D}
\docType{package}
\title{
  Interpolation and extrapolation with beta diversity for three dimensions of biodiversity
}
\description{
  This package \code{iNEXT.beta3D} (iNterpolation and EXTrapolation with beta diversity for three dimensions of biodiversity)
  is a sequel to \code{iNEXT}. Here the three dimensions (3D) of diversity include taxonomic diversity (TD),
  phylogenetic diversity (PD) and functional diversity (FD). An online version "iNEXT.beta3D Online"
  (https://chao.shinyapps.io/iNEXT_beta3D/) is also available for users without an R background.\cr
  
  A unified framework based on Hill numbers (for TD) and their generalizations (Hill-Chao numbers, for PD and FD) 
  is adopted to quantify 3D. In this framework, TD quantifies the effective number of species, PD quantifies the 
  effective total branch length or total evolutionary history, mean-PD (PD divided by 
  tree depth) quantifies the effective number of lineages, and FD quantifies the effective number of virtual functional
  groups (or functional "species"). Thus, TD, mean-PD, and FD are all in the same units of species/lineage equivalents
  and can be meaningfully compared; see Chao et al. (2021) for a review of the unified framework. \cr
  
  For each of the three dimensions, \code{iNEXT.beta3D} focuses on the multiplicative diversity decomposition 
  (gamma = alpha x beta) of orders q = 0, 1 and 2 based on sampling data. Beta diversity quantifies the extent 
  of among-assemblage differentiation, or the changes in species/lineages/functional-groups composition and
  their abundance among assemblages. \code{iNEXT.beta3D} features standardized 3D estimates with a common sample size
  (for alpha and gamma diversity) or sample coverage (for alpha, beta and gamma diversity). \code{iNEXT.beta3D} also 
  features standardized estimates of four classes of coverage-based dissimilarity measures. Based on the 
  rarefaction and extrapolation (R/E) method for Hill numbers (TD) for q = 0, 1 and 2, Chao et al. (2023b) 
  developed the pertinent R/E theory for taxonomic beta diversity with applications to real-world spatial,
  temporal and spatio-temporal data. An application to Gentry's global forest data along with a concise 
  description of the theory is provided in Chao et al. (2023a). The extension to phylogenetic and functional
  beta diversity is generally parallel.\cr
  
  \code{iNEXT.beta3D} also features two types of R/E sampling curves:
  
  *	Sample-size-based (or size-based) R/E curves: This type of sampling curve plots standardized 3D gamma and 
  alpha diversity with respect to sample size. Note that the size-based beta diversity is not a statistically 
  valid measure (Chao et al. 2023b) and thus the corresponding sampling curve is not provided. 
  
  *	Coverage-based R/E sampling curves: This type of sampling curve plots standardized 3D gamma, alpha, and
  beta diversity as well as four classes of dissimilarity measures with respect to sample coverage 
  (an objective measure of sample completeness).\cr
  
  Sufficient data are needed to run \code{iNEXT.beta3D}. If your data comprise only a few species and their 
  abundances/phylogenies/traits, it is probable that the data lack sufficient information to run \code{iNEXT.beta3D}.  
  
  
  This package contains three main functions: \cr
  
  1. \code{iNEXTbeta3D} computes standardized 3D estimates with a common sample size (for alpha and gamma diversity)
  or sample coverage (for alpha, beta and gamma diversity) for default sample sizes or coverage values. This function 
  also computes coverage-based standardized 3D estimates of four classes of dissimilarity measures for default coverage 
  values. In addition, this function also computes standardized 3D 
  estimates with a particular vector of user-specified sample sizes or coverage values. \cr
  
  2. \code{ggiNEXTbeta3D} visualizes the output from the function \code{iNEXTbeta3D}. \cr
  
  3. \code{DataInfobeta3D} provides basic data information for (1) the reference sample in each assemblage, 
  (2) the gamma reference sample in the pooled assemblage, and (3) the alpha reference sample in the joint assemblage.
}
\author{
Anne Chao, Yu-Ti Lee, Keng-Lei Lin, Po-Yen Chuang \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
Chao, A., Chiu, C.-H., Hu, K.-H., and Zeleny, D. (2023a). Revisiting Alwyn H. Gentry's forest transect data: 
  a statistical sampling-model-based approach. Japanese Journal of Statistics and Data Science, 6, 861-884. 
https://doi.org/10.1007/s42081-023-00214-1 \cr\cr
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M. and Magurran, A. E. (2021).
Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity 
and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023b). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: The iNEXT.beta3D standardization. Ecological Monographs e1588. https://doi.org/10.1002/ecm.1588 \cr\cr
}
\keyword{ package }
