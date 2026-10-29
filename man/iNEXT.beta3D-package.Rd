\name{iNEXT.beta3D-package}
\alias{iNEXT.beta3D-package}
\alias{iNEXT.beta3D}
\docType{package}
\title{
  Interpolation and extrapolation with beta diversity for three dimensions of biodiversity
}


\description{
  This package \code{iNEXT.beta3D} (iNterpolation and EXTrapolation with beta diversity for three dimensions of biodiversity)
  is a sequel to \code{iNEXT}. The three dimensions (3D) of biodiversity include taxonomic diversity (TD),
  phylogenetic diversity (PD), and functional diversity (FD). An online version, "iNEXT.beta3D Online"
  (https://chao.shinyapps.io/iNEXTbeta3D/), is also available for users without an R background.\cr

  A unified framework based on Hill numbers (for TD) and their generalizations, Hill-Chao
  numbers (for PD and FD), is adopted in this package. TD quantifies the effective number of
  species, meanPD (PD divided by reference time) quantifies the effective number of lineages,
  and FD quantifies the effective number of virtual functional groups (or functional
  "species"). Thus, TD, meanPD, and FD are all expressed in comparable units of
  species/lineage equivalents, allowing meaningful comparisons across the three
  dimensions of biodiversity. See Chao et al. (2021) for a review of this unified framework.\cr

  For each of the three dimensions, \code{iNEXT.beta3D} focuses on the multiplicative
  decomposition of diversity (beta = gamma/alpha, expressed in units of effective number of
  assemblages) for orders q = 0, 1, and 2, based on sampling data. The resulting 3D beta
  diversity quantifies the extent of among-assemblage differentiation, reflecting changes in
  species, lineages, or functional-groups composition and their abundances across
  assemblages.\cr

  \code{iNEXT.beta3D} provides standardized 3D estimates at a common sample size (for alpha
  and gamma diversity only, because sample-size-based beta diversity is not generally valid)
  or at a common sample coverage (for alpha, beta, and gamma diversity).
  \code{iNEXT.beta3D} also provides standardized estimates for four classes of
  coverage-based dissimilarity measures.\cr

  Based on rarefaction and extrapolation (R/E) theory for Hill-numbers-based TD of orders
  q = 0, 1, and 2, Chao et al. (2023b) developed the corresponding R/E theory for taxonomic
  beta diversity, with applications to real-world spatial, temporal, and spatio-temporal data.
  An application to Gentry's global forest data together with a concise description of the
  theory is provided in Chao et al. (2023a). The extension to phylogenetic and functional
  beta diversity follows a generally parallel framework.
  \code{iNEXT.beta3D} features two types of R/E sampling curves:

  \itemize{
    \item Sample-size-based (or size-based) R/E sampling curves: These curves plot standardized 3D
    gamma and alpha diversity as functions of sample size. Note that size-based beta
    diversity is not generally a statistically valid measure (Chao et al. 2023b);
    therefore, the corresponding beta diversity sampling curve is not provided.

    \item Sample-coverage-based (or coverage-based) R/E sampling curves: These curves
    plot standardized 3D gamma, alpha, and beta diversity, as well as four classes of
    dissimilarity measures, as functions of sample coverage, an objective measure of
    sample completeness.
  }

  Note that sufficient data are required to run \code{iNEXT.beta3D}. If a dataset contains
  only a few species, along with their abundances, phylogenies, or traits, it may not
  contain sufficient information for reliable \code{iNEXT.beta3D} analyses.\cr

  This package contains three main functions:

  \enumerate{
    \item \code{iNEXTbeta3D()} computes standardized 3D estimates at common sample sizes
    (for gamma and alpha diversity) or common sample coverage levels (for gamma, alpha,
    and beta diversity), as well as four classes of coverage-based dissimilarity measures.
    Users can obtain estimates at default or user-specified sample sizes or coverage levels.

    \item \code{ggiNEXTbeta3D()} visualizes the output from \code{iNEXTbeta3D()}.

    \item \code{DataInfobeta3D()} provides basic information about (1) the reference sample
    from each individual assemblage, (2) the gamma reference sample from the pooled
    assemblage, and (3) the alpha reference sample from the joint assemblage.
  }
}


\author{
Anne Chao, Yu-Ti Lee, Keng-Lei Lin, Po-Yen Chuang \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{

Chao, A., Chiu, C.-H., Hu, K.-H., and Zelený, D. (2023a).
Revisiting Alwyn H. Gentry's forest transect data: a statistical sampling-model-based approach.
Japanese Journal of Statistics and Data Science, 6, 861-884.
https://doi.org/10.1007/s42081-023-00214-1 \cr\cr

Chao, A., Chiu, C.-H., Villéger, S., Sun, I.-F., Thorn, S., Lin, Y.-C., Chiang, J. M., and Sherwin, W. B. (2019).
An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures.
Ecological Monographs, 89, e01343. \cr\cr

Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M., and Magurran, A. E. (2021).
Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity
and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr

Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zelený, D., Colwell, R. K., and Magurran, A. E. (2023b).
Rarefaction and extrapolation with beta diversity under a framework of Hill numbers:
The iNEXT.beta3D standardization. Ecological Monographs e1588.
https://doi.org/10.1002/ecm.1588
}
\keyword{ package }
