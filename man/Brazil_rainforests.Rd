\name{Brazil_rainforests}
\alias{Brazil_rainforests}
\docType{data}

\title{
  Brazil rainforests data (datatype = "abundance" for TD, PD, and FD)
}

\description{
  This dataset includes tree species abundance counts in two rainforests (Marim and Rebio2)
  in Brazil, collected by Magnago et al. (2014, 2015, 2017). Within each forest,
  there are two assemblages/habitats (Edge and Interior). The data are slightly different
  from those used in Chao et al. (2023) because some species are excluded due to lack of
  phylogenetic information. \cr
}

\usage{data(Brazil_rainforests)}

\format{
  \code{Brazil_rainforests} is a list containing two species-by-assemblage matrices.
  The two columns represent the Edge and Interior habitats. \cr

  A list of 2 matrices: \cr
  $ Marim: a matrix with 119 species (rows) and 2 assemblages (columns),
  \code{"Edge"} and \code{"Interior"}. \cr
  $ Rebio2: a matrix with 118 species (rows) and 2 assemblages (columns),
  \code{"Edge"} and \code{"Interior"}. \cr
}

\references{

  Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M. and Magurran, A. E. (2021).
  Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity
  and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.\cr\cr

  Magnago, L. F. S., Edwards, D. P., Edwards, F. A., Magrach, A., Martins, S. V., & Laurance, W. F. (2014).
  Functional attributes change but functional richness is unchanged after fragmentation of Brazilian Atlantic forests.
  Journal of Ecology, 102, 475-485. \cr\cr

  Magnago, L. F. S., Magrach, A., Barlow, J., Schaefer, C. E. G. R., Laurance, W. F., Martins, S. V., & Edwards, D. P. (2017).
  Do fragment size and edge effects predict carbon stocks in trees and lianas in tropical forests?
  Functional Ecology, 31, 542-552. \cr\cr

  Magnago, L. F. S., Rocha, M. F., Meyer, L., Martins, S. V., & Meira-Neto, J. A. A. (2015).
  Microclimatic conditions at forest edges have significant impacts on vegetation structure in large Atlantic forest fragments.
  Biodiversity and Conservation, 24, 2305-2318. \cr\cr
}

\keyword{datasets}