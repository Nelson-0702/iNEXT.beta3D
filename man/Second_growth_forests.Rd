\name{Second_growth_forests}
\alias{Second_growth_forests}
\docType{data}
\title{
  Second-growth forests data (datatype = "incidence_raw" for taxonomic diversity)
}
\description{
  This dataset includes tree incidence data in 100 subplots (each with 0.01 ha) collected from two 
  second-growth forests, namely Cuatro Rios (CR) and Juan Enriquez (JE) in Costa Rica.
  Each 1-ha forest was divided into 100 subplots (each with 0.01 ha) and only species' incidence records 
  in each subplot were used to compute the incidence frequency for a species 
  (i.e., the number of subplots in which that species occurred); see Chazdon et al. (2021, 2022)
  for sampling details and pertinent analyses. The original time series data cover 2005 to 2017; here only data from three years
  (2005, 2011, and 2017) are used for assessing temporal beta diversity among the three years within each forest. 
  

}
\usage{data(Second_growth_forests)}
\format{
  Second_growth_forests is a list with two forests. The input format for each forest is a list with three species-by-sampling-unit matrices ("Year_2005", "Year_2011", and "Year_2017"). Each matrix records the species as 0 (non-detection) or 1 (detection) in each sampling unit. \cr
  A list of 2  \cr
  $ CR 2005 vs. 2011 vs. 2017: A list of 3 \cr
  .. ..$ Year_2005 (152 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2011 (152 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2017 (152 (species) x 100 (quadrats)) \cr 
  $ JE 2005 vs. 2011 vs. 2017: A list of 3 \cr
  .. ..$ Year_2005 (108 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2011 (108 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2017 (108 (species) x 100 (quadrats)) \cr 
}
\references{
  Chazdon, R. (2021). Tree abundance in eight 1-ha tropical forest plots in northeastern Costa Rica from 1997-2017, https://doi.org/10.5061/dryad.ncjsxksvr \cr\cr
  Chazdon, R. L., N. Norden, R. K. Colwell, and A. Chao. (2022). Monitoring recovery of tree diversity during tropical forest restoration: lessons from long-term trajectories of natural regeneration. Philosophical Transactions of the Royal Society B, 378: 20210069. \cr\cr
}

\keyword{datasets}
