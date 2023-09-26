[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/virtualspecies)](https://cran.r-project.org/package=virtualspecies)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/virtualspecies)](http://r-pkg.org/pkg/virtualspecies)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/virtualspecies)](http://r-pkg.org/pkg/virtualspecies)
[![Research software impact](http://depsy.org/api/package/cran/virtualspecies/badge.svg)](http://depsy.org/package/r/virtualspecies)

==============

virtualspecies
==============

The virtualspecies R package


This package provides a user-friendly framework for generating virtual species 
distribution, integrating the existing methodological approaches.
The package includes (1) generating virtual species’ suitability from a spatial 
set of environmental conditions, with two different approaches; (2) converting 
the environmental suitability into presence-absence with a probabilistic 
approach; (3) introducing dispersal limitations in the realised virtual species 
distributions and (4) sampling occurrences with different biases in the sampling 
procedure.

There is a full online tutorial available here:
http://borisleroy.com/virtualspecies/


# Update 1.6 : shift from `raster` to `terra`

`virtualspecies` has been updated to version 1.6 in September 2023. The major
change in this update is that now this package relies upon the R package 
`terra` rather than `raster`; in addition, internally it uses functions from
the package `sf` rather than `sp`. This update was required since multiple 
spatial packages are going to be retired from R, and future developments and
improvements of spatial packages will rely on `terra` and `sf`. 

**This update may cause some breaking changes**. I maintained the changes as 
minimal as possible so that previous code (version <= 1.5) will still be
compatible with the new code (version >= 1.6). However, objects generated
in versions <= 1.5 will likely not be compatible with the package version 1.6.
They can still be converted manually to become compatible, using functions 
from the package `terra` such as `rast()`. If you find yourself
in trouble, email me about that. In addition, the methods in functions
`limitDistribution()` and `sampleOccurrences()` have changed slightly because
they now use the package `rnaturaleath` (formerly they used `rworldmap`).
