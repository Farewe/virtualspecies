#' Limit a virtual species distribution to a defined area
#' 
#' @description
#' This function is designed to limit species distributions to a subsample of
#' their total distribution range. It will thus generate a species which is not
#' at the equilibrium with its environment (i.e., which did not occupy the full
#' range of suitable environmental conditions).
#' 
#' This function basically takes any type of raster and will limit values above
#' 0 to areas where the species is allowed to disperse.
#' @param x a \code{rasterLayer} object composed of 0, 1 and NA, or the output list from 
#' \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}} 
#' or \code{\link{generateRandomSp}}
#' @param geographical.limit \code{"country"}, \code{"region"}, 
#' \code{"continent"}, \code{"polygon"}, \code{"raster"} or \code{"extent"}. The method used
#' to limit the distribution range: see details.
#' @param area \code{NULL}, a character string, a \code{polygon}, a \code{raster} or an \code{extent} object.
#' The area in which the distribution range will be limited: see details. If \code{NULL}
#' and \code{geographical.limit = "extent"}, then you will be asked to draw an
#' extent on the map.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the resulting limited
#' distribution will be plotted.
#' @details
#' \href{http://borisleroy.com/virtualspecies_tutorial/08-dispersallimitation.html}{Online tutorial for this function}
#' 
#' 
#' \bold{How the function works:}
#' 
#' The function will remove occurrences of the species outside the chosen area:
#' \itemize{
#' \item{NA are kept unchanged}
#' \item{0 are kept unchanged}
#' \item{values > 0 within the limits of \code{area} are kept unchanged}
#' \item{values > 0 outside the limits of \code{area} are set to 0}
#' }
#' 
#' 
#' \bold{How to define the area in which the range is limited:}
#' 
#' You can choose to limit the distribution range of the species to:
#' \enumerate{
#' \item{a particular country, region or continent (assuming your raster has
#' the WGS84 projection): 
#' 
#' Set the argument
#' \code{geographical.limit} to \code{"country"}, \code{"region"} or
#' \code{"continent"}, and provide the name(s) of the associated countries,
#' regions or continents to \code{area} (see examples). 
#' 
#' List of possible \code{area} names:
#' \itemize{
#' \item{Countries: type \code{levels(getMap()@@data$SOVEREIGNT)} in the console}
#' \item{Regions: "Africa", "Antarctica", "Asia", "Australia", "Europe", 
#' "North America", "South America"}
#' \item{Continents: "Africa", "Antarctica", "Australia", "Eurasia", 
#' "North America", "South America"}}
#' }
#' \item{a polygon:
#' 
#' Set \code{geographical.limit} to \code{"polygon"}, and provide your
#' polygon to \code{area}.
#' }
#' \item{a raster:
#' 
#' Set \code{geographical.limit} to \code{"raster"}, and provide your
#' raster to \code{area}. Your raster values should be 1 (suitable area),
#' 0 (unsuitable area) or NA (outside your mask).
#' }
#' \item{an extent object:
#' 
#' Set \code{geographical.limit} to \code{"extent"}, and either provide your
#' extent object to \code{area}, or leave it \code{NULL} to draw an extent on
#' the map.}
#' } 
#' @return
#' a \code{list} containing 7 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the virtual species distribution, as a Raster object containing the
#' environmental suitability)}
#' \item{\code{PA.conversion}: the parameters used to convert the suitability into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a Raster object containing 0 (absence) / 1 (presence) / NA}
#' \item{\code{geographical.limit}: the method used to
#' limit the distribution and the area in which the distribution is restricted}
#' \item{\code{occupied.area}: the area occupied by the virtual species as a
#' Raster of presence-absence}
#' }
#' The structure of the virtualspecies object can be seen using \code{str()}
#' @export
#' @import raster
#' @importFrom utils installed.packages
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' # Create an example stack with six environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100),
#'              raster(a * logisticFun(1:100, alpha = 10, beta = 70)),
#'              raster(t(a)),
#'              raster(exp(a)),
#'              raster(log(a)))
#' names(env) <- paste("Var", 1:6, sep = "")   
#' 
#' # More than 6 variables: by default a PCA approach will be used
#' sp <- generateRandomSp(env)
#' 
#' # limiting the distribution to a specific extent
#' limit <- extent(0.5, 0.7, 0.6, 0.8)
#' 
#' limitDistribution(sp, area = limit)
#' 
#' 
#' # Example of a raster of habitat patches
#' habitat.raster <- setValues(sp$pa.raster, 
#'                           sample(c(0, 1), size = ncell(sp$pa.raster), replace = TRUE))
#' 
#' plot(habitat.raster) # 1 = suitable habitat; 0 = unsuitable habitat
#' sp <- limitDistribution(sp, geographical.limit = "raster", area = habitat.raster)
#' par(mfrow = c(2, 1))
#' plot(sp$pa.raster)
#' plot(sp$occupied.area) # Species could not occur in many cells because
#' # habitat patches were unsuitable



limitDistribution <- function(x,
                              geographical.limit = "extent",
                              area = NULL,
                              plot = TRUE)
{
  
  
  if(inherits(x, "virtualspecies"))
  {
    if(inherits(x$pa.raster, "RasterLayer"))
    {
      sp.raster <- x$pa.raster
      results <- x
    } else stop("x must be:\n- a raster layer object\nor\n- the output list from function convertToPA(), generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  } else if (inherits(x, "RasterLayer"))
  {
    sp.raster <- x
    results <- list(sp.raster)
  } else stop("x must be:\n- a raster layer object\nor\n- the output list from function convertToPA(), generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  
  
  if(length(geographical.limit) > 1)
  {
    stop('Only one dispersal limit method can be applied at a time')
  }
  
  if (!(geographical.limit %in% c("country", "region", "extent", "polygon", "continent", "raster")))
  {
    stop('Argument geographical.limit must be one of : country, region, continent, polygon, extent, raster')
  }

  
  if (geographical.limit %in% c("country", "region", "continent"))
  {
    if(!("rworldmap" %in% rownames(installed.packages())))
    {
      stop('You need to install the package "rworldmap" in order to use geographical.limit = "region" or "country" or "continent"')
    }
    worldmap <- rworldmap::getMap()
    
    if(geographical.limit == "country")
    {
      if (any(!(area %in% levels(worldmap@data$SOVEREIGNT))))
      {
        stop("area name(s) must be correctly spelled. Type 'levels(getMap()@data$SOVEREIGNT)' to obtain valid names.")
      }    
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    } else if (geographical.limit == "region")
    {
      if (any(!(area %in% levels(worldmap@data$REGION))))
      {
        stop(paste("region name(s) must be correctly spelled, according to one of the following : ", 
                   paste(levels(worldmap@data$REGION), collapse = ", "), sep = "\n"))
      } 
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    } else if (geographical.limit == "continent")
    {
      if (any(!(area %in% levels(worldmap@data$continent))))
      {
        stop(paste("region name(s) must be correctly spelled, according to one of the following : ", 
                   paste(levels(worldmap@data$continent), collapse = ", "), sep = "\n"))
      } 
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    }
  } else if (geographical.limit == "polygon") # Projections are not checked here. Perhaps we should add projection check between raster & polygon in the future?
    # This is especially important given that randomPoints weights samplings by the cell area (because cells closer to
    # the equator are larger)
  {
    if(!(inherits(area, c("SpatialPolygons", "SpatialPolygonsDataFrame"))))
    {
      stop("If you choose geographical.limit = 'polygon', please provide a polygon of class SpatialPolygons or SpatialPolygonsDataFrame to argument area")
    }
    warning("Polygon projection is not checked. Please make sure you have the same projections between your polygon and your presence-absence raster")
    results$geographical.limit <- list(method = geographical.limit,
                                       area = area)
  } else if(geographical.limit == "extent")
  {
    if(!(inherits(area, "Extent")))
    {
      message("No object of class extent (or wrong object) provided: click twice on the map to draw the extent in which presence points will be sampled")
      plot(sp.raster)
      area <- drawExtent(show = TRUE)      
    }
  } else if(geographical.limit == "raster")
  {
    if (!(inherits(area, "RasterLayer")))
    {
      stop("If you choose to limit the distribution with a raster, please provide the raster to argument 'area'")    
    }
    if (!area@data@haveminmax)
    {
      area <- setMinMax(area)
    }
    if(area@data@max > 1)
    {
      warning("The raster used to limit species distribution has values greater than 1. This will likely result in something very different to a distribution raster.")
    }
    if(area@data@min < 0)
    {
      warning("The raster used to limit species distribution has values lower than 0. This will likely result in something very different to a distribution raster.")
    }
    results$geographical.limit <- list(method = geographical.limit,
                                       area = area)
  }

  geographical.limit.raster <- sp.raster   
  
  if(geographical.limit == "country")
  {
    geographical.limit.raster1 <- rasterize(worldmap[which(worldmap@data$SOVEREIGNT %in% area), ],
                                            geographical.limit.raster, 
                                            field = 1,
                                            background = 0,
                                            silent = TRUE)
    geographical.limit.raster <- geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "region")
  {
    geographical.limit.raster1 <- rasterize(worldmap[which(worldmap@data$REGION %in% area), ],
                              geographical.limit.raster, 
                              field = 1,
                              background = 0,
                              silent = TRUE)
    geographical.limit.raster <- geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "continent")
  {
    geographical.limit.raster1 <- rasterize(worldmap[which(worldmap@data$continent %in% area), ],
                                            geographical.limit.raster,
                                            field = 1,
                                            background = 0,
                                            silent = TRUE)
    geographical.limit.raster <- geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "extent")
  {
    geographical.limit.raster <- geographical.limit.raster * rasterize(area, 
                                                                       sp.raster, 
                                                                       field = 1, 
                                                                       background = 0)
    results$geographical.limit <- list(method = geographical.limit,
                                       extent = area)
  } else if(geographical.limit == "polygon")
  {
    geographical.limit.raster1 <- rasterize(area,
                                            geographical.limit.raster, 
                                            field = 1,
                                            background = 0,
                                            silent = TRUE)
    geographical.limit.raster <- geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "raster")
  {
    geographical.limit.raster <- geographical.limit.raster * area
  }
  
  if(plot)
  {
    plot(geographical.limit.raster)
  }
  results$occupied.area <- geographical.limit.raster   
  return(results)
}