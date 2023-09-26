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
#' @param x a \code{SpatRaster} object composed of 0, 1 and NA, or the output 
#' list from 
#' \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}} 
#' or \code{\link{generateRandomSp}}
#' @param geographical.limit \code{"country"}, \code{"region"}, 
#' \code{"continent"}, \code{"polygon"}, \code{"raster"} or \code{"extent"}.
#'  The method used
#' to limit the distribution range: see details.
#' @param area \code{NULL}, a character string, a \code{polygon}, a 
#' \code{raster} or an \code{extent} object.
#' The area in which the distribution range will be limited: see details. 
#' If \code{NULL}
#' and \code{geographical.limit = "extent"}, then you will be asked to draw an
#' extent on the map.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the resulting 
#' limited
#' distribution will be plotted.
#' @details
#' \href{http://borisleroy.com/virtualspecies_tutorial/08-dispersallimitation.html}{
#' Online tutorial for this function}
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
#' \item{Countries: type 
#' \code{unique(rnaturalearth::ne_countries(returnclass ='sf')$sovereignt)} 
#' in the console}
#' \item{Regions: "Africa", "Antarctica", "Asia", "Oceania", "Europe", 
#' "Americas"}
#' \item{Continents: "Africa", "Antarctica", "Asia", "Europe", 
#' "North America", "Oceania", "South America"}}
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
#' @import terra
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' # Create an example stack with six environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- c(rast(a * dnorm(1:100, 50, sd = 25)),
#'          rast(a * 1:100),
#'          rast(a * logisticFun(1:100, alpha = 10, beta = 70)),
#'          rast(t(a)),
#'          rast(exp(a)),
#'          rast(log(a)))
#' names(env) <- paste("Var", 1:6, sep = "")   
#' 
#' # More than 6 variables: by default a PCA approach will be used
#' sp <- generateRandomSp(env)
#' 
#' # limiting the distribution to a specific extent
#' limit <- ext(1, 50, 1, 50)
#' 
#' limitDistribution(sp, area = limit)
#' 
#' 
#' # Example of a raster of habitat patches
#' habitat.raster <- setValues(sp$pa.raster, 
#'                             sample(c(0, 1), size = ncell(sp$pa.raster), 
#'                             replace = TRUE))
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
    if(inherits(x$pa.raster, "SpatRaster"))
    {
      sp.raster <- x$pa.raster
      results <- x
    } else stop("x must be:\n- a raster layer object\nor\n- the output list", 
                " from function convertToPA(), generateSpFromFun(), ", 
                "generateSpFromPCA() or generateRandomSp()")
  } else if (inherits(x, "RasterLayer"))
  {
    sp.raster <- rast(x)
    results <- list(wrap(sp.raster))
  } else if (inherits(x, "SpatRaster"))
  {
    sp.raster <- x
    results <- list(wrap(sp.raster))
  } else stop("x must be:\n- a raster layer object\nor\n- the output list", 
              " from function convertToPA(), generateSpFromFun(),", 
              " generateSpFromPCA() or generateRandomSp()")
  
  
  if(length(geographical.limit) > 1)
  {
    stop('Only one dispersal limit method can be applied at a time')
  }
  
  if (!(geographical.limit %in% c("country", "region", "extent", "polygon", 
                                  "continent", "raster")))
  {
    stop('Argument geographical.limit must be one of : country, region,', 
         ' continent, polygon, extent, raster')
  }

  
  if (geographical.limit %in% c("country", "region", "continent"))
  {
    if(!("rnaturalearth" %in% 
         rownames(utils::installed.packages())))
    {
      stop('You need to install the package "rnaturalearth".')
    }
    worldmap <- rnaturalearth::ne_countries(returnclass = "sf")
    
    if(geographical.limit == "country")
    {
      if (any(!(area %in% worldmap$sovereignt)))
      {
        stop("area name(s) must be correctly spelled. Type ", 
             "unique(rnaturalearth::ne_countries(returnclass =",
             "'sf')$sovereignt) to obtain valid names.")
      }    
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    } else if (geographical.limit == "region")
    {
      if (any(!(area %in% worldmap$region_un)))
      {
        stop(paste("region name(s) must be correctly spelled, according to", 
                   " one of the following : ", 
                   paste(unique(worldmap$region_un), collapse = ", "), 
                   sep = "\n"))
      } 
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    } else if (geographical.limit == "continent")
    {
      if (any(!(area %in% worldmap$continent)))
      {
        stop(paste("region name(s) must be correctly spelled,",
                   "according to one of the following : ", 
                   paste(unique(worldmap$continent), collapse = ", "), 
                   sep = "\n"))
      } 
      results$geographical.limit <- list(method = geographical.limit,
                                         area = area)
    }
  } else if (geographical.limit == "polygon") 
  {
    if(is.null(area)) {
      message("No object of class SpatVector or sf provided to area (or wrong",
              "  object class). A window with a map ",
              "will open, click on the map to draw the polygon of the area", 
              ".\n Once finished, press ",
              "escape to close the polygon.")
      if("RStudioGD" %in% names(grDevices::dev.list())) {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      plot(sp.raster)
      area <- draw(x = "polygon")
    } else if(!(inherits(area, c("SpatVector", "sf"))))
    {
      stop("If you choose geographical.limit = 'polygon', please provide a",
           " polygon of class ",
           "sf or SpatVector to argument area. You can also set",
           " area = NULL to draw the polygon manually.")
    }
    results$geographical.limit <- list(method = geographical.limit,
                                       area = area)
  } else if(geographical.limit == "extent")
  {
    if(!(inherits(area, "SpatExtent"))) {
      
      message("No object of class SpatExtent provided (or wrong object class).", 
              " A window with a map ",
              "will open, click on the map to draw the extent of the area", 
              ".\n Once finished, press ",
              "escape to close the polygon.")
      if("RStudioGD" %in% names(grDevices::dev.list())) {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      plot(sp.raster)
      area <- vect(draw())
    } else {
      area <- vect(area)
    }

  } else if(geographical.limit == "raster")
  {
    if(inherits(area, "RasterLayer")) {
      message("area was a raster object, converting to terra...")
      area <- rast(area)
    }
    if (!(inherits(area, "SpatRaster")))
    {
      stop("If you choose to limit the distribution with a raster, please", 
           " provide the raster to argument 'area'")    
    }

    if(global(area, "max", na.rm = TRUE)[1, 1] > 1)
    {
      warning("The raster used to limit species distribution has values", 
              " greater than 1. This will likely result in something very", 
              " different to a distribution raster.")
    }
    if(global(area, "min", na.rm = TRUE)[1, 1] < 0)
    {
      warning("The raster used to limit species distribution has values lower", 
              " than 0. This will likely result in something very different to", 
              " a distribution raster.")
    }
    results$geographical.limit <- list(method = geographical.limit,
                                       area = area)
  }

  geographical.limit.raster <- sp.raster   
  
  if(geographical.limit == "country")
  {
    geographical.limit.raster1 <- rasterize(
      worldmap[which(worldmap$sovereignt %in% area), ],
      geographical.limit.raster, 
      field = 1,
      background = 0,
      silent = TRUE)
    geographical.limit.raster <- 
      geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "region")
  {
    geographical.limit.raster1 <- 
      rasterize(worldmap[which(worldmap$region_un %in% area), ],
                geographical.limit.raster, 
                field = 1,
                background = 0,
                silent = TRUE)
    geographical.limit.raster <-
      geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "continent")
  {
    geographical.limit.raster1 <- rasterize(
      worldmap[which(worldmap$continent %in% area), ],
      geographical.limit.raster,
      field = 1,
      background = 0,
      silent = TRUE)
    geographical.limit.raster <- 
      geographical.limit.raster * geographical.limit.raster1
  } else if(geographical.limit == "extent")
  {
    geographical.limit.raster <-
      geographical.limit.raster * rasterize(area, 
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
    geographical.limit.raster <- 
      geographical.limit.raster * geographical.limit.raster1
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