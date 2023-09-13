#' Synchronise NA values among layers of a stack
#' 
#' @description
#' This function ensures that cells containing NAs are the same among all the
#' layers of a raster stack, i.e.that for any given pixel of the stack, if one layer has a NA, then 
#' all layers should be set to NA for that pixel.
#' @details
#' This function can do that in two different ways; if your computer has enough RAM a fast way will be
#' used; otherwise a slower but memory-safe way will be used.
#' @param x a raster stack object which needs to be synchronised.
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @importFrom terra mem_info mask app
#' @importFrom utils capture.output
#' @examples
#' # Creation of a stack with different NAs across layers
#' m <- matrix(nr = 10, nc = 10, 1:100)
#' r1 <- rast(m)
#' r2 <- rast(m)
#' r1[sample(1:ncell(r1), 20)] <- NA
#' r2[sample(1:ncell(r2), 20)] <- NA
#' s <- c(r1, r2)
#' 
#' 
#' # Effect of the synchroniseNA() function
#' plot(s) # Not yet synchronised
#' s <- synchroniseNA(s)
#' plot(s) # Synchronised
#' 
synchroniseNA <- function(x)
{
  if(inherits(x, "SpatRaster")) {
    # The gain in performance is only about 7% so this code is currently
    # commented.
    # opt <- terra:::spatOptions()
    # opt$ncopies = 2
    # mem <- x@ptr$mem_needs(opt)
    # 
    # if(mem[1] < prod(mem[2:3])) {
    #   val <- values(x)
    #   if(any(is.na(val)))
    #   {
    #     NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
    #   }
    #   val[NA.pos, ] <- NA
    #   x <- setValues(x, val)
    #   return(x)
    # } else {
      x <- terra::mask(x, terra::app(x, fun = sum))
      return(x)
    # }
  } else if(inherits(x, "RasterStack")){
    if(canProcessInMemory(x, n = 2))
    {
      val <- getValues(x)
      if(any(is.na(val)))
      {
        NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
      }
      val[NA.pos, ] <- NA
      x <- setValues(x, val)
      return(x)
    } else
    {
      x <- mask(x, calc(x, fun = sum))
      return(x)
    }
  }
}