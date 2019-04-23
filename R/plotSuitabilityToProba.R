#' Visualise the function that was used to transform environmental suitability into
#' probability of occurrence 
#' 
#' This function plots the relationships between the environmental suitability
#' and the probability of occurrence, which is used to generate the presence-
#' absence distribution. 
#' It requires the output from \code{\link{convertToPA}}.
#' 
#' @param sp the output from \code{\link{convertToPA}}.
#' @param add \code{TRUE} or \code{FALSE}. If \code{TRUE}, the relationship
#' will be added to the currently active graph.
#' @param ... further arguments to be passed to \code{plot}. See 
#' \code{\link[graphics]{plot}} and \code{\link[graphics]{par}}.
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100))
#' names(env) <- c("variable1", "variable2")
#'
#' parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 1e-04, 
#'                                              sd = 1e-04),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' # Generation of the virtual species
#' sp1 <- generateSpFromFun(env, parameters)
#' sp1
#' 
#' 
#' # Converting to presence-absence, probablistic method, logistic conversion
#' # A species with a low prevalence:
#' 
#' sp1.lowprev <- convertToPA(sp1, species.prevalence = 0.1)
#' plotSuitabilityToProba(sp1.lowprev)
#' 
#' # A species with a high prevalence:
#' 
#' sp1.highprev <- convertToPA(sp1, species.prevalence = 0.9)
#' plotSuitabilityToProba(sp1.lowprev)
#' 
#' # Converting to presence-absence, probablistic method, linear conversion
#' # A species with a low prevalence:
#' 
#' sp1.lowprev <- convertToPA(sp1, species.prevalence = 0.1,
#'                            prob.method = "linear")
#' plotSuitabilityToProba(sp1.highprev)
#' 
#' # A species with a high prevalence:
#' 
#' sp1.highprev <- convertToPA(sp1, species.prevalence = 0.9,
#'                            prob.method = "linear")
#' plotSuitabilityToProba(sp1.highprev)
#' 

plotSuitabilityToProba <- function(sp, add = FALSE, ...)
{
  if("virtualspecies" %in% class(sp))
  {
    if(!("PA.conversion" %in% names(sp)))
    {
      stop("sp does not seem to be a valid object: provide the output from convertToPA()")
    }
  } else
  {
    stop("sp does not seem to be a valid object: provide the output from convertToPA()")
  }
  method <- sp$PA.conversion[1]
  
  x <- seq(minValue(sp$suitab.raster),
           maxValue(sp$suitab.raster), length = 1000)
  
  if(method == "probability")
  {
    if(sp$PA.conversion[2] == "logistic")
    {
      y <- logisticFun(x, alpha = as.numeric(sp$PA.conversion["alpha"]),
                       beta =  as.numeric(sp$PA.conversion["beta"]))
    } else
    {
      y <- as.numeric(sp$PA.conversion["a"]) * x + 
        as.numeric(sp$PA.conversion["b"])
      
      y[y < 0] <- 0
      y[y > 1] <- 1
    }
  } else if (method == "threshold")
  {
    y <- x
    y[y >= as.numeric(sp$PA.conversion["cutoff"])] <- 1
    y[y < as.numeric(sp$PA.conversion["cutoff"])] <- 0
  }

  if(!add)
  {
    plot(x = x, y = y, type = "l", bty = "l", las = 1, cex.axis = .8, 
         ylim = c(0, 1), xlab = "Environmental suitability",
         ylab = "Probability of ocurrence",
         ...)
  } else
  {
    lines(x = x, y = y)
  }
}
