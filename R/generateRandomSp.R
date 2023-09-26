#' Generate a random virtual species distribution from environmental variables
#' 
#' @description
#' This function generates randomly a virtual species distribution.
#' 
#' @param raster.stack a SpatRaster object, in which each layer represent an
#'  environmental 
#' variable.
#' @param approach \code{"automatic"}, \code{"random"}, \code{"response"}
#' or \code{"pca"}. This parameters defines how species will be generated. 
#' \code{"automatic"}: If less than 6 variables in \code{raster.stack}, a 
#' response approach will be used, otherwise a pca approach will be used.
#' \code{"random"}: the approach will be randomly picked. Otherwise choose
#' \code{"response"} or \code{"pca"}. See details.
#' @param rescale \code{TRUE} or \code{FALSE}. If \code{TRUE}, the final 
#' probability of presence is rescaled between 0 and 1.
#' @param convert.to.PA \code{TRUE} or \code{FALSE}. If \code{TRUE}, the 
#' virtual species distribution will also be converted into Presence-Absence.
#' @param relations [response approach] a vector containing the possible types 
#' of response function.
#' The implemented type of relations are \code{"gaussian"}, \code{"linear"},
#' \code{"logistic"} and \code{"quadratic"}.
#' @param rescale.each.response \code{TRUE} or \code{FALSE}. If \code{TRUE}, 
#' the individual responses to
#' each environmental variable are rescaled between 0 and 1
#' @param realistic.sp [response approach] \code{TRUE} or \code{FALSE}. If 
#' \code{TRUE}, the function will try to define responses that can form a viable
#' species. If \code{FALSE}, the responses will be randomly generated
#' (may result in environmental conditions that do not co-exist).
#' @param species.type [response approach] \code{"additive"} or 
#' \code{"multiplicative"}. Defines 
#' how the final probability of presence is calculated: if \code{"additive"},
#'  responses to each
#' variable are summed; if \code{"multiplicative"}, responses are multiplied.
#' See \code{\link{generateSpFromFun}}
#' @param niche.breadth [pca approach] \code{"any"}, \code{"narrow"} or 
#' \code{"wide"}. This parameter
#' defines how tolerant is the species regarding environmental conditions by 
#' adjusting
#' the standard deviations of the gaussian functions. See 
#' \code{\link{generateSpFromPCA}}
#' @param sample.points [pca approach] \code{TRUE} of \code{FALSE}. If you have 
#' a large
#' raster file then use this parameter to sample a number of points equal to
#' \code{nb.points}.
#' @param nb.points [pca approach] a numeric value. Only useful if 
#' \code{sample.points = TRUE}.
#' The number of sampled points from the raster, to perform the PCA. A too small
#' value may not be representative of the environmental conditions in your 
#' raster.
#' @param PA.method \code{"threshold"} or \code{"probability"}. If 
#' \code{"threshold"}, then occurrence probabilities are simply converted into
#' presence-absence according to the threshold \code{beta}. If 
#' \code{"probability"}, then
#' probabilities are converted according to a logistic function of threshold 
#' \code{beta} and slope \code{alpha}.
#' @param beta \code{"random"}, a numeric value in the range of your 
#' probabilities or \code{NULL}. This is the threshold of conversion into
#' presence-absence (= the inflexion point if \code{PA.method = "probability"}).
#' If \code{"random"}, a numeric value will be randomly generated within 
#' the range
#' of probabilities of occurrence. See \code{\link{convertToPA}}
#' @param alpha \code{NULL} or a negative numeric value. Only useful if 
#' \code{PA.method = "probability"}. The value of \code{alpha} will
#' shape the logistic function transforming occurrences into presence-absences.
#' See \code{\link{logisticFun}} and examples therein for the choice of 
#' \code{alpha}
#' @param adjust.alpha \code{TRUE} or \code{FALSE}. Only useful if 
#' \code{rescale = FALSE}. If  \code{adjust.alpha = TRUE}, then the value 
#' of \code{alpha} will
#' be adjusted to an appropriate value  for the range of suitabilities.  
#' @param species.prevalence \code{NULL} or a numeric value between 0 and 1.
#' The species prevalence is the proportion of sites actually occupied by the
#' species. See \code{\link{convertToPA}}
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the generated 
#' virtual species will be plotted.
#' @details
#' \href{http://borisleroy.com/virtualspecies_tutorial/05-randomspecies.html}{
#' Online tutorial for this function}
#' 
#' 
#' 
#' This function generate random virtual species, either using a PCA 
#' approach, or using
#' a response approach. In case of a response approach, only four response 
#' functions are
#' currently used: gaussian, linear, logistic and quadratic functions.
#' 
#' Note that in case of numerous predictor variables, the "response" 
#' approach will
#' not work well because it will often generate contradicting response functions 
#' (e.g., mean annual temperature optimum at degrees C and temperature 
#' of the coldest month at
#' 10 degrees C). In these case, it is advised to use the PCA approach 
#' (by default, a PCA approach
#' will be used if there are more than 6 predictor variables).
#' 
#' If \code{rescale.each.response = TRUE}, then the probability response to each
#' variable will be normalised between 0 and 1 according to the following 
#' formula:
#' P.rescaled = (P - min(P)) / (max(P) - min (P)). Similarly, if 
#' \code{rescale = TRUE},
#' the final environmental suitability will be rescaled between 0 and 1 
#' with the same formula.
#' 
#' By default, the function will perform a probabilistic conversion into 
#' presence-
#' absence, with a randomly chosen beta threshold. If you want to customise the 
#' conversion parameters, you have to define \bold{two} of the three 
#' following parameters:
#' \itemize{
#' \item{\code{beta}: the 'threshold' of the logistic function (i.e. the 
#' inflexion point)}
#' \item{\code{alpha}: the slope of the logistic function}
#' \item{\code{species.prevalence}: the proportion of sites in which the species
#' occur}
#' }
#' 
#' If you provide \code{beta} and \code{alpha}, the \code{species.prevalence}
#' is calculated immediately calculated after conversion into presence-absence.
#' 
#' As explained in \code{\link{convertToPA}}, if you choose choose a precise
#' \code{species.prevalence}, it may not be possible to reach this particular 
#' value because of the availability of environmental conditions. Several
#' runs may be necessary to reach the desired \code{species.prevalence}.
#' @import terra
#' @export
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @return a \code{list} with 3 to 5 elements (depending if the conversion 
#' to presence-absence was performed):
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, 
#' \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the 
#' species}
#' \item{\code{suitab.raster}: the virtual species distribution, as a 
#' SpatRaster object containing the
#' environmental suitability)}
#' \item{\code{PA.conversion}: the parameters used to convert the suitability 
#' into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a SpatRaster object 
#' containing 0 (absence) / 1 (presence) / NA}
#' }
#' The structure of the virtualspecies can object be seen using \code{str()}
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
#' generateRandomSp(env)
#' 
#' # Manually choosing a response approach: this may fail because it is hard
#' # to find a realistic species with six distinct responses to six variables
#' \donttest{
#' generateRandomSp(env, approach = "response")
#' }
#' 
#' # Randomly choosing the approach
#' generateRandomSp(env, approach = "random")
#' 
#' 

generateRandomSp <- function(raster.stack, 
                             approach = "automatic",
                             rescale = TRUE,
                             convert.to.PA = TRUE,
                             relations = c("gaussian", "linear", "logistic", "quadratic"),
                             rescale.each.response = TRUE,
                             realistic.sp = TRUE,
                             species.type = "multiplicative",
                             niche.breadth = "any",
                             sample.points = FALSE, 
                             nb.points = 10000,
                             PA.method = "probability",
                             alpha = -.1,
                             adjust.alpha = TRUE,
                             beta = "random",
                             species.prevalence = NULL,
                             plot = TRUE)
{
  if(inherits(raster.stack, "Raster")) {
    raster.stack <- rast(raster.stack)
  }
  if(!(inherits(raster.stack, "SpatRaster")))
  {
    stop("raster.stack must be a SpatRaster object")
  }

  
  if(approach == "automatic")
  {
    if(nlyr(raster.stack) <= 5)
    {
      approach <- "response"
    } else
    {
      approach <- "pca"
    }
  } else if (approach == "random")
  {
    approach <- sample(c("response", "pca"), 1)
  } else if(!(approach %in% c("response", "pca")))
  {
    stop("Argument approach was misspecified. Either choose 'automatic', ",
         "'random', 'response' or 'pca'.")
  }
  
  var.names <- names(raster.stack)
  
  if(approach == "pca")
  {
    results <- generateSpFromPCA(raster.stack,
                                 rescale = rescale,
                                 niche.breadth = niche.breadth,
                                 sample.points = sample.points, 
                                 nb.points = nb.points,
                                 plot = FALSE)
  } else if (approach == "response")
  {
    parameters <- list()
    message(" - Determining species' response to predictor variables\n")
    if(any(!(relations %in% c("gaussian", "linear", "logistic", "quadratic"))))
    {
      stop(paste("Wrong relation type specified: pick among '", 
                 paste(c("gaussian", "linear", "logistic", "quadratic"), 
                       collapse = " "), "'",
                 collapse = " "))
    }
    valid.cells <- setValues(raster.stack[[1]], 1)
    var.order <- sample(var.names, length(var.names), replace = F)
    for (i in 1:length(var.order))
    {
      
      cur.var <- var.order[i]
      cur.rast <- raster.stack[[cur.var]]
      if(realistic.sp) cur.rast <- cur.rast * valid.cells # Cur.rast is 
      # here restricted to current suitable conds
      
      type <- sample(relations, 1)
      
      min_ <- global(cur.rast, "min", na.rm = TRUE)[1, 1]
      max_ <- global(cur.rast, "max", na.rm = TRUE)[1, 1]
      
      
      if (type == "gaussian")
      {
        parameters[[cur.var]] <- list(
          fun = 'dnorm',
          args = c(mean = sample(seq(min_,
                                     max_, 
                                     length = 100000), 1),
                   sd = sample(seq(0, 
                                   (max_ - min_), 
                                   length = 100000), 1))
        )
      } else if (type == "linear")
      { # At the moment this is not really useful because the rescale will transforme the results in either 0:1 or 1:0, regardless of the slope
        # To be improved later
        parameters[[cur.var]] <- list(
          fun = 'linearFun',
          args = c(a = sample(seq(-1, 1, length = 100), 1),
                   b = sample(seq(min_, 
                                  max_, 
                                  length = 100000), 1))
        )
      } else if (type == "logistic")
      {
        beta.t <- sample(seq(min_,
                             max_,
                             length = 1000000), 1)
        alpha.t <-  sample(c(seq((max_ - min_)/1000,
                                 (max_ - min_)/100, length = 10),
                             seq((max_ - min_)/100,
                                 (max_ - min_)/10, length = 100),
                             seq((max_ - min_)/10,
                                 (max_ - min_)*10, length = 10)), size = 1)
        if(realistic.sp == TRUE)
        {
          if(beta.t > max_)
          {
            alpha.t <- alpha.t
          } else if (beta.t < min_)
          {
            alpha.t <- -alpha.t
          } else
          {
            alpha.t <- sample(c(alpha.t, -alpha.t), 1)
          }
        }
        
        parameters[[cur.var]] <- list(fun = 'logisticFun',
                                      args = c(alpha = alpha.t,
                                               beta = beta.t)
        )
      } else if (type == "quadratic")
      {
        max.point <- sample(seq(min_,
                                max_, 
                                length = 1000), 1)
        a <- sample(seq(-.01, -20, length = 10000), 1)
        b <- - max.point * 2 * a
        parameters[[cur.var]] <- list(fun = 'quadraticFun',
                                      args = c(a = a,
                                               b = b,
                                               c = 0)
        )
        
      }
      
      # Restricting values to suitable conditions
      tmp.rast <- app(raster.stack[[cur.var]], fun = function(x)
      {
        do.call(match.fun(parameters[[cur.var]]$fun), 
                args = c(list(x), parameters[[cur.var]]$args))
      }
      )
      tmp.rast <- (tmp.rast - global(tmp.rast, "min", na.rm = TRUE)[1, 1]) /
        (global(tmp.rast, "max", na.rm = TRUE)[1, 1] - 
           global(tmp.rast, "min", na.rm = TRUE)[1, 1])
      valid.cells <- valid.cells * (tmp.rast > 0.05)
    }
    message(" - Calculating species suitability\n")
    results <- generateSpFromFun(raster.stack, 
                                 parameters, 
                                 rescale = rescale,
                                 species.type = species.type, 
                                 plot = FALSE,
                                 rescale.each.response = rescale.each.response)
  }
  
  
  
  
  
  if(convert.to.PA == TRUE) {
    message(" - Converting into Presence - Absence\n")
    
    # Need to adjust alpha to appropriate scale if rescale = FALSE
    if(rescale == FALSE) {
      if(adjust.alpha)
      {
        alpha <- diff(c(global(results$suitab.raster,
                               min, na.rm = TRUE)[1, 1],
                        global(results$suitab.raster,
                               max, na.rm = TRUE)[1, 1])) * alpha
      }
      results <- convertToPA(results, 
                             PA.method = PA.method,
                             alpha = alpha,
                             beta = beta,
                             species.prevalence = species.prevalence,
                             plot = FALSE)
    
      if(plot) plot(results)
    } else {
      
      results <- convertToPA(results, 
                             PA.method = PA.method,
                             alpha = alpha,
                             beta = beta,
                             species.prevalence = species.prevalence,
                             plot = FALSE)
      
      if(plot) plot(results)
    }
  } else {
    if(plot) plot(results)
  }
  
  return(results)
}
