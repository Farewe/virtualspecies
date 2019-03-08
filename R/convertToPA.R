#' Convert a virtual species distribution (or a suitability raster) into presence-absence
#' 
#' This functions converts the probabilities of presence from the output of
#'  \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}}, \code{\link{generateRandomSp}}
#' or a suitability raster into
#' a presence-absence raster. The conversion can be threshold-based, or based
#' on a probability of conversion (see details).
#' 
#' @param x a suitability raster, or the output from functions 
#' \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}} 
#' or \code{\link{generateRandomSp}}
#' @param PA.method \code{"threshold"} or \code{"probability"}. If 
#' \code{"threshold"}, then occurrence probabilities are simply converted into
#' presence-absence according to the threshold \code{beta}. If \code{"probability"}, then
#' probabilities are converted according to a logistic function of threshold 
#' \code{beta} and slope \code{alpha}.
#' @param proba.method \code{"logistic"} or \code{"linear"}. Defines how 
#' the initial environmental suitability is translated into probabilities of 
#' presence/absence.
#' @param beta \code{"random"}, a numeric value in the range of your 
#' probabilities or \code{NULL}. This is the threshold of conversion into 
#' presence-absence (if \code{PA.method = "probability"} and 
#' \code{proba.method = "logistic"}, then beta = the inflexion point).
#' If \code{"random"}, a numeric value will be randomly generated within the range
#' of \code{x}.
#' @param alpha \code{NULL} or a negative numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "logistic"}. 
#' The value of \code{alpha} will
#' shape the logistic function transforming occurrences into presence-absences.
#' See \code{\link{logisticFun}} and examples therein for the choice of \code{alpha}
#' @param a \code{NULL} or a  numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "linear"}. 
#' Slope of the linear conversion of environmental suitability.
#' @param b \code{NULL} or a  numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "linear"}. 
#' Intercept of the linear conversion of environmental suitability.
#' @param species.prevalence \code{NULL} or a numeric value between 0 and 1.
#' The species prevalence is the proportion of sites actually occupied by the
#' species.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, maps of probabilities
#' of occurrence and presence-absence will be plotted.
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @references
#' Meynard C.N. & Kaplan D.M. 2013. Using virtual species to study species 
#' distributions and model performance. 
#' \emph{Journal of Biogeography} \bold{40}:1-8
#' 
#' Meynard C.N. & Kaplan D.M. 2011. The effect of a gradual response to the 
#' environment on species distribution model performance.
#' \emph{Ecography} \bold{35}:499-509
#' 
#' @details 
#' The conversion of probabilities of occurrence into presence - absence is
#' usually performed by selecting a threshold above which presence always occurs,
#' and never below. However, this approach may be too much unrealistic because
#' species may sometime be present in areas with a low probability of occurrence,
#' or be absent from areas with a high probability of occurrence. In addition,
#' when using a threshold you erase the previously generated response shapes: 
#' it all becomes threshold. Thus, this threshold approach should be avoided.
#' 
#'  
#' Hence, a more
#' realistic conversion consists in converting probabilities into presence -
#' absence with a probability function (see references). Such a probability 
#' conversion can be performed with two different methods here:
#' \enumerate{
#' \item{Using a logistic transformation of  environmental suitability
#' (see \code{\link{logisticFun}}). This method is designed to simulate 
#' a sharp transition from highly suitable sites to poorly suitable sites. This
#' method alters the shape of the environmental suitability.}
#' \item{Using a linear transformation of environmental suitability. This method
#' will use environmental suitability as the probability of drawing presences 
#' in a Bernoulli trial. This method does not alter significantly the shape of
#' the environmental suitability, but it may alter suitability values
#' depending on the chosen parameters.}}
#' 
#' --- note ---
#' 
#' Regardless of the method you choose to generate your presence-absence
#' distribution, the initial environmental suitability will likely be altered:
#' \enumerate{
#' \item{Threshold method: in this case, the suitability response to the 
#' environment is completely lost, which is why it is not recommended.}
#' \item{Probabilistic method, logistic conversion: the shape of the
#' functional response to the environment is altered}
#' \item{Probabilistic method, linear conversion: if the coefficients of the
#' linear regression are not slope = 1 and intercept = 0, then the range of
#' values of the suitability distribution is altered}}
#' 
#' Therefore, the function stores the modified environmental suitability in the 
#' output list, under the object \code{suitab.raster}. The original suitability 
#' (before logistic or linear conversion) will be stored in the output list as 
#' \code{suitab.before.PA.conversion}. 
#' 
#' --------------------------------------------------------------------------
#' 
#' \bold{PROBABILISTC CONVERSION - LOGISTIC METHOD}
#' 
#' To perform the logistic transformation of environmental suitability
#' you have to define two of the
#' three following parameters:
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
#' On the other hand, if you provide \code{species.prevalence} and either
#' \code{beta} or \code{alpha}, the function will try to determine \code{alpha}
#' (if you provided \code{beta}) or \code{beta} (if you provided \code{alpha}).
#' 
#' The relationship between species prevalence, alpha and beta is dependent
#' on the available range of environmental conditions (see Meynard and Kaplan,
#' 2011 and especially the Supporting Information). As a consequence, the 
#' desired species prevalence may not be available for the defined \code{alpha} 
#' or \code{beta}. In these conditions, the function will retain the \code{alpha} or
#' \code{beta} which provides the closest prevalence to your \code{species.prevalence},
#' but you may also provide another value of \code{alpha} or \code{beta} to obtain
#' a closer prevalence. 
#'  
#' --------------------------------------------------------------------------
#' 
#' \bold{PROBABILISTIC CONVERSION - LINEAR METHOD }
#' 
#' To perform the linear transformation of environmental suitability
#' you have to define *one* of the following:
#' \itemize{
#' \item{nothing - in which case your input environmental suitability will be
#' used as the probability of occurrence for the Bernoulli trial (it is 
#' equivalent to defining a slope \code{a} of 1 and intercept \code{b} of  0.)}
#' \item{the coefficients of the linear regression: slope \code{a} and 
#' intercept \code{b}. The transformed environmental suitability will
#' be used as the probability of occurrence for the Bernoulli trial.}
#' \item{\code{species.prevalence}: the proportion of sites in which the 
#' species occur. In this case, the function will try to find coefficients
#' of a linear regression which results in the requested 
#' \code{species.prevalence}.}
#' } 
#' 
#' Using a linear transformation may result in suitabilities above 1 or below 0,
#' which is not possible. Therefore, all suitabilities above 1 are converted to 
#' 1, and all suitabilities below 0 are converted to 0 before running the 
#' Bernoulli trial. 
#' 
#' 
#'  --------------------------------------------------------------------------
#'  
#' In all cases, the \code{species.prevalence} indicated in the output is the
#' prevalence measured on the output presence-absence map.
#' 
#' @note
#' The approximation of \code{alpha} or \code{beta} to the chosen 
#' \code{species.prevalence} may take time if you work on very large rasters.
#' @return
#' a \code{list} containing 5 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the virtual species distribution, as a Raster object containing the
#' environmental suitability)}
#' \item{\code{PA.conversion}: the parameters used to convert the suitability into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a Raster object containing 0 (absence) / 1 (presence) / NA}
#' }
#' The structure of the virtualspecies object can be seen using \code{str()}
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100))
#' names(env) <- c("variable1", "variable2")
#' 
#' # Creation of the parameter list
#' parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 0.00012,
#'                                             sd = 0.0001),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' sp1 <- generateSpFromFun(env, parameters, plot = FALSE)
#' 
#' # Conversion into presence-absence with a threshold-based approach
#' convertToPA(sp1, PA.method = "threshold", beta = 0.2,  plot = TRUE)
#' convertToPA(sp1, PA.method = "threshold", beta = "random", plot = TRUE)
#' 
#' # Conversion into presence-absence with a probability approach using logistic
#' # method
#' convertToPA(sp1, PA.method = "probability", beta = 0.4, 
#'               alpha = -0.05, plot = TRUE)
#' convertToPA(sp1, PA.method = "probability", beta = "random", 
#'               alpha = -0.1, plot = TRUE)
#'               
#' # Conversion into presence-absence with a probability approach using linear 
#' # method
#' convertToPA(sp1, PA.method = "probability", prob.method = "linear",
#'             a = 1, b = 0, plot = TRUE)         
#'               
#'               
#' # Conversion into presence-absence by choosing the prevalence
#' # Threshold method
#' convertToPA(sp1, PA.method = "threshold",
#'               species.prevalence = 0.3, plot = TRUE)
#' # Logistic method, with alpha provided              
#' convertToPA(sp1, PA.method = "probability", alpha = -0.1, 
#'               species.prevalence = 0.2, plot = TRUE)        
#' # Logistic method, with beta provided              
#' convertToPA(sp1, PA.method = "probability", beta = 0.5, 
#'               species.prevalence = 0.2, alpha = NULL, 
#'               plot = TRUE)    
#' # Linear method
#' convertToPA(sp1, PA.method = "probability", prob.method = "linear",
#'             species.prevalence = 0.2,
#'             plot = TRUE)              
#' convertToPA(sp1, PA.method = "probability", prob.method = "linear",
#'             species.prevalence = 0.5,
#'             plot = TRUE) 
#' convertToPA(sp1, PA.method = "probability", prob.method = "linear",
#'             species.prevalence = 0.8,
#'             plot = TRUE)                
#'  
#' # Plot the output Presence-Absence raster only             
#' sp1 <- convertToPA(sp1, plot = FALSE)
#' plot(sp1$pa.raster)
                    

convertToPA <- function(x, 
                        PA.method = "probability",
                        prob.method = "logistic",
                        beta = "random",
                        alpha = -.05,
                        a = NULL,
                        b = NULL,
                        species.prevalence = NULL,
                        plot = TRUE)

{
  if("virtualspecies" %in% class(x))
  {
    if("RasterLayer" %in% class(x$suitab.raster))
    {
      sp.raster <- x$suitab.raster
    } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions
               generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  } else if ("RasterLayer" %in% class(x))
  {
    sp.raster <- x
  } else stop("x must be:\n- a raster layer object\nor\n- the output list from functions
               generateSpFromFun(), generateSpFromPCA() or generateRandomSp()")
  
  if(any(is.na(maxValue(sp.raster))))
  {
    sp.raster <- setMinMax(sp.raster)
  }
  
  if(PA.method == "threshold")
  {
    if(is.numeric(beta))
    {
      if(is.numeric(species.prevalence))
      {
        warning("Both beta and species.prevalence were provided. beta will be
                ignored.")
        beta <- NULL
      } else if(beta < sp.raster@data@min) 
      {
        warning("beta is lower than all values in your suitability raster. The
                species will most likely be present everywhere")
      } else if(beta > sp.raster@data@max)
      {
        warning("beta is higher than all values in your suitability raster. 
                The species will most likely be absent everywhere")
      }
    } else if(beta == "random")
    {
      if(is.numeric(species.prevalence))
      {
        beta <- NULL 
      } else
      {
        beta <- sample(seq(sp.raster@data@min, 
                           sp.raster@data@max, length = 1000), 1)
      }
    } else if(is.null(beta))
    {
      if(is.null(species.prevalence))
      {
        stop("Either provide beta or species.prevalence when choosing
             PA.method = 'threshold'")
      }
      } else
      {
        stop("beta must either be 'random', a numeric value (preferably within the range of
             your data or NULL")
      }
      } else if(PA.method == "probability")
      {
        if(prob.method == "logistic")
        {
          if(length(c(alpha, beta, species.prevalence)) <= 1)
          {
            if(!is.null(species.prevalence))
            {
              warning("Neither alpha nor beta were provided. As a consequence, alpha
                      will be determined to a random value, and beta will be adjusted 
                      automatically to the desired species prevalence.")
              alpha <- -sample(c(seq((sp.raster@data@max - sp.raster@data@min)/1000,
                                     (sp.raster@data@max - sp.raster@data@min)/100, length = 10),
                                 seq((sp.raster@data@max - sp.raster@data@min)/100,
                                     (sp.raster@data@max - sp.raster@data@min)/10, length = 100),
                                 seq((sp.raster@data@max - sp.raster@data@min)/10,
                                     (sp.raster@data@max - sp.raster@data@min)*10, length = 10)), size = 1)
            } else
            {
              stop("If you choose PA.method = 'probability', you must provide two of the
                   three following parameters: beta, alpha and species.prevalence.")
            }
            } else if(length(c(alpha, beta, species.prevalence)) > 2)
          {
              if(beta != "random")
              {
                stop("You should not provide the three parameters beta, alpha and 
                     species.prevalence. Set beta to 'random' if you want to
                     specify species.prevalence.")
              }
              beta <- NULL
              } 
          # Checking the arguments. At this stage only two of them should be not NULL
          if(!is.null(beta))
          {
            if(is.numeric(beta))
            {
              if(beta < sp.raster@data@min) 
              {
                warning("beta is lower than all values in your suitability raster. The
                        species will most likely be present everywhere")
              } else if(beta > sp.raster@data@max)
              {
                warning("beta is higher than all values in your suitability raster. 
                        The species will most likely be absent everywhere")
              }
            } else if(beta == "random")
            {
                beta <- sample(seq(sp.raster@data@min, 
                                   sp.raster@data@max, length = 1000), 1)
            } else
            {
                stop("beta must either be 'random', a numeric value (preferably within the range of
                     your data) or NULL")
            }
          }
          
          if(!is.null(species.prevalence))
          {
            if(is.numeric(species.prevalence))
            {
              if(!(species.prevalence >= 0 & 
                   species.prevalence <= 1))
              {
                stop("species.prevalence must be a numeric value between 0 and 1.")
              }
            } else 
            {
              stop("species.prevalence must either be a numeric value between 0 and 1
                   or NULL")
            }
          }
          
          if(!is.null(alpha))
          {
            if(!is.numeric(alpha))
            {
              stop("Please provide a numeric value to alpha")
            } else if(alpha > 0)
            {
              warning("alpha was provided > 0. 
                      This means that low probabilities will be converted to presences, 
                      and high probabilities will be converted to absences.
                      If this is not what was intended, provide a negative alpha.")
            }
            }
          
          if(!is.null(species.prevalence))
          {
            if(!is.null(beta))
            {
              message("   --- Determing alpha automatically according to beta and species.prevalence\n\n")
            } else
            {
              message("   --- Determing beta automatically according to alpha and species.prevalence\n\n")
            }
          } else
          {
            message("   --- Determing species.prevalence automatically according to alpha and beta\n\n")
          }
        } else if(prob.method == "linear")
        {
          if(length(c(a, b, species.prevalence)) <= 1)
          {
            if(!is.null(species.prevalence))
            {
              if(is.numeric(species.prevalence))
              {
                if(!(species.prevalence >= 0 & 
                     species.prevalence <= 1))
                {
                  stop("species.prevalence must be a numeric value between 0 and 1.")
                }
              } else 
              {
                stop("species.prevalence must either be a numeric value between 0 and 1
                     or NULL")
              }
            } else
            {
              if(length(c(a, b, species.prevalence)) < 1)
              {
                message("   --- No target prevalence provided; setting the linear
                      transformation to slope 1 and intercept 0.")
                a <- 1
                b <- 0
              } else if(!is.null(a))
              {
                stop("Only the slope (a) of the linear transformation was 
                     provided, please also provide the intercept (b).")
              } else if(!is.null(b))
              {
                stop("Only the intercept (b) of the linear transformation was 
                     provided; please also provide the slope (a).")
              }
            }
          } else if(!is.null(species.prevalence) & 
                    (!is.null(a) | !is.null(b)))
          {
            stop("You should either provide species.prevalence or both a and b,
                 but not a combination of species.prevalence, a and b.")
          }
          
          if(!is.null(b))
          {
            if(!is.numeric(b))
            {
               stop("b must be a numeric value or NULL")
            }
          }
          if(!is.null(a))
          {
            if(!is.numeric(a))
            {
              stop("a must be a numeric value or NULL")
            }
          }
          
          if(!is.null(species.prevalence))
          {
            message("   --- Searching for a linear transformation of 
                    environmental suitability that fits the chosen 
                    species.prevalence.")
          } else
          {
            message(paste0("   --- Determing species prevalence automatically based
                           on a linear transformation of environmental suitability of
                           slope a = ", a,
                           " and intercept b = ", b))
          }
        } else
        {
          stop("prob.method must be either 'logistic' or 'linear'")
        }
        
      }
  
  if (PA.method == "probability")
  {
    if(prob.method == "logistic")
    {
      if(!is.null(species.prevalence))
      {
        if(!is.null(beta))
        {
          alpha.test <- NULL
          for (alpha in c((sp.raster@data@max - sp.raster@data@min)/1000, (sp.raster@data@max - sp.raster@data@min) * 10))
          {
            if(alpha > 0) alpha <- -alpha
            
            
            
            PA.raster <- calc(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- calc(PA.raster, fun = function(x)
            {
              sapply(x, FUN = function(y)
              {
                if(is.na(y))
                { NA } else
                {
                  sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
                }
              }
              )
            }
            )
            freqs <- freq(PA.raster)
            if(any(is.na(freqs[, 1])))
            {
              freqs <- freqs[-which(is.na(freqs[, 1])), ]
            }
            alpha.test <- rbind(alpha.test, c(alpha, 
                                              ifelse(nrow(freqs) == 2,
                                                     freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                     ifelse(freqs[, "value"] == 0,
                                                            0, 1))))
          }
          epsilon <- species.prevalence - alpha.test[, 2]
          if(all(epsilon > 0))
          {
            warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen beta and available environmental conditions (see details).
                          The closest possible estimate of prevalence was", round(alpha.test[2, 2], 2),
                          "\nPerhaps you can try a lower beta value."))
            alpha <- alpha.test[2, 1]
          } else if (all(epsilon < 0))
          {
            warning(paste("Warning, the desired species prevalence cannot be obtained, because of the chosen beta and available environmental conditions (see details).
                          The closest possible estimate of prevalence was", round(alpha.test[1, 2], 2),
                          "\nPerhaps you can try a higher beta value."))
            alpha <- alpha.test[1, 1]
          } else 
          {
            while (all(abs(epsilon) > 0.01))
            {
              alpha <- (alpha.test[which(epsilon == max(epsilon[epsilon < 0])), 1] + 
                          alpha.test[which(epsilon == min(epsilon[epsilon > 0])), 1]) / 2
              PA.raster <- calc(sp.raster, fun = function(x)
              {
                logisticFun(x, beta = beta, alpha = alpha)
              })
              
              PA.raster <- calc(PA.raster, fun = function(x)
              {
                sapply(x, FUN = function(y)
                {
                  if(is.na(y))
                  { NA } else
                  {
                    sample(x = c(0, 1),  size = 1, prob = c(1 - y, y))
                  }
                }
                )
              }
              )
              freqs <- freq(PA.raster)
              if(any(is.na(freqs[, 1])))
              {
                freqs <- freqs[-which(is.na(freqs[, 1])), ]
              }
              alpha.test <- rbind(alpha.test, c(alpha, 
                                                ifelse(nrow(freqs) == 2,
                                                       freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                       ifelse(freqs[, "value"] == 0,
                                                              0, 1))))
              
              epsilon <- species.prevalence - alpha.test[, 2]
            }
          }
        } else
        {
          beta.test <- NULL
          for (beta in c(sp.raster@data@min, sp.raster@data@max))
          {
            transformed.suitability <- calc(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- calc(transformed.suitability, fun=function (x)
            {
              sapply(x, FUN = function(y)
              {
                if(is.na(y))
                { NA } else
                {
                  rbinom(n = 1, size = 1, prob = y)
                }
              }
              )
            })
            # plot(PA.raster)
            freqs <- freq(PA.raster)
            if(any(is.na(freqs[, 1])))
            {
              freqs <- freqs[-which(is.na(freqs[, 1])), , drop = FALSE]
            }
            beta.test <- rbind(beta.test, c(beta, 
                                            ifelse(nrow(freqs) == 2,
                                                   freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                   ifelse(freqs[, "value"] == 0,
                                                          0, 1))))
          }
          epsilon <- data.frame(epsi = species.prevalence - beta.test[, 2], 
                                prevalence = beta.test[, 2])
          if(all(epsilon$epsi > 0))
          {
            warning(paste("Warning, the desired species prevalence may not be obtained, because of the chosen alpha and available environmental conditions (see details).
                          The closest possible estimate of prevalence was", round(beta.test[1, 2], 3),
                          "\nPerhaps you can try an alpha value closer to 0."))
            beta <- beta.test[1, 1]
          } else if (all(epsilon$epsi < 0))
          {
            warning(paste("Warning, the desired species prevalence may be obtained, because of the chosen beta and available environmental conditions (see details).
                          The closest possible estimate of prevalence was", round(beta.test[2, 2], 3),
                          "\nPerhaps you can try an alpha value closer to 0."))
            beta <- beta.test[2, 1]
          } else 
          {
            while (all(apply(epsilon, 1, 
                             function(x) ifelse(abs(x[1]) > 0.001,
                                                TRUE,
                                                ifelse(x[2] == 0, TRUE, FALSE)))))
            {
              beta <- (beta.test[which(epsilon$epsi == max(epsilon$epsi[epsilon$epsi < 0])), 1][1] + 
                         beta.test[which(epsilon$epsi == min(epsilon$epsi[epsilon$epsi > 0])), 1][1]) / 2
              transformed.suitability <- calc(sp.raster, fun = function(x)
              {
                logisticFun(x, beta = beta, alpha = alpha)
              })
              
              PA.raster <- calc(transformed.suitability, fun=function (x)
              {
                sapply(x, FUN = function(y)
                {
                  if(is.na(y))
                  { NA } else
                  {
                    rbinom(n = 1, size = 1, prob = y)
                  }
                }
                )
              })
              freqs <- freq(PA.raster)
              if(any(is.na(freqs[, 1])))
              {
                freqs <- freqs[-which(is.na(freqs[, 1])), , drop = FALSE]
              }
              beta.test <- rbind(beta.test, c(beta, 
                                              ifelse(nrow(freqs) == 2,
                                                     freqs[freqs[, "value"] == 1, "count"] / sum(freqs[, "count"]),
                                                     ifelse(freqs[, "value"] == 0,
                                                            0, 1))))
              epsilon <- data.frame(epsi = species.prevalence - beta.test[, 2], 
                                    prevalence = beta.test[, 2])
            }
          }
        }
      } 
      transformed.suitability <- calc(sp.raster, fun = function(x)
      {
        logisticFun(x, beta = beta, alpha = alpha)
      })
      
      PA.raster <- calc(transformed.suitability, fun=function (x)
      {
        sapply(x, FUN = function(y)
        {
          if(is.na(y))
          { NA } else
          {
            rbinom(n = 1, size = 1, prob = y)
          }
        }
        )
      })
      if(PA.raster@data@max == 0) # Necessary to generate species with very low prevalence
      {                           # Without this step, rasters with only zeros can be generated
        while(PA.raster@data@max == 0)
        {
          transformed.suitability <- calc(sp.raster, fun = function(x)
          {
            logisticFun(x, beta = beta, alpha = alpha)
          })
          
          PA.raster <- calc(transformed.suitability, fun=function (x)
          {
            sapply(x, FUN = function(y)
            {
              if(is.na(y))
              { NA } else
              {
                rbinom(n = 1, size = 1, prob = y)
              }
            }
            )
          })
        }
      }
    } else if(prob.method == "linear")
    {
      if(!is.null(species.prevalence))
      {
        tmp <- .findLinearConversion(sp.raster, 
                                     target.prevalence = species.prevalence)
        
        a <- tmp$a
        b <- tmp$b
        transformed.suitability <- tmp$transformed.suitability
        PA.raster <- tmp$distribution
      } else 
      {
        transformed.suitability <- .transf(sp.raster, c(a, b))
        PA.raster <- calc(transformed.suitability, fun=function (x)
        {
          sapply(x, FUN = function(y)
          {
            if(is.na(y))
            { NA } else
            {
              rbinom(n = 1, size = 1, prob = y)
            }
          }
          )
        })
      }
    } 
  } else if (PA.method == "threshold")
  {
    if(!is.null(species.prevalence))
    {
      beta <- quantile(sp.raster, 1 - species.prevalence)
      names(beta) <- NULL
    }
    message("    - Threshold of conversion into PA:", round(beta, 3), "; method = ", PA.method, "\n")
    
    PA.raster <- reclassify(sp.raster,
                            c(-Inf, beta, 0,
                              beta, +Inf, 1))
  } else {stop("Wrong PA.method entered (either 'probability' or 'threshold')")}
  
  species.prevalence <- round(freq(PA.raster)[2, 2] / sum(freq(PA.raster)[1:2, 2]), 3)
  names(species.prevalence) <- NULL
  
  if("virtualspecies" %in% class(x))
  {
    if(PA.method == "threshold")
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          cutoff = beta,
                          species.prevalence = species.prevalence)
    } else if(prob.method == "logistic")
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          probabilistic.method = prob.method,
                          alpha = alpha,
                          beta = beta,
                          species.prevalence = species.prevalence)
    } else if(prob.method == "linear")
    {
      names(a) <- NULL
      names(b) <- NULL
      x$PA.conversion = c(conversion.method = PA.method,
                          probabilistic.method = prob.method,
                          a = a,
                          b = b,
                          species.prevalence = species.prevalence)
    }
    x$suitab.before.PA.conversion <- x$suitab.raster
    x$suitab.raster <- transformed.suitability
    x$pa.raster <- PA.raster
    results <- x
    if(plot) plot(stack(results$suitab.raster, results$pa.raster), main = c("Suitability", "Presence-absence"))
  } else if ("RasterLayer" %in% class(x))
  {
    if(PA.method == "threshold")
    {
      PA.conversion = c(cutoff = beta,
                        conversion.method = PA.method, 
                        species.prevalence = species.prevalence)
    } else if(prob.method == "logistic")
    {
      PA.conversion = c(conversion.method = PA.method, 
                        probabilistic.method = prob.method,
                        alpha = alpha,
                        beta = beta,
                        species.prevalence = species.prevalence)
    } else if(prob.method == "linear")
    {
      PA.conversion = c(conversion.method = PA.method, 
                        probabilistic.method = prob.method,
                        a = a,
                        b = b,
                        species.prevalence = species.prevalence)
    }
    results <- list(suitab.before.PA.conversion = x,
                    suitab.raster = transformed.suitability,
                    PA.conversion = PA.conversion,
                    pa.raster = PA.raster)
    if(plot) plot(stack(results$suitab.raster, results$pa.raster), main = c("Suitability", "Presence-absence"))
  }
  return(results)
}


.plotTransf <- function(coefs, suit.raster, ...)
{
  x <- seq(suit.raster@data@min,
           suit.raster@data@max, length = 100)
  y <- coefs[1] * x + coefs[2]
  
  y[y < 0] <- 0
  y[y > 1] <- 1
  
  plot(y ~ x, type = "l", bty = "l", las = 1, cex.axis = .8, ...)
}

.transf <- function(x, coefs)
{
  x <- x * coefs[1] + coefs[2]
  x[Which(x < 0)] <- 0
  x[Which(x > 1)] <- 1
  return(x)
}


.findLinearConversion <- function(suit.raster, target.prevalence,
                                 m = FALSE,
                                 plot.conv = FALSE)
{
  
  max.suit <- cellStats(suit.raster, stat = "range")[2]
  min.suit <- cellStats(suit.raster, stat = "range")[1]
  mean.suit <- cellStats(suit.raster, stat = "mean")
  
  slopemax <- (1 - target.prevalence) / (max.suit - mean.suit)
  slopemin <- target.prevalence / (mean.suit - min.suit)
  
  bmax <- slopemax * mean.suit - target.prevalence
  bmin <- slopemin * mean.suit - target.prevalence
  
  new.suit.max = .transf(suit.raster, coefs = c(slopemax, bmax))
  distr.max <- calc(new.suit.max, fun=function (x)
  {
    sapply(x, FUN = function(y)
    {
      if(is.na(y))
      { NA } else
      {
        rbinom(n = 1, size = 1, prob = y)
      }
    }
    )
  })
  prev.max <- cellStats(distr.max, stat = "mean")
  new.suit.min = .transf(suit.raster, coefs = c(slopemin, bmin))
  distr.min <- calc(new.suit.min, fun=function (x)
  {
    sapply(x, FUN = function(y)
    {
      if(is.na(y))
      { NA } else
      {
        rbinom(n = 1, size = 1, prob = y)
      }
    }
    )
  })
  prev.min <- cellStats(distr.min, stat = "mean")
  
  if(target.prevalence > max(c(prev.max, prev.min)))
  {
    if(max(c(prev.max, prev.min)) == prev.min)
    {
      params <-   c(slope = slopemin, 
                    intercept = bmin)
      prev0 <- prev.min
      distr0 <- distr.min
      raster0 <- new.suit.min
    } else
    {
      params <-   c(slope = slopemax, 
                    intercept = bmax)
      prev0 <- prev.max
      distr0 <- distr.max
      raster0 <- new.suit.max
    }
  } else if (target.prevalence < min(c(prev.max, prev.min)))
  {
    if(min(c(prev.max, prev.min)) == prev.min)
    {
      params <-   c(slope = slopemin, 
                    intercept = bmin)
      prev0 <- prev.min
      distr0 <- distr.min
      raster0 <- new.suit.min
    } else
    {
      params <-   c(slope = slopemax, 
                    intercept = bmax)
      prev0 <- prev.max
      distr0 <- distr.max
      raster0 <- new.suit.max
    }
  } else
  {
    pos <- (target.prevalence - min(c(prev.max, prev.min))) /
      (max(c(prev.max, prev.min)) - min(c(prev.max, prev.min))) 
    slopepos <- pos * (max(c(slopemax, slopemin)) - min(c(slopemax, slopemin))) +
      min(c(slopemax, slopemin))
    bpos <- slopepos * mean.suit - target.prevalence
    params <- c(slope = slopepos, 
                intercept = bpos)
    raster0 <- .transf(suit.raster, coefs = params)
    distr0 <- calc(raster0, fun=function (x)
    {
      sapply(x, FUN = function(y)
      {
        if(is.na(y))
        { NA } else
        {
          rbinom(n = 1, size = 1, prob = y)
        }
      }
      )
    })
    prev0 <- cellStats(distr0, stat = "mean") 
  }
  
  
  if(m) cat("prev0 = ", prev0, "\n")
  slop <- params[1]
  
  b0 <- params[2]
  if(abs(prev0 - target.prevalence) > 0.01)
  {
    if(m) cat("Determining b1...\n")
    if(prev0 < target.prevalence)
    {
      b1 <- .99 # Just below 1 to avoid a reponse fully equal to 1
    } else
    {
      b1 <- - cellStats(suit.raster, stat = "max") * slop + 0.01
      # Just above minimum value to avoid a reponse fully equal to 0
    }
    
    raster1 = .transf(suit.raster, coefs = c(slop, b1))
    distr1 = calc(raster1, fun=function (x)
    {
      sapply(x, FUN = function(y)
      {
        if(is.na(y))
        { NA } else
        {
          rbinom(n = 1, size = 1, prob = y)
        }
      }
      )
    })
    prev1 <- cellStats(distr1, stat = "mean")
    if(m) cat("prev1 = ", prev1, "\n")
    if(abs(prev1 - target.prevalence) > 0.01)
    {
      if(m)  cat("Finding a better b...\n")
      b0.5 <- median(c(b0, b1))
      raster0.5 = .transf(suit.raster, coefs = c(slop, b0.5))
      distr0.5 = calc(raster0.5, fun=function (x)
      {
        sapply(x, FUN = function(y)
        {
          if(is.na(y))
          { NA } else
          {
            rbinom(n = 1, size = 1, prob = y)
          }
        }
        )
      })
      prev0.5 <- cellStats(distr0.5, stat = "mean")
      while(abs(prev0.5 - target.prevalence) > 0.01)
      {
        if(target.prevalence >= min(c(prev0, prev0.5)) & 
           target.prevalence <= max(c(prev0, prev0.5)))
        {
          b1 <- b0.5
        } else
        {
          b0 <- b0.5
        }
        b0.5 <- median(c(b0, b1))
        raster0.5 = .transf(suit.raster, coefs = c(slop, b0.5))
        distr0.5 = calc(raster0.5, fun=function (x)
        {
          sapply(x, FUN = function(y)
          {
            if(is.na(y))
            { NA } else
            {
              rbinom(n = 1, size = 1, prob = y)
            }
          }
          )
        })
        prev0.5 <- cellStats(distr0.5, stat = "mean")
        if(m) cat("b0.5 = ", b0.5, " - prevalence = ", prev0.5, "\n")
      }
      rasterfinal <- raster0.5
      bfinal <- b0.5
      prevfinal <- prev0.5
      distrfinal <- distr0.5
      if(m) cat("Best b is equal to b0.5\n")
    } else
    {
      if(m) cat("Best b is equal to b1\n")
      rasterfinal <- raster1
      bfinal <- b1
      prevfinal <- prev1
      distrfinal <- distr1
    }
  } else
  {
    if(m) cat("Best b is equal to b0\n")
    bfinal <- b0
    rasterfinal <- raster0
    prevfinal <- prev0
    distrfinal <- distr0
  }
  if(plot.conv)
  {
    par(mfrow = c(2, 2))
    plot(suit.raster, main = "initial raster")
    plot(rasterfinal,
         main = "transformed raster")
    .plotTransf(c(params[1], bfinal),
                suit.raster, main = paste0("transformation function\ntarget prevalence = ", target.prevalence))
    plot(distrfinal, main = paste0("presence absence raster\nobserved prevalence = ", prevfinal))
  }
  return(list(a = slop,
              b = bfinal,
              prevalence = prevfinal,
              transformed.suitability = rasterfinal,
              distribution = distrfinal))
}

