#' Convert a virtual species distribution (or a suitability raster) 
#' into presence-absence
#' 
#' This functions converts the probabilities of presence from the output of
#'  \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}}, 
#'  \code{\link{generateRandomSp}}
#' or a suitability raster into
#' a presence-absence raster. The conversion can be threshold-based, or based
#' on a probability of conversion (see details).
#' 
#' @param x the output from functions 
#' \code{\link{generateSpFromFun}}, \code{\link{generateSpFromPCA}} 
#' or \code{\link{generateRandomSp}}, or a suitability SpatRaster
#' @param PA.method \code{"threshold"} or \code{"probability"}. If 
#' \code{"threshold"}, then occurrence probabilities are simply converted into
#' presence-absence according to the threshold \code{beta}. If 
#' \code{"probability"}, then
#' probabilities are converted according to a logistic function of threshold 
#' \code{beta} and slope \code{alpha}.
#' @param prob.method \code{"logistic"} or \code{"linear"}. Defines how 
#' the initial environmental suitability is translated into probabilities of 
#' presence/absence.
#' @param beta \code{"random"}, a numeric value in the range of your 
#' probabilities or \code{NULL}. This is the threshold of conversion into 
#' presence-absence (if \code{PA.method = "probability"} and 
#' \code{prob.method = "logistic"}, then beta = the inflexion point).
#' If \code{"random"}, a numeric value will be randomly generated within the range
#' of \code{x}.
#' @param alpha \code{NULL} or a negative numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "logistic"}. 
#' The value of \code{alpha} will
#' shape the logistic function transforming occurrences into presence-absences.
#' See \code{\link{logisticFun}} and examples therein for the choice of 
#' \code{alpha}
#' @param a \code{NULL} or a  numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "linear"}. 
#' Slope of the linear conversion of environmental suitability.
#' @param b \code{NULL} or a  numeric value. Only useful if 
#' \code{PA.method = "probability"} and  \code{proba.method = "linear"}. 
#' Intercept of the linear conversion of environmental suitability.
#' @param species.prevalence \code{NULL} or a numeric value between 0 and 1.
#' The species prevalence is the proportion of sites actually occupied by the
#' species.
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, maps of 
#' probabilities
#' of occurrence and presence-absence will be plotted.
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, D.M. Kaplan, C. Bellard & F. Courchamp
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
#' \href{http://borisleroy.com/virtualspecies_tutorial/04-presenceabsence.html}{
#' Online tutorial for this function}
#' 
#' 
#' 
#' The conversion of environmental suitability into presence - absence used to 
#' be
#' performed by selecting a threshold above which presence always occurs,
#' and never below. However, this approach may is unrealistic because
#' species may sometime be present in areas with a low probability of 
#' occurrence,
#' or be absent from areas with a high probability of occurrence. In addition,
#' when using a threshold you erase the previously generated response shapes: 
#' it all becomes threshold. Thus, this threshold approach should be avoided.
#' 
#'  
#' A more
#' realistic conversion consists in converting environmental suitability into
#' presence -
#' absence with a probability function (see references). Such a probability 
#' conversion can be performed with two different methods here:
#' \enumerate{
#' \item{Using a logistic transformation of  environmental suitability
#' (see \code{\link{logisticFun}}). A logistic function on the other hand, will 
#' ensure that the simulated probability is within the 0-1 range and allow easy 
#' control of species prevalence. However, the 
#' logistic function will also flatten out the relationship at the extreme 
#' suitability values, and narrow or broaden the intermediate probability values
#' depending on the slope of the logistic curve }
#' \item{Using a linear transformation of environmental suitability. A linear 
#' transformation will preserve the shape of the originally simulated 
#' occurrence-environment relationships, uniformly increasing or decreasing the
#' probabilities of occurrence across the landscape.}}
#' 
#' --- note ---
#' 
#' If the Virtual Species study aims at comparing simulated and predicted 
#' probability values, it is important to recover the correct simulated 
#' probability instead of directly using the initial suitability function. 
#' 
#' Therefore, the function stores the probability of occurrence in the 
#' output list, under the object \code{probability.of.occurrence}.
#' The initial suitability function (before logistic or linear conversion)
#' will still be stored in the output list as \code{suitab.raster}. 
#' 
#' --------------------------------------------------------------------------
#' 
#' \bold{PROBABILISTIC CONVERSION - LOGISTIC METHOD}
#' 
#' To perform the logistic transformation of environmental suitability
#' you have to define two of the
#' three following parameters:
#' \itemize{
#' \item{\code{beta}: the 'threshold' of the logistic function (i.e. the 
#' inflexion point. It should normaly be in the range of values of your 
#' environmental suitability.)}
#' \item{\code{alpha}: the slope of the logistic function. It should generally
#' be in value equal to something like 1/20 or 1/10 of your environmental 
#' suitability range}
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
#' or \code{beta}. In these conditions, the function will retain the 
#' \code{alpha} or
#' \code{beta} which provides the closest prevalence to your 
#' \code{species.prevalence},
#' but you may also provide another value of \code{alpha} or \code{beta} to 
#' obtain
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
#' \code{species.prevalence}
#' (see below).}
#' } 
#' 
#' Method used to find coefficients of a linear regression which results in the
#' requested \code{species.prevalence}:
#' 
#' \enumerate{
#' \item{The simplest linear transformation of habitat suitability would
#' be to just multiply the raw suitability by a constant. For example, if the 
#' raw average suitability in the area is 0.04, it means an expected prevalence
#' of 0.40. To to go from this expected prevalence of 0.04 to an expected
#' prevalence of 0.4, we can just multiply the raw suitability by 10. It is the
#' default choice, unless it results in probabilities superior to 1 or raw
#' suitability have values below 0, in which case the function proceeds to
#'  method 2.}
#' \item{If it does not work, then we look at the line that passes through 
#' (min suitability, 0) and (mean suitability, desired prevalence). For this 
#' line, we only need to ensure that the maximum probability of occurence is 
#' lower than 1. Otherwise, the function proceeds to method 3.}
#' \item{If method 2 fails, then we test the line going through (mean 
#' suitability, desired prevalence) and (max suitability, 1). If the minimum
#' probability resulting from this line is greater than 0, then this method is 
#' correct.
#' }
#' }
#' 
#' One of these 3 lines should always work. In fact, one of the last two has to 
#' work, and it does not hurt to try the first one which is simpler.
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
#' a \code{list} containing 6 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, 
#' \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the 
#' species}
#' \item{\code{suitab.raster}: the environmental suitability of your virtual 
#' species, as a Raster object }
#' \item{\code{probability.of.occurrence}: the probability of occurrence of your 
#' species, based on the chosen transformation of environmental suitability,
#' as a Raster object }
#' \item{\code{PA.conversion}: the parameters used to convert the suitability 
#' into presence-absence}
#' \item{\code{pa.raster}: the presence-absence map, as a Raster object 
#' containing 0 (absence) / 1 (presence) / NA}
#' }
#' The structure of the virtualspecies object can be seen using \code{str()}
#' 
#' 
#' @importFrom stats median rbinom runif
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- c(rast(a * dnorm(1:100, 50, sd = 25)),
#'          rast(a * 1:100))
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
  if(inherits(x, "virtualspecies"))
  {
    if(x$approach == "bca") {
      approach <- "bca"
      sp.raster <- x$suitab.raster.current
    } else if(inherits(x$suitab.raster, "SpatRaster"))
    {
      approach <- "not bca"
      sp.raster <- x$suitab.raster
    } else stop("x must be:\n- a raster layer object from package terra", 
                " \nor\n- the output list", 
                " from functions generateSpFromFun(), generateSpFromPCA(), ", 
                "generateRandomSp() or generateSpFromBCA()")
  } else if (inherits(x, "RasterLayer"))
  {
    approach <- "not bca"
    sp.raster <- rast(x)
  } else if (inherits(x, "SpatRaster")) {
    approach <- "not bca"
    sp.raster <- x
  } else stop("x must be:\n- a raster layer object from package terra", 
              " \nor\n- the output list", 
              " from functions generateSpFromFun(), generateSpFromPCA(), ", 
              "generateRandomSp() or generateSpFromBCA()")
  min_p <- global(sp.raster, min, na.rm = TRUE)[1, 1]
  max_p <- global(sp.raster, max, na.rm = TRUE)[1, 1]
  

  if(PA.method == "threshold")
  {
    if(is.numeric(beta))
    {
      if(is.numeric(species.prevalence))
      {
        warning("Both beta and species.prevalence were provided. beta will be
                ignored.")
        beta <- NULL
      } else if(beta < min_p) 
      {
        warning("beta is lower than all values in your suitability raster. The
                species will most likely be present everywhere")
      } else if(beta > max_p)
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
        beta <- sample(seq(min_p, 
                           max_p, length = 1000), 1)
        
        message("   --- Generating a random value of beta for the threshold conversion\n\n")
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
          alpha <- -sample(c(seq((max_p - min_p)/1000,
                                 (max_p - min_p)/100, length = 10),
                             seq((max_p - min_p)/100,
                                 (max_p - min_p)/10, length = 100),
                             seq((max_p - min_p)/10,
                                 (max_p - min_p)*10, length = 10)), size = 1)
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
          if(beta < min_p) 
          {
            warning("beta is lower than all values in your suitability raster. The
                        species will most likely be present everywhere")
          } else if(beta > max_p)
          {
            warning("beta is higher than all values in your suitability raster. 
                        The species will most likely be absent everywhere")
          }
        } else if(beta == "random")
        {
          beta <- sample(seq(min_p, 
                             max_p, length = 1000), 1)
          message("   --- Generating a random value of beta for the logistic conversion\n\n")
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
                    species.prevalence.\n")
      } else
      {
        message(paste0("   --- Determing species prevalence automatically based
                           on a linear transformation of environmental suitability of
                           slope a = ", a,
                       " and intercept b = ", b, "\n"))
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
          for (alpha in c((max_p - min_p)/1000, (max_p - min_p) * 10))
          {
            if(alpha > 0) alpha <- -alpha
            
            
            
            proba.of.occurrence <- app(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
            
            alpha.test <- rbind(alpha.test, c(alpha, 
                                              global(PA.raster,
                                                     stat = "mean", 
                                                     na.rm = TRUE)[1, 1]))
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
              alpha <- (alpha.test[which(epsilon == 
                                           max(epsilon[epsilon < 0])), 1] + 
                          alpha.test[which(epsilon == 
                                             min(epsilon[epsilon > 0])), 1]) / 2
              proba.of.occurrence <- app(sp.raster, fun = function(x)
              {
                logisticFun(x, beta = beta, alpha = alpha)
              })
              PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
              
              alpha.test <- rbind(alpha.test, c(alpha, 
                                                global(PA.raster,
                                                       stat = "mean",
                                                       na.rm = TRUE)[1, 1]))
              
              epsilon <- species.prevalence - alpha.test[, 2]
            }
          }
        } else
        {
          beta.test <- NULL
          # We define the upper and lower boundaries for beta.
          # We choose to be able to define beta values beyond the boundaries of
          # our probability of occurrence raster, so we can have a larger range 
          # of prevalence
          for (beta in c(min_p - 
                         diff(c(min_p,
                                max_p)) / 2,
                         max_p + 
                         diff(c(min_p,
                                max_p)) / 2))
          {
            proba.of.occurrence <- app(sp.raster, fun = function(x)
            {
              logisticFun(x, beta = beta, alpha = alpha)
            })
            
            PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
            
            beta.test <- rbind(beta.test, 
                               c(beta, 
                                 global(PA.raster, stat = "mean", 
                                        na.rm = TRUE)[1, 1]))
          }
          epsilon <- data.frame(epsi = species.prevalence - beta.test[, 2], 
                                prevalence = beta.test[, 2])
          if(all(epsilon$epsi > 0))
          {
            warning(paste("Warning, the desired species prevalence may not", 
                          " be obtained, because of the chosen alpha and ", 
                          "available environmental conditions (see details). ", 
                          "The closest possible estimate of prevalence was", 
                          round(beta.test[1, 2], 3),
                          "\nPerhaps you can try an alpha value closer to 0."))
            beta <- beta.test[1, 1]
          } else if (all(epsilon$epsi < 0))
          {
            warning(paste("Warning, the desired species prevalence may be ", 
            "obtained, because of the chosen beta and available environmental ",
            "conditions (see details). ", 
            "The closest possible estimate of prevalence was", 
            round(beta.test[2, 2], 3),
            "Perhaps you can try an alpha value closer to 0."))
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
              proba.of.occurrence <- app(sp.raster, fun = function(x)
              {
                logisticFun(x, beta = beta, alpha = alpha)
              })
              
              PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
              
              beta.test <- rbind(beta.test,
                                 c(beta, 
                                   global(PA.raster, stat = "mean",
                                          na.rm = TRUE)[1, 1]))
              epsilon <- data.frame(epsi = species.prevalence - beta.test[, 2], 
                                    prevalence = beta.test[, 2])
            }
          }
        }
      } 
      proba.of.occurrence <- app(sp.raster, fun = function(x)
      {
        logisticFun(x, beta = beta, alpha = alpha)
      })
      
      PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
      
      if(global(PA.raster, stat = "max", 
                na.rm = TRUE)[1, 1] == 0) # Necessary to generate 
        # species with very low prevalence
        # Without this step, rasters with only zeros can be generated
      {                           
        while(global(PA.raster, stat = "max", na.rm = TRUE)[1, 1] == 0)
        {
          proba.of.occurrence <- app(sp.raster, fun = function(x)
          {
            logisticFun(x, beta = beta, alpha = alpha)
          })
          
          PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
        }
      }
      
      if(approach == "bca") {
        proba.of.occurrence.future <- app(x$suitab.raster.future,
                                          fun = function(x)
        {
          logisticFun(x, beta = beta, alpha = alpha)
        })
        
        PA.raster.future <- .quickBernoulliTrial(proba.of.occurrence.future)
      }
    } else if(prob.method == "linear")
    {
      if(!is.null(species.prevalence))
      {
        tmp <- .findLinearConversion(sp.raster, 
                                     target.prevalence = species.prevalence)
        
        a <- tmp$a
        b <- tmp$b
        proba.of.occurrence <- unwrap(tmp$probability.of.occurrence)
        PA.raster <- unwrap(tmp$distribution)
        if(approach == "bca") {
          proba.of.occurrence.future <- .transf(x$suitab.raster.future, c(a, b))
          PA.raster.future <- .quickBernoulliTrial(proba.of.occurrence.future)
        }
      } else 
      {
        proba.of.occurrence <- .transf(sp.raster, c(a, b))
        PA.raster <- .quickBernoulliTrial(proba.of.occurrence)
        if(approach == "bca") {
          proba.of.occurrence.future <- .transf(x$suitab.raster.future, c(a, b))
          PA.raster.future <- .quickBernoulliTrial(proba.of.occurrence.future)
        }
      }
    } else if(prob.method == "truncated linear")
    { # This one is not documented in the package, I have only used it for
      # self testings. It works if anyone is interested enough to read these 
      # lines ;)
      tmp <- .findTruncatedLinearConversion(
        sp.raster, 
        target.prevalence = species.prevalence)
      
      a <- tmp$a
      b <- tmp$b
      proba.of.occurrence <- unwrap(tmp$probability.of.occurrence)
      PA.raster <- unwrap(tmp$distribution)
      if(approach == "bca") {
        proba.of.occurrence.future <- .transf(x$suitab.raster.future, c(a, b))
        PA.raster.future <- .quickBernoulliTrial(proba.of.occurrence.future)
      }
    }
  } else if (PA.method == "threshold")
  {
    if(!is.null(species.prevalence))
    {
      beta <- global(sp.raster, 
                     fun = quantile,
                     probs = 1 - species.prevalence,
                     na.rm = TRUE)[1, 1]
    }
    
    PA.raster <- proba.of.occurrence <-  
      classify(sp.raster,
               matrix(data = c(-Inf, beta, 0,
                               beta, +Inf, 1), 
                      nc = 3, nr = 2, byrow = TRUE))
    if(approach == "bca") {
      PA.raster.future <- proba.of.occurrence.future <-
        classify(x$suitab.raster.future,
                 matrix(data = c(-Inf, beta, 0,
                                 beta, +Inf, 1), 
                        nc = 3, nr = 2, byrow = TRUE))
        
    }
    
  } else {
    stop("Wrong PA.method entered (either 'probability' or 'threshold')")
  }
  
  species.prevalence <- global(PA.raster, stat = "mean", na.rm = TRUE)[1, 1]
  
  if(inherits(x, "virtualspecies"))
  {
    if(PA.method == "threshold")
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          cutoff = beta,
                          species.prevalence = species.prevalence)
      message(paste0("   Threshold conversion finished:
              \n- cutoff = ", beta, 
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    } else if(prob.method == "logistic")
    {
      x$PA.conversion = c(conversion.method = PA.method,
                          probabilistic.method = prob.method,
                          alpha = alpha,
                          beta = beta,
                          species.prevalence = species.prevalence)
      message(paste0("   Logistic conversion finished:
              \n- beta = ", beta, 
                     "\n- alpha = ", alpha,
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    } else if(prob.method == "linear")
    {
      names(a) <- NULL
      names(b) <- NULL
      x$PA.conversion = c(conversion.method = PA.method,
                          probabilistic.method = prob.method,
                          a = a,
                          b = b,
                          species.prevalence = species.prevalence)
      message(paste0("   Linear conversion finished:
              \n- slope (a) = ", a, 
                     "\n- intercept (b) = ", b,
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    }
    if(approach == "bca") {
      x$probability.of.occurrence.current <- wrap(proba.of.occurrence)
      x$pa.raster.current <- wrap(PA.raster)
      x$probability.of.occurrence.future <- wrap(proba.of.occurrence.future)
      x$pa.raster.future <- wrap(PA.raster.future)
    } else {
      x$probability.of.occurrence <- wrap(proba.of.occurrence)
      x$pa.raster <- wrap(PA.raster)
    }

    
    
    results <- x
    if(plot) {
      if(approach == "bca") {
        plot(c(results$suitab.raster.current,
               results$probability.of.occurrence.current, 
               results$pa.raster.current,
               results$suitab.raster.future,
               results$probability.of.occurrence.future, 
               results$pa.raster.future),
             main = c("Current suitability",
                      "Current prob of occ", 
                      "Current pres-abs",
                      "Future suitability",
                      "Future prob of occ", 
                      "Future pres-abs"),
             col = rev(viridis::magma(10)))
      } else {
        plot(c(results$suitab.raster,
               results$probability.of.occurrence, 
               results$pa.raster),
             main = c("Environmental suitability",
                      "Probability of occurrence", 
                      "Presence-absence"),
             col = rev(viridis::magma(10)))
      }

    }
  } else if (inherits(x, "SpatRaster"))
  {
    if(PA.method == "threshold")
    {
      PA.conversion = c(cutoff = beta,
                        conversion.method = PA.method, 
                        species.prevalence = species.prevalence)
      message(paste0("   Threshold conversion finished:
                     \n- cutoff = ", beta, 
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    } else if(prob.method == "logistic")
    {
      PA.conversion = c(conversion.method = PA.method, 
                        probabilistic.method = prob.method,
                        alpha = alpha,
                        beta = beta,
                        species.prevalence = species.prevalence)
      
      message(paste0("   Logistic conversion finished:
                     \n- beta = ", beta, 
                     "\n- alpha = ", alpha,
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    } else if(prob.method == "linear")
    {
      PA.conversion = c(conversion.method = PA.method, 
                        probabilistic.method = prob.method,
                        a = a,
                        b = b,
                        species.prevalence = species.prevalence)
      
      message(paste0("   Linear conversion finished:
                     \n- slope (a) = ", a, 
                     "\n- intercept (b) = ", b,
                     "\n- species prevalence =", species.prevalence, "\n\n"))
    }

    results <- list(suitab.raster = wrap(x),
                    probability.of.occurrence = wrap(proba.of.occurrence),
                    PA.conversion = PA.conversion,
                    pa.raster = wrap(PA.raster))
    
    if(plot) plot(c(results$suitab.raster,
                        results$probability.of.occurrence,
                        results$pa.raster), 
                  main = c("Suitability",
                           "Probability of ocurrence", 
                           "Presence-absence"),
                  col = rev(viridis::magma(10)))
  }
  return(results)
}

.quickBernoulliTrial <- function(prob.raster, ...)
{
  # Raster of same dimentions than the  probability raster
  random.numbers <- rast(x = prob.raster)
  # Generate random numbers between 0 and 1 from uniform distribution
  random.numbers <- setValues (random.numbers, runif(ncell(prob.raster), 0, 1)) 
  # Attribute presence or absence on the basis of whether the random number
  # is above or below the probability of occurrence
  pa.raster <- prob.raster > random.numbers
}



.transf <- function(x, coefs)
{
  x <- x * coefs[1] + coefs[2]
  minval <- global(x, min, na.rm = TRUE)[1, 1]
  maxval <- global(x, max, na.rm = TRUE)[1, 1]
    
  if(minval < 0 | maxval > 1)
  {
    if(minval < 0 & maxval > 1)
    {
      message("The linear transformation resulted in probability values below 0
and above 1, so these were respectively truncated to 0 and 1\n")
    } else if(minval < 0 )
    {
      message("The linear transformation resulted in probability values below 0
so these were truncated to 0\n")
    } else
    {
      message("The linear transformation resulted in probability values above 1
so these were truncated to 1\n")
    }
  }
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x)
}

.findLinearConversion = function(suit.raster,
                                 target.prevalence)
{
  suit.max <- global(suit.raster, max, na.rm = TRUE)[1, 1]
  suit.mean <- global(suit.raster, stat = "mean", na.rm = TRUE)[1, 1]
  suit.min <- global(suit.raster, na.rm = TRUE)[1, 1]
  
  xs = c(suit.min, suit.max)
  ys = c(0, 1)
  
  # Only include (0,0) case if suitability >= 0
  if (suit.min >= 0) 
  {
    xs = c(0, xs)
    ys = c(0, ys) 
  }
  
  AB = .abcoefs(suit.mean,
                target.prevalence,
                xs,
                ys)
  
  ymn = .lab(suit.min, AB$a, AB$b)
  ymx = .lab(suit.max, AB$a, AB$b)
  
  # Round to avoid very small floating point calculation errors
  ymn = round(ymn, 6)
  ymx = round(ymx, 6)
  
  I = min(which(ymn >= 0 & ymx <= 1)) # Find first one that works
  
  
  #Calculate the resulting prevalence:
  new.suit <- AB$a[I] * suit.raster + AB$b[I]
  distr = .quickBernoulliTrial(new.suit)
  prev = global(distr, stat = "mean", na.rm = TRUE)[1, 1]
  
  return(list(a = AB$a[I],
              b = AB$b[I],
              prevalence = prev,
              probability.of.occurrence = wrap(new.suit),
              distribution = wrap(distr)))
}

# Get line coefficients from two points
.abcoefs = function(x1, y1, x2, y2)
{
  list(b = y1 - x1 * (y1 - y2) / (x1 - x2), 
       a = (y1 - y2) / (x1 - x2)) 
}

# Function for a line with intercept (b) and slope (a)
.lab = function(x, b, a)
{
  a * x + b 
}


.findTruncatedLinearConversion <- function(suit.raster, target.prevalence,
                                 m = FALSE,
                                 plot.conv = FALSE)
{
  
  max.suit <- global(suit.raster, stat = "max", na.rm = TRUE)[1, 1]
  min.suit <- global(suit.raster, stat = "min", na.rm = TRUE)[1, 1]
  mean.suit <- global(suit.raster, stat = "mean", na.rm = TRUE)[1, 1]
  
  slopemax <- (1 - target.prevalence) / (max.suit - mean.suit)
  slopemin <- target.prevalence / (mean.suit - min.suit)
  
  bmax <- slopemax * mean.suit - target.prevalence
  bmin <- slopemin * mean.suit - target.prevalence
  
  new.suit.max = .transf(suit.raster, coefs = c(slopemax, bmax))
  distr.max <- .quickBernoulliTrial(new.suit.max)
  prev.max <- global(distr.max, stat = "mean", na.rm = TRUE)[1, 1]
  new.suit.min = .transf(suit.raster, coefs = c(slopemin, bmin))
  distr.min <- .quickBernoulliTrial(new.suit.min)
  prev.min <- global(distr.min, stat = "mean", na.rm = TRUE)[1, 1]
  
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
    slopepos <- pos * (max(c(slopemax, slopemin)) -
                         min(c(slopemax, slopemin))) +
      min(c(slopemax, slopemin))
    bpos <- slopepos * mean.suit - target.prevalence
    params <- c(slope = slopepos, 
                intercept = bpos)
    raster0 <- .transf(suit.raster, coefs = params)
    distr0 <- .quickBernoulliTrial(raster0)
    prev0 <- global(distr0, stat = "mean", na.rm = TRUE)[1, 1]
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
      b1 <- - global(suit.raster, stat = "max",
                     na.rm = TRUE)[1, 1] * slop + 0.01
      # Just above minimum value to avoid a reponse fully equal to 0
    }
    
    raster1 = .transf(suit.raster, coefs = c(slop, b1))
    distr1 = .quickBernoulliTrial(raster1)
    prev1 <- global(distr1, stat = "mean", na.rm = TRUE)[1, 1]
    if(m) cat("prev1 = ", prev1, "\n")
    if(abs(prev1 - target.prevalence) > 0.01)
    {
      if(m)  cat("Finding a better b...\n")
      b0.5 <- median(c(b0, b1))
      raster0.5 = .transf(suit.raster, coefs = c(slop, b0.5))
      distr0.5 = .quickBernoulliTrial(raster0.5)
      prev0.5 <- global(distr0.5, stat = "mean", na.rm = TRUE)[1, 1]
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
        b0.5 <- stats::median(c(b0, b1))
        raster0.5 = .transf(suit.raster, coefs = c(slop, b0.5))
        distr0.5 = .quickBernoulliTrial(raster0.5)
        prev0.5 <- global(distr0.5, stat = "mean", na.rm = TRUE)[1, 1]
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

  return(list(a = slop,
              b = bfinal,
              prevalence = prevfinal,
              probability.of.occurrence = wrap(rasterfinal),
              distribution = wrap(distrfinal)))
}

