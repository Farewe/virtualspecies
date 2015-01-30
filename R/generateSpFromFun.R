#' Generate a virtual species distributions with responses to environmental variables
#' 
#' This function generates a virtual species distribution from a RasterStack of environmental
#' variables and a defined set of responses to each environmental parameter.
#' 
#' @param raster.stack a RasterStack object, in which each layer represent an environmental 
#' variable.
#' @param parameters a list containing the functions of response of the species to 
#' environmental variables with their parameters. See details.
#' @param rescale \code{TRUE} or \code{FALSE}. If \code{TRUE}, the final probability of presence is rescaled between 0 and 1.
#' @param species.type \code{"additive"}, \code{"multiplicative"} or \code{"mixed"}. Defines 
#' how the final probability of presence is calculated: if \code{"additive"}, responses to each
#' variable are summed; if \code{"multiplicative"}, responses are multiplicated; if \code{"mixed"}
#' responses are both summed and multiplicated depending on argument \code{formula}
#' @param formula \code{NULL} to create a random formula to calculate the final probability 
#' of presence, or a character string of the form: \code{"layername1 + layername2 *
#' layername3 * etc."} to manually define it. Only used if \code{species.type} is set to
#' \code{"mixed"}
#' @param rescale.each.response \code{TRUE} or \code{FALSE}. If \code{TRUE}, the individual responses to
#' each environmental variable are rescaled between 0 and 1 (see details).
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the generated virtual species will be plotted.
#' @return a \code{list} with 3 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the raster containing the environmental suitability of the virtual species}
#' }
#' The structure of the virtualspecies object can be seen using str()
#' @seealso \code{\link{generateSpFromPCA}} to generate a virtual species with a PCA approach
#' @details
#' This functions proceeds into several steps:
#' \enumerate{
#' \item{The response to each environmental variable is calculated with the functions provided
#' in \code{parameters}. This results in a probability of presence for each variable.}
#' \item{If \code{rescale.each.response} is \code{TRUE}, each probability of presence is rescaled between 0 and 1.}
#' \item{The final probability of presence is calculated according to the chosen \code{species.type}.}
#' \item{If \code{rescale} is \code{TRUE}, the final probability of presence is rescaled between 0 and 1,
#' with the formula (val - min) / (max - min).}
#' }
#' The RasterStack containing environmental variables must have consistent names, 
#' because they will be checked with the \code{parameters}. For example, they can be named
#' var1, var2, etc. Names can be checked and set with \code{names(my.stack)}.
#' 
#' The \code{parameters} have to be carefully created, otherwise the function will not work:
#' \itemize{
#' \item{Either see \code{\link{formatFunctions}} to easily create your list of parameters}
#' \item{Or create a \code{list} defined according to the following template: \cr
#' \code{list(
#'            var1 = list(fun = 'fun1', args = list(arg1 = ..., arg2 = ..., etc.)),
#'            var2 = list(fun = 'fun2', args = list(arg1 = ..., arg2 = ..., etc.)))}\cr
#' It is important to keep the same names in the parameters as in the stack of environmental
#' variables. Similarly, argument names must be identical to argument names in the associated 
#' function (e.g., if you use \code{fun = 'dnorm'}, then args should look like \code{list(mean = 0, sd = 1)}).
#' 
#' See the example section below for more examples.}}
#'            
#' 
#' Any response function that can be applied to the environmental variables can
#' be chosen here. Several functions are proposed in this package:
#' \code{\link{linearFun}}, \code{\link{logisticFun}} and \code{\link{quadraticFun}}.
#' Another classical example is the normal distribution: \code{\link[stats]{dnorm}}.
#' Ther users can also create and use their own functions.
#' 
#'   
#' If \code{rescale.each.response = TRUE}, then the probability response to each
#' variable will be normalised between 0 and 1 according to the following formula:
#' P.rescaled = (P - min(P)) / (max(P) - min (P))
#' This rescaling has a strong impact on response functions, so users may prefer to
#' use \code{rescale.each.response = FALSE} and apply their own rescaling within
#' their response functions.
#' 
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100))
#' names(env) <- c("variable1", "variable2")
#' plot(env) # Illustration of the variables
#' 
#' # Easy creation of the parameter list:
#' # see in real time the shape of the response functions
#' parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 1e-04, 
#'                                              sd = 1e-04),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#'                               
#' # If you provide env, then you can see the shape of response functions:
#' parameters <- formatFunctions(x = env,
#'                               variable1 = c(fun = 'dnorm', mean = 1e-04, 
#'                                              sd = 1e-04),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' 
#' # Generation of the virtual species
#' sp1 <- generateSpFromFun(env, parameters)
#' sp1
#' par(mfrow = c(1, 1))
#' plot(sp1)
#' 
#' 
#' # Manual creation of the parameter list
#' # Note that the variable names are the same as above
#' parameters <- list(variable1 = list(fun = 'dnorm',
#'                                     args = list(mean = 0.00012,
#'                                                 sd = 0.0001)),
#'                    variable2 = list(fun = 'linearFun',
#'                                     args = list(a = 1, b = 0)))
#' # Generation of the virtual species
#' sp1 <- generateSpFromFun(env, parameters, plot = TRUE)
#' sp1
#' plot(sp1)


generateSpFromFun <- function(raster.stack, parameters, 
                           rescale = TRUE, 
                           species.type = "multiplicative", formula = NULL, rescale.each.response = TRUE,
                           plot = FALSE)
{
  approach <- "response"
  if(!(is(raster.stack, "Raster")))
  {
    stop("raster.stack must be a raster stack object")
  }
  if(any(is.na(maxValue(raster.stack))))
  {
    raster.stack <- setMinMax(raster.stack)
  }
  n.l <- nlayers(raster.stack)
  if(n.l != length(parameters)) 
  {stop("Provide as many layers in raster.stack as functions on parameters")}
  if(any(!(names(parameters) %in% names(raster.stack)) |
           !(names(raster.stack) %in% names(parameters))))
     {stop("Layer names and names of parameters must be identical")}
  # Checking the structure and consistency of parameters
  for (i in 1:length(parameters))
  {
    if(any(!(c("fun", "args") %in% names(parameters[[i]]))))
    {stop("The structure of parameters does not seem correct. 
          Please provide function and arguments for variable '",
          names(parameters)[i], "'. See help(generateSpFromFun) for more details.",
          sep = "")}
    test <- tryCatch(match.fun(parameters[[i]]$fun), error = function(c) "error")
    if(class(test) != "function")
    {
      stop(paste("The function ", parameters[[i]]$fun, " does not exist, please verify spelling.", sep = ""))
    }
    if(any(!(names(parameters[[i]]$args) %in% names(formals(fun = test)))))
    {
      stop(paste("Arguments of variable '", names(parameters)[i], "' (", 
                 paste(names(parameters[[i]]$args), collapse = ", "), 
                 ") do not match arguments of the associated function\n
                 List of possible arguments for this function: ",
                 paste(names(formals(fun = test)), collapse = ", "), sep = ""))
    }
    rm(test)
  }
  suitab.raster <- stack(sapply(names(raster.stack), FUN = function(y)
  {
    calc(raster.stack[[y]], fun = function(x)
    {
      do.call(match.fun(parameters[[y]]$fun), args = c(list(x), parameters[[y]]$args))
    }
    )
  }))
  
  for (var in names(raster.stack))
  {
    parameters[[var]]$min <- raster.stack[[var]]@data@min
    parameters[[var]]$max <- raster.stack[[var]]@data@max
  }
  
  if(rescale.each.response)
  {
    suitab.raster <- stack(sapply(names(suitab.raster), function(y)
      {
        (suitab.raster[[y]] - suitab.raster[[y]]@data@min) / (suitab.raster[[y]]@data@max - suitab.raster[[y]]@data@min)
      }))
  }

  if(species.type == "multiplicative")
  {
    suitab.raster <- raster::overlay(suitab.raster, fun = prod)
  } else if (species.type == "additive")
  {
    suitab.raster <- raster::overlay(suitab.raster, fun = sum)
  } else if (species.type == "mixed")
  {
    layer.names <- sample(names(suitab.raster))
    if(!is.null(formula))
    {
      if (length(strsplit(formula, " ")[[1]]) != (nlayers(suitab.raster) * 2 - 1))
      {stop("The entered formula isn't correct. Check that the layer names are correct, and that operators are correctly spaced:\n
             formula = 'layername1 + layername2 * layername3 * etc.'")}
      else if (any(!(strsplit(formula, " ")[[1]][!(strsplit(formula, " ")[[1]] %in% c("+", "*", "/", "-"))] %in% names(raster.stack))))
      {
        stop("The entered formula isn't correct. The layer names in the formula do not seem to correspond to layer names in raster.stack:\n
             formula = 'layername1 + layername2 * layername3 * etc.'")
      }
    } else
    {
      formula <- layer.names[1]
      for (i in layer.names[2:nlayers(suitab.raster)])
      {
        formula <- paste(formula, sample(c("*", "+"), 1), i)
      }
    }

    operators <- strsplit(formula, " ")[[1]]
    operators <- operators[(1:length(operators)) %% 2 == 0]
    id <- c(1:nlayers(suitab.raster), 1:(nlayers(suitab.raster) - 1) + 0.5)
    
    mixed.fun <- NULL
    eval(parse(text = paste("mixed.fun <- function(",
                            paste("x", 1:nlayers(suitab.raster), sep = "", collapse = ", "),
                            ") {",
                            paste(c(paste("x", 1:nlayers(suitab.raster), sep = ""), operators)[order(id)], collapse = " "),
                            "}"
                            )))
    
    suitab.raster <- overlay(suitab.raster, fun = mixed.fun)
    print(formula)
  } else
  {
    stop("species.type must be either 'multiplicative', 'additive' or 'mixed'")
  }
  if(rescale)
  {
    suitab.raster <- (suitab.raster - suitab.raster@data@min) / (suitab.raster@data@max - suitab.raster@data@min)
  }
    
  if(species.type == "mixed")
  {
    results <- list(approach = approach,
                    
                    details = list(variables = names(parameters),
                                   sp.type = c(sp.type = species.type,
                                               formula = formula),
                                   rescale.each.response = rescale.each.response,
                                   rescale = rescale,
                                   parameters = parameters),
                    suitab.raster = suitab.raster
    )
  } else
  {
    results <- list(approach = approach,
                    details = list(variables = names(parameters),
                                   sp.type = species.type,
                                   rescale.each.response = rescale.each.response,
                                   rescale = rescale,
                                   parameters = parameters),
                    suitab.raster = suitab.raster
    )
  }
  if(plot)
  {
    plot(results$suitab.raster, main = "Environmental suitability of the virtual species")
  }
  
  class(results) <- append(class(results), "virtualspecies")

  return(results)
}
