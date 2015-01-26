% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/generateSpFromFun.R
\name{generateSpFromFun}
\alias{generateSpFromFun}
\title{Generate a virtual species distributions with responses to environmental variables}
\usage{
generateSpFromFun(raster.stack, parameters, rescale = TRUE,
  species.type = "multiplicative", formula = NULL,
  rescale.each.response = TRUE, plot = FALSE)
}
\arguments{
\item{raster.stack}{a RasterStack object, in which each layer represent an environmental
variable.}

\item{parameters}{a list containing the functions of response of the species to
environmental variables with their parameters. See details.}

\item{rescale}{\code{TRUE} or \code{FALSE}. If \code{TRUE}, the final probability of presence is rescaled between 0 and 1.}

\item{species.type}{\code{"additive"}, \code{"multiplicative"} or \code{"mixed"}. Defines
how the final probability of presence is calculated: if \code{"additive"}, responses to each
variable are summed; if \code{"multiplicative"}, responses are multiplicated; if \code{"mixed"}
responses are both summed and multiplicated depending on argument \code{formula}}

\item{formula}{\code{NULL} to create a random formula to calculate the final probability
of presence, or a character string of the form: \code{"layername1 + layername2 *
layername3 * etc."} to manually define it. Only used if \code{species.type} is set to
\code{"mixed"}}

\item{rescale.each.response}{\code{TRUE} or \code{FALSE}. If \code{TRUE}, the individual responses to
each environmental variable are rescaled between 0 and 1 (see details).}

\item{plot}{\code{TRUE} or \code{FALSE}. If \code{TRUE}, the generated virtual species will be plotted.}
}
\value{
a \code{list} with 3 elements:
\itemize{
\item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
\item{\code{details}: the details and parameters used to generate the species}
\item{\code{suitab.raster}: the virtual species distribution, as a Raster object containing the
environmental suitability)}
}
#' The structure of the virtualspecies object can be seen using str()
}
\description{
This function generates a virtual species distribution from a RasterStack of environmental
variables and a defined set of responses to each environmental parameter.
}
\details{
This functions proceeds into several steps:
\enumerate{
\item{The response to each environmental variable is calculated with the functions provided
in \code{parameters}. This results in a probability of presence for each variable.}
\item{If \code{rescale.each.response} is \code{TRUE}, each probability of presence is rescaled between 0 and 1.}
\item{The final probability of presence is calculated according to the chosen \code{species.type}.}
\item{If \code{rescale} is \code{TRUE}, the final probability of presence is rescaled between 0 and 1.}
}
The RasterStack containing environmental variables must have consistent names,
because they will be checked with the \code{parameters}. For example, they can be named
var1, var2, etc. Names can be checked and set with \code{names(my.stack)}.

The \code{parameters} have to be carefully created, otherwise the function will not work:
\itemize{
\item{Either see \code{\link{formatFunctions}} to easily create your list of parameters}
\item{Or create a \code{list} defined according to the following template: \cr
\code{list(
           var1 = list(fun = 'fun1', args = list(arg1 = ..., arg2 = ..., etc.)),
           var2 = list(fun = 'fun2', args = list(arg1 = ..., arg2 = ..., etc.)))}\cr
It is important to keep the same names in the parameters as in the stack of environmental
variables. Similarly, argument names must be identical to argument names in the associated
function (e.g., if you use \code{fun = 'dnorm'}, then args should look like \code{list(mean = 0, sd = 1)}).

See the example section below for more examples.}}


Any response function that can be applied to the environmental variables can
be chosen here. Several functions are proposed in this package:
\code{\link{linearFun}}, \code{\link{logisticFun}} and \code{\link{quadraticFun}}.
Another classical example is the normal distribution: \code{\link[stats]{dnorm}}.
Ther users can also create and use their own functions.


If \code{rescale.each.response = TRUE}, then the probability response to each
variable will be normalised between 0 and 1 according to the following formula:
P.rescaled = (P - min(P)) / (max(P) - min (P))
This rescaling has a strong impact on response functions, so users may prefer to
use \code{rescale.each.response = FALSE} and apply their own rescaling within
their response functions.
}
\examples{
# Create an example stack with two environmental variables
a <- matrix(rep(dnorm(1:100, 50, sd = 25)),
            nrow = 100, ncol = 100, byrow = TRUE)
env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
             raster(a * 1:100))
names(env) <- c("variable1", "variable2")
plot(env) # Illustration of the variables

# Easy creation of the parameter list:
# see in real time the shape of the response functions
parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 1e-04,
                                             sd = 1e-04),
                              variable2 = c(fun = 'linearFun', a = 1, b = 0))

# If you provide env, then you can see the shape of response functions:
parameters <- formatFunctions(x = env,
                              variable1 = c(fun = 'dnorm', mean = 1e-04,
                                             sd = 1e-04),
                              variable2 = c(fun = 'linearFun', a = 1, b = 0))

# Generation of the virtual species
sp1 <- generateSpFromFun(env, parameters)
sp1
par(mfrow = c(1, 1))
plot(sp1)


# Manual creation of the parameter list
# Note that the variable names are the same as above
parameters <- list(variable1 = list(fun = 'dnorm',
                                    args = list(mean = 0.00012,
                                                sd = 0.0001)),
                   variable2 = list(fun = 'linearFun',
                                    args = list(a = 1, b = 0)))
# Generation of the virtual species
sp1 <- generateSpFromFun(env, parameters, plot = TRUE)
sp1
plot(sp1)
}
\author{
Boris Leroy \email{leroy.boris@gmail.com}

with help from C. N. Meynard, C. Bellard & F. Courchamp
}
\seealso{
\code{\link{generateSpFromPCA}} to generate a virtual species with a PCA approach
}
