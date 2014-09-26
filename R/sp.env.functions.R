#' Linear function
#' @description A simple linear function of the form
#' \deqn{ax+b}{a*x+b}
#' @param x a numeric value or vector
#' @param a a numeric value or vector
#' @param b a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{logisticFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' x <- 1:100
#' y <- linearFun(x, a = 0.5, b = 0)
#' plot(y ~ x, type = "l")
linearFun <- function(x, a, b) {a * x + b}

#' Logistic function
#' 
#' @description A simple logistic function of the form
#' \deqn{\frac{1}{{1 + e^{\frac{x - \beta}{\alpha}}}}}{
#' 1 / (1 + exp((x - \beta)/\alpha))}
#' @param x a numeric value or vector
#' @param alpha a numeric value or vector
#' @param beta a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @details
#' The value of \code{beta} determines the 'threshold' of the logistic curve 
#' (i.e. the inflexion point).
#' 
#' The value of \code{alpha} determines the slope of the curve (see examples):
#' \itemize{
#' \item{\code{alpha} very close to 0 will result in a threshold-like response.}
#' \item{Values of \code{alpha} with the same order of magnitude as the range of
#' \code{x} (e.g., the range of\code{x} / 10) will result in a 
#' logistic function.} 
#' \item{\code{alpha} very far from 0 will result in a linear function.}
#' }
#' @examples
#' x <- 1:100
#' y <- logisticFun(x, alpha = -10, b = 50)
#' plot(y ~ x, type = "l")
#' 
#' # The effect of alpha:
#' y1 <- logisticFun(x, alpha = -0.01, b = 50)
#' y2 <- logisticFun(x, alpha = -10, b = 50)
#' y3 <- logisticFun(x, alpha = -1000, b = 50)
#' 
#' par(mfrow = c(1, 3))
#' plot(y1 ~ x, type = "l", main = expression(alpha %->% 0))
#' plot(y2 ~ x, type = "l", main = expression(alpha %~~% range(x)/10))
#' plot(y3 ~ x, type = "l", main = expression(alpha %->% infinity))
logisticFun <- function(x, alpha, beta) {1 / (1 + exp((x - beta)/alpha))}

#' Quadratic function
#' 
#' @description A simple quadratic function of the form
#' \deqn{ax^2+bx+c}{
#' a*x^2+b*x+c}
#' @param x a numeric value or vector
#' @param a a numeric value or vector
#' @param b a numeric value or vector
#' @param c a numeric value or vector
#' @return a numeric value or vector resulting from the function
#' @export
#' @seealso \code{\link{linearFun}}, \code{\link{quadraticFun}}
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#'
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @examples
#' x <- 1:100
#' y <- quadraticFun(x, a = 2, b = 2, c = 3)
#' plot(y ~ x, type = "l")
quadraticFun <- function(x, a, b, c) {a * x^2 + b * x + c}


# Functions useful for the PCA approach

.f <- function(x, co) x %*% co

.pca.coordinates <- function(x, pca, na.rm)
{
  x <- sweep(x, 2L, pca$cent, check.margin=FALSE)
  x <- sweep(x, 2L, pca$norm, "/", check.margin=FALSE)
  x1 <- apply(x, 1, .f, co = pca$c1[, 1])
  x2 <- apply(x, 1, .f, co = pca$c1[, 2])
  return(cbind(x1, x2))
}

.prob.gaussian <- function(x, means, sds)
{
  dnorm(x[1], mean = means[1], sd = sds[1]) * dnorm(x[2], mean = means[2], sd = sds[2])
}


.thermalFun <- function(Pmax, Tb, To, rho, sigma)
{
  Pmax * exp(-exp(rho * (Tb - To) - 6) - sigma * (Tb - To)^2)
}
