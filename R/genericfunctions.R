#' @export
#' @method print virtualspecies
print.virtualspecies <- function(x, ...)
{
  cat(paste("Virtual species generated from", 
            length(x$details$variables), 
            "variables:\n",
            paste(x$details$variables, collapse = ", ")))
  cat("\n\n- Approach used:")
  if(x$approach == "response")
  {
    cat(" Responses to each variable")

    cat("\n- Response functions:")
    sapply(x$details$variables, FUN = function(y)
    {
      cat("\n   .", y, 
          "  [min=", x$details$parameters[[y]]$min, 
          "; max=", x$details$parameters[[y]]$max,
          "] : ", 
          x$details$parameters[[y]]$fun,
          "   (", 
          paste(names(x$details$parameters[[y]]$args), 
                x$details$parameters[[y]]$args, sep = '=', collapse = "; "),
          ")", sep = "")
    })
    
    if (x$details$rescale.each.response)
    {
      cat("\n- Each response function was rescaled between 0 and 1")
    } else
    {
      cat("\n- Response functions were not rescaled between 0 and 1")
    }
    
    cat("\n- Environmental suitability formula = ", x$details$formula, sep = "")
    
    if (x$details$rescale)
    {
      cat("\n- Environmental suitability was rescaled between 0 and 1\n")
    } else
    {
      cat("\n- Environmental suitability was not rescaled between 0 and 1\n")
    }
    
  } else if(x$approach == "pca")
  {
    cat(" Response to axes of a PCA")
    cat("\n- Axes: ", 
        paste(x$details$axes, collapse = ", "),
        "; ", round(sum(x$details$pca$eig[x$details$axes])/sum(x$details$pca$eig) * 100, 2),
        "% explained by these axes")
    cat("\n- Responses to axes:")
    sapply(1:length(x$details$axes), function(y)
    {
      cat("\n   .Axis ", x$details$axes[y],
          "  [min=", round(min(x$details$pca$li[, x$details$axes[y]]), 2),
          "; max=", round(max(x$details$pca$li[, x$details$axes[y]]), 2),
          "] : dnorm (mean=", x$details$means[y], "; sd=", x$details$sds[y],
          ")", sep = "")
    })
    if (x$details$rescale)
    {
      cat("\n- Environmental suitability was rescaled between 0 and 1\n")
    } else
    {
      cat("\n- Environmental suitability was not rescaled between 0 and 1\n")
    }
  } else if(x$approach == "bca"){
    cat(" Response to axes of a BCA")
    cat("\n- Axes: ", 
        paste(x$details$axes, collapse = " & "),
        "; ", round(sum(x$details$bca$eig[x$details$axes])/sum(x$details$bca$eig) * 100, 2),
        "% explained by these axes")
    cat("\n- Responses to axes:")
    sapply(1:length(x$details$axes), function(y)
    {
      cat("\n   .Axis ", x$details$axes[y],
          "  [min=", round(min(x$details$bca$ls[, x$details$axes[y]]), 2),
          "; max=", round(max(x$details$bca$ls[, x$details$axes[y]]), 2),
          "] : dnorm (mean=", x$details$means[y], "; sd=", x$details$sds[y],
          ")", sep = "")
    })
    if (x$details$rescale)
    {
      cat("\n- Environmental suitability was rescaled between 0 and 1")
    } else
    {
      cat("\n- Environmental suitability was not rescaled between 0 and 1")
    }
  }
  if(!is.null(x$PA.conversion))
  {
    cat("\n- Converted into presence-absence:")
    cat("\n   .Method =", x$PA.conversion["conversion.method"])
    if(x$PA.conversion["conversion.method"] == "probability")
    {
      cat("\n   .alpha (slope)           =", x$PA.conversion["alpha"])
      cat("\n   .beta  (inflexion point) =", x$PA.conversion["beta"])
      cat("\n   .species prevalence      =", x$PA.conversion["species.prevalence"])
    } else if(x$PA.conversion["conversion.method"] == "threshold")
    {
      cat("\n   .threshold           =", x$PA.conversion["cutoff"])
      cat("\n   .species prevalence  =", x$PA.conversion["species.prevalence"], "\n")
    }
  }
  if(!is.null(x$occupied.area))
  {
    if(!is.null(x$geographical.limit))
    {
      cat("\n- Distribution bias introduced:")
      cat("\n   .method used :", x$geographical.limit$method)
      if(x$geographical.limit$method %in% c("country", "region", "continent"))
      {
        cat("\n   .area(s)     :", x$geographical.limit$area, "\n")
      } else if(x$geographical.limit$method == "extent")
      {
        cat("\n   .extent      : [Xmin; Xmax] = [", 
            x$geographical.limit$extent@xmin, "; ",
            x$geographical.limit$extent@xmax, "] - [Ymin; Ymax] = [",
            x$geographical.limit$extent@ymin, "; ",
            x$geographical.limit$extent@ymax, "]", "\n", sep = "")
      } else if(x$geographical.limit$method == "polygon")
      {
        cat("\n   .polygon    : Object of class ", class(x$geographical.limit$area), "\n", sep = "")
      }
    }
  }
}

#' @export
#' @method str virtualspecies
str.virtualspecies <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}

#' @export
#' @method plot virtualspecies
plot.virtualspecies <- function(x, ...)
{
  y <- raster::stack(x$suitab.raster)
  names(y) <- "Suitability.raster"
  if(!is.null(x$pa.raster))
  {
    y <- stack(y,
               x$pa.raster)
    names(y)[[nlayers(y)]] <- "Presence.absence.raster"
  }
  if(!is.null(x$occupied.area))
  {
    y <- stack(y,
               x$occupied.area)
    names(y)[[nlayers(y)]] <- "Occupied.area.raster"
  }
  x <- y
  plot(x, ...)
}

#' @export
#' @method print VSSampledPoints
print.VSSampledPoints <- function(x, ...)
{
  # Next line is to ensure retrocompatibility with earlier versions of
  # virtualspecies where no print function was designed for VSSampledPoints
  if(!is.list(x$detection.probability))
  {
    print(x)
  } else
  {
    cat(paste("Occurrence points sampled from a virtual species"))
    cat(paste("\n\n- Type:", x$type))
    cat(paste("\n- Number of points:", nrow(x$sample.points)))
    if(length(x$bias))
    {
      cat("\n- Sampling bias: ")
      cat(paste("\n   .Bias type:", 
                x$bias$bias))
      cat(paste("\n   .Bias strength:",
                x$bias$bias.strength))
    } else
    {
      cat("\n- No sampling bias")
    }
    cat(paste0("\n- Detection probability: "))
    cat(paste0("\n   .Probability: ", x$detection.probability$detection.probability))
    cat(paste0("\n   .Corrected by suitability: ", x$detection.probability$correct.by.suitability))
    cat(paste0("\n- Probability of identification error (false positive): ", x$error.probability))
    if(length(x$sample.prevalence))
    {
      cat(paste0("\n- Sample prevalence: "))
      cat(paste0("\n   .True:", x$sample.prevalence["true.sample.prevalence"]))
      cat(paste0("\n   .Observed:", x$sample.prevalence["observed.sample.prevalence"]))
    }
    cat(paste0("\n- Multiple samples can occur in a single cell: ", 
               ifelse(x$replacement, "Yes", "No")))
    cat("\n\n")
    if(nrow(x$sample.points) > 10)
    {
      cat("First 10 lines: \n")
      print(x$sample.points[1:10, ])
      cat(paste0("... ", nrow(x$sample.points) - 10, " more lines."))
    } else
    {
      print(x$sample.points)
    }
  }
}

#' @export
#' @method str VSSampledPoints
str.VSSampledPoints <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}