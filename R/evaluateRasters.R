# #' Compare two rasters to evaluate prediction accuracy
# #' 
# #' This function compare two rasters using continuous and/or binary metrics. 
# #' If needed, a binary conversion using thresholds will be made to turn suitability rasters into presence/absence rasters.
# #' 
# #' 
# #' @keywords internal
# #' @param obs a RasterLayer object, representing the Observed distribution used as a reference.
# #' @param pred a RasterLayer object, representing the Predicted distribution from a model output.
# #' @param metrics a vector containing a selection of following available metrics : \code{"RMSE"}, \code{"AUC"}, \code{"TSS"}, \code{"KAPPA"}, \code{"JACCARD"}, \code{"SORENSEN"} or \code{"all"} . 
# #' See details.
# #' @param binary \code{TRUE} or \code{FALSE}. Should a binary transformation be used for suitability rasters.
# #' @param thresholds a vector containing two numeric values. Will be used for binary tranformation of suitability rasters. 
# #' 
# #' @details
# #' This function can compute one continuous metric:
# #' \itemize{
# #' \item{\code{RMSE} : Root Main Squared Error (Caruana et Niculescu-Mizil, 2004). 
# #' The mean error between two suitability rasters.
# #' }}
# #' 
# #' This function can also compute 5 binary metrics:
# #' \itemize{
# #' \item{\code{AUC} : Area Under the Receiver Operating Characteristic (Mason et Graham, 2002).
# #' This metric compare a suitability prediction to a binary observation without the need of a threshold conversion.
# #' }
# #' \item{\code{TSS} : True Skill Statistic (Alouche et al., 2006). \code{TSS = specificity + sensitivity -1} 
# #' \itemize{
# #' \item{\code{Sensitivity} (Fielding and Bell, 1997). Proportion of presences which are correctly identified as such. 
# #' }
# #' \item{\code{Specificity} (Fielding and Bell, 1997). Proportion of absences which are correctly identified as such.
# #' }
# #' }}
# #' \item{\code{KAPPA} : Cohen's Kappa (Cohen, 1960). The inter-rater agreement for qualitative items.
# #' }
# #' \item{\code{JACCARD} : Jaccard similarity coefficient (Jaccard, 1908). Similarity between finite sample sets.
# #' \itemize{
# #' \item{\code{OPR} : Over Prediction Rate (Barbosa et al., 2013). Proportion of false presences predicted.
# #' }
# #' \item{\code{UTP} : Unpredicted True Presences (Fielding and Bell, 1997). Proportion of true presences missed.
# #' }
# #' }}
# #' \item{\code{SORENSEN} : Sorensen-Dice index, also known as F-measure (Daskallaki et al., 2006). Similarity between finite sample sets.
# #' }
# #' 
# #' \item{\code{all} : Compute all of the above binary metrics.
# #' }
# #' }
# #' @author
# #' Robin Delsol \email{robin.delsol@@wanadoo.fr}
# #' 
# #' Maintainer : Boris Leroy \email{leroy.boris@@gmail.com}
# #' 
# #' @return \code{data.frame} with 2 columns, "Metric" and "Value"
# #' 
# #' @examples
# #' # Create an example observed suitability raster
# #' a <- matrix(rep(dnorm(1:100, 40, sd = 25)), 
# #'             nrow = 100, ncol = 100, byrow = TRUE)
# #' obs <- raster(a* dnorm(1:100, 50, sd = 25))
# #' 
# #' # Create an example of a bad predicted suitability raster
# #' b1 <- matrix(rep(dnorm(1:100, 75, sd = 40)), 
# #'             nrow = 100, ncol = 100, byrow = TRUE)
# #' pred1 <- raster(b1* dnorm(1:100, 50, sd = 25))
# #' 
# #' # Create an example of a good predicted suitability raster
# #' b2 <- matrix(rep(dnorm(1:100, 45 , sd = 27)), 
# #'             nrow = 100, ncol = 100, byrow = TRUE)
# #' pred2 <- raster(b2* dnorm(1:100, 50, sd = 25))
# #' 
# #' par(mfrow=c(1,3))
# #' plot(obs, main="Reference raster")
# #' plot(pred1, main='"Bad" prediction')
# #' plot(pred2, main='"Good" prediction')
# #' 
# #' evaluateRasters(obs, pred1, binary=TRUE)
# #' evaluateRasters(obs, pred2, binary=TRUE, metrics="all")
# 
# 
# 
# evaluateRasters <- function(obs, pred, metrics = "all", binary = FALSE, thresholds = NULL)
# {
#   if(!(is(obs, "RasterLayer")))        {
#     stop("Obs must be a RasterLayer object.\n")
#   }
#   if(!(is(pred, "RasterLayer")))       {
#     stop("Pred must be a RasterLayer object.\n")
#   } 
#   if(!(ncell(obs) == ncell(pred)) )    {
#     stop("Obs and pred must have the same number of cells.\n")
#   }
#   if(!(canProcessInMemory(obs, n = 2))){
#     stop("Sorry, your computer does not have enough memory to extract all the values from the rasters.\n")
#   }
#   if ("all"  %in% metrics)             {
#     metrics=c("RMSE", "AUC", "TSS", "KAPPA", "JACCARD", "SORENSEN" )
#   }
#   # Condition min max a tester
#   type.obs  <- (obs*obs==obs)   @data@min
#   type.pred <- (pred*pred==pred)@data@min
#   
#   if(type.obs==0 & type.pred==0){ convert = 'both'}
#   if(type.obs==0 & type.pred==1){ convert = "obs" }
#   if(type.obs==1 & type.pred==0){ convert = "pred"}
#   if(type.obs==1 & type.pred==1){ convert = "none"; binary=T}
#   
#   evals <- matrix(ncol=2, nrow=0)
#   
#   
#   if (convert %in% c("pred", "obs") & binary==F){stop("Nothing computed\n")
#   }
#   
#   if (convert=="both"){
#     if ("RMSE"  %in% metrics){
#       obs.values  <- getValues(obs)
#       pred.values <- getValues(pred)
#       RMSE <- sqrt(mean((obs.values - pred.values )^2, na.rm=T))
#       evals <- rbind(evals, c("RMSE", RMSE))
#     }  
#   }
#   
#   if (binary==T){
#     
#     if (convert=="none") { 
#       
#       obs.values  <- getValues(obs)
#       pred.values <- getValues(pred)
#       
#       if ("AUC"   %in% metrics){
#         message(" -- Warning: Pred must be a suitability raster in order to have an AUC evaluation.\n")}
#     }    
#     if (convert=="pred") {
#       obs.values  <- getValues(obs)
#       pred.values <- getValues(pred)
#       
#       if ("AUC"   %in% metrics){
#         # roc1 <- pROC::roc(obs.values, pred.values) # annotated until the function is released
#         # AUC <- as.numeric(pROC::auc(roc1)) 
#         evals <- rbind(evals, c("AUC", AUC)) 
#       }  
#       
#       if (is.null(thresholds)){
#         thresholds  <- 0.6*pred@data@max
#         message(" -- Warning: You did not choose any threshold. Default values were be used for 'pred' binary conversion.\n")
#       }
#       
#       if (thresholds[1]>pred@data@max | thresholds[1]<pred@data@min ){
#         stop("Selected threshold is not in the good range of values.\n")
#       }
#       
#       if (length(thresholds)>1){
#         thresholds <- thresholds[1]
#         message(' -- Warning: thresholds length > 1. Only the first value was used for "pred" binary conversion.\n')
#       }
#       pred.values <- getValues(reclassify(pred, c(-Inf, thresholds[1], 0, thresholds[1], +Inf, 1) ))
#     }    
#     if (convert=="obs")  {
#       
#       if (is.null(thresholds)){
#         thresholds <- 0.6*obs@data@max
#         message(" -- Warning: You did not choose any threshold. Default values were be used for 'obs' binary conversion.\n")
#       }
#       
#       if (thresholds[1]>obs@data@max | thresholds[1]<obs@data@min ){
#         stop("Selected threshold is not in the good range of values.\n")
#       }
#       
#       if (length(thresholds)>1){
#         thresholds <- thresholds[1]
#         message(' -- Warning: thresholds length > 1. Only the first value was used for "obs" binary conversion.\n')
#       }
#       obs.values  <- getValues(reclassify(obs, c(-Inf, thresholds[1], 0, thresholds[1], +Inf, 1) ))
#       pred.values <- getValues(pred)
#       if ("AUC"   %in% metrics){
#         message(" -- Warning: Pred must be a suitability raster in order to have an AUC evaluation.\n")}
#     }    
#     if (convert=="both") {
#       
#       if (is.null(thresholds)){
#         thresholds=c(0.6*obs@data@max, 0.6*pred@data@max)
#         message(" -- Warning: You did not choose any threshold. Default values will be used for binary conversions.\n")
#       }
#       
#       if (length(thresholds)>2){
#         thresholds <- thresholds[1:2]
#         message(' -- Warning: thresholds length > 2. Only the first two values were used for binary conversions.\n')
#       }
#       if (length(thresholds)<2){
#         thresholds <- rep(thresholds,2)
#         message(' -- Warning: thresholds length = 1. The same values was used for both binary conversions.\n')
#       }
#       
#       if (thresholds[2]>pred@data@max | thresholds[2]<pred@data@min | thresholds[1]>obs@data@max | thresholds[1]<obs@data@min ){
#         stop("Selected threshold are not in the good range of values.\n")
#       }
#       
#       obs.values  <- getValues(reclassify(obs, c(-Inf, thresholds[1], 0, thresholds[1], +Inf, 1) ))
#       
#       if ("AUC"   %in% metrics){
#         pred.values <- getValues(pred)
#         # roc1 <- pROC::roc(obs.values, pred.values)
#         # AUC <- as.numeric(pROC::auc(roc1)) 
#         evals <- rbind(evals, c("AUC", AUC)) 
#       }
#       
#       pred.values <- getValues(reclassify(pred, c(-Inf, thresholds[2], 0, thresholds[2], +Inf, 1) ))
#     }
#     
#     diff <- obs.values - 2 * pred.values
#     TP  <-  sum(diff == -1, na.rm=T)
#     TN  <-  sum(diff ==  0, na.rm=T) 
#     FP  <-  sum(diff == -2, na.rm=T)
#     FN  <-  sum(diff ==  1, na.rm=T) 
#     
#     if ("TSS"     %in% metrics) {
#       sensitivity  <-  TP/(TP+FN)
#       specificity <-   TN/(TN+FP)
#       TSS <- (sensitivity + specificity -1)
#       evals <- rbind(evals, 
#                      c("Specificity", specificity),
#                      c("Sensitivity", sensitivity),
#                      c("TSS", TSS))
#     }
#     if ("KAPPA"   %in% metrics) {
#       x = TP - (TP+FP)*(TP+FN) / (TP+TN+FP+FN)
#       Kappa = x / ( x + FP * (1/2) + FN*(1/2) )
#       evals <- rbind(evals, 
#                      c("Kappa", Kappa))}   
#     if ("JACCARD" %in% metrics) {
#       Jaccard  =    TP / (  TP+FN+FP) 
#       OPR  =  FP / ( TP+FP )
#       UTP  =  FN / ( TP+FN ) 
#       evals <- rbind(evals,
#                      c("Jaccard", Jaccard), 
#                      c("OPR", OPR),
#                      c("UTP", UTP))}
#     if ("SORENSEN" %in% metrics){
#       Sorensen =  2*TP / (2*TP+FN+FP) 
#       evals <- rbind(evals, 
#                      c("Sorensen", Sorensen))
#     } 
#   } 
#   
#   evals[,2] <- round(as.numeric(evals[,2]),5)
#   evals <- as.data.frame(evals)
#   colnames(evals) <- c("Metric","Value")
#   
#   return(evals)
# }
