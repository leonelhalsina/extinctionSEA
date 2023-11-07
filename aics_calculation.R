#rm(list=ls(all=TRUE))
library(AICcmodavg)


ready_files<-list.files(pattern = "lemad_ext_")
models_to_compare<-ready_files

Aics<-NULL

loglikelihoods<-NULL
free_parameters <- NULL
for(i in 1:length(models_to_compare)){
  cat("loading:",paste0(models_to_compare[[i]]),"\n")
  mod <- readRDS(paste0(models_to_compare[[i]]))
  ll <- mod$model_ml
  free_pars <- mod$number_free_pars
  free_parameters<-c(free_parameters,free_pars)
  loglikelihoods <- c(loglikelihoods,ll)
  Aics <- c(Aics,AICcCustom(ll,K=free_pars,second.ord = F)) 
}

ICbweights2 <- function(IC){
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}


calculated_weights<-ICbweights2(Aics)

table_calculated_weights <- data.frame(models_to_compare,loglikelihoods,free_parameters,calculated_weights)

class(table_calculated_weights)
table_calculated_weights<-table_calculated_weights[order(table_calculated_weights$calculated_weights,decreasing = TRUE),]



write.csv(table_calculated_weights,"table_aic.csv")

best_model <- readRDS(paste0(table_calculated_weights[1,1]))



best_model$estimated_rates
