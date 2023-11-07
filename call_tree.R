do_things <- function(i){
library(secsse)
library(ape)
library(geiger)
library(BioGeoBEARS)
library(ape)
library(phytools)
library(stringr)
library(deSolve)
library(apTreeshape)
library(phytools)
library(Matrix)
library(lemad)
library(DDD)
library(phylocanvas)
library(phylobase)
library(nodiv)
library(doParallel)
library(foreach)
library(doMC)
source("function_divalike.R")
II <- i
num_max_multiregion <- 3
DEC_events_lemad <- TRUE


gryllo_tree <- read.tree("CardioBEAST.newick")
gryllo_distribution <- read.csv("distribution.csv")



rownames(gryllo_distribution) <- gryllo_distribution[,1]

name.check(gryllo_tree, gryllo_distribution, data.names=NULL)
gryllo_distribution[,1] == gryllo_tree$tip.label




all_areas <- c("A","B","C","D","E","F","G","H","I","J","K","L")


## Homogenize locations over studies###



uniform_coding <- NULL 

for(i in 1:nrow(gryllo_distribution)){
  uniform_coding_onespp <- NULL 
  
  take_this_spp <- as.character(gryllo_distribution[,2][i])
  
  for(j in 1:nchar(take_this_spp)){
    one_location <-  strsplit(take_this_spp,split="")[[1]][j]
    
    if(one_location ==  "F" || one_location ==  "E"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"a") # Borneo
    }
    if(one_location ==  "H" || one_location ==  "J"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"b") # Sulawesi
    }
    if(one_location ==  "C"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"c") # Sumatra
    }
    if(one_location ==  "G" || one_location ==  "I"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"d") # Java
    }
    if(one_location ==  "A" || one_location ==  "B" ){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"e") # Main
    }
    if(one_location ==  "D"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"f") # Philippines
    }
    if(one_location ==  "K"  ||one_location ==  "L"  ){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"g") # New Guinea 
    }
  }
  uniform_coding <- c(uniform_coding,uniform_coding_onespp)
  
}

uniform_coding[which(uniform_coding=="gg")] <- "g"
as.character(gryllo_distribution[,2])

gryllo_distribution <- as.data.frame(cbind(gryllo_distribution,uniform_coding=uniform_coding))
all_areas_uniform_coding <- c("a","b","c","d","e","f","g")

####

#divalike_analysis <- do_divalike(areas = all_areas_uniform_coding,
# #                                phylotree_recons = gryllo_tree,
#                                 traitdata = as.character(gryllo_distribution[,3]),
#                                 num_max_multiregion = num_max_multiregion,
#                                 do_plot = FALSE )
#saveRDS(divalike_analysis,file="divalike.RDS")

#### The lemad part 
phylotree_recons <- gryllo_tree

# species locations should be in the same order than tree tips
species_presence <- gryllo_distribution[match(phylotree_recons$tip.label,gryllo_distribution[,1]),]
species_presence <- as.character(species_presence[,3])

# we need to order the letters for each species:
for(ik in 1:length(species_presence)){
  string_to_sort <- species_presence[ik]
  species_presence[ik] <- paste(sort(unlist(strsplit(string_to_sort, ""))), collapse = "")
  
}
bd_estimates <- bd_ML(branching.times(phylotree_recons))
initial_lambda <- c(bd_estimates$lambda0/2,bd_estimates$lambda0/2)
initial_disperextirpation <- bd_estimates$lambda0/5

all_try_this_extinction <- c(0,0.0057,0.057,0.57)
try_this_extinction <- all_try_this_extinction[II]
if(is.na(try_this_extinction)){
  try_this_extinction <- "free"
}

cat("________doing this extinction: ",try_this_extinction,"_________\n")
cat("________max number of regions ",num_max_multiregion,"_________\n")
  starting_time <- Sys.time()
output <- lemad_analysis(
  phylotree_recons,
  species_presence,
  areas = all_areas_uniform_coding,
  num_max_multiregion = num_max_multiregion,
  DEC_events = DEC_events_lemad,
  missing_spp_areas = NULL,
  lineage_extinction = try_this_extinction,
  initial_lambda = initial_lambda,
  initial_disperextirpation,
  run_parallel = TRUE,
  use_fortran_code = TRUE)
  ending_time <- Sys.time()
  cat("time:",ending_time - starting_time)
  ending_time <- Sys.time()
  cat("time:",ending_time - starting_time)
  print(starting_time)
  print(ending_time)

saveRDS(output,file=paste0("lemad_ext_",try_this_extinction,".RDS"))
}

