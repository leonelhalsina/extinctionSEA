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

II <- i
num_max_multiregion <- 3
DEC_events_lemad <- TRUE
pseuduvaria_tree <- read.nexus("A10_36_finaltree.out.tre")

outgroups_to_remove <- c("Orophea_enterocarpa",
                         "Orophea_celebica",
                         "Alphonsea_kinabaluensis",
                         "Mitrephora_keithii",
                         "Polyalthia_suberosa",
                         "Haplostichanthus_longirostris",
                         "Miliusa_campanulata",
                         "Miliusa_horsfieldii",
                         "Polyalthia_korinti",
                         "Sapranthus_viridiflorus",
                         "Monocarpia_euneura",
                         "Polyalthia_coffeoides",
                         "Polyalthia_longifolia")


pseuduvaria_tree <- drop.tip(pseuduvaria_tree,outgroups_to_remove)

pseuduvaria_distribution <- read.csv("distribution_pseu.csv")
pseuduvaria_distribution <- pseuduvaria_distribution[-c(3,4),] # outgroups

all_areas <- c("A","B","C","D","E","F","G","H")




## Homogenize locations over studies###



uniform_coding <- NULL 

for(i in 1:nrow(pseuduvaria_distribution)){
  uniform_coding_onespp <- NULL 
  
  take_this_spp <- as.character(pseuduvaria_distribution[,2][i])
  
  for(j in 1:nchar(take_this_spp)){
    one_location <-  strsplit(take_this_spp,split="")[[1]][j]
    
    if(one_location ==  "C"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"a") # Borneo
    }
    if(one_location ==  "E"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"b") # Sulawesi
    }
    # if(one_location==  ""){
    #   uniform_coding_onespp <- paste0(uniform_coding_onespp,"c") # Sumatra
    # }
    if(one_location ==  "B"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"d") # Java
    }
      if(one_location ==  "A"){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"e") # Main
      }
    # this also includes Sumatra
    
    if(one_location ==  "D" ){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"f") # Philippines
    }
    if(one_location ==  "F" || one_location ==  "G"   ){
      uniform_coding_onespp <- paste0(uniform_coding_onespp,"g") # New Guinea 
    }
  }
  uniform_coding <- c(uniform_coding,uniform_coding_onespp)
  
}


pseuduvaria_distribution <- as.data.frame(cbind(pseuduvaria_distribution,uniform_coding))

all_areas_uniform_coding <- c("a","b","d","e","f","g")


#### The lemad part 
phylotree_recons <- pseuduvaria_tree

# species locations should be in the same order than tree tips
species_presence <- pseuduvaria_distribution[match(phylotree_recons$tip.label,pseuduvaria_distribution[,1]),]
species_presence <- as.character(species_presence[,3])

# we need to order the letters for each species:
for(ik in 1:length(species_presence)){
  string_to_sort <- species_presence[ik]
  species_presence[ik] <- paste(sort(unlist(strsplit(string_to_sort, ""))), collapse = "")
  
}
bd_estimates <- bd_ML(branching.times(phylotree_recons))
initial_lambda <- c(bd_estimates$lambda0/2,bd_estimates$lambda0/2)
initial_disperextirpation <- bd_estimates$lambda0/5


all_try_this_extinction <- c(0,0.025,0.25,2.5)
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
  condition_on_origin = NULL,
  DEC_events = DEC_events_lemad,
  missing_spp_areas = NULL,
  lineage_extinction = try_this_extinction,
  initial_lambda = initial_lambda,
  initial_disperextirpation)

ending_time <- Sys.time()
cat("time:",ending_time - starting_time)
print(starting_time)
print(ending_time)

saveRDS(output,file=paste0("lemad_ext_",try_this_extinction,".RDS"))

}
