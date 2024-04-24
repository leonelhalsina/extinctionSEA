do_things <- function(i){
#library(secsse)
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
orchid_tree <-
  read.nexus("Paph-ITS+trnL+F+atpB-com-new-exclude gap.trees")

# remove outgroups
orchid_tree <- drop.tip(orchid_tree,c("Phragmipedium_longifolium","Phragmipedium_pearcei" ))
## Revell's way of re-scaling
scale_to_this <- 7.9
orchid_tree$edge.length <-
  orchid_tree$edge.length/max(nodeHeights(orchid_tree)[,2]) * scale_to_this
##

orchid_distribution <- read.csv("distrib.csv")
orchid_distribution <- orchid_distribution[-c(79,80),] # remove outgroups

rownames(orchid_distribution) <- orchid_distribution[, 1]

name.check(orchid_tree, orchid_distribution, data.names = NULL)


orchid_distribution <-
  orchid_distribution[order(match(orchid_distribution[, 1], orchid_tree$tip.label)), ]

rownames(orchid_distribution) == orchid_tree$tip.label
table(orchid_distribution[, 2])

all_areas <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")

## Homogenize locations over studies###



uniform_coding <- NULL

for (i in 1:nrow(orchid_distribution)) {
  uniform_coding_onespp <- NULL
  
  take_this_spp <- as.character(orchid_distribution[, 2][i])
  
  for (j in 1:nchar(take_this_spp)) {
    one_location <-  strsplit(take_this_spp, split = "")[[1]][j]
    
    if (one_location ==  "B" || one_location ==  "G" ||  one_location ==  "H") {
      uniform_coding_onespp <- paste0(uniform_coding_onespp, "a") # Borneo
    }

    # if (one_location ==  "B") {
    #   uniform_coding_onespp <- paste0(uniform_coding_onespp, "c") # Sumatra
    # }
    # 
    # if (one_location ==  "B") {
    #   uniform_coding_onespp <- paste0(uniform_coding_onespp, "d") # Java
    # }
    
    
    # but also, SUmatra and Java are considered a single region by the authors
    # I am using BORNEO as a code because it is the largest island 
    if (one_location ==  "C") {
      uniform_coding_onespp <-
        paste0(uniform_coding_onespp, "b") # Sulawesi
    }
    
    if (one_location ==  "A") {
      uniform_coding_onespp <- paste0(uniform_coding_onespp, "e") # Main
    }
    if (one_location ==  "E") {
      uniform_coding_onespp <-
        paste0(uniform_coding_onespp, "f") # Philippines
    }
    if (one_location ==  "D" || one_location ==  "F") {
      uniform_coding_onespp <-
        paste0(uniform_coding_onespp, "g") # New Guinea
    }

    # for some letters like G, H, F, I had a look at their figure 5 to sort it out
  }
  uniform_coding <- c(uniform_coding, uniform_coding_onespp)
  
}
uniform_coding[which(uniform_coding == "aa")] <- "a"
uniform_coding[which(uniform_coding == "abag")] <- "abg"

orchid_distribution <-
  as.data.frame(cbind(orchid_distribution, uniform_coding))


all_areas_uniform_coding <- c("a","b","e","f","g")
##

phylotree_recons <- orchid_tree

# species locations should be in the same order than tree tips
species_presence <- orchid_distribution[match(phylotree_recons$tip.label,orchid_distribution[,1]),]
species_presence <- as.character(species_presence[,3])

# we need to order the letters for each species:
for(ik in 1:length(species_presence)){
  string_to_sort <- species_presence[ik]
  species_presence[ik] <- paste(sort(unlist(strsplit(string_to_sort, ""))), collapse = "")
  
}
bd_estimates <- bd_ML(branching.times(phylotree_recons))
initial_lambda <- c(bd_estimates$lambda0/2,bd_estimates$lambda0/2)
initial_disperextirpation <- bd_estimates$lambda0/5

all_try_this_extinction <- c(0,0.054,0.54,5.4)
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
  initial_disperextirpation,
  )
ending_time <- Sys.time()
cat("time:",ending_time - starting_time)
print(starting_time)
print(ending_time)

saveRDS(output,file=paste0("lemad_ext_",try_this_extinction,".RDS"))

}

