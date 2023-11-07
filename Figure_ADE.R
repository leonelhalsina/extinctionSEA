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


best_model <- readRDS("lemad_ext_2regions0.0057.RDS")
num_max_multiregion <- 3
output <- best_model
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
probabilities_at_node <- output$ancestral_states

all_states <- give_me_states_combination(all_areas_uniform_coding,num_max_multiregion)
highest_prob_state_lemad <- NULL
for(ij in 1:nrow(probabilities_at_node)){
  
  id_order_state_highest <- which(max(probabilities_at_node[ij,])==probabilities_at_node[ij,])[1]
  # cat(id_order_state_highest,"\n")
  highest_prob_state_lemad <- c(highest_prob_state_lemad,all_states[id_order_state_highest])
}


tip_state <- species_presence




##### trying using the pie argument, to make easy the plotting of multi-area lineages

tip_and_internal_nodes <- c(tip_state,highest_prob_state_lemad)
matrix_colors <- matrix(0,nrow = length(tip_and_internal_nodes),
                        ncol = length(all_areas_uniform_coding))
colnames(matrix_colors) <- all_areas_uniform_coding
for(ij in 1:length(tip_and_internal_nodes)){
      this_tip_state <- tip_and_internal_nodes[ij]
    for(IJ in 1:nchar(this_tip_state)){
    number_to_accomodate <- 1/nchar(this_tip_state)
    
    find_this_state <- substr(this_tip_state,IJ,IJ)
   this_column <- which(colnames(matrix_colors) == find_this_state)
   matrix_colors[ij,this_column] <- number_to_accomodate
  }
 
}

color_borneo <- "red"
color_sulawesi <- "dimgray"
color_sumatra <- "blue"
color_java <- "darkgoldenrod"
color_main <- "seagreen3"
color_phili <- "black"
color_papua <- "darkred"


colors_for_pie <-  c(color_borneo,color_sulawesi,color_sumatra,color_java,color_main,color_phili,color_papua)

plot.phylo(main=best_model$estimated_rates[4],phylotree_recons,show.tip.label = FALSE,cex=0.6,label.offset=0.7)
nodes_matrix_colors <- matrix_colors[(length(tip_state)+1):nrow(matrix_colors),]
nodelabels(pie = nodes_matrix_colors , 
           piecol=colors_for_pie ,
           cex=0.6)
tips_matrix_colors <- matrix_colors[1:length(tip_state),]
tiplabels(pie = tips_matrix_colors, 
           piecol= colors_for_pie,
           cex=0.4)

axisPhylo()


