# ---------------------- motifint_counter ----------------------

# Author: Kate P Maia
# Checked: 12/2023

# motifint_counter: SAME PURPOSE OF groupint_counter FOR SUBGRAPHS. 
# Counts the number of interactions inside the network's largest component, and inside network subgraphs (all combined since motifs 6 have connectance = 1 and links can take part in multiple motifs). 
# input: prf (motif profile of a given network, see Script 7) and M_GC (interaction matrix of the network's largest component).
# output: a vector with the number of interactions inside the network's largest component (All) and inside all motifs (Mot).

# --------------------------------------------------------------

motifint_counter <- function(prf, M_GC){ 
  
  # only interactions in prfint (RC, not RR/CC), reverse would double count interactions (RC and not CR)
  colnames(prf) <- rep("name", 4)
  prfint <- rbind(prf[,c(1,3)], prf[,c(1,4)], prf[,c(2,3)], prf[,c(2,4)])
  motifint <- distinct(prfint)
  if (any(motifint != prfint[!duplicated(prfint),])) {"ERROR IN DUPLICATED"}
  motifint <- c(sum(M_GC), nrow(motifint))
  names(motifint) <- c("All", "Mot")
  return(motifint)
  
}