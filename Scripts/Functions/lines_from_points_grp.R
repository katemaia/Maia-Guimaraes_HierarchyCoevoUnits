# -------------------- lines_from_points_grp -------------------

# Author: Kate P Maia
# Checked: 12/2023

# Creates a dataframe of lines for a dataframe of points, i.e., a dataframe with data to draw lines connecting points in the input dataframe.
# input: points dataframe with columns Code (equivalent to IntType), eff_measure (direct and indirect effects), m (0.1, 0.5, 0.9), group (Subg, Mod, Sct) and value (proportion of effect per interaction).
# output: lines dataframe with columns IntType, eff_measure, x, xend, y, yend which are the starting and ending points of each line to be drawn. 

# --------------------------------------------------------------

lines_from_points_grp <- function(points) {
  
  line1_start <- points %>% filter(group == "Subg")
  line1_end <- points %>% filter(group == "Mod")
  if (any(line1_start[, 1:3] != line1_end[, 1:3])) {print("ERROR IN LINE 1")}
  
  line2_start <- points %>% filter(group == "Mod")
  line2_end <- points %>% filter(group == "Sct")
  if (any(line2_start[, 1:3] != line2_end[, 1:3])) {print("ERROR IN LINE 2")}
  
  line1 <- cbind(line1_start, line1_end[, -c(1:3)])
  line2 <- cbind(line2_start, line2_end[, -c(1:3)])
  
  lines <- rbind(line1, line2)
  colnames(lines)[4:7] <- c("x", "y", "xend", "yend")
  
  return(lines)
}
