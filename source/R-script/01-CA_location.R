### calculate TDA characters from pdb data
##
#
library(bio3d)
library(TDA)
library(readr)
library(data.table)


### search the path and number of input data
f.list  <-  list.files("./", full.names = TRUE, recursive = TRUE)
pdb.id  <-  grep(".pdb", f.list)
f.list  <-  f.list[pdb.id]
sep.id  <-  grep("Separate_", f.list)
sep.id  <-  grep("CA-data", f.list)

f.list  <-  f.list[sep.id]


for (ugeid in 1:length(f.list)) {
  ### read data and calculate
  f.name    <-  f.list[ugeid]
  pdb.conf  <-  read.pdb2(f.name)
  print(paste0("Now :  ",  f.name))
  
  num.atom  <-  nrow(pdb.conf$atom)
  conf.mat  <-  matrix(NA, nrow = num.atom, ncol = 3)
  
  conf.mat[, 1]  <-  pdb.conf$atom$x
  conf.mat[, 2]  <-  pdb.conf$atom$y
  conf.mat[, 3]  <-  pdb.conf$atom$z
  
  
  if (file.exists("Pickup_Atoms.csv") == FALSE) {
    pu.id   <-  seq(1, nrow(conf.mat))
    write.table(pu.id, file = "Pickup_Atoms.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  pu.id     <-  read_table("Pickup_Atoms.csv", col_names = FALSE)
  
  conf.mat  <-  conf.mat[c(pu.id$X1), ]
  
  #maxs      <-  quantile(dist(conf.mat), probs=0.5)
  maxs      <-  20
  diag.cir  <-  ripsDiag(conf.mat, maxdimension = 2, maxscale = maxs, 
                         library = "Dionysus", location = TRUE, printProgress = TRUE)
  
  
  
  out.num    <-  sapply(strsplit(f.name, "\\.p"), function(x){x[1]})
  dir.name   <-  sapply(strsplit(out.num, "Sep"), function(x){x[1]})
  out.num    <-  sapply(strsplit(out.num, "_"),   function(x){x[2]})
  bd.plots   <-  paste(dir.name, "Separate_", out.num, ".csv", sep = "")
  bd.data    <-  data.frame("dimension"  =  diag.cir$diagram[, 1],
                            "Birth"      =  diag.cir$diagram[, 2],
                            "Death"      =  diag.cir$diagram[, 3])
  write_csv(bd.data, file = bd.plots)
  
  
  all_points       <-  apply(conf.mat, 1, paste, collapse="_")
  id               <-  which(diag.cir[["diagram"]][, 1] >= 1)
  
  location_points  <-  vector("list",length(id))
  
  for(i in 1:length(id)){
    xx                    <-  diag.cir[["cycleLocation"]][[id[i]]]
    location_points[[i]]  <-  sort(match(unique(as.vector(apply(xx, c(1,2), paste, collapse="_"))), all_points))
  }
  
  mrd.out    <-  data.frame(diag.cir[["diagram"]][id,],location=sapply(location_points,paste,collapse=";"))
  
  loc.data   <-  data.frame("Dimension"  =  mrd.out[, 1],
                            "Birth"      =  mrd.out[, 2],
                            "Death"      =  mrd.out[, 3],
                            "Location"   =  mrd.out[, 4])
  
  loc.name   <-  paste(dir.name, "Location-Separate_", out.num, ".csv", sep = "")
  
  write_csv(loc.data, file = loc.name)
  
  
  
  
  out.name   <-  paste(dir.name, "Persistent_homology-Separate_", out.num, ".pdf", sep = "")
  
  pdf(out.name, width = 960/96, height = 600/96)
    par(las = 1)
    plot(x = NULL,   type = "n", xlab = "Birth",  ylab = "Death", xlim = c(0, 20), ylim = c(0, 20),
         xaxs = "i", yaxs = "i", xaxt = "n",      yaxt = "n",     bty = "n", 
         main = "Persistent_homology")
    axis(side = 1, at = seq(0, 20, 5), tck = 1.0, lty = "blank", 
         labels = expression(0, 5, 10, 15, 20))
    axis(side = 2, at = seq(0, 20, 5), tck = 1.0, lty = "blank", 
         labels = expression(0, 5, 10, 15, 20))
  
  
    lines(x = 1:20,   y = 1:20)
    points(x = diag.cir$diagram[which(diag.cir$diagram[, 1] == 0), 2],   
           diag.cir$diagram[which(diag.cir$diagram[, 1] == 0), 3], pch = 16, col = "black")
    points(x = diag.cir$diagram[which(diag.cir$diagram[, 1] == 1), 2],   
           diag.cir$diagram[which(diag.cir$diagram[, 1] == 1), 3], pch = 16, col = "red")
    points(x = diag.cir$diagram[which(diag.cir$diagram[, 1] == 2), 2],   
           diag.cir$diagram[which(diag.cir$diagram[, 1] == 2), 3], pch = 16, col = "blue")
  
    legend("bottomright", legend = c("0D", "1D", "2D"), col = c("black", "red", "blue"), pch = 16)
    box(bty = "l")
  dev.off()

}  

