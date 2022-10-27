### Estimating Deviation for each points
##
#
# load library
library(data.table)
library(tidyverse)
library(kernlab)


target.dir  <-  list.dirs("./")
id          <-  grep("CA-data", target.dir)
target.dir  <-  target.dir[id]

np.target   <-  sapply(strsplit(target.dir, "/CA"), function(x){x[1]})
np.target   <-  sapply(strsplit(np.target,  "//"),  function(x){x[2]})


for (n.target in np.target) {
  # load file list
  mut.list  <-  list.files(paste0("./", n.target), full.names = TRUE, recursive = TRUE)
  mut.list  <-  mut.list[grep("Location", mut.list)]
  
  mut.location   <-  mut.list[grep("Location-Separate", mut.list)]
  mut.locatref   <-  mut.list[grep("Ref_location", mut.list)]
  
  mut.location   <-  mut.location[order(nchar(mut.location))]
  mut.locatref   <-  mut.locatref[order(nchar(mut.locatref))]
  
  num.location  <-  length(mut.location)
  num.locatref  <-  length(mut.locatref)
  
  
  rm(mut.list)
  
  
  
  ### Coordinate pickup
  ##
  ref.loc  <-  mut.locatref[grep("_1st", mut.locatref)][grep("detect/", mut.locatref[grep("_1st", mut.locatref)])]
  
  ref.num  <-  fread(ref.loc)
  num.loc  <-  nrow(ref.num)
  num.tms  <-  ncol(ref.num)
  
  tc.coord       <-  as.data.frame(matrix(data = NA, nrow = (num.loc * num.tms), ncol = 3))
  
  
  colnames(tc.coord)  <-  c("point No.", "birth", "death")
  
  for (i in 1:num.loc) {
    sp  <-  1 + (i - 1) * num.tms
    ep  <-  i * num.tms
    tc.coord[sp:ep, 1]  <-  i
  }
  
  
  for (j in 1:num.loc) {
    for (i in 1:num.tms) {
      t.num    <-  i
      t.sch    <-  paste(t.num, ".csv", sep = "")
      t.dat    <-  mut.location[grep(t.sch, mut.location)][1]
      t.coord  <-  fread(t.dat)
      
      ID  <-  as.numeric(ref.num[j, i:i])
      
      if (is.na(ID) == FALSE) {
        tc.coord[(i + (j - 1) * num.tms), 2]  <-  t.coord[ID, 2]
        tc.coord[(i + (j - 1) * num.tms), 3]  <-  t.coord[ID, 3]
      }
      
    }
    

    if (j >= 1000) {
      j.num  <-  as.character(j)
    } else if (j >= 100) {
      j.num  <-  paste("0",   as.character(j), sep = "")
    } else if (j >= 10) {
      j.num  <-  paste("00",  as.character(j), sep = "")
    } else {
      j.num  <-  paste("000", as.character(j), sep = "")
    }
    
    out       <-  subset(tc.coord, tc.coord[, 1] == j)
    out.name  <-  paste("./", n.target, "/Pu-loc_2/", j.num, "-point_coutcrse-coordinate.csv", sep = "")
    write.csv(subset(out, is.na(out[, 2]) == FALSE), out.name)
    
  }
  
}

