###  For converting from decalt coordinate to polar coordinate
##
#
library(data.table)
#library(tidyverse)
library(dplyr)

target.dir  <-  list.dirs("./")
id          <-  grep("CA-data", target.dir)
target.dir  <-  target.dir[id]

n.target    <-  sapply(strsplit(target.dir, "/CA"), function(x){x[1]})
n.target    <-  sapply(strsplit(n.target,   "//"),  function(x){x[2]})


for (t.type in n.target) {
  t.path         <-  paste("./", t.type, "/Pu-loc_2/", sep = "")
  pu.coord.list  <-  list.files(t.path, full.names = TRUE)
  
  file.num       <-  length(pu.coord.list)
  
  
  for (i in 1:file.num) {
    i.file  <-  fread(pu.coord.list[i])
    
    i.file                                           %>%
      mutate(r     = sqrt(birth^2 + death^2))        %>%
      mutate(Theta = atan((death - birth) / birth))  %>%
      mutate(V1    = V1 - .$V1[1] + 1)               %>%
      select(V1, r, Theta)                           -> o.file
    
    colnames(o.file)[1]  <-  "Step"
    
    o.name  <-  paste("./", t.type, "/Polar-loc_2/", formatC(i,width=4,flag="0"), "-point_polar-coordinate.csv", sep = "")
    
    write.csv(o.file, file = o.name, quote = FALSE, row.names = FALSE)
  }
}


