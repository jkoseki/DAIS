### detection of the position
##
#
#library(tidyr)
library(data.table)
library(stringr)

#set your  contl
n.contl  <-  "WT"


case.dir  <-  list.dirs("./")
id        <-  grep("CA-data", case.dir)
case.dir  <-  case.dir[id]
id        <-  grep(paste0("//", n.contl, "/CA-data"), case.dir)
case.dir  <-  case.dir[-id]

np.case   <-  sapply(strsplit(case.dir, "/CA"), function(x){x[1]})
np.case   <-  sapply(strsplit(np.case,  "//"),  function(x){x[2]})


for (n.case in np.case) {
  file.list.mut     <-  list.files(paste0("./", n.case, "/CA-data/"), full.names = TRUE, recursive = TRUE)
  id                <-  grep("Location", file.list.mut)
  file.list.mut     <-  file.list.mut[id]
  num.files         <-  length(file.list.mut)
  
  
  ## separation source for SHIROKANE
  f.data      <-  fread(file.list.mut[1])
  loc.dim     <-  f.data$Dimension
  loc.number  <-  str_split(f.data$Location[], pattern = ";")
  
  num.f.data    <-  nrow(f.data)
  
  dim.1.range <-  which(loc.dim == 1)
  dim.2.range <-  which(loc.dim == 2)
  num.dim.1   <-  length(dim.1.range)
  
  common.loc  <-  data.frame("XXX" = c(1:num.f.data))
  colnames(common.loc)  <-  n.case
  
  
  file.list.wt       <-  list.files(paste0("./", n.contl, "/CA-data/"), full.names = TRUE, recursive = TRUE)
  id                 <-  grep("Location", file.list.wt)
  file.list.wt       <-  file.list.wt[id]
  
  s.data             <-  fread(file.list.wt[1])
  s.loc.dim          <-  s.data$Dimension
  s.loc.num          <-  str_split(s.data$Location[], pattern = ";") 
  
  ref.loc            <-  data.frame("n.contl" = c(1:num.f.data))
  colnames(ref.loc)  <-  n.contl
  
  ref.loc[]   <- NA 
  ref.length  <- mapply(length, s.loc.num)
  com.points  <- NULL
  
  for (k in 1:num.f.data) {
    com.points        <-  mapply(intersect, loc.number[k], s.loc.num)
    num.com           <-  mapply(length, com.points[])
    
    if (ref.length >= length(loc.number[[k]])){
      max.length  <-  ref.length
    } else {
      max.length  <-  length(loc.number[[k]])
    }
    
    com.rate          <-  num.com / max.length
    
    if (loc.dim[k] == 1) {
      id                <-  which.max(com.rate[dim.1.range])
    } else {
      id                <-  which.max(com.rate[dim.2.range]) + num.dim.1
    }

    if ((length(id) != 0)) {    
      if (com.rate[id] > 0.75) {
        ref.loc[k, 1]     <-  id
      } else {
        ref.loc[k, 1]     <-  NA
      }
    }
  }
  
  common.loc            <-  cbind(common.loc, ref.loc)
  
  cnum.contl            <-  which(colnames(common.loc) == n.contl)
  id                    <-  which(is.na(common.loc[, cnum.contl]) == FALSE)
  
  common.loc            <-  common.loc[id, ]
  
  out.name  <-  paste("./", n.case, "-", n.contl, "_Common-location_ref-0001.csv", sep = "")
  
  write.csv(common.loc, file = out.name, row.names = FALSE, quote = FALSE)
}



