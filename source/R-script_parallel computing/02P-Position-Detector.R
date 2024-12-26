### detection of the position for SHIROKANE
##
#
library(tidyr)
library(data.table)
library(stringr)

val         <-  commandArgs(trailingOnly = TRUE)

target.dir  <-  list.dirs("./")
id          <-  grep("CA-data", target.dir)
target.dir  <-  target.dir[id]

t.dir <- target.dir[as.integer(val[1])]
print(paste0("now:", t.dir))
file.list.wt     <-  list.files(t.dir, full.names = TRUE, recursive = TRUE)
id               <-  grep("Location", file.list.wt)
file.list.wt     <-  file.list.wt[id]
num.files        <-  length(file.list.wt)

t.dir_loc        <-  sapply(strsplit(t.dir, "/CA"), function(x){x[1]})

f.data      <-  fread(file.list.wt[1])
loc.dim     <-  f.data$Dimension
loc.number  <-  str_split(f.data$Location[], pattern = ";")

num.f.data    <-  nrow(f.data)

dim.1.range <-  which(loc.dim == 1)
dim.2.range <-  which(loc.dim == 2)
num.dim.1   <-  length(dim.1.range)

nowref  <-  "1st"

common.loc  <-  data.frame("XXX" = c(1:num.f.data))
colnames(common.loc)  <-  nowref


for (j in 2:num.files) {
  s.data        <-  fread(file.list.wt[j])
  s.loc.dim     <-  s.data$Dimension
  s.loc.num     <-  str_split(s.data$Location[], pattern = ";") 
  
  if (j == 2) {
    ref.loc     <-  data.frame("2nd" = c(1:num.f.data))
  } else if (j == 3) {
    ref.loc  <-  data.frame("3rd" = c(1:num.f.data))
  } else if (j != 1) {
    ref.loc  <-  data.frame("XXX" = c(1:num.f.data))
    
    colnames(ref.loc)  <-  paste(j, "th", sep = "")
  }
  
  ref.loc[]   <- NA 
  ref.length  <- mapply(length, s.loc.num)
  com.points  <- NULL
  
  for (k in 1:num.f.data) {
    com.points        <-  mapply(intersect, loc.number[k], s.loc.num)
    num.com           <-  mapply(length, com.points[])
    
    if (ref.length[1] >= length(loc.number[[k]])){
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
    
    if (length(com.rate[id]) != 0) {      
      if (com.rate[id] > 0.75) {
        ref.loc[k, 1]     <-  id
      } else {
        ref.loc[k, 1]     <-  NA
      }
    }
    
  }
  
  common.loc          <-  cbind(common.loc, ref.loc)
}


out.name  <-  paste(t.dir_loc, "/Location_detect/Ref_location_", nowref, ".csv", sep = "")

write.csv(common.loc, file = out.name, row.names = FALSE, quote = FALSE)

