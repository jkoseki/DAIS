### make a change point matrix
##
#
library(data.table)
library(tidyverse)
library(bio3d)


#set your  contl
n.contl  <-  "WT"


cp.list    <-  list.files("./")
cp.list    <-  cp.list[grep("Change-point.csv", cp.list)]

case.list  <-  sapply(strsplit(cp.list, "-"), function(x){x[1]})


for (i in case.list) {
  t.pdb  <-  read.pdb2(paste0("./", i, "/Separate_00001.pdb"))
  ca.id  <-  which(t.pdb$atom$elety == "CA")
  t.pdb  <-  t.pdb$atom[ca.id, ]
  r.cnt  <-  as.data.frame(matrix(0, nrow = 1, ncol = length(paste0(c(t.pdb$resid), c(t.pdb$resno)))))

  colnames(r.cnt)  <-  paste0(c(t.pdb$resid), c(t.pdb$resno))

  
  cp.dat  <-  fread(paste0("./", i, "-Change-point.csv"))
  
  if (nrow(cp.dat) != 0) {
    for (j in 1:nrow(cp.dat)) {
      pu.id  <-  unlist(strsplit(cp.dat$`PDB Number`[j], ";"))
      
      for (k in 1:length(pu.id)) {
         res.id  <-  which(pu.id[k] == t.pdb$resno)
         
         r.cnt[, res.id]  <-  r.cnt[, res.id] + 1
      }
    }
  }
  
  out.file  <-  paste0(i, "-", n.contl, "_change-point_res-counts.csv")
  write.csv(r.cnt, out.file, row.names = FALSE, quote = FALSE)
}



cp.list    <-  list.files("./")
cp.list    <-  cp.list[grep("change-point_res-counts.csv", cp.list)]

case.list  <-  sapply(strsplit(cp.list, "-"), function(x){x[1]})
pre.mtr    <-  NA


for (i in 1:length(cp.list)) {
  t.cprc  <-  fread(cp.list[i])
  t.cprc  <-  t(t.cprc)
  r.num   <-  function(str) as.numeric(regmatches(str, regexpr("[0-9]+", str, perl = TRUE)))
  t.cprc  <-  data.frame("res_num"  = sapply(row.names(t.cprc), r.num),
                         "target"   = t.cprc[, 1])
  colnames(t.cprc)[2] <- c(case.list[i])
    
  if (i != 1) {
    t.cprc  <-  full_join(pre.mtr, t.cprc, by = "res_num")
  } else {
    first.r.name  <-  row.names(t.cprc)
  }
  
  pre.mtr <- t.cprc
}

row.names(t.cprc)  <-  first.r.name

if (length(cp.list) != 1) {
  t.cprc             <-  t.cprc[, -1]
}

out.file  <-  "Summary_change-point_res-counts.csv"
write.csv(t.cprc, out.file, quote = FALSE)
