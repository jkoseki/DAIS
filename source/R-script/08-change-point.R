### make a change point data
##
#
library(data.table)
library(tidyverse)
library(bio3d)


f.list   <-  list.files("./Pickup-Variation-Data/")
f.id     <-  grep("_vs_", f.list)
f.lists  <-  f.list[f.id]
f.lists  <-  sapply(strsplit(f.list, "\\.pdf"), function(x){x[1]})


#set your  contl
n.contl  <-  "WT"

case.dir  <-  list.dirs("./")
id        <-  grep("CA-data", case.dir)
case.dir  <-  case.dir[id]
id        <-  grep(paste0("//", n.contl, "/CA-data"), case.dir)
case.dir  <-  case.dir[-id]

np.case   <-  sapply(strsplit(case.dir, "/CA"), function(x){x[1]})
np.case   <-  sapply(strsplit(np.case,  "//"),  function(x){x[2]})

ncase       <-  length(np.case)


for (nc in 1:ncase) {
  case.name  <-  np.case[nc]

  f.list  <-  f.lists[grep(case.name, f.lists)]



  case    <-  sapply(strsplit(f.list, "_vs_"), function(x){x[1]})
  contl   <-  sapply(strsplit(f.list, "_vs_"), function(x){x[2]})

  case.l  <-  length(unlist(strsplit(case[1],  "-")))
  contl.l <-  length(unlist(strsplit(contl[1], "-")))


  case.num    <-  sapply(strsplit(case,  "-"), function(x){x[case.l]})

  contl.name  <-  paste(c(unlist(strsplit(contl[1],  "-"))[1:(contl.l - 1)]), collapse="-")
  contl.num   <-  sapply(strsplit(contl, "-"), function(x){x[contl.l]})

  c.point     <-  data.frame("Dim"        =  NA,
                             "#"          =  NA,
                             "Case"       =  as.numeric(case.num),
                             "Contl"      =  as.numeric(contl.num),
                             "Common"     =  NA,
                             "PDB Number" =  NA)

  colnames(c.point)  <-  c("Dim", "#", case.name, contl.name, "Common", "PDB Number")

  loc.name    <-  list.files("./")
  loc.name.id <-  grep("Common-location_ref-0001.csv", loc.name)
  loc.name    <-  loc.name[loc.name.id]
  loc.name.id <-  grep(case.name, loc.name)
  loc.name    <-  loc.name[loc.name.id]
  
  corresp.loc <-  fread(loc.name)
  pu.atom     <-  fread("./Pickup_Atoms.csv")
  loc.case    <-  fread(paste0("./", case.name,  "/CA-data/Location-Separate_00001.csv"))
  loc.contl   <-  fread(paste0("./", contl.name, "/CA-data/Location-Separate_00001.csv"))



  ref.prot    <-  read.pdb2(paste0("./", case.name,  "/Separate_00001.pdb"))
  ref.prot$atom %>%
    select("resid", "resno", "elety") -> ref.seq

  ref.seq  <-  subset(ref.seq, ref.seq$elety == "CA")



  for (i in 1:nrow(c.point)) {
    c.point[i, 2]  <-  which(corresp.loc[, 1] == c.point[i, 3])
    c.point[i, 1]  <-  loc.case$Dimension[c.point[i, 2]]  
  
    pu.loc.case    <-  unlist(strsplit(loc.case$Location[ c.point[i, 3]], ";"))
    pu.loc.contl   <-  unlist(strsplit(loc.contl$Location[c.point[i, 4]], ";"))
  
    common.loc     <-  c(unlist(intersect(pu.loc.case, pu.loc.contl)))
    c.point[i, 5]  <-  paste0(common.loc, collapse=";")

    pu.cor.loc     <-  NULL

    for (j in 1:length(common.loc)) {
      id           <-  unlist(pu.atom[as.numeric(common.loc[j])])
      pu.cor.loc   <-  c(pu.cor.loc, ref.seq$resno[id])
    }
  
    c.point[i, 6]  <-  paste0(pu.cor.loc, collapse=";")
  }



  out.file  <-  paste0("./", case.name, "-Change-point.csv")
  write.csv(c.point, out.file, row.names = FALSE, quote = FALSE)
}


