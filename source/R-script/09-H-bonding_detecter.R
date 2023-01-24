### detective H-bonding
##
#
library(data.table)
library(tidyverse)
library(bio3d)
library(stringr)
library(dplyr)

#set your  contl
n.contl  <-  "WT"

case.dir  <-  list.dirs("./")
id        <-  grep("CA-data", case.dir)
case.dir  <-  case.dir[id]
id        <-  grep(paste0("//", n.contl, "/CA-data"), case.dir)
case.dir  <-  case.dir[-id]

np.case   <-  sapply(strsplit(case.dir, "/CA"), function(x){x[1]})
np.case   <-  sapply(strsplit(np.case,  "//"),  function(x){x[2]})


pdb.wt   <-  list.files(paste0("./", n.contl, "/"), full.names = TRUE, recursive = FALSE)
wt.id    <-  grep(".pdb", pdb.wt)
pdb.wt   <-  pdb.wt[wt.id]


#==========
for (case.name in np.case) {
  print(case.name)
  pdb.mut    <-  list.files(paste0("./", case.name, "/"), full.names = TRUE, recursive = FALSE)
  mut.id     <-  grep(".pdb", pdb.mut)
  pdb.mut    <-  pdb.mut[mut.id]
  
  
  open.file  <-  paste0("./", case.name, "-Change-point.csv")
  target.ca  <-  read.csv(open.file)
  
  target.ca                     %>%
    select(Dim, PDB.Number)     %>%
    filter(Dim == 1 | Dim == 2) -> target.ca
  
  
  
  tc.h.bonding.wt   <-  NULL
  tc.h.bonding.mut  <-  NULL
  
  
  for (i in 1:length(pdb.wt)) {
    wt.data   <-  read.pdb2(pdb.wt[i])
    mut.data  <-  read.pdb2(pdb.mut[i])
    
    wt.data$atom                            %>%
      select(elety, resid, resno, x, y, z)  -> wt.data
    
    mut.data$atom                           %>%
      select(elety, resid, resno, x, y, z)  -> mut.data
    
    
    for (j in 1:nrow(target.ca)) {
      ca.list  <-  strsplit(as.character(target.ca$PDB.Number[j]), ";")
      num.ca   <-  length(ca.list[[1]])
      pu.res   <-  as.numeric(ca.list[[1]])
      pu.res   <-  pu.res[order(pu.res)]
      
      wt.data                 %>%
        filter(elety == "CA") -> wt.ca
      mut.data                %>%
        filter(elety == "CA") -> mut.ca
      
      ca.dis.wt   <-  as.matrix(dist(data.frame(wt.ca$x,  wt.ca$y,  wt.ca$z )))
      ca.dis.mut  <-  as.matrix(dist(data.frame(mut.ca$x, mut.ca$y, mut.ca$z)))
      
      
      
      id.wt   <-  NULL
      id.mut  <-  NULL
      
      for (k in pu.res) {
        id.wt   <-  unique(c(id.wt,  wt.ca [which(ca.dis.wt [which(wt.ca$resno  == k), ] <= 10), ]$resno))
        id.mut  <-  unique(c(id.mut, mut.ca[which(ca.dis.mut[which(mut.ca$resno == k), ] <= 10), ]$resno))
      }
      
      pu.wt    <-  id.wt [order(id.wt)]
      pu.mut   <-  id.mut[order(id.mut)]
      
      
      
      id.wt   <-  NULL
      id.mut  <-  NULL
      
      for (k in pu.wt) {
        id.wt   <-  c(id.wt,  which(wt.data$resno  == k))
      }
      for (k in pu.mut) {
        id.mut  <-  c(id.mut, which(mut.data$resno == k))
      }
      
      
      wt.data2   <-  wt.data[id.wt, ]
      mut.data2  <-  mut.data[id.mut, ]
      
      id.wt      <-  which(str_sub(wt.data2$elety,  1 , 1) == "N" |
                             str_sub(wt.data2$elety,  1 , 1) == "O" |
                             str_sub(wt.data2$elety,  1 , 1) == "H" )
      
      id.mut     <-  which(str_sub(mut.data2$elety,  1 , 1) == "N" |
                             str_sub(mut.data2$elety,  1 , 1) == "O" |
                             str_sub(mut.data2$elety,  1 , 1) == "H" )
      
      wt.data2   <-  wt.data2[id.wt, ]
      mut.data2  <-  mut.data2[id.mut, ]
      
      rownames(wt.data2)   <-  paste(wt.data2$resno,  wt.data2$elety,  sep = "-")
      rownames(mut.data2)  <-  paste(mut.data2$resno, mut.data2$elety, sep = "-")
      
      
      dis.wt    <-  as.matrix(dist(data.frame(wt.data2$x,  wt.data2$y,  wt.data2$z )))
      dis.mut   <-  as.matrix(dist(data.frame(mut.data2$x, mut.data2$y, mut.data2$z)))
      
      id.wt     <-  which(dis.wt  >= 1.5 & dis.wt  <= 2.2, arr.ind = TRUE)
      id.mut    <-  which(dis.mut >= 1.5 & dis.mut <= 2.2, arr.ind = TRUE)
      
      
      id.wt     <-  id.wt [which(wt.data2$resno [id.wt [, 1]] != wt.data2$resno [id.wt [, 2]] &
                                   abs(wt.data2$resno [id.wt [, 1]] - wt.data2$resno [id.wt [, 2]]) != 1 &
                                   str_sub(wt.data2$elety [id.wt [, 1]], 1, 1) 
                                 != str_sub(wt.data2$elety [id.wt [, 2]], 1, 1) &
                                   str_sub(wt.data2$elety [id.wt [, 1]], 1, 1) == "H"), ]
      id.mut    <-  id.mut[which(mut.data2$resno[id.mut[, 1]] != mut.data2$resno[id.mut[, 2]] &
                                   abs(mut.data2$resno[id.mut[, 1]] - mut.data2$resno[id.mut[, 2]]) != 1 &
                                   str_sub(mut.data2$elety [id.mut [, 1]], 1, 1) 
                                 != str_sub(mut.data2$elety[id.mut[, 2]], 1, 1) &
                                   str_sub(mut.data2$elety[id.mut[, 1]], 1, 1) == "H"), ]
      
      
      ang.wt  <-  data.frame(matrix(FALSE, nrow = nrow(id.wt),  ncol = 1))
      ang.mut <-  data.frame(matrix(FALSE, nrow = nrow(id.mut), ncol = 1))
      colnames(ang.wt)   <-  c("T_or_F")
      colnames(ang.mut)  <-  c("T_or_F")
      
      for (k in 1:nrow(id.wt)) {
        detC     <-  as.numeric(names(dis.wt[id.wt[k, 1], order(dis.wt[id.wt[k, 1], ])][1:3]))
        Vec_1_X  <-  wt.data2[detC[2], ]$x - wt.data2[detC[1], ]$x
        Vec_1_Y  <-  wt.data2[detC[2], ]$y - wt.data2[detC[1], ]$y
        Vec_1_Z  <-  wt.data2[detC[2], ]$z - wt.data2[detC[1], ]$z
        Vec_2_X  <-  wt.data2[detC[3], ]$x - wt.data2[detC[1], ]$x
        Vec_2_Y  <-  wt.data2[detC[3], ]$y - wt.data2[detC[1], ]$y
        Vec_2_Z  <-  wt.data2[detC[3], ]$z - wt.data2[detC[1], ]$z
        
        Theta    <-  ((Vec_1_X * Vec_2_X) + (Vec_1_Y * Vec_2_Y) + (Vec_1_Z * Vec_2_Z)) /
          sqrt((Vec_1_X * Vec_1_X) + (Vec_1_Y * Vec_1_Y) + (Vec_1_Z * Vec_1_Z)) /
          sqrt((Vec_2_X * Vec_2_X) + (Vec_2_Y * Vec_2_Y) + (Vec_2_Z * Vec_2_Z))
        
        if (180 * acos(Theta) / pi > 160) {
          ang.wt[k, 1]  <-  TRUE
        } 
      }
      
      for (k in 1:nrow(id.mut)) {
        detC     <-  as.numeric(names(dis.mut[id.mut[k, 1], order(dis.mut[id.mut[k, 1], ])][1:3]))
        Vec_1_X  <-  mut.data2[detC[2], ]$x - mut.data2[detC[1], ]$x
        Vec_1_Y  <-  mut.data2[detC[2], ]$y - mut.data2[detC[1], ]$y
        Vec_1_Z  <-  mut.data2[detC[2], ]$z - mut.data2[detC[1], ]$z
        Vec_2_X  <-  mut.data2[detC[3], ]$x - mut.data2[detC[1], ]$x
        Vec_2_Y  <-  mut.data2[detC[3], ]$y - mut.data2[detC[1], ]$y
        Vec_2_Z  <-  mut.data2[detC[3], ]$z - mut.data2[detC[1], ]$z
        
        Theta    <-  ((Vec_1_X * Vec_2_X) + (Vec_1_Y * Vec_2_Y) + (Vec_1_Z * Vec_2_Z)) /
          sqrt((Vec_1_X * Vec_1_X) + (Vec_1_Y * Vec_1_Y) + (Vec_1_Z * Vec_1_Z)) /
          sqrt((Vec_2_X * Vec_2_X) + (Vec_2_Y * Vec_2_Y) + (Vec_2_Z * Vec_2_Z))
        
        if (180 * acos(Theta) / pi > 160) {
          ang.mut[k, 1]  <-  TRUE
        } 
      }
      
      
      id.wt   <-  subset(id.wt,  ang.wt [, 1] == TRUE)
      id.mut  <-  subset(id.mut, ang.mut [, 1] == TRUE)
      
      
      id.wt[, 1]   <-  wt.data2$resno [id.wt [, 1]]
      id.wt[, 2]   <-  wt.data2$resno [id.wt [, 2]]
      id.mut[, 1]  <-  mut.data2$resno[id.mut[, 1]]
      id.mut[, 2]  <-  mut.data2$resno[id.mut[, 2]]
      
      if (j != 1) {
        h.bonding.wt   <-  rbind(h.bonding.wt,  id.wt)
        h.bonding.mut  <-  rbind(h.bonding.mut, id.mut)
      } else {
        h.bonding.wt   <-  id.wt
        h.bonding.mut  <-  id.mut
      }
      
      
      
    }
    
    colnames(h.bonding.wt)   <-  c("H-Donor", "H-Acceptor")
    colnames(h.bonding.mut)  <-  c("H-Donor", "H-Acceptor")
    
    if (i != 1) {
      tc.h.bonding.wt   <-  c(tc.h.bonding.wt,  list(unique(h.bonding.wt [order(h.bonding.wt [, 1]), ])))
      tc.h.bonding.mut  <-  c(tc.h.bonding.mut, list(unique(h.bonding.mut[order(h.bonding.mut[, 1]), ])))
    } else {
      tc.h.bonding.wt   <-  list(unique(h.bonding.wt [order(h.bonding.wt [, 1]), ]))
      tc.h.bonding.mut  <-  list(unique(h.bonding.mut[order(h.bonding.mut[, 1]), ]))
    }
    
  }
  
  
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_wt.rds")
  saveRDS(tc.h.bonding.wt,  out.file)
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_mut.rds")
  saveRDS(tc.h.bonding.mut, out.file)
  
  
  
  
  sum.tc.h.bonding.wt   <-  NULL
  sum.tc.h.bonding.mut  <-  NULL
  sep.row               <-  data.frame("", "")
  colnames(sep.row)     <-  c("H-Donor", "H-Acceptor")
  
  for (i in 1:500) {
    nrwt   <-  nrow(tc.h.bonding.wt[[i]])
    nrmut  <-  nrow(tc.h.bonding.mut[[i]])
    
    if (is.null(nrwt) != TRUE) { 
      row.names(tc.h.bonding.wt[[i]])   <-  paste0("Step-", formatC(i, width = 3, flag = 0), "-", 1:nrwt)
    } else {
      tc.h.bonding.wt[[i]]              <-  rbind(tc.h.bonding.wt[[i]], tc.h.bonding.wt[[i]])
      row.names(tc.h.bonding.wt[[i]])   <-  paste0("Step-", formatC(i, width = 3, flag = 0), "-", 1:2)
    }
    if (is.null(nrmut) != TRUE) {
      row.names(tc.h.bonding.mut[[i]])  <-  paste0("Step-", formatC(i, width = 3, flag = 0), "-", 1:nrmut)
    } else {
      tc.h.bonding.mut[[i]]             <-  rbind(tc.h.bonding.mut[[i]], tc.h.bonding.mut[[i]])
      row.names(tc.h.bonding.mut[[i]])  <-  paste0("Step-", formatC(i, width = 3, flag = 0), "-", 1:1)
    }
    
    
    row.names(sep.row)    <-  paste0("Step-", formatC(i, width = 3, flag = 0), "-END")

    if (is.null(nrwt) != TRUE) {     
      sum.tc.h.bonding.wt   <-  rbind(sum.tc.h.bonding.wt,  tc.h.bonding.wt [[i]], sep.row)
    } else {
      sum.tc.h.bonding.wt   <-  rbind(sum.tc.h.bonding.wt,  tc.h.bonding.wt [[i]][1, ], sep.row)
    }
    if (is.null(nrmut) != TRUE) {
      sum.tc.h.bonding.mut  <-  rbind(sum.tc.h.bonding.mut, tc.h.bonding.mut[[i]], sep.row)
    } else {
      sum.tc.h.bonding.mut  <-  rbind(sum.tc.h.bonding.mut, tc.h.bonding.mut[[i]][1, ], sep.row)
    }
  }
  
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_wt.csv")
  write.csv(sum.tc.h.bonding.wt,  out.file)
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_mut.csv")
  write.csv(sum.tc.h.bonding.mut, out.file)
  
  
  
  
  
  ### Analyses
  ##
  #
  sum.tc.h.bonding.wt                                                                     %>%
    mutate(Det_name  =  paste0(formatC(as.numeric(.$`H-Donor`),    width = 4, flag = 0),
                               ";", 
                               formatC(as.numeric(.$`H-Acceptor`), width = 4, flag = 0))) ->  sum.tc.h.bonding.wt2
  
  sum.tc.h.bonding.mut                                                                     %>%
    mutate(Det_name  =  paste0(formatC(as.numeric(.$`H-Donor`),    width = 4, flag = 0),
                               ";", 
                               formatC(as.numeric(.$`H-Acceptor`), width = 4, flag = 0))) ->  sum.tc.h.bonding.mut2
  
  rownames(sum.tc.h.bonding.wt2)   <-  rownames(sum.tc.h.bonding.wt)
  rownames(sum.tc.h.bonding.mut2)  <-  rownames(sum.tc.h.bonding.mut)
  
  
  u.det.name.wt   <-  unique(sum.tc.h.bonding.wt2$Det_name)
  u.det.name.mut  <-  unique(sum.tc.h.bonding.mut2$Det_name)
  
  idna            <-  grep("NA", u.det.name.wt)
  u.det.name.wt   <-  u.det.name.wt[-idna]
  idna            <-  grep("NA", u.det.name.mut)
  u.det.name.mut  <-  u.det.name.mut[-idna]
  
  
  
  dif.wt          <-  setdiff(u.det.name.wt,   u.det.name.mut)
  dif.mut         <-  setdiff(u.det.name.mut,  u.det.name.wt)
  com.bonding     <-  intersect(u.det.name.wt, u.det.name.mut)
  
  
  out.file  <-  paste0(case.name, "-WT-only_Hydrogen-bonding.csv")
  write.csv(dif.wt,      out.file, quote = FALSE, row.names = FALSE)
  out.file  <-  paste0(case.name, "-Mut-only_Hydrogen-bonding.csv")
  write.csv(dif.mut,     out.file, quote = FALSE, row.names = FALSE)
  out.file  <-  paste0(case.name, "-Common_Hydrogen-bonding.csv")
  write.csv(com.bonding, out.file, quote = FALSE, row.names = FALSE)
  
  
  tc.count.wt     <-  data.frame("WT_H-bonding"     =  u.det.name.wt)
  tc.count.mut    <-  data.frame("Mut_H-bonding"  =  u.det.name.mut)
  
  
  id.end.wt       <-  c(0, grep("END", rownames(sum.tc.h.bonding.wt2)))
  id.end.mut      <-  c(0, grep("END", rownames(sum.tc.h.bonding.mut2)))
  
  
  for (i in 1:500) {
    count.wt   <-  data.frame(matrix(0, nrow = length(u.det.name.wt),  ncol = 1))
    count.mut  <-  data.frame(matrix(0, nrow = length(u.det.name.mut), ncol = 1))
    
    for (j in 1:length(u.det.name.wt)) {
      num.wt  <-  length(grep(u.det.name.wt[j],  
                              sum.tc.h.bonding.wt2 $Det_name[(id.end.wt[i]  + 1):(id.end.wt[i + 1]  - 1)]))
      
      if (num.wt == 1) count.wt[j, 1] <- 1
    }
    
    for (j in 1:length(u.det.name.mut)) {
      num.mut <-  length(grep(u.det.name.mut[j],  
                              sum.tc.h.bonding.mut2 $Det_name[(id.end.mut[i]  + 1):(id.end.mut[i + 1]  - 1)]))
      
      if (num.mut == 1) count.mut[j, 1] <- 1
    }
    
    colnames(count.wt)   <-  paste0("Count-", formatC(i, width = 3, flag = 0))
    colnames(count.mut)  <-  paste0("Count-", formatC(i, width = 3, flag = 0))
    
    tc.count.wt          <-  cbind(tc.count.wt , count.wt)
    tc.count.mut         <-  cbind(tc.count.mut, count.mut)
  }
  
  
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_wt-Count.csv")
  write.csv(tc.count.wt,  out.file,  row.names = FALSE, quote = FALSE)
  out.file  <-  paste0(case.name, "-Time-course_H-bonding_mut-Count.csv")
  write.csv(tc.count.mut, out.file, row.names = FALSE, quote = FALSE)
  
  
  common.hb      <-  intersect(tc.count.wt[, 1], tc.count.mut[, 1])
  wt.only.hb     <-  setdiff(tc.count.wt[, 1], tc.count.mut[, 1])
  mut.only.hb    <-  setdiff(tc.count.mut[, 1], tc.count.wt[, 1])
  
  
  summarized.hb  <-  data.frame("  Type "  =  "Common",
                                "H_bonds"  =  common.hb,
                                "001_100"  =  NA,
                                "101_200"  =  NA,
                                "201_300"  =  NA,
                                "301_400"  =  NA,
                                "401_500"  =  NA,
                                "bd_Rate"  =  NA)
  
  summarized.hb  <-  rbind(summarized.hb, data.frame("  Type "  =  "WT_only",
                                                     "H_bonds"  =  wt.only.hb,
                                                     "001_100"  =  NA,
                                                     "101_200"  =  NA,
                                                     "201_300"  =  NA,
                                                     "301_400"  =  NA,
                                                     "401_500"  =  NA,
                                                     "bd_Rate"  =  NA))
  
  summarized.hb  <-  rbind(summarized.hb, data.frame("  Type "  =  "Mut_only",
                                                     "H_bonds"  =  mut.only.hb,
                                                     "001_100"  =  NA,
                                                     "101_200"  =  NA,
                                                     "201_300"  =  NA,
                                                     "301_400"  =  NA,
                                                     "401_500"  =  NA,
                                                     "bd_Rate"  =  NA))
  
  for (i in 1:nrow(summarized.hb)) {
    id.wt      <-  which(tc.count.wt[,  1] == as.character(summarized.hb[i, 2]))
    id.mut     <-  which(tc.count.mut[, 1] == as.character(summarized.hb[i, 2]))
    
    if ((length(id.wt) != 0) && (length(id.mut) != 0)) {
      summarized.hb[i, 3]  <-  sum(tc.count.mut[id.mut,   2:101]) - sum(tc.count.wt[id.wt,   2:101])
      summarized.hb[i, 4]  <-  sum(tc.count.mut[id.mut, 102:201]) - sum(tc.count.wt[id.wt, 102:201])
      summarized.hb[i, 5]  <-  sum(tc.count.mut[id.mut, 202:301]) - sum(tc.count.wt[id.wt, 202:301])
      summarized.hb[i, 6]  <-  sum(tc.count.mut[id.mut, 302:401]) - sum(tc.count.wt[id.wt, 302:401])
      summarized.hb[i, 7]  <-  sum(tc.count.mut[id.mut, 402:501]) - sum(tc.count.wt[id.wt, 402:501])
    } else if (length(id.wt) == 0) {
      summarized.hb[i, 3]  <-  sum(tc.count.mut[id.mut,   2:101])
      summarized.hb[i, 4]  <-  sum(tc.count.mut[id.mut, 102:201])
      summarized.hb[i, 5]  <-  sum(tc.count.mut[id.mut, 202:301])
      summarized.hb[i, 6]  <-  sum(tc.count.mut[id.mut, 302:401])
      summarized.hb[i, 7]  <-  sum(tc.count.mut[id.mut, 402:501])
    } else if (length(id.mut) == 0) {
      summarized.hb[i, 3]  <-  0 - sum(tc.count.wt[id.wt,   2:101])
      summarized.hb[i, 4]  <-  0 - sum(tc.count.wt[id.wt, 102:201])
      summarized.hb[i, 5]  <-  0 - sum(tc.count.wt[id.wt, 202:301])
      summarized.hb[i, 6]  <-  0 - sum(tc.count.wt[id.wt, 302:401])
      summarized.hb[i, 7]  <-  0 - sum(tc.count.wt[id.wt, 402:501])
    }
    
    summarized.hb[i, 8]  <-  sum(summarized.hb[i, 3:7] / 500)
    
  }
  
  
  summarized.hb  <-  summarized.hb[order(-abs(summarized.hb[, 8])), ]
  
  
  out.file  <-  paste0(case.name, "-H-bonding_time-course-variation.csv")
  write.csv(summarized.hb, out.file, row.names = FALSE, quote = FALSE)
}
