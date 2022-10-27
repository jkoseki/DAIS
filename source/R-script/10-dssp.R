### detective 
##
#
library(data.table)
library(tidyverse)

dssp.head   <-  c("A", "B", "C") 
write.table(dssp.head, "DSSP-header", quote = FALSE, row.names = FALSE, col.names = FALSE)

system('cat DSSP-header ./WT/Separate_00500.pdb > DSSP-input.pdb')
system('rm DSSP-header')

system('mkdssp -i DSSP-input.pdb > DSSP-out.txt')
system('rm DSSP-input.pdb')

system('grep " X " DSSP-out.txt > 2nd-struct.out')
system('rm DSSP-out.txt')

system('cut -b 6-10 2nd-struct.out > res_num.txt')
system('cut -b 12   2nd-struct.out > chain.txt')
system('cut -b 14   2nd-struct.out > amino_acid.txt')
system('cut -b 17   2nd-struct.out > 2nd_struct.txt')
system('sed -i "s/ /-/g" 2nd_struct.txt')
system('rm 2nd-struct.out')

sec.struct  <-  data.frame("res_num"    = fread('res_num.txt',    header = FALSE),
                           "chain"      = fread('chain.txt',      header = FALSE),
                           "amino_acid" = fread('amino_acid.txt', header = FALSE),
                           "2nd-struct" = fread('2nd_struct.txt', header = FALSE))

colnames(sec.struct)  <-  c("res_num", "chain", "amino_acid", "2nd-struct")

system('rm res_num.txt chain.txt amino_acid.txt 2nd_struct.txt')


sec.struct  %>%
  mutate(group = .$'2nd-struct')  ->  sec.struct

sec.struct$'group'[which(sec.struct$'group'  ==  'B')]  <-  '-'
sec.struct$'group'[which(sec.struct$'group'  ==  'E')]  <-  'B'
sec.struct$'group'[which(sec.struct$'group'  ==  'G')]  <-  'H'
sec.struct$'group'[which(sec.struct$'group'  ==  'I')]  <-  'H'
sec.struct$'group'[which(sec.struct$'group'  ==  'T')]  <-  '-'

write.table(sec.struct, file = "2nd-struct.csv", quote = FALSE, col.names = FALSE, sep = ",")


###  Analysis
##
#
#Intrinsically disordered proteins
idp        <-  sec.struct$res_num[which(sec.struct$'2nd-struct'  ==  '-')]

c.point     <-  list.files("./")
id          <-  grep("Change-point.csv", c.point)
c.point     <-  c.point[id]

cp.num      <-  length(c.point)


for (j in 1:cp.num) {
   mut.name   <-  sapply(strsplit(c.point[j], "-"), function(x){x[1]})

   st.change  <-  fread(paste0("./", c.point[j]))
   st.change  <-  st.change[grep(";", st.change$'PDB Number'), ]

   loc.num    <-  str_split(st.change$'PDB Number'[], pattern = ";")
   num.st.c   <-  length(loc.num)
   det.st.c   <-  data.frame("Determinant"  =  matrix(NA, num.st.c, 1))

   for (i in 1:num.st.c) {
      det.comm  <-  intersect(idp, loc.num[[i]])
      
      if (length(det.comm) / length(loc.num[[i]]) < 0.2) {
         det.st.c[i, 1]  <- 'PU'
      } else {
      det.st.c[i, 1]  <- 'exc'
      }
   }


   pu.change.point  <-  subset(st.change, det.st.c[, 1]  ==  'PU')
   out.name         <-  paste0("./PU_", mut.name, "_Change-point.csv")
   write.csv(pu.change.point,  file = out.name, quote = FALSE)



   cp.ca.all  <-  NULL

   for (i in 1:num.st.c) {
      cp.ca.all  <-  c(cp.ca.all, as.numeric(loc.num[[i]]))
   }


   ### H-bond refilter
   ##
   input.name <-  paste0("./", mut.name, "-H-bonding_time-course-variation.csv")
   h.bonding  <-  fread(input.name)

   loc.num_h  <-  str_split(h.bonding$'H_bonds'[], pattern = ";")
   num.h.bd   <-  length(loc.num_h)
   det.h.bd   <-  data.frame("Determinant"  =  matrix(NA, num.h.bd))


   for (i in 1:num.h.bd) {
      det.comm  <-  intersect(idp, as.numeric(loc.num_h[[i]]))
   
      if (length(det.comm)  ==  0  &&  h.bonding$bd_Rate[i]  >  0.6) { 
         id1  <-  which(sec.struct$res_num  ==  as.numeric(loc.num_h[[i]][1]))
         id2  <-  which(sec.struct$res_num  ==  as.numeric(loc.num_h[[i]][2]))
      
         if (length(unique(sec.struct$'group'[id1:id2])) > 1) {
            ca1  <-  as.numeric(loc.num_h[[i]][1])
            ca2  <-  as.numeric(loc.num_h[[i]][2])
            id.comp  <-  c((ca1 -2), (ca1 -1), ca1, (ca1 +1), (ca1 +2), 
                           (ca2 -2), (ca2 -1), ca2, (ca2 +1), (ca2 +2))

            if (length(intersect(id.comp, cp.ca.all)) != 0) {
               det.h.bd[i, 1]  <- 'PU'
            } else {
               det.h.bd[i, 1]  <- 'exc'
            }

         } else {
            det.h.bd[i, 1]  <- 'exc'
         }
      } else {
         det.h.bd[i, 1]  <- 'exc'
      }
   }

   pu.h.bonding.tc  <-  subset(h.bonding, det.h.bd[, 1]  ==  'PU')

   out.name         <-  paste0("./PU_", mut.name, 
                               "_H-binding_time-course-variation.csv")
   write.csv(pu.h.bonding.tc,  file = out.name, quote = FALSE)

}


