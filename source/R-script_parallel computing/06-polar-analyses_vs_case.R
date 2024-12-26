###  For plot r
##
#
library(data.table)
library(openxlsx)
library(dplyr)


#set your  contl
n.contl  <-  "WT"

t.type    <-  2
d.name.1  <-  "Polar-loc_2"
x.max     <-  500
w.type    <-  "0001-0500_"


case.dir  <-  list.dirs("./")
id        <-  grep("CA-data", case.dir)
case.dir  <-  case.dir[id]
id        <-  grep(paste0("//", n.contl, "/CA-data"), case.dir)
case.dir  <-  case.dir[-id]

np.case   <-  sapply(strsplit(case.dir, "/CA"), function(x){x[1]})
np.case   <-  sapply(strsplit(np.case,  "//"),  function(x){x[2]})


for (n.case in np.case) {
  n.com.loc  <-  paste0("./", n.case, "-", n.contl, "_Common-location_ref-0001.csv")  
  cor.num    <-  fread(n.com.loc)
  
  
  
  t.dir  <-  paste("./", n.contl, "/", d.name.1, sep = "")
  f.list.1.contl   <-  list.files(t.dir, full.names = TRUE)
  
  t.dir  <-  paste("./", n.case, "/", d.name.1, sep = "")
  f.list.1.case  <-  list.files(t.dir, full.names = TRUE)
  
  
  pu.num  <- nrow(cor.num)
  
  
  for (i in 1:pu.num) {
    id.contl   <-  cor.num[i, ..n.contl]
    id.case    <-  cor.num[i, ..n.case]
    
    contl.pol  <-  fread(f.list.1.contl[as.numeric(id.contl)])
    case.pol   <-  fread(f.list.1.case[as.numeric(id.case)])
    
    if (t.type == 2) {
      contl.pol                 %>%
        filter(Step <= 500)  -> contl.pol
      
      case.pol                  %>%
        filter(Step <= 500)  -> case.pol
    }
    
    y.min    <-  as.integer(min(min(contl.pol$r), min(case.pol$r)) -3)
    y.max    <-  as.integer(max(max(contl.pol$r), max(case.pol$r)) +1)
    
    
    out.name <-  paste("./r-plot/",
                       w.type, formatC(i, width = 3, flag = "0"),
                       "_", n.case, "-", formatC(as.numeric(id.case), width = 3, flag = "0"),
                       "_vs_", n.contl, "-", formatC(as.numeric(id.contl),  width = 3, flag = "0"),
                       ".pdf", sep = "")
    
    pdf(out.name, width = 960/96, height = 600/96)
    par(las = 1)
    plot(x = NULL,   type = "n", xlab = "Calculation Step",  ylab = "Distance", 
         xlim = c(0, x.max), ylim = c(y.min, y.max),
         xaxs = "i", yaxs = "i", xaxt = "n",      yaxt = "n",     bty = "n", 
         main = "Time cource Distance")
    axis(side = 1, at = seq(0, x.max, 500), tck = 1.0, 
         lty = "blank", labels = seq(0, x.max, 500))
    axis(side = 2, at = seq(y.min, y.max, 1), tck = 1.0, 
         lty = "blank", labels = seq(y.min, y.max, 1))
    
    
    points(contl.pol$Step, contl.pol$r, pch = 16, col = "blue")
    points(case.pol$Step,  case.pol$r,  pch = 16, col = "red")
    
    legend("bottomright", legend = c(n.case, n.contl), 
           col = c("red", "blue"), pch = 16)
    box(bty = "l")
    dev.off()
    
  }
  
}






