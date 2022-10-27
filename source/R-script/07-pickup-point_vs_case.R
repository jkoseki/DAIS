### pickup point
##
#
library(data.table)
library(earth)
library(Rtsne)
library(mclust)

#set your  contl
n.contl  <-  "WT"


my_outlier <- function(x){
  medy <- median(x)
  mady <- median(abs(x - medy)) * 0.67449
  maxy <- medy + 5 * mady
  miny <- medy - 5 * mady
  x[which(x>maxy)] <- NA
  x[which(x<miny)] <- NA
  return(x)
}


case.dir  <-  list.dirs("./")
id        <-  grep("CA-data", case.dir)
case.dir  <-  case.dir[id]
id        <-  grep(paste0("//", n.contl, "/CA-data"), case.dir)
case.dir  <-  case.dir[-id]

np.case   <-  sapply(strsplit(case.dir, "/CA"), function(x){x[1]})
np.case   <-  sapply(strsplit(np.case,  "//"),  function(x){x[2]})


#for (n.case in np.case) {
for (nc in 1:704) {
  n.case      <-  np.case[nc]
  print(n.case)
  anno        <-  read.csv(paste0("./", n.case, "-", n.contl, "_Common-location_ref-0001.csv"))
  cnum.contl  <-  which(colnames(anno) == n.contl)
  cnum.case   <-  which(colnames(anno) == n.case)
  contl_id    <-  formatC(anno[, cnum.contl], width=4, flag="0")
  case_id     <-  formatC(anno[, cnum.case],  width=4, flag="0")
  case        <-  paste0("./", n.case,  "/Polar-loc_2/", case_id,  "-point_polar-coordinate.csv")
  contl       <-  paste0("./", n.contl, "/Polar-loc_2/", contl_id, "-point_polar-coordinate.csv")
  
  n   <-  length(case)
  r   <-  100
  y1  <-  z1 <- y2 <- z2 <- matrix(0,n,r)
  
  for(i in 1:n){
    cat(i,"/",n,"\n")
    x <- fread(contl[i])
    id <- which(x$Step <=500)
    x <- as.matrix(x[id,])
    
    if (nrow(x) > 1) {
      x[,2] <- x[,2] / x[1,2]
      x[,2] <- my_outlier(x[,2])
      x[,3] <- x[,3] - x[1,3]    # adding koseki 2021_0415
      x[,3] <- my_outlier(x[,3])
      df <- data.frame(x=x[,1],y=x[,2],z=x[,3])
      id <- which(!is.na(df$y))

      earth.obj <- earth(y ~ x, data=df[id,])
      fitted.y1 <- predict(earth.obj,data.frame(x=seq(from=1,to=500,length=r)))
      id <- which(!is.na(df$z))
      loess.obj <- earth(z ~ x, data=df[id,]) # 25% smoothing span
      fitted.z1 <- predict(earth.obj,data.frame(x=seq(from=1,to=500,length=r)))
      y1[i,] <- fitted.y1
      z1[i,] <- fitted.z1
    }
          
    x <- fread(case[i])
    id <- which(x$Step <=500)
    x <- as.matrix(x[id,])
        
    if (nrow(x) > 1) {
      x[,2] <- x[,2] / x[1,2]
      x[,2] <- my_outlier(x[,2])
      x[,3] <- x[,3] - x[1,3]    # adding koseki 2021_0415
      x[,3] <- my_outlier(x[,3])
      df <- data.frame(x=x[,1],y=x[,2],z=x[,3])
      id <- which(!is.na(df$y))
      
      earth.obj <- earth(y ~ x, data=df[id,])
      fitted.y2 <- predict(earth.obj,data.frame(x=seq(from=1,to=500,length=r)))
      id <- which(!is.na(df$z))
      earth.obj <- earth(z ~ x, data=df[id,]) # 25% smoothing span
      fitted.z2 <- predict(earth.obj,data.frame(x=seq(from=1,to=500,length=r)))
      y2[i,] <- fitted.y2
      z2[i,] <- fitted.z2
    }
  }
  
  X <- y2 - y1
  rownames(X) <- 1:n
  
  heatmapCol <- function (data, col, lim) {
    nrcol <- length(col)
    data.range <- range(data)
    if (diff(data.range) == 0) stop("data has range 0")
    if(missing(lim)) lim <- min(abs(data.range))*0.7
    nrcol <- length(col)
    reps1 <- ceiling(nrcol * (-lim - data.range[1])/(2 * lim))
    reps2 <- ceiling(nrcol * (data.range[2] - lim)/(2 * lim))
    col1 <- c(rep(col[1], reps1), col, rep(col[nrcol], reps2))
    return(col1)
  }
  
  library(gplots)
  library(pheatmap)
  
  col  <- heatmapCol(X, col=colorpanel(100, high = "red", low = "blue", mid = "white"))
  id   <- which(abs(X[, 1]) <= 0.15)
  pheatmap(X[id, ], fontsize_row = 5, col = col, cluster_col = FALSE)
  
  
  anno[sort.list(apply(X,1,sum),decreasing=TRUE)[1:20],]
  anno[sort.list(apply(X,1,sum),decreasing=FALSE)[1:20],]
  
  Tgnum1 <- sort.list(apply(X,1,sum),decreasing=TRUE)[1:20]
  Tgnum2 <- sort.list(apply(X,1,sum),decreasing=FALSE)[1:20]
  
  
  for (i in Tgnum1) {
    if (abs(y1[i, 1] - y2[i, 1]) < 0.05) {
      cnum.contl  <-  which(colnames(anno) == n.contl)
      cnum.case   <-  which(colnames(anno) == n.case)
      nout  <-  paste0("./Pickup-Variation-Data/", n.case, "-", formatC(anno[i, cnum.case], width = 3, flag = "0"), 
                       "_vs_", n.contl, "-", formatC(anno[i, cnum.contl], width = 3, flag = "0"), ".pdf")
      pdf(file = nout, width = 11.69, height = 8.27)
      matplot(data.frame(y1[i, ], y2[i, ]), pch = 1:2, main = paste0(i, "-plot"))
      dev.off()
      print(paste0('01-point-', i))
    }  
  }
  
  for (i in Tgnum2) {
    if (abs(y1[i, 1] - y2[i, 1]) < 0.05) {
      cnum.contl  <-  which(colnames(anno) == n.contl)
      cnum.case   <-  which(colnames(anno) == n.case)
      nout  <-  paste0("./Pickup-Variation-Data/", n.case, "-", formatC(anno[i, cnum.case], width = 3, flag = "0"), 
                       "_vs_", n.contl, "-", formatC(anno[i, cnum.contl], width = 3, flag = "0"), ".pdf")
      pdf(file = nout, width = 11.69, height = 8.27)
      matplot(data.frame(y1[i, ], y2[i, ]), pch = 1:2, main = paste0(i, "-plot"))
      dev.off()
      print(paste0('02-point-', i))
    }
  }
  
}
