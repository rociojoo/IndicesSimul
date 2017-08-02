# setwd('/home/rocio/mnt/data/data2/salidas')

library('pdist')
pdist.pdist <- pdist::pdist
library('polyclip','colorspace')
library(biwavelet)
library(geoR)

getEllipse <- function (X_, Y_, D_, npts=300) {  ## computing ellipses for jPPA
  x0 = mean(X_)
  y0= mean(Y_)
  a <- D_/2
  c2 <- (diff(X_)^2+diff(Y_)^2)/2
  if (a^2-c2 < 0){
    print('Error: max dist is too small regarding the distance between points. Change it. ')
  }else{
    b <- sqrt(a^2-c2)
    phi <- atan( diff(Y_)/diff(X_))

    theta <- seq(0, 2 * pi, length = npts)
    xh <-  a * cos(theta)
    yh <-  b * sin(theta)
    co <- cos(phi)
    si <- sin(phi)
    x  <- x0 + co*xh - si*yh
    y  <- y0 + si*xh + co*yh
    # x  <- x0 + xh
    # y  <- y0 + yh

    return(list(x = x, y = y))
  }
}

jPPA <- function(A,B,Dmax,factor.color=1,tam.cel=3){ 
  # in contrast to the jppa function in the wildlifeTG package, there is no need to have ltraj objects, nor spatial polygons
  # by gridding the space to compute the areas, the function efficiently computes the intersection and union areas without encountering problems related to disjoint sets
  XA <- A[,1]
  YA <- A[,2]
  XB <- B[,1]
  YB <- B[,2]
  testA <- lapply(1:(length(XA)-1), function(i){ getEllipse(X_=XA[c(i,i+1)],Y_=YA[c(i,i+1)],D_=Dmax)})
  testB <- lapply(1:(length(XB)-1), function(i){ getEllipse(X_=XB[c(i,i+1)],Y_=YB[c(i,i+1)],D_=Dmax)})
  range.x <- range(do.call( rbind, testA)[,1],do.call( rbind, testB)[,1],na.rm = TRUE)
  range.y <- range(do.call( rbind, testA)[,2],do.call( rbind, testB)[,2],na.rm = TRUE)
  xgrid=seq(from=range.x[1],to=range.x[2],by=tam.cel)
  ygrid=seq(from=range.y[1],to=range.y[2],by=tam.cel)
  grilla.inter <- grilla.union <- matrix(0,ncol=length(xgrid),
                                         nrow=length(ygrid),byrow = FALSE)
  test.2 <- sapply(1:length(testA),function(i){
    UnionAB <- polyclip(testA[[i]], testB[[i]], op="union", fillA = 'positive', fillB = 'positive', x0 = mean(c(XA[i],XB[i])), y0 = mean(c(YA[i],YB[i])), eps=1e-5 )
    InterAB <- polyclip(testA[[i]], testB[[i]], op="intersection", fillA = 'positive', fillB = 'positive', x0 = mean(c(XA[i],XB[i])), y0 = mean(c(YA[i],YB[i])), eps=1e-5 )
    if (length(InterAB) > 0){
      poly.inter.sub <- polygrid(xgrid=seq(from=range.x[1],to=range.x[2],by=tam.cel),
                                 ygrid=seq(from=range.y[1],to=range.y[2],by=tam.cel),
                                 borders = InterAB[[1]], vec=TRUE)
      toto <- matrix(as.numeric(poly.inter.sub$vec.inout),ncol=length(xgrid),
                     nrow=length(ygrid),byrow = FALSE)
      grilla.inter <<- toto + grilla.inter
    }
    test <- sapply(1:length(UnionAB),function(j){
      poly.union.sub <- polygrid(xgrid=seq(from=range.x[1],to=range.x[2],by=tam.cel),
                                 ygrid=seq(from=range.y[1],to=range.y[2],by=tam.cel),
                                 borders = UnionAB[[j]], vec=TRUE)
      toto <- matrix(as.numeric(poly.union.sub$vec.inout),ncol=length(xgrid),
                     nrow=length(ygrid),byrow = FALSE)
      grilla.union <<- toto + grilla.union
    })

  })
  grilla.den <- sum(grilla.union>0)
  grilla.num <- sum(grilla.inter>0)

  jPPA <- grilla.num/grilla.den

  return(jPPA)
}

gamma.s <- function(ind.s,epsilon,cono.x,WC.edge){ # necessary for WCavglss and WCmaxlss indices
  L <- length(cono.x)
  beta <- alpha <- rep(0,L+1)
  for (j in 2:(L+1)){
    if (sum(WC.edge[ind.s,] > epsilon,na.rm = TRUE) > 0){
      if (is.na(WC.edge[ind.s,j-1]) == FALSE && (WC.edge[ind.s,j-1]> epsilon) == TRUE){
        alpha[j] <- max(alpha[j-1],1+beta[j-1])
        beta[j] <- 1 + beta[j-1]
      }else{
        alpha[j] <- alpha[j-1]
        beta[j] <- 0
      }
    }else{
      beta[j] <- alpha[j] <- 0
    }
  }
  return(alpha[L+1])
}

Cs.function <- function(A,B,distances){
  do <- mean(distances)
  distances.de <- as.matrix(pdist.pdist(X = A,Y = B))
  de <- mean(distances.de)
  Cs <- (de-do)/(de+do)
  p <- binom.test(x=sum(distances<=de),n=length(distances),p=0.5,alternative="greater")$p.value
  Cs.2 <- (1-p)/(1-(0.5)^nrow(A))
  return(list(Cs=Cs,Cs.2=Cs.2))
}

CSEM.function <- function(distances,longueur,delta){
  m <- 1
  PM <- 1
  P <- NULL
  while (PM[m] != 0 & m < longueur){
    seq1 <- 1:(longueur-m)
    seq2 <- seq1+m
    pairs <- as.matrix(cbind(seq1,seq2))
    P[[m]] <- as.numeric(sapply(1:length(seq1),function(x){max(distances[pairs[x,1]:pairs[x,2]])}) < delta)
    PM <- c(PM,sum(P[[m]]))
    (m <- m + 1)
  }
  PM <- PM[-1]

  if (PM[length(P)] == 0){
    CSEM <- (length(P)-1)/(longueur-1)
  }else{
    CSEM <- (length(P))/(longueur-1)
  }
  return(CSEM)
}

DI.function <- function(A,B,delta.DI,longueur,distance1,distance2){ # adapted from the wildlifeDI package without requiring ltraj objects
  DId0 <- 1 - (abs(distance1-distance2)/(distance1+distance2))^delta.DI
  DId0[which(distance1+distance2 == 0)] <- 1
  DI0.d <- mean(DId0,na.rm=TRUE)
  x1 <- A[-1, ]
  x2 <- A[-longueur, ]
  dx <- c(x1[,1] - x2[,1])
  dy <- c(x1[,2] - x2[,2])
  abs.angle.A <- atan2(dy, dx)
  x1 <- B[-1, ]
  x2 <- B[-longueur, ]
  dx <- c(x1[,1] - x2[,1])
  dy <- c(x1[,2] - x2[,2])
  abs.angle.B <- atan2(dy, dx)
  DItheta <- cos(abs.angle.A - abs.angle.B)
  DItheta[which(is.na(abs.angle.A)== TRUE & is.na(abs.angle.B) == TRUE)] <- 1
  DItheta[which(is.na(DItheta))] <- 0
  DI.theta <- mean(DItheta,na.rm=TRUE)

  # DI
  DI0 <- mean(DItheta*DId0,na.rm = TRUE)
  return(list(DI=DI0,DI.d=DI0.d,DI.theta=DI.theta))
}

# We basically modify the WTC function in the biwavelet package to avoid fitting an arima (which is not always possible)
wtc.modif <- function (d1, d2, pad = TRUE, dj = 1/12, s0 = 2 * dt, J1 = NULL,
                       max.scale = NULL, mother = "morlet", param = -1, lag1 = NULL,
                       sig.level = 0.95, sig.test = 0, nrands = 300, quiet = FALSE){
  mother <- match.arg(tolower(mother), MOTHERS)
  checked <- check.data(y = d1, x1 = d2)
  xaxis <- d1[, 1]
  dt <- checked$y$dt
  t <- checked$y$t
  n <- checked$y$n.obs
  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt
    }
    J1 <- round(log2(max.scale/s0)/dj)
  }
  wt1 <- wt(d = d1, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale,
            mother = mother, param = param, sig.level = sig.level,
            sig.test = sig.test, lag1 = lag1[1],do.sig = FALSE)
  wt2 <- wt(d = d2, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale,
            mother = mother, param = param, sig.level = sig.level,
            sig.test = sig.test, lag1 = lag1[2],do.sig = FALSE)
  d1.sigma <- sd(d1[, 2], na.rm = T)
  d2.sigma <- sd(d2[, 2], na.rm = T)
  s.inv <- 1/t(wt1$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt1$wave))
  smooth.wt1 <- smooth.wavelet(s.inv * (abs(wt1$wave)^2), dt,
                               dj, wt1$scale)
  smooth.wt2 <- smooth.wavelet(s.inv * (abs(wt2$wave)^2), dt,
                               dj, wt2$scale)
  coi <- pmin(wt1$coi, wt2$coi, na.rm = T)
  CW <- wt1$wave * Conj(wt2$wave)
  CW.corr <- (wt1$wave * Conj(wt2$wave) * max(wt1$period))/matrix(rep(wt1$period,
                                                                      length(t)), nrow = NROW(wt1$period))
  power <- abs(CW)^2
  power.corr <- (abs(CW)^2 * max.scale)/matrix(rep(wt1$period,
                                                   length(t)), nrow = NROW(wt1$period))
  smooth.CW <- smooth.wavelet(s.inv * (CW), dt, dj, wt1$scale)
  rsq <- abs(smooth.CW)^2/(smooth.wt1 * smooth.wt2)
  phase <- atan2(Im(CW), Re(CW))
  if (nrands > 0) {
    signif <- wtc.sig(nrands = nrands, lag1 = lag1, dt = dt,
                      n, pad = pad, dj = dj, J1 = J1, s0 = s0, max.scale = max.scale,
                      mother = mother, sig.level = sig.level, quiet = quiet)
  }else {
    signif <- NA
  }
  results <- list(coi = coi, wave = CW, wave.corr = CW.corr,
                  power = power, power.corr = power.corr, rsq = rsq, phase = phase,
                  period = wt1$period, scale = wt1$scale, dt = dt, t = t,
                  xaxis = xaxis, s0 = s0, dj = dj, d1.sigma = d1.sigma,
                  d2.sigma = d2.sigma, mother = mother, type = "wtc", signif = signif)
  class(results) <- "biwavelet"
  return(results)
}

WC.function <- function(distance1,distance2,sl,su,epsilon){
  wtc.t1t2 <- wtc.modif(cbind(1:length(distance1),distance1),cbind(1:length(distance2),distance2), nrands = 0,sig.level = 0.95,quiet = TRUE)
  cono.x <- 1:length(wtc.t1t2$coi)
  matriz.ind.edge <- matrix(NA,nrow=length(wtc.t1t2$scale),ncol=length(cono.x))
  test <- sapply(1:length(wtc.t1t2$coi),function(x){
    matriz.ind.edge[which(wtc.t1t2$scale < wtc.t1t2$coi[x]),x] <<- 1
  })
  WC.edge <- wtc.t1t2$rsq*matriz.ind.edge
  WC.edge.scale <- apply(WC.edge,MARGIN = 1,FUN=mean,na.rm=TRUE)
  if (is.null(sl)){
    sl = wtc.t1t2$scale[1]
  }
  if (is.null(su)){
    su = max(wtc.t1t2$scale)
  }
  ind.scales <- which(wtc.t1t2$scale >= sl & wtc.t1t2$scale <= su)
  WC.avg <- mean(WC.edge.scale[ind.scales],na.rm=TRUE)
  WC.max.scale <- apply(WC.edge,MARGIN = 1,FUN=max,na.rm=TRUE)
  WC.max <- max(WC.max.scale[ind.scales])
  gamma.ns <- sapply(ind.scales,function(x){
    gamma.s(x,epsilon,cono.x,WC.edge)/sum(is.na(WC.edge[x,])==FALSE)
  })
  WC.avg.lss <- mean(gamma.ns,na.rm=TRUE)
  WC.max.lss <- max(gamma.ns,na.rm=TRUE)
  return(list(WC.avg=WC.avg,WC.max=WC.max,WC.avg.lss=WC.avg.lss,WC.max.lss=WC.max.lss))
}

Sigma=  1 
Sd.Sigma.1 = 1
Sd.Sigma.2 = 1 
Amount.Scale = 0.05
Amount.Translate = seq(from=0,to=1,by=0.05) # that is omega
Nstep <- 100 
NSim <- 1000
delta <- 5
Dmax <- 20*1.852 # 20 knots for data in meters (vmax*delta.t)
delta.DI <- 1
# parameters for WC
sl=2
su=6
epsilon=0.5


Cs1 <- Cs2 <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
jPPA.index <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
Lon <- Lat <- Lonlat <- Speed <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
CSEM <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
Prox1 <- Prox2 <- Prox3 <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
DI.stat <- DId <- DItheta <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))
WC.avg <- WC.max <- WC.avg.lss <- WC.max.lss <- matrix(NA, nrow=NSim,ncol=length(Amount.Translate))

for (iw in 1:length(Amount.Translate)){
  for (i6 in 1:NSim){
    set.seed(i6+50000)
    A  <- cbind(cumsum(rnorm(Nstep, mean=0, sd=Sigma)),cumsum(rnorm(Nstep, mean=0, sd=Sigma)))
    probas <- rbinom(n=Nstep,size=100,prob = Amount.Translate[iw])/100
    B  <- probas*jitter(A,amount = Amount.Scale) + (1-probas)*cbind(cumsum((rnorm(Nstep, mean=0, sd=Sigma))),cumsum(rnorm(Nstep, mean=0, sd=Sigma)))
    distances <- diag(as.matrix(pdist.pdist(X = A,Y = B)))
    longueur <- length(A[,1])
    distance1 <- sqrt((diff(A[,1]))^2 + (diff(A[,2]))^2)
    distance2 <- sqrt((diff(B[,1]))^2 + (diff(B[,2]))^2)
    Prox1[i6,iw] <- sum(as.numeric(distances < 1))/longueur
    Prox2[i6,iw] <- sum(as.numeric(distances < 3))/longueur
    Prox3[i6,iw] <- sum(as.numeric(distances < 5))/longueur
    Cs.res <- Cs.function(A=A,B=B,distances=distances)
    Cs1[i6,iw] <- Cs.res$Cs
    Cs2[i6,iw] <- Cs.res$Cs.2
    jPPA.index[i6,iw] <- jPPA(A = A,B=B,Dmax=Dmax)
    Lon[i6,iw] = cor(A[,1],B[,1])
    Lat[i6,iw] = cor(A[,2],B[,2])
    Lonlat[i6,iw] = (Lon[i6,iw]+Lat[i6,iw])/2
    Speed[i6,iw] = cor(distance1,distance2,use='pairwise.complete.obs')
    CSEM[i6,iw] <- CSEM.function(distances=distances,longueur = longueur,delta=delta)
    DI <- DI.function(A=A,B=B,delta.DI = delta.DI,longueur = longueur,distance1 = distance1,distance2 = distance2)
    DId[i6,iw] <- DI$DI.d
    DItheta[i6,iw] <- DI$DI.theta
    DI.stat[i6,iw] <- DI$DI
    WC <- WC.function(distance1 = distance1,distance2 = distance2,sl=sl,su=su,epsilon=epsilon)
    WC.avg[i6,iw] <- WC$WC.avg
    WC.max[i6,iw] <- WC$WC.max
    WC.avg.lss[i6,iw] <- WC$WC.avg.lss
    WC.max.lss[i6,iw] <- WC$WC.max.lss
    
  }
  
}

