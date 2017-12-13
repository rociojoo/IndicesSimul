
library('pdist')
pdist.pdist <- pdist::pdist
library('polyclip','colorspace')
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

Cs.function <- function(A,B,distances){
  do <- mean(distances)
  distances.de <- as.matrix(pdist.pdist(X = A,Y = B))
  de <- mean(distances.de)
  Cs <- (de-do)/(de+do)
  Cs.2 <- (1-binom.test(x=sum(distances<=de),n=length(distances),p=0.5,alternative="greater")$p.value) 
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

# An example for indices derivation:
# first, simulating a dyad:
Nstep <- 100
A  <- cbind(cumsum(rnorm(Nstep, mean=0, sd=1)),cumsum(rnorm(Nstep, mean=0, sd=1)))
probas <- rbinom(n=Nstep,size=100,prob = 0.5)/100
B  <- probas*jitter(A,amount = 0.05) + (1-probas)*cbind(cumsum((rnorm(Nstep, mean=0, sd=1))),
                                                        cumsum(rnorm(Nstep, mean=0, sd=1)))

Dmax <- 10
delta.DI <- 1
distances <- diag(as.matrix(pdist.pdist(X = A,Y = B)))
longueur <- length(A[,1])
distance1 <- sqrt((diff(A[,1]))^2 + (diff(A[,2]))^2)
distance2 <- sqrt((diff(B[,1]))^2 + (diff(B[,2]))^2)
Prox1 <- sum(as.numeric(distances < 1))/longueur
Prox2 <- sum(as.numeric(distances < 3))/longueur
Prox3 <- sum(as.numeric(distances < 5))/longueur
Cs.res <- Cs.function(A=A,B=B,distances=distances) # I used to compute another statistic (Cs2) in the same function
Cs1 <- Cs.res$Cs
jPPA.index <- jPPA(A = A,B=B,Dmax=Dmax)
Lon = cor(A[,1],B[,1])
Lat = cor(A[,2],B[,2])
Lonlat = (Lon+Lat)/2
Speed = cor(distance1,distance2,use='pairwise.complete.obs')
CSEM1 <- CSEM.function(distances=distances,longueur = longueur,delta=1)
CSEM2 <- CSEM.function(distances=distances,longueur = longueur,delta=2)
CSEM3 <- CSEM.function(distances=distances,longueur = longueur,delta=3)
DI <- DI.function(A=A,B=B,delta.DI = delta.DI,longueur = longueur,distance1 = distance1,distance2 = distance2)
DId <- DI$DI.d
DItheta <- DI$DI.theta
DI.stat <- DI$DI


