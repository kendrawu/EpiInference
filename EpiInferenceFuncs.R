fn.one <- function (X0=10,R0=1.8,Tg=2.5,k=25,N=1000,noreals=20,maxdays=60) {
  fn.fig1(R0=1.8,Tg=2.5, k=25,N=1000,noreals=20,maxdays=60,X0=10)
  fn.fig1(R0=3.0,Tg=6.25,k=25,N=1000,noreals=20,maxdays=60,X0=10)
  fn.fig2(R0=1.8,Tg=2.5, k=10,N=1000,noreals=10, maxdays=60,X0=50)
  fn.fig2(R0=1.8,Tg=2.5, k=10,N=1000,noreals=20, maxdays=60,X0=50)
  fn.fig2(R0=1.8,Tg=2.5, k=10,N=1000,noreals=100,maxdays=60,X0=50)
  fn.fig2(R0=1.8,Tg=2.5, k=25,N=1000,noreals=10, maxdays=60,X0=50)
  fn.fig2(R0=1.8,Tg=2.5, k=25,N=1000,noreals=20, maxdays=60,X0=50)
  fn.fig2(R0=1.8,Tg=2.5, k=25,N=1000,noreals=100,maxdays=60,X0=50)
  fn.fig2(R0=3.0,Tg=6.25,k=10,N=1000,noreals=10, maxdays=60,X0=100)
  fn.fig2(R0=3.0,Tg=6.25,k=10,N=1000,noreals=20, maxdays=60,X0=100)
  fn.fig2(R0=3.0,Tg=6.25,k=10,N=1000,noreals=100,maxdays=60,X0=100)
  fn.fig2(R0=3.0,Tg=6.25,k=25,N=1000,noreals=10, maxdays=60,X0=100)
  fn.fig2(R0=3.0,Tg=6.25,k=25,N=1000,noreals=20, maxdays=60,X0=100)
  fn.fig2(R0=3.0,Tg=6.25,k=25,N=1000,noreals=100,maxdays=60,X0=100)
  fn.R0resid1(R0set=1.8,Tgset=2.5,seed=50,dlike=50,maxdays=60,k=25,notests=50)
  fn.Tgresid1(R0set=1.8,Tgset=2.5,seed=50,dlike=50,maxdays=60,k=10,notests=50)
  fn.R0resid2(Tgset=2.5,seed=50,dlike=50,maxdays=60,k=10,notests=10)
  fn.Tgresid2(R0set=1.8,seed=50,dlike=50,maxdays=60,k=10,notests=10)
}


## Functions for Exp242p1script.r##

fn.epicurve.wj <- function(y,theta,N,k,i) {
  
  ## Eqn 1: Generate secondary cases w_j from one infector of day i
  
  ## Compute initial conditions
  R0 <- theta[1]
  Tg <- theta[2]
  S <- N - sum(y) #number of susceptibles in the population
  w_j <- 0 #initial cumulative incidence
  
  # Constrains index of infectees
  tmax <- length(y) #max observed days
  #infector index begins, i
  i.max <- min(i+k,tmax) #infector index ends
  
  if (S>0) {
    #w_j <- w_j + y[i] * R0 
    w_j <- w_j + rpois(1,R0)
    if ((sum(y)+w_j)>N) w_j <- S
  } else w_j <- 0
  
  #print(noquote(paste(i,S,y_j)))
  
  ## Output
  return(w_j)
} 


fn.epicurve.y.oneday <- function(y,theta,N,k,i) {
  
  ## Eqn 1: Generate secondary cases according to expected cases from fn.epicurve.wj (eqn 2) and 
  ## binomial distribution and distribute them along the time course according to Poisson distribution
  
  ## Compute initial conditions
  R0 <- theta[1]
  Tg <- theta[2]
  tmax <- length(y)
  
  S <- N-sum(y)
  y_i <- y[i]
  p <- 1
  
  while (p <= y_i) {
    
    # Generate secondary cases of _each_ infector on day i
    # Each generates different number of secondary cases and has different serial interval
    #w_j <- fn.epicurve.wj(y,theta,N,k,i)
    w_j <- rpois(1,R0)
    
    # Constrain j, infectee index
    j <- i + rpois(1,Tg)
    j.max <- min(i+k,tmax)
    while (j>j.max) j <- i + rpois(1,Tg)
    
    # Put secondary cases onto y_j
    if (S>0 & w_j>0) {
      Z_j <- rbinom(1,S,w_j/N) # eqn 1
      if ((sum(y)+Z_j)>N) Z_j <- S
      y[j] <- y[j] + Z_j
    }
    
    p <- p+1
  }
  
  ## Output
  return(y)
  
}


fn.epicurve.y.oneday.new <- function(y,theta,N,k,i) {
  
  ## Eqn 1: Generate secondary cases according to expected cases from fn.epicurve.wj (eqn 2) and 
  ## binomial distribution and distribute them along the time course according to Poisson distribution
  
  ## Compute initial conditions
  R0 <- theta[1]
  Tg <- theta[2]
  tmax <- length(y)
  
  S <- N-sum(y)
  p <- 1
  
  while (p <= y[i]) {
    
    # Generate secondary cases of _each_ infector on day i
    # Each generates different number of secondary cases and has different serial interval
    #w_j <- fn.epicurve.wj(y,theta,N,k,i)
    w_j <- rpois(1,R0)
    
    # Put secondary cases onto y_j
    if (S>0 & w_j>0) {
      Z_j <- rbinom(1,S,w_j/N) # eqn 1
      q <- 1
      
      while (q <= Z_j) {
        
        # Constrain j, infectee index
        j <- i + rpois(1,Tg)
        j.max <- min(i+k,tmax)
        while (j>j.max) j <- i + rpois(1,Tg)
        
        y[j] <- y[j] + 1
        q <- q+1
      }
    }  
    p <- p+1
  }
  
  ## Output
  return(y)
  
}

fn.epicurve <- function(y0,theta,N,k,tmax) {
  
  ## Eqn 2: Put all the y_j along the time course, where t={1,...,j,...tmax}
  
  y <- vector(mode="numeric", length=tmax); y[]<-0
  
  ## plant seeds
  if (sum(y0)<1) stop ('seed has to be at least 1')
  seed.length <- length(y0)
  #t.begin <- seed.length+1
  t.begin <- 1 #infector index begins
  for (q in 1:seed.length) y[q] <- y0[q]
  
  ## Generate an epidemic from t.begin to tmax
  for (i in t.begin: tmax) {
    y <- fn.epicurve.y.oneday.new(y,theta,N,k,i)
    #print(noquote('Y is:')); print(y)
  }
  
  ## Output
  return(y)
}


## -------------------------------------------------------------------------------------------- ##
# Simulate some epidemics and make a plot (Fig 1)
fn.fig1 <- function(
  R0,
  Tg,
  N,
  noreals,
  maxdays,
  k,
  X0,
  filenamestem="Exp242Fig1",
  ymax=100,
  nolines=5) {
  
  # Make a matrix for the number of runs
  output <- matrix(nrow=noreals,ncol=maxdays)
  
  # Run the simulation the correct number of times
  i <- 0
  if (nolines>noreals) nolines <- noreals
  while (i <= noreals) {
    theta <- c(R0,Tg)
    #output[i,] <- kp1.sim.all.days(seed=seed,maxdays=maxdays,N=N,R0=R0,Tg=Tg,k=k)
    output[i,] <- fn.epicurve(y0=X0,theta,N=N,k=k,tmax=maxdays)
    i <- i + 1
  }  
  
  # Calculate the mean and the 95% CIs
  vecMean <- apply(output,2,mean)
  vecLower <- apply(output,2,quantile,probs=c(0.05))
  vecUpper <- apply(output,2,quantile,probs=c(0.95))
  
  # Return a plot
  pdf(paste(filenamestem,"_R0",R0,"_Tg",Tg,".pdf",sep=""))
  plot(1:2,type="n",xlim=c(0,maxdays),ylim=c(0,ymax),axes=FALSE,
       xlab="Days of epidemic",
       ylab="Incidence")
  axis(1)
  axis(2)
  for (j in sample(1:noreals,nolines)) points(output[j,],type="l",col="blue",lwd=0.5)
  points(1:maxdays,vecMean,type="l",col="black",lwd=2)
  points(1:maxdays,vecLower,type="l",col="red",lwd=1,lty=2)
  points(1:maxdays,vecUpper,type="l",col="red",lwd=1,lty=2)
  legend("topright",c("Examples","Mean","5% percentile","95% percentile"),
         col=c("blue","black","red","red"), lty=c(1,1,2,2), cex=0.6)
  dev.off()
  
}

# Generate all the required runs for a 4 part figure
fn.fig2 <- function(
  R0,
  Tg,
  k,
  filenamestem="Exp242fig2",
  noreals,
  N,
  X0,
  centiles=c(0,1),
  maxdays=maxdays,
  nolines=8
) {
  
  # Define data array
  estimates <- array(dim=c(maxdays-1,nrow=noreals,4))
  
  for (i  in 1:noreals) {
    print(noquote(paste('iteration:',i)))
    #test <- kp1.sim.all.days(seed=10,maxdays=maxdays,N=1000,R0=R0,Tg=Tg,k=10)
    theta <- c(R0,Tg)
    test <- fn.epicurve(y0=X0,theta,N=N,k=k,tmax=maxdays)
    for (j in 2:maxdays) {
      startvec <- c(runif(1,min=0.5,max=5),runif(1,min=1.01,max=5)) 
      x <- optim(
        startvec,
        kp1.like, #fn.like,
        y=test,d=j,N=N,k=k,R0=-1,Tg=-1,
        method="L-BFGS-B",
        lower=c(0.5,1.001),upper=c(10.0,10.0),
        control=c(fnscale=-1)
      )
      estimates[j-1,i,] <- x$par
      
      # Example optims for 
      y <- optim(
        c(runif(1,min=0.5,max=5)),
        kp1.like, #fn.like,
        y=test,d=j,N=N,k=k,R0=-1,Tg=Tg,
        method="Brent",
        lower=0.5,upper=10,
        control=c(fnscale=-1)
      )
      estimates[j-1,i,3] <- y$par
      
      # Example optims for 
      z <- optim(
        c(runif(1,min=1.001,max=8)),
        kp1.like, #fn.like,
        y=test,d=j,N=N,k=k,R0=R0,Tg=-1,
        method="Brent",
        lower=1.01,upper=10,
        control=c(fnscale=-1)
      )
      estimates[j-1,i,4] <- z$par
      
    }
  }
  
  # Univariate R0
  if (nolines>noreals) nolines <- noreals
  est.tmp <- t(estimates[,,3])
  vecMean <- apply(est.tmp,2,mean)
  vecMedian <- apply(est.tmp,2,quantile,probs=c(0.5))
  vecLower <- apply(est.tmp,2,quantile,probs=c(0.05))
  vecUpper <- apply(est.tmp,2,quantile,probs=c(0.95))
  print(noquote('Univarate R0: Mean, Median, 5%CI, 95%CI'))
  print(noquote(paste(vecMean[maxdays-1], vecMedian[maxdays-1], vecLower[maxdays-1], vecUpper[maxdays-1])))
  
  pdf(paste(filenamestem,"A_R0_",R0,"_",Tg,"_rep",noreals,".pdf",sep=""))
  plot(1:2,type="n",xlim=c(0,maxdays),ylim=c(0,10),axes=FALSE,
       main="(a) Univariate R0 estimates",
       xlab="Days of epidemic",
       ylab="Estimated R0")
  axis(1)
  axis(2)
  #for (i in 1:noreals) points(2:maxdays,estimates[,i,3],type="l",col="grey",lwd=0.5)
  for (i in sample(1:noreals,nolines)) points(2:maxdays,estimates[,i,3],type="l",col="grey",lwd=0.5)
  lines(vecMean,col="black")
  lines(vecMedian,col="red")
  lines(vecLower,col="green",lty=2)
  lines(vecUpper,col="green",lty=2)
  abline(h=R0,col="dark blue",lwd=2)
  legend("topright",c("True","Examples","Mean","Median","5% percentile","95% percentile"),
         col=c("dark blue","grey","black","red","green","green"), lty=c(1,1,1,1,2,2), cex=0.6)
  dev.off()
  
  # Univariate Tg
  est.tmp <- t(estimates[,,4])
  vecMean <- apply(est.tmp,2,mean)
  vecMedian <- apply(est.tmp,2,quantile,probs=c(0.5))
  vecLower <- apply(est.tmp,2,quantile,probs=c(0.05))
  vecUpper <- apply(est.tmp,2,quantile,probs=c(0.95))
  print(noquote('Univarate Tg: Mean, Median, 5%CI, 95%CI'))
  print(noquote(paste(vecMean[maxdays-1], vecMedian[maxdays-1], vecLower[maxdays-1], vecUpper[maxdays-1])))
  
  pdf(paste(filenamestem,"B_Tg_",R0,"_",Tg,"_rep",noreals,".pdf",sep=""))
  plot(1:2,type="n",xlim=c(0,maxdays),ylim=c(0,10),axes=FALSE,
       main="(b) Univariate Tg estimates",
       xlab="Days of epidemic",
       ylab="Estimated Tg")
  axis(1)
  axis(2)
  #for (i in 1:noreals) points(2:maxdays,estimates[,i,4],type="l",col="grey",lwd=0.5)
  for (i in sample(1:noreals,nolines)) points(2:maxdays,estimates[,i,4],type="l",col="grey",lwd=0.5)
  lines(vecMean,col="black")
  lines(vecMedian,col="red")
  lines(vecLower,col="green",lty=2)
  lines(vecUpper,col="green",lty=2)
  abline(h=Tg,col="dark blue",lwd=2)
  legend("topright",c("True","Examples","Mean","Median","5% percentile","95% percentile"),
         col=c("dark blue","grey","black","red","green","green"), lty=c(1,1,1,1,2,2), cex=0.6)
  dev.off()
  
  # Bivariate R0
  est.tmp <- t(estimates[,,1])
  vecMean <- apply(est.tmp,2,mean)
  vecMedian <- apply(est.tmp,2,quantile,probs=c(0.5))
  vecLower <- apply(est.tmp,2,quantile,probs=c(0.05))
  vecUpper <- apply(est.tmp,2,quantile,probs=c(0.95))
  print(noquote('Bivarate R0: Mean, Median, 5%CI, 95%CI'))
  print(noquote(paste(vecMean[maxdays-1], vecMedian[maxdays-1], vecLower[maxdays-1], vecUpper[maxdays-1])))
  
  pdf(paste(filenamestem,"C_R0_",R0,"_",Tg,"_rep",noreals,".pdf",sep=""))
  plot(1:2,type="n",xlim=c(0,maxdays),ylim=c(0,10),axes=FALSE,
       main="(c) Bivariate R0 estimates",
       xlab="Days of epidemic",
       ylab="Estimated R0")
  axis(1)
  axis(2)
  #for (i in 1:noreals) points(2:maxdays,estimates[,i,1],type="l",col="grey",lwd=0.5)
  for (i in sample(1:noreals,nolines)) points(2:maxdays,estimates[,i,1],type="l",col="grey",lwd=0.5)
  lines(vecMean,col="black")
  lines(vecMedian,col="red")
  lines(vecLower,col="green",lty=2)
  lines(vecUpper,col="green",lty=2)
  abline(h=R0,col="dark blue",lwd=2)
  legend("topright",c("True","Examples","Mean","Median","5% percentile","95% percentile"),
         col=c("dark blue","grey","black","red","green","green"), lty=c(1,1,1,1,2,2), cex=0.6)
  dev.off()
  
  # Bivariate Tg
  est.tmp <- t(estimates[,,2])
  vecMean <- apply(est.tmp,2,mean)
  vecMedian <- apply(est.tmp,2,quantile,probs=c(0.5))
  vecLower <- apply(est.tmp,2,quantile,probs=c(0.05))
  vecUpper <- apply(est.tmp,2,quantile,probs=c(0.95))
  print(noquote('Bivarate Tg: Mean, Median, 5%CI, 95%CI'))
  print(noquote(paste(vecMean[maxdays-1], vecMedian[maxdays-1], vecLower[maxdays-1], vecUpper[maxdays-1])))
  
  pdf(paste(filenamestem,"D_Tg_",R0,"_",Tg,"_rep",noreals,".pdf",sep=""))
  plot(1:2,type="n",xlim=c(0,maxdays),ylim=c(0,10),axes=FALSE,
       main="(d) Bivariate Tg estimates",
       xlab="Days of epidemic",
       ylab="Estimated Tg")
  axis(1)
  axis(2)
  #for (i in 1:noreals) points(2:maxdays,estimates[,i,2],type="l",col="grey",lwd=0.5)
  for (i in sample(1:noreals,nolines)) points(2:maxdays,estimates[,i,2],type="l",col="grey",lwd=0.5)
  lines(vecMean,col="black")
  lines(vecMedian,col="red")
  lines(vecLower,col="green",lty=2)
  lines(vecUpper,col="green",lty=2)
  abline(h=Tg,col="dark blue",lwd=2)
  legend("topright",c("True","Examples","Mean","Median","5% percentile","95% percentile"),
         col=c("dark blue","grey","black","red","green","green"), lty=c(1,1,1,1,2,2), cex=0.6)
  dev.off()
  
}

# kp1.debug.1()
fn.R0resid1 <- function(R0set, Tgset, seed, dlike, maxdays, k, notests) {
  
  xvals <- seq(min(0.1,R0set-1.5),R0set+1.5,0.1)
  yvals <- vector(mode="numeric",length=length(xvals))
  resid <- vector(mode="numeric",length=notests)
  #j <- 0
  
  for (j in 1:notests) {
    print(noquote(paste('iteration:',j)))
    #test <- kp1.sim.all.days(maxdays=90,N=1000,R0=R0set,Tg=Tgset,k=10)
    theta <- c(R0set,Tgset)
    test <- fn.epicurve(y0=seed,theta,N=1000,k=k,tmax=maxdays+1)
    # plot(test,type="l")
    
    for (i in 1:length(xvals)) {
      yvals[i] <- kp1.log.like.upto.day(test,dlike,N=1000,R0=xvals[i],Tg=Tgset,k=k)
    }
    
    resid[j] <- R0set-xvals[match(max(yvals),yvals)]
    #j <- j + 1
  }
  
  ## Output
  resid.tab <- table(resid)
  resid.prop <- prop.table(resid.tab)*100 #Convert frequency to percentage
  
  pdf(paste(filenamestem="Exp242Fig4","A_R0",".pdf",sep=""))
  plot(resid.prop,xlim=c(-1,1),ylim=c(0,50),axes=FALSE,
       xlab="R0 residual",
       ylab="Occurrences (%)")
  axis(1)
  axis(2)
  dev.off()
}


fn.Tgresid1 <- function(R0set, Tgset, seed, dlike, maxdays, k, notests) {
  
  xvals <- seq(max(1.1,Tgset-1.5),Tgset+1.5,0.1)
  yvals <- vector(mode="numeric",length=length(xvals))
  resid <- vector(mode="numeric",length=notests)
  #j <- 1
  
  for (j in 1:notests) {
    print(noquote(paste('iteration:',j)))
    #test <- kp1.sim.all.days(maxdays=90,N=1000,R0=R0set,Tg=Tgset,k=10)
    theta <- c(R0set,Tgset)
    test <- fn.epicurve(y0=seed,theta,N=1000,k=k,tmax=maxdays+1)
    
    for (i in 1:length(xvals)) {
      yvals[i] <- kp1.log.like.upto.day(test,dlike,N=1000,R0=R0set,Tg=xvals[i],k=k)
    }

    resid[j] <- Tgset-xvals[match(max(yvals),yvals)]
    #j <- j + 1
  }
  

  ## Output
  resid.tab <- table(resid)
  resid.prop <- prop.table(resid.tab)*100 #Convert frequency to percentage

  pdf(paste(filenamestem="Exp242Fig4","B_Tg",".pdf",sep=""))
  plot(resid.prop,xlim=c(-3,3),ylim=c(0,50),axes=FALSE,
       xlab="Tg residual",
       ylab="Occurrences (%)")
  axis(1)
  axis(2)
  dev.off()
}


fn.R0resid2 <- function(Tgset,seed,dlike,maxdays,k,notests) {
  
  r0.range <- seq(1,3,0.1)
  resid <- matrix(0, nrow=length(r0.range), ncol=notests)
  resid.med <- vector(mode="numeric",length=length(r0.range))
  resid.lower <- vector(mode="numeric",length=length(r0.range))
  resid.upper <- vector(mode="numeric",length=length(r0.range))
  
  
  for (p in 1:length(r0.range)) {
    
    for (j in 1:notests) {
      print(noquote(paste('iteration:',p,j)))
      # test <- kp1.sim.all.days(maxdays=90,N=1000,R0=R0set,Tg=Tgset,k=10)
      theta <- c(R0=r0.range[p],Tgset)
      test <- fn.epicurve(y0=seed,theta,N=1000,k=k,tmax=maxdays+1)
      
      # Set range of trial r0 (regardless of R0set=r0.range[p])
      xvals <- seq(max(0.1,r0.range[p]-1.5),r0.range[p]+1.5,0.1)
      yvals <- vector(mode="numeric",length=length(xvals))
      
      for (i in 1:length(xvals)) {
        yvals[i] <- kp1.log.like.upto.day(test,dlike,N=1000,R0=xvals[i],Tg=Tgset,k=k)
      }
      
      resid[p,j] <- r0.range[p]-xvals[match(max(yvals),yvals)]
    }  
    
    resid.vec <- resid[p,]
    resid.med[p] <- as.numeric(quantile(resid.vec,prob=0.5,na.rm=FALSE))
    resid.lower[p] <- as.numeric(quantile(resid.vec,prob=0.05,na.rm=FALSE))
    resid.upper[p] <- as.numeric(quantile(resid.vec,prob=0.95,na.rm=FALSE))
  }
  
  list(med=resid.med, lower=resid.lower, upper=resid.upper)
  pdf(paste(filenamestem="Exp242Fig4","C_R0",".pdf",sep=""))
  plot(r0.range,resid.med, type="l",ylim=c(-1,1), col="black",
       xlab="Trial R0",
       ylab="R0 residual")
  lines(r0.range,resid.lower,type="l",col="green",lty=2)
  lines(r0.range,resid.upper,type="l",col="green",lty=2)
  abline(h=0,col="dark blue",lwd=2)
  legend("topright",c("Median","5% percentile","95% percentile"),
         col=c("black","green","green"), lty=c(1,2,2), cex=0.6)
  dev.off()
}


fn.Tgresid2 <- function(R0set,seed,dlike,maxdays,k,notests) {
  
  Tg.range <- seq(1.1,10,0.25)
  resid <- matrix(0, nrow=length(Tg.range), ncol=notests)
  resid.med <- vector(mode="numeric",length=length(Tg.range))
  resid.lower <- vector(mode="numeric",length=length(Tg.range))
  resid.upper <- vector(mode="numeric",length=length(Tg.range))
  
  
  for (p in 1:length(Tg.range)) {
    
    for (j in 1:notests) {
      print(noquote(paste('iteration:',p,j)))
      #test <- kp1.sim.all.days(maxdays=90,N=1000,R0=R0set,Tg=Tgset,k=10)
      theta <- c(R0=R0set,Tg=Tg.range[p])
      #test <- 1
      #if (max(test)<20) test <- fn.epicurve(y0=seed,theta,N=1000,k=k,tmax=maxdays+1)
      test <- fn.epicurve(y0=seed,theta,N=1000,k=k,tmax=maxdays)
      
      # Set range of trial Tg (regardless of Tgset=Tg.range[p])
      xvals <- seq(max(1.1,Tg.range[p]-2),Tg.range[p]+2,0.1)
      yvals <- vector(mode="numeric",length=length(xvals))
      
      for (i in 1:length(xvals)) {
        yvals[i] <- kp1.log.like.upto.day(test,dlike,N=1000,R0=R0set,Tg=xvals[i],k=k)
      }
      
      resid[p,j] <- Tg.range[p]-xvals[match(max(yvals),yvals)]
    }  
    
    resid.vec <- resid[p,]
    resid.med[p] <- as.numeric(quantile(resid.vec,prob=0.5,na.rm=FALSE))
    resid.lower[p] <- as.numeric(quantile(resid.vec,prob=0.05,na.rm=FALSE))
    resid.upper[p] <- as.numeric(quantile(resid.vec,prob=0.95,na.rm=FALSE))
  }
  
  list(med=resid.med, lower=resid.lower, upper=resid.upper)
  pdf(paste(filenamestem="Exp242Fig4","D_Tg",".pdf",sep=""))
  plot(Tg.range,resid.med, type="l",ylim=c(-3,3), col="black",
       xlab="Trial Tg",
       ylab="Tg residual")
  lines(Tg.range,resid.lower,type="l",col="green",lty=2)
  lines(Tg.range,resid.upper,type="l",col="green",lty=2)
  abline(h=0,col="dark blue",lwd=2)
  legend("topright",c("Median","5% percentile","95% percentile"),
         col=c("black","green","green"), lty=c(1,2,2), cex=0.6)
  dev.off()
}

## -------------------------------------------------------------------------------------------- ##
kp1.sim.one.day <- function(d, y, N, R0, Tg, k) {
  rtn <- y
  exppic <- kp1.calc.exp.inf(d,y,N,R0,Tg,k)
  S_d <- N - sum(y[1:(d-1)])
  if (exppic > 0) {
    rtn[d] <- rtn[d] + rbinom(1,S_d,exppic/N)
  }
  rtn
}

# Compute number of expected potential infectious contacts of day 'currentinf' that may result in cases with onset on day d
kp1.calc.exp.inf <- function(d, y, N, R0, Tg, k) {
  
  maxdays <- length(y)
  if ((d > maxdays) | (d < 2)) stop("Problem in kp1.calc.exp.inf. d is out of range ", d)
  if (Tg <= 1) stop("Problem in kp1.calc.exp.inf: value of Tg is too small ", Tg)
  S_d <- N - sum(y[1:(d-1)]) #number of susceptible at the start of day d
  S <- S_d
  
  # Construct data structure
  yj.d <- vector(mode="numeric", length=maxdays); yj.d[]<-0
  I.sum <- vector(mode="numeric", length=maxdays); I.sum[]<-0
  
  # Construct probability distributions
  prob.t <- dpois(1:k, Tg); prob.t <- prob.t/sum(prob.t)
  
  # Compute i.begin (for y[i])
  i.begin <- d-k
  
  # Compute i.end (for y[i])
  i.end <- d-1
  if (i.end>maxdays) i.end<-maxdays
  
  # Compute expected secondary case time-series retrospectively from day d to day d-k
  # that is, i.begin does not always start from 1, but i.end is always d-1  
  if ((d - k) < 1) i.begin <- 1 else
    i.begin <- d - k
  
  for (i in i.begin: i.end) {
    
    if (S>0) yj.d[i] <- y[i] * R0 * prob.t[d-i]
    
  }
  
  exp.y <- sum(yj.d)
  
  return(exp.y)
}

kp1.calc.exp.inf.tmp <- function(d, y, N, R0, Tg, k) {
  maxdays <- length(y)
  if ((d > maxdays) | (d < 2)) stop("problem in kp1.sim.one.day")
  if (Tg <= 1) stop("value of Tg too small in kp1.sim.one.day")
  exppic <- 0
  S_d <- N - sum(y[1:(d-1)]) #number of susceptible at the start of day d
  if (S_d > 0) {
    if (d - k < 2) {
      currentinf <- 1 #i in eqn
    } else {
      currentinf <- d - k
    }
    while (currentinf < d) {
      #exppic <- exppic + y[currentinf]*R0*dpois(d-currentinf-1,Tg-1)
      #exppic <- exppic + y[currentinf]*R0*dpois(d-currentinf-1,Tg-0.5)
      exppic <- exppic + y[currentinf]*R0*dpois(d-currentinf,Tg)
      currentinf <- currentinf + 1
      #print(exppic)
    }
  }
  exppic
}

# plot(kp1.sim.all.days(N=1000,R0=1.8,Tg=3,maxdays=75,seed=10,k=12),type="l")
kp1.sim.all.days <- function(seed, maxdays, N, R0, Tg, k) {
  y <- vector(mode="numeric",length=maxdays)
  y[1] <- seed
  for (d in 2:(maxdays-1)) {
    y <- kp1.sim.one.day(d,y,N,R0,Tg,k)
  }
  y
}

# Inference function
# kp1.log.like.one.day(y=c(1,1,0,0),d=2,N=1000,R0=1.8,Tg=2,k=6)
kp1.log.like.one.day <- function(y,d,N,R0,Tg,k) {
  rtn <- 0
  S_d <- N - sum(y[1:(d-1)])
  if (S_d > 0) {
    exppic <- kp1.calc.exp.inf(d,y,N,R0,Tg,k)
    rtn <- rtn + dbinom(y[d],S_d,exppic/N,log=TRUE)
  } 
  rtn
}


kp1.log.like.upto.day <- function(y,d,...) {
  rtn <- 0
  if ((d < 2) | (d > length(y))) stop("d out of range in kp1.log.like.upto.day")
  for (i in 2:d) {
    rtn <- rtn + kp1.log.like.one.day(y,i,...)
  }
  rtn
}


# debugging line above here
kp1.like <- function(vecPar,y,d,N,k,R0,Tg) {
  R0tmp <- R0
  Tgtmp <- Tg
  if (R0 < 0 & Tg > 0 & length(vecPar)==1) {
    R0tmp <- vecPar[1]
  } else if (R0 > 0 & Tg < 0 & length(vecPar)==1) {
    Tgtmp <- vecPar[1]    	
  } else if (R0 < 0 & Tg < 0 & length(vecPar)==2) {
    R0tmp <- vecPar[1]
    Tgtmp <- vecPar[2]
  } else stop("problem here kp1.like")
  lnlike <- kp1.log.like.upto.day(y,d,R0=R0tmp,Tg=Tgtmp,N=N,k=k)
  lnlike
}


## ------------------------------------------------------------------------------------------- ##
