#Reading data
#------------



#To input the leukaemia data:

 leukaemia<-read.table("http://www.mas.ncl.ac.uk/~nrh3/dagstat/leukaemia.dat",header=T)
#

#To input the kidney data:

 kidney<-read.table("http://www.mas.ncl.ac.uk/~nrh3/dagstat/kidney.dat",header=T)


#To input the lung cancer data:

 lung<-read.table("http://www.mas.ncl.ac.uk/~nrh3/dagstat/lung.dat",header=T)


#To input the retinopathy data:

 retinopathy<-read.table("http://www.mas.ncl.ac.uk/~nrh3/dagstat/retinopathy.dat",header=T)



#To input the ships data:

 ships<-read.table("http://www.mas.ncl.ac.uk/~nrh3/dagstat/ships.dat",header=T)

#Alternative if data already loaded
#leukaemia<-read.table("leukaemia.dat",header=T)
#kidney<-read.table("kidney.dat",header=T)
#lung<-read.table("lung.dat",header=T)
#retinopathy<-read.table("retinopathy.dat",header=T)
#ships<-read.table("ships.dat",header=T)



#Routines for practical sessions follow....


####################################################################
stateprobs=function(fit1,fit2,newcase,addplot=FALSE,lty=1){
n=dim(newcase)[1]
b1=basehaz(fit1,centered=FALSE)
b2=basehaz(fit2,centered=FALSE)
time=c(0,b1$time)
plot(time,time*0,pch=" ",ylim=c(0,1),xlab="Time t",ylab="P(t)")
for(it in 1:n){
A1=(b1$hazard)*exp(sum(fit1$coefficients*newcase[it,]))
A2=(b2$hazard)*exp(sum(fit2$coefficients*newcase[it,]))
a1=A1-c(0,A1[-length(A1)])
a2=A2-c(0,A2[-length(A2)]) 
a0=1-a1-a2
a0=c(1,a0)
a1=c(0,a1)
a2=c(0,a2)
time=c(0,b1$time)
p01=cumsum(cumprod(a0)*a1)
p02=cumsum(cumprod(a0)*a2)
p00=1-p01-p02
lines(time,p00,type="s",col=2,lty=it)
lines(time,p01,type="s",col=3,lty=it)
lines(time,p02,type="s",col=4,lty=it)
}
legend(0.8*max(time),0.9,c("OK","Failed","Dead"),lty=c(1,1,1),col=c(2,3,4),bty="n")
if(n>1) {legend(0.2*max(time),0.9,as.character(1:n),lty=1:n,bty="n")}
}


#################################################################
simAhat=function(n=100,shape=1,rancens=TRUE,cenrate=1)
{ 
  #Weibull survival with cumulative hazard t^shape
  #Exponential censoring (if any) with rate cenrate
  
  t=(-log(runif(n)))^(1/shape)
  ic=rep(1,n)
  if(rancens==TRUE){
    c1=(-log(runif(n))/cenrate)
    i=c1<t
    t[i]=c1[i]
    ic[i]=0
  }
  i=order(t)
  time=t[i]
  cens=ic[i]
  cat(100*sum(1-cens)/n, "percent censored\n")
  Y=n:1
  dN=cens
  alpha=dN/Y
  Ahat=cumsum(alpha)
  VarAhat=cumsum(dN/Y^2)
  sdA=sqrt(VarAhat)
  upper=Ahat+2*sdA
  lower=Ahat-2*sdA
  tm=seq(0,max(time),length=100)
  A=tm^shape
  yl=range(c(0,lower,upper,A))
  plot(time,Ahat,col=4,type="s",xlab="Time",ylab="Events",ylim=yl)
  lines(tm,A,col=2)
  lines(time,lower,col=4,lty=2,type="s")
  lines(time,upper,col=4,lty=2,type="s") 
}

#######################################################################

simfrail=function(n=1000,shape=1,b=1,usefrail=TRUE,frailvar=0.5,cenrate=1,rancens=TRUE)
{ 
  #Weibull survival with cumulative baseline hazard t^shape, a single N(0,1) covariate and gamma frailty
  #Exponential censoring (if any) with rate cenrate
  z=rgamma(n,1/frailvar,1/frailvar)
  if(usefrail==FALSE) z=rep(1,length=n)
  x=rnorm(n)
  ebx=exp(b*x)
  zebx=z*ebx
  t=(-log(runif(n))/zebx)^(1/shape)
  ic=rep(1,n)
  if(rancens==TRUE){
    c1=(-log(runif(n))/cenrate)
    i=c1<t
    t[i]=c1[i]
    ic[i]=0
  }
  i=order(t)
  time=t[i]
  cens=ic[i]
  x=x[i]
  cat(100*sum(1-cens)/n, "percent censored\n")
  data.frame(time=time,cens=cens, x=x, id=1:n)
}
###################################################################
plaallung=function(cols=3:8, maxtime=18, mpl = c(3, 3)){
time=lung$time
cens=lung$cens
p=length(cols)
X0=as.matrix(lung[,cols])
X0=cbind(1,X0)
n=length(time)
ii=order(time)
time=time[ii]
cens=cens[ii]
X0=X0[ii,]
tms=unique(time)
m=length(tms)
events=risk=matrix(0,n,m)
for(it in 1:m){
events[,it][(time==tms[it])&(cens==1)]=1
risk[,it][time>=tms[it]]=1
}
ii=tms<=maxtime
tms=tms[ii]
events=events[,ii]
risk=risk[,ii]
m=dim(events)[2]
acum=vcum=NULL
for(it in 1:m){
 ii=(1:n)[risk[,it]==1]
 dN=events[ii,it]
 X=as.matrix(X0[ii,])
 X1=solve(t(X)%*%X)%*%t(X)
 ah=X1%*%matrix(dN,ncol=1)
 pp=as.vector(X%*%ah)
 vh=(X1*X1)%*%(pp*(1-pp)) 
 vh=as.vector(vh)
 ij=vh<0
 vh[ij]=0
 acum=rbind(acum,as.vector(ah))
 vcum=rbind(vcum,vh)
 }

acum=apply(acum,2,cumsum)
vcum=apply(vcum,2,cumsum)
upper=acum+2*sqrt(vcum)
lower=acum-2*sqrt(vcum)

lbs=as.vector(dimnames(X0)[[2]])
lbs[1]="baseline"
 par(mfrow = mpl)
  for(it in 1:(p+1)) {
    yl=range(c(lower[,it],upper[,it]))
    plot(tms, acum[,it], type = "s", xlab = 
         "Time", ylab = lbs[it], ylim=yl)
    lines(tms,upper[,it],type="s",lty=2)
    lines(tms,lower[,it],type="s",lty=2)
    abline(h=0,lty=3)
  }

  par(mfrow=c(1,1))
}
###################################################################
plaalships=function( maxtime=300, mpl = c(2, 3),useowner=FALSE){
time1=ships$time1
time2=ships$time2
cens=ships$cens
type=ships$type
dwt=ships$dwt
speed=ships$speed
container=tanker=speed*0
container[type==2]=1
tanker[type==3]=1
owner=ships$owner-1
p=4
X0=cbind(1,container,tanker,dwt,speed)
if(useowner==TRUE){
p=5
X0=cbind(1,container,tanker,dwt,speed,owner)}

n=length(time2)
tms=sort(unique(time2))
m=length(tms)
events=risk=matrix(0,n,m)
for(it in 1:m){
events[,it][(time2==tms[it])&(cens==1)]=1
risk[,it][(time1<tms[it])&(time2>=tms[it])]=1
}
ii=tms<=maxtime
tms=tms[ii]
events=events[,ii]
risk=risk[,ii]
m=dim(events)[2]
acum=vcum=NULL
for(it in 1:m){
 ii=(1:n)[risk[,it]==1]
 dN=events[ii,it]
 X=as.matrix(X0[ii,])
 X1=solve(t(X)%*%X)%*%t(X)
 ah=X1%*%matrix(dN,ncol=1)
 pp=as.vector(X%*%ah)
 vh=(X1*X1)%*%(pp*(1-pp)) 
 vh=as.vector(vh)
 ij=vh<0
 vh[ij]=0
 acum=rbind(acum,as.vector(ah))
 vcum=rbind(vcum,vh)
 }

acum=apply(acum,2,cumsum)
vcum=apply(vcum,2,cumsum)
upper=acum+2*sqrt(vcum)
lower=acum-2*sqrt(vcum)

lbs=c("baseline","container","tanker","dwt","speed","owner")
 par(mfrow = mpl)
  for(it in 1:(p+1)) {
    yl=range(c(lower[,it],upper[,it]))
    plot(tms, acum[,it], type = "s", xlab = 
         "Time", ylab = lbs[it], ylim=yl)
    lines(tms,upper[,it],type="s",lty=2)
    lines(tms,lower[,it],type="s",lty=2)
    abline(h=0,lty=3)
  }

  par(mfrow=c(1,1))
}

###############################################################
simcrdat=function(n=1000,b=c(1,1),alpha00=1,alpha01=1,shape=2,censrate=1,xi=0.5,nu=0.5,usegamma=FALSE, shared=TRUE,userandom=FALSE, tmax=5){
#Shared frailty model, common margins, exponential censoring

#Default is positive stable z with parameter nu (between 0 and 1). Taking nu=1 means no frailty.
#Under the positive stable distribution and independent censoring, the marginals are prop hazards with regression parameters beta*nu.  So to compensate we redefine new beta as old beta divided by nu.  That way we expect the estimates to be close to the betas we put in.

#We allow an option for time1 and time2 to have shared z (so dependent censoring) or separate z (so we can investigate the marginals).

  

if(usegamma==FALSE){

if(nu==1){z1=z2=rep(1,length=n)
}else{
b=b/nu
z1=getposstable(n=n,nu=nu)
z2=getposstable(n=n,nu=nu)
}
}

#Alternative in case of need is gamma.  Note that even under independent censoring the marginals are not PH
if(usegamma==TRUE){
if(xi==0){z=rep(1,length=n)}else{
z1=rgamma(n,shape=1/xi,scale=xi)
z2=rgamma(n,shape=1/xi,scale=xi)
}
}


if(shared==TRUE){
z2=z1
}


x1=rnorm(n)
x2=x1*0
x2[runif(n)>0.5]=1
bx=b[1]*x1+b[2]*x2
ebx=exp(bx)
h1=alpha00*z1*ebx
h2=alpha01*z2*ebx
t1=(-log(runif(n))/h1)^(1/shape)
t2=(-log(runif(n))/h2)^(1/shape)
tcens=(-log(runif(n))/censrate)
if(userandom==FALSE) tcens=rep(tmax,length=n) #Set tmax=Inf if needed
tcens[tcens>tmax]=tmax
time=tcens


delta=rep(0,length=n)
i1=(t1<t2)&(t1<tcens)
delta[i1]=1
time[i1]=t1[i1]
i2=(t2<t1)&(t2<tcens)
delta[i2]=2
time[i2]=t2[i2]
cens1=cens2=delta*0
cens1[delta==1]=1
cens2[delta==2]=1
#We will keep the latent event times (randomly censored) as a bivariate survival example
cens11=cens21=rep(1,length=n)
i11=t1>tcens
t1[i11]=tcens[i11]
cens11[i11]=0
i21=t2>tcens
t2[i21]=tcens[i21]
cens21[i21]=0



#The survfit type="mstate" option requires some censoring. Hence the last survival time is always censored here. 
ii=time==max(time)
delta[ii]=cens1[ii]=cens2[ii]=0
#Put the data in time order
jj=order(time)
time=time[jj]
delta=delta[jj]
cens1=cens1[jj]
cens2=cens2[jj]
x1=x1[jj]
x2=x2[jj]
#Keeping the underlying times
time1=t1[jj]
time2=t2[jj]
cens11=cens11[jj]
cens21=cens21[jj]

list(time=time, cens=delta,cens1=cens1,cens2=cens2, x1=x1,x2=x2,time1=time1, time2=time2,cens11=cens11,cens21=cens21)
}

##################################################################################
#This is a wrapper to simulate positive stable frailties
getposstable=function(n=500,nu=0.5){
z=rlaptrans(n=n,ltpdf=ltposstable,nu=nu)
z
}


##################################################################################
ltposstable=function(s,nu){
#nu should be between 0 and 1
#s any positive number
exp(-s^nu)
}
############################################################################################
#Martin Ridout code, as described in Ridout, M.S. (2009) Generating random numbers from a distribution specified by its Laplace transform. Statistics and Computing, 19, 439-450. 


#======================================================================================
rlaptrans <- function(n, ltpdf, ..., tol=1e-7, x0=1, xinc=2, m=11, L=1, A=19, nburn=38)
#======================================================================================
{

    #---------------------------------------------------------
    # Function for generating a random sample of size n from a 
    # distribution, given the Laplace transform of its p.d.f.
    #---------------------------------------------------------

    maxiter = 500

      # -----------------------------------------------------
      # Derived quantities that need only be calculated once,
      # including the binomial coefficients
      # -----------------------------------------------------
    nterms = nburn + m*L
    seqbtL = seq(nburn,nterms,L)
    y = pi * (1i) * seq(1:nterms) / L
    expy = exp(y)
    A2L = 0.5 * A / L
    expxt = exp(A2L) / L
    coef = choose(m,c(0:m)) / 2^m


      # --------------------------------------------------
      # Generate sorted uniform random numbers. xrand will
      # store the corresponding x values
      # --------------------------------------------------
    u = sort(runif(n), method="qu")
    xrand = u

      #------------------------------------------------------------
      # Begin by finding an x-value that can act as an upper bound
      # throughout. This will be stored in upplim. Its value is
      # based on the maximum value in u. We also use the first
      # value calculated (along with its pdf and cdf) as a starting
      # value for finding the solution to F(x) = u_min. (This is
      # used only once, so doesn't need to be a good starting value
      #------------------------------------------------------------
    t = x0/xinc
    cdf = 0   
    kount0 = 0
    set1st = FALSE
    while (kount0 < maxiter & cdf < u[n]) {
        t = xinc * t
        kount0 = kount0 + 1
        x = A2L / t
        z = x + y/t
        ltx = ltpdf(x, ...)
        ltzexpy = ltpdf(z, ...) * expy
        par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
        par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
        pdf = expxt * sum(coef * par.sum[seqbtL]) / t
        cdf = expxt * sum(coef * par.sum2[seqbtL]) / t
        if (!set1st & cdf > u[1]) {
            cdf1 = cdf
            pdf1 = pdf
            t1 = t
            set1st = TRUE
        }
    }
    if (kount0 >= maxiter) {
       stop('Cannot locate upper quantile')
    }
    upplim = t

      #--------------------------------
      # Now use modified Newton-Raphson
      #--------------------------------

    lower = 0
    t = t1
    cdf = cdf1
    pdf = pdf1
    kount = numeric(n)

    maxiter = 1000

    for (j in 1:n) {

          #-------------------------------
          # Initial bracketing of solution
          #-------------------------------
        upper = upplim

        kount[j] = 0
        while (kount[j] < maxiter & abs(u[j]-cdf) > tol) {
            kount[j] = kount[j] + 1

              #-----------------------------------------------
              # Update t. Try Newton-Raphson approach. If this 
              # goes outside the bounds, use midpoint instead
              #-----------------------------------------------
            t = t - (cdf-u[j])/pdf 
            if (t < lower | t > upper) {
               t = 0.5 * (lower + upper)
            }

              #----------------------------------------------------
              # Calculate the cdf and pdf at the updated value of t
              #----------------------------------------------------
            x = A2L / t
            z = x + y/t
            ltx = ltpdf(x, ...)
            ltzexpy = ltpdf(z, ...) * expy
            par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
            par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
            pdf = expxt * sum(coef * par.sum[seqbtL]) / t
            cdf = expxt * sum(coef * par.sum2[seqbtL]) / t

              #------------------
              # Update the bounds 
              #------------------
            if (cdf <= u[j]) {
                lower = t}
              else {
                upper = t}
        }
        if (kount[j] >= maxiter) {
           warning('Desired accuracy not achieved for F(x)=u')
        }
        xrand[j] = t
        lower = t
    }

    if (n > 1) {
       rsample <- sample(xrand) }
     else {
       rsample <- xrand} 
    rsample
}






