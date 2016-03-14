# own code

library(survival)
source("someR.R")

# exercise 1

time <- c(1, 3, 4, 5, 6, 8, 10, 13)
cens <- c(1, 0, 1, 1, 0, 1,  1,  1)

dN <- cens
N <- cumsum(dN)
Y <- length(time):1

alphahat <- 1/Y
alphahat[cens == 0] <- 0
Ahat <- cumsum(alphahat)

lambdahat <- Y*alphahat
VarAhat <- cumsum(1/Y^2)

ex1 <- data.frame(time=time, cens=cens, N=N, Y=Y, alphahat=alphahat, Ahat=Ahat, VarAhat=VarAhat)

plot(c(0,time), c(0,Ahat), type="s", xlab="t", ylab=expression(hat(A)(t)))
AhatCIl <- Ahat - qnorm(0.975)*VarAhat
AhatCIu <- Ahat + qnorm(0.975)*VarAhat
lines(time, AhatCIl, type="s", lty=2)
lines(time, AhatCIu, type="s", lty=2)

KMS <- cumprod(c(1, (1-dN/Y)))
plot(c(0,time), KMS, type="s", xlab="t")
lines(c(0,time), c(1, exp(-Ahat)), type="s", col=2)

# exercise 2
tem1 <- basehaz(coxph(Surv(time, cens) ~ 1, data=leukaemia))
Ahat1 <- tem1$hazard
time1 <- tem1$time

tem2 <- survfit(Surv(time, cens) ~ 1, data=leukaemia)
Ahat2 <- cumsum(tem2$n.event/tem2$n.risk)
time2 <- tem2$time
VarAhat2 <- cumsum(tem2$n.event/(tem2$n.risk)^2)

plot(time1, Ahat1, type="s")
plot(tem2$time, tem2$surv, type="s")
plot(tem2)

# ....

# exercise 3
head(ships)
hist(ships$speed)
hist(table(ships$id), breaks=(1:9)-0.5)

tem2 <- survfit(Surv(time1, time2, cens) ~ 1, data=ships)
Ahat2 <- cumsum(tem2$n.event[-1]/tem2$n.risk[-1])
time2 <- tem2$time[-1]
plot(time2, Ahat2, type="s")
plot(time2, exp(-Ahat2), type="s")
