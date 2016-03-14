# own code

library(survival)
source("someR.R")

# Session 1
# ---------
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

# ....

# Session 3
# ---------

# exercise 1
plaallung()
# baseline approx 0 -> baseline almost zero hazard (for age=0, male, activity 0 etc.)
# age approx 0 -> no effect
# female constantly negative (tapes off slightly) -> hazard lower for women
# activity constant increase (borderline sig) -> hazard higher when less activity
# anorexia increasing until ca 8, then flat -> increases hazard before 8, later not
#, hoarseness, metastases approx 0

# exercise 2
tau <- 8
fit0 <- coxph(Surv(time, cens) ~ age + female + activity + anorexia + hoarseness + metastases, data = lung)
lung1 <- lung
lung2 <- lung
i1 <- lung$time > tau
i2 <- lung$time <= tau
lung1$cens[i1] <- 0
lung2$cens[i2] <- 0
fit1 <- coxph(Surv(time, cens) ~ age + female + activity + anorexia + hoarseness + metastases, data = lung1)
fit2 <- coxph(Surv(time, cens) ~ age + female + activity + anorexia + hoarseness + metastases, data = lung2)
fit0$loglik
fit1$loglik
fit2$loglik
fit1$loglik + fit2$loglik
coef(fit1)
coef(fit2)

# exercise 3
plaalships()
# baseline almost linear -> constant positive hazard
# type: bulker (ref category), container, tanker
# container decreasing until ca 150, then constant until ca 220, then increasing or zero
# -> lower hazard until 150, then no effect or increased hazard
# tanker (0 until ca 80, then) decreasing to ca 150, then increasing
# -> lower hazard until 150, then increased hazard 
# -> similar to container (merge container & tanker?)
# dwt decreasing until ca 150, then increasing
# speed slightly decreasing -> bordeline negative effect (lower hazard for faster ships)

# exercise 4
tau <- 150
fit0 <- coxph(Surv(time1, time2, cens) ~ as.factor(type) + dwt + speed, data = ships)
ships1 <- ships
ships2 <- ships
i1 <- ships$time2 > tau
i2 <- ships$time2 <= tau
ships1$cens[i1] <- 0
ships2$cens[i2] <- 0
fit1 <- coxph(Surv(time1, time2, cens) ~ as.factor(type) + dwt + speed, data = ships1)
fit2 <- coxph(Surv(time1, time2, cens) ~ as.factor(type) + dwt + speed, data = ships2)
fit0$loglik
fit1$loglik
fit2$loglik
fit1$loglik + fit2$loglik
coef(fit1)
coef(fit2)
fit1
fit2

# exercise 5
plaalships(useowner=TRUE)
# ....


# Session 4
# ---------

# exercise 1
table(kidney$txnum)
cif0 <- survfit(Surv(time, cens, type="mstate") ~ 1, data=kidney)
plot(cif0$time, cif0$prev[,1], type="s")
lines(cif0$time, cif0$prev[,2], type="s", col=2)
legend("bottomright", legend=c("failure", "death"), col=c(1,2), lty=c(1,1))

cif1 <- survfit(Surv(time, cens, type="mstate") ~ 1, data=subset(kidney, txnum==1))
cif2 <- survfit(Surv(time, cens, type="mstate") ~ 1, data=subset(kidney, txnum==2))
cif3 <- survfit(Surv(time, cens, type="mstate") ~ 1, data=subset(kidney, txnum>=3))

par(mfrow=c(1,3))
plot(cif1$time, cif1$prev[,1], type="s", ylim=c(0,1))
lines(cif1$time, cif1$prev[,2], type="s", col=2)
legend("topleft", legend=c("failure", "death"), col=c(1,2), lty=c(1,1))

plot(cif2$time, cif2$prev[,1], type="s", ylim=c(0,1))
lines(cif2$time, cif2$prev[,2], type="s", col=2)
legend("topleft", legend=c("failure", "death"), col=c(1,2), lty=c(1,1))

plot(cif3$time, cif3$prev[,1], type="s", ylim=c(0,1))
lines(cif3$time, cif3$prev[,2], type="s", col=2)
legend("topleft", legend=c("failure", "death"), col=c(1,2), lty=c(1,1))
par(mfrow=c(1,1))

# exercise 2

coxfail <- coxph(Surv(time, censfail) ~ rage + capd + drmm + bmm, data=kidney)
coxfail
coxdeath <- coxph(Surv(time, censdeath) ~ rage + capd + drmm + bmm, data=kidney)
coxdeath
#