# VER Ch.19: Modelling Survival Data

# either download data from the web:
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/calf_pneu.rdata"))
#unlink(temp)  # remove the temporary file

# or use the unzipped archive
load("../../data/ver2_data_R/calf_pneu.rdata")

dim(calf_pneu)
str(calf_pneu)
library(Hmisc)
calf_pneu <- upData(calf_pneu, labels = c(calf = 'Calf id',
                                  stock = 'Stocking method',
                             days = 'Time to onset of pneumonia or censoring',
                                  pn = 'Pneumonia'),
                  levels = list(stock = list('batch' = 0, 'continuous' = 1)))


# Ex. 19.2 Actuarial life table
library(KMsurv)
interval <- seq(from = 30, to = 165, by = 15)
interval <- floor(calf_pneu$days / 15)
interval.censor <- data.frame(interval, calf_pneu$pn)
library(nlme)
pneumonia <- gsummary(interval.censor, sum, groups = interval)
total <- gsummary(interval.censor, length, groups = interval)
rm(interval.censor)
lt.data <- cbind(pneumonia[, 1:2], total[, 2])
length <- length(lt.data$interval)
lt.data[length + 1, ]$interval <- NA
nevent <- lt.data[, 2]
nlost <- lt.data[, 3] - lt.data[, 2]
(life.table <- lifetab(lt.data$interval, 24, nlost, nevent))


# Ex. 19.2 Kaplan-Meier survivor function
library(survival)
km.sf <- survfit(Surv(days, pn == 1) ~ 1, data = calf_pneu)
summary(km.sf)
plot(km.sf, xlab = "time (days)", ylab = "cumulative survival probability",
     conf.int = TRUE)
# Nelson-Aalen estimate of the cumulative hazard obtained by transforming the
# Fleming-Harrington estimate of survival
fh.sf <- survfit(Surv(days, pn == 1) ~ 1, data = calf_pneu, type = "fleming")
plot(fh.sf$time, -log(fh.sf$surv), xlab = "time (days)",
     ylab = "cumulative hazard",
     main = "Nelson-Aalen cumulative hazard estimate", ylim = c(0, 1.5),
     type = "l")
lines(fh.sf$time, -log(fh.sf$upper), lty = 3)
lines(fh.sf$time, -log(fh.sf$lower), lty = 3)
# It's more accurate to draw the estimated cumulative hazard as a step function.
# R provides a function to conveniently draw step functions:
plot(stepfun(fh.sf$time, c(0, -log(fh.sf$surv))), do.points = FALSE, 
     xlab = "time (days)", ylab = "cumulative hazard",
     main = "", ylim = c(0, 1.5))
lines(stepfun(fh.sf$time, c(0, -log(fh.sf$upper))), lty = 5, do.points = FALSE)
lines(stepfun(fh.sf$time, c(0, -log(fh.sf$lower))), lty = 5, do.points = FALSE)
# N-A estimate can also be computed directly from the contents of a survfit
# object.
nelson.aalen <- cumsum(fh.sf$n.event / fh.sf$n.risk)
t <- fh.sf$time
points(t, nelson.aalen, pch = "+")


# Ex. 19.3 Log-rank test
survdiff(Surv(days, pn == 1) ~ stock, data = calf_pneu, rho = 0)
# rho is optional. Default is rho = 0 which gives the log-rank test.
# Rho = 1 gives the "Peto & Peto modiﬁcation of the Gehan-Wilcoxon test".
# Rho larger than zero gives greater weight to the ﬁrst part of the survival curves.
# Rho smaller than zero gives weight to the later part of the survival curves.
survdiff(Surv(days, pn == 1) ~ stock, data = calf_pneu, rho = 1)


# Ex. 19.4 Comparing survivor functions
(km.stock <- survfit(Surv(days, pn == 1) ~ stock, data = calf_pneu))
plot(km.stock, conf.int = FALSE, col = c("blue4", "darkorange"),
     xlab = "time (days)", ylab = "cumulative survival probability")
legend("bottomleft", inset = .05, c("batch", "continuous"),
       text.col = c("blue4", "darkorange"))

km.df <- data.frame(
   time    = km.stock$time,
   n.risk  = km.stock$n.risk,
   n.event = km.stock$n.event,
   surv    = km.stock$surv,
   strata  = gsub("stock=", "", summary(km.stock, censored = T)$strata),
   upper   = km.stock$upper,
   lower   = km.stock$lower
) 
zeros <- data.frame(time = 0, surv = 1, strata = gsub("stock=", "",
                                          levels(summary(km.stock)$strata)), 
                    upper = 1, lower = 1)
library(plyr)
km.df <- rbind.fill(zeros, km.df)
km.df$strata <- ordered(km.df$strata, levels = c("batch", "continuous"))
library(ggplot2)
ggplot(km.df, aes(time, surv, colour = strata)) + 
  geom_step(size = 0.6) + xlim(0, 150) + ylim(0, 1) + 
  xlab("time (days)") + ylab("cumulative survival probability") +
  labs(colour = "stock")


# Ex. 19.5 Cox proportional hazards model
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/pgtrial.rdata"))
#unlink(temp)
# or use the unzipped archive
load("../../data/ver2_data_R/pgtrial.rdata")

dim(pgtrial)
str(pgtrial)
pgtrial <- upData(pgtrial, labels = c(herd = 'Herd id', cow = 'Cow id',
                             tx = 'Treatment', lact = 'Lactation number',
                             thin = 'Body condition', dar = 'Days at risk',
                             preg = 'Pregnant or censored'),
                  levels = list(thin = list('normal' = 0, 'thin' = 1),
                    preg = list('censored' = 0, 'pregnant' = 1)))
pgtrial$herd <- as.factor(pgtrial$herd)

coxph.mod <- coxph(Surv(dar, preg == 'pregnant') ~ herd + tx + lact + thin,
                   data = pgtrial, ties = 'breslow')
(coxph.sum <- summary(coxph.mod))
# R gives several options to control ties in case several events occurred at the
# same time: the Efron method (default in R), Breslow method (default in software
# like SAS or Stata), and the exact method. Breslow is the simplest and adequate
# if not too many ties in the dataset. Efron is closer to the exact approximation.


# Ex. 19.7 Stratified Cox proportional hazards model
scoxph.mod <- coxph(Surv(dar, preg == 'pregnant') ~ tx + tx*herd + lact +
                    thin + strata(herd), data = pgtrial, method = 'breslow')
summary(scoxph.mod)


# Evaluating assumption of proportional hazards
# fig. 19.15 Log cumulative hazard plot
coxph.mod2 <- coxph(Surv(dar, preg == 'pregnant') ~ tx, data = pgtrial,
                   ties = 'breslow')
pgtrial2 <- with(pgtrial, data.frame(tx = c(0, 1)))
tfit.add <- survfit(coxph.mod2, newdata = pgtrial2)
df1 <- data.frame(
    time    = tfit.add[1, ]$time,
    n.risk  = tfit.add[1, ]$n.risk,
    n.event = tfit.add[1, ]$n.event,
    surv    = tfit.add[1, ]$surv,
    strata  = "0",
    upper   = tfit.add[1, ]$upper,
    lower   = tfit.add[1, ]$lower,
    log.surv = log(-log(tfit.add[1, ]$surv))
 ) 
df2 <- data.frame(
    time    = tfit.add[2, ]$time,
    n.risk  = tfit.add[2, ]$n.risk,
    n.event = tfit.add[2, ]$n.event,
    surv    = tfit.add[2, ]$surv,
    strata  = "1",
    upper   = tfit.add[2, ]$upper,
    lower   = tfit.add[2, ]$lower,
    log.surv = log(-log(tfit.add[2, ]$surv))
 )
dfpar.add <- rbind(df1, df2)
zeros <- data.frame(time = 0, surv = 1, strata = c(1, 2), 
                     upper = 1, lower = 1)
dfpar.add <- rbind.fill(zeros, dfpar.add)
dfpar.add$strata <- factor(dfpar.add$strata, labels = c("No tx", "Tx"))
ggplot(dfpar.add, aes(log(time), log.surv, colour = strata)) + 
  geom_step(size = 0.6) + 
  scale_color_manual("Tx", values = c('blue4', 'darkorange')) + 
  xlab("ln(time)") + ylab("Log-log survival")

# fig. 19.16 Kaplan-Meier Cox plot
tfit.km <- survfit(Surv(dar, preg == 'pregnant') ~ tx, data = pgtrial)
df3.km <- data.frame(
    time    = tfit.km$time,
    n.risk  = tfit.km$n.risk,
    n.event = tfit.km$n.event,
    surv    = tfit.km$surv,
    strata  = gsub("tx=", "", summary(tfit.km, censored = T)$strata),
    upper   = tfit.km$upper,
    lower   = tfit.km$lower
 ) 
zeros <- data.frame(time = 0, surv = 1, strata = gsub("tx=", "", 
                                          levels(summary(tfit.km)$strata)), 
                     upper = 1, lower = 1)
df3.km <- rbind.fill(df3.km, zeros)
df3.km$cat <- with(df3.km, ifelse(strata == "0", "No tx, observed",
                                  "Tx, observed"))
dfpar.add$cat <- with(dfpar.add, ifelse(strata == "No tx", "No tx, expected",
                                        "Tx, expected"))
dfpar.obs <- rbind.fill(dfpar.add, df3.km)
ggplot(dfpar.obs, aes(time, surv, colour = cat)) + 
   geom_step(size = 0.6) + 
   scale_color_manual("", values = c('blue1', 'blue4', 'darkorange1',
                            'darkorange4')) + 
   xlab("time") + ylab("survival probability")

# Ex. 19.11 Assessing the proportional hazards assumption - Schoenfeld residuals
(schoen <- cox.zph(coxph.mod))
plot(schoen, var = 4)

#dfbeta <- residuals(coxph.mod, type = 'dfbeta')
#par(mfrow = c(2, 2))
#for (j in 1:3) {
#  plot(dfbeta[ , j], ylab = names(coef(coxph.mod))[j])
#  abline(h = 0, lty = 2)
#}


# Ex. 19.13 Evaluating overall fit of a model
# The default residuals of coxph in R are the martingale residuals, not the
# Cox-Snell, but it can be computed:
cox.snell <- (as.numeric(pgtrial$preg) - 1) -
    resid(coxph.mod, type = "martingale")
# Then using Nelson-Aalen method to estimate the cumulative hazard function
# for residuals
coxph.res <- survfit(coxph(Surv(cox.snell, pgtrial$preg == 'pregnant') ~ 1,
                           method = 'breslow'), type = 'aalen')
plot(coxph.res$time, -log(coxph.res$surv), type = 's',
     xlab = 'Modified Cox-Snell residuals', ylab = 'Cumulative hazard')
abline(0, 1, col = 'red', lty = 2)
# alternatively:
coxph.res2 <- survfit(Surv(cox.snell, pgtrial$preg == 'pregnant') ~ 1)
Htilde <- cumsum(coxph.res2$n.event / coxph.res$n.risk)
plot(coxph.res2$time, Htilde, type = 's', col = 'blue')
abline(0, 1, col = 'red', lty = 2)

# GOF (Gronnesby and Borgan omnibus gof)
library(gof)
cumres(coxph.mod)

# Concordance
# Harrell's c index computes the proportion of all pairs of subjects in which
# the model correctly predicts the sequence of events. It ranges from 0 to 1
# with 0.5 for random predictions and 1 for a perfectly discriminating model.
# It is obtained from the Somer’s Dxy rank correlation:
library(rms)
fit.cph <- cph(Surv(dar, preg == 'pregnant') ~ herd + tx + lact + thin,
               data = pgtrial, x = TRUE, y = TRUE, surv = TRUE)
# Get the Dxy
v <- validate(fit.cph, dxy = TRUE, B = 100)
Dxy <- v[rownames(v) == "Dxy", colnames(v) == "index.corrected"]
# c-statistic according to the Dxy = 2*(c - 0.5)
(Dxy / 2) + 0.5


# Ex. 19.14 Evaluating functional form of predictors
lact.mod <- coxph(Surv(dar, preg == 'pregnant') ~ lact, data = pgtrial,
                  ties = 'breslow')
lact.res <- resid(lact.mod, type = "martingale")
# fig. 19.19 Plot of martingale residuals vs lactation number
plot(pgtrial$lact, lact.res, xlab = 'lactation', ylab = 'Martingale residuals')
lines(lowess(pgtrial$lact, lact.res, iter = 0))
# adding quadratic term
lact.mod <- update(lact.mod, . ~ . + I(lact^2))
lact.res <- resid(lact.mod, type = "martingale")
plot(pgtrial$lact, lact.res, xlab = 'lactation', ylab = 'Martingale residuals')
lines(lowess(pgtrial$lact, lact.res, iter = 0))


# Checking for outliers
# deviance residuals
dev.res <- resid(coxph.mod, type = "deviance")
# fig. 19.20 Deviance residuals
plot(pgtrial$dar, dev.res, xlab = 'time (days)', ylab = 'deviance residuals')
# list of extreme residuals
cbind(dev.res, pgtrial)[abs(dev.res) > 2, ]

# Detecting influential points
# score residuals
score.res <- resid(coxph.mod, type = "score")
# score residuals for tx
# fig. 19.21 Score residuals
plot(pgtrial$dar, score.res[ , 3], xlab = 'time (days)',
     ylab = 'score residuals')
text(pgtrial$dar, score.res[ , 3], rownames(pgtrial), cex = 0.6, pos = 4)
cbind(score.res[ , 3], pgtrial)[abs(score.res[ , 3]) > 2, ]

# influential observations
dfbeta <- resid(coxph.mod, type = "dfbeta")
# dfbeta residuals for tx
plot(pgtrial$dar, dfbeta[ , 3], xlab = 'time (days)',
     ylab = 'scaled score residual')
text(pgtrial$dar, dfbeta[ , 3], rownames(pgtrial), cex = 0.6, pos = 4)
# with standardized dfbeta
dfbetas <- resid(coxph.mod, type = "dfbetas")
plot(pgtrial$dar, dfbetas[ , 3], xlab = 'time (days)',
     ylab = 'standardized score residuals')
text(pgtrial$dar, dfbetas[ , 3], rownames(pgtrial), cex = 0.6, pos = 4)


# Ex. 19.15 Exponential regression
exp.mod <- survreg(Surv(dar, preg == 'pregnant') ~ herd + tx + (lact - 1) +
                   thin, data = pgtrial, dist = "exponential")
summary(exp.mod)
# note that R outputs the parameter estimates of the AFT form of the exponential
# model. Multiply the estimated coefficients by minus one to get estimates
# consistent with the PH parameterization of the model. So for tx, the estimated
# hazard ratio is exp(0.2178) = 1.24 (at any given point in time, a treated cow is
# 1.24 times more likely to conceive that a non-treated cow). The corresponding
# accelerating factor for an exponential model is just the reciprocal of the hazard
# ratio, exp(-0.2178) = 0.80: treating a cow accelerates the time to conception
# by a factor of 0.80.


# Ex. 19.16 Weibull model
library(car)
pgtrial$parity <- recode(pgtrial$lact, "1 = 1; 2:hi = '2+'")
weib.mod <- survreg(Surv(dar, preg == 'pregnant') ~ herd + tx + parity +
                   thin, data = pgtrial, dist = "weibull")
summary(weib.mod)
# the shape parameter is the reciprocal of what is called by R the scale
# parameter. The shape parameter is then 1/1.15 = 0.869.


# Ex. 19.17 Piecewise constant exponential regression model
# piecewise exponential in time intervals of 40 units
# see blog post for explanation (http://denishaine.wordpress.com/2013/07/05/veterinary-epidemiologic-research-modelling-survival-data-parametric-and-frailty-models/)
interval.width <- 40 
# function to compute time breaks given the exit time = dar
cow.breaks <- function(dar) unique(c(seq(0, dar, by = interval.width),
                                      dar))
# list of each subject's time periods
the.breaks <- lapply(unique(pgtrial$cow), function(id){
  cow.breaks(max(pgtrial$dar[pgtrial$cow == id]))
})
# the expanded period of observation:
start <- lapply(the.breaks, function(x) x[-length(x)])  # for left time points 
stop <-  lapply(the.breaks, function(x) x[-1])  # for right time points
# needed to complete the long version of the pgtrial data set
count.per.cow <- sapply(start, length)
index <- tapply(pgtrial$cow, pgtrial$cow, length)
index <- cumsum(index)  # index of last observation for each cow

event <- rep(0, sum(count.per.cow))
event[cumsum(count.per.cow)] <- pgtrial$preg[index]

# bring all of this together to create the expanded dataset
pw.pgtrial <- data.frame(
    cow = rep(pgtrial$cow[index], count.per.cow),
    dar = rep(pgtrial$dar[index], count.per.cow),
    herd = rep(pgtrial$herd[index], count.per.cow),
    tx = rep(pgtrial$tx[index], count.per.cow),
    lact = rep(pgtrial$lact[index], count.per.cow),
    thin = rep(pgtrial$thin[index], count.per.cow),
    start = unlist(start),
    stop = unlist(stop),
    event = event
  )
# create time variable which indicates the period of observation (offset in
# the Poisson model)
pw.pgtrial$time <- pw.pgtrial$stop - pw.pgtrial$start  # length of observation
# create a factor for each interval, allowing to have a different rate for each
# period
pw.pgtrial$interval <- factor(pw.pgtrial$start)
pw.pgtrial[100:110, ]
# Poisson model
pw.model <- glm(event ~ offset(log(time)) + interval + herd + tx + lact +
                thin, data = pw.pgtrial, family = "poisson")
summary(pw.model)


# Ex. 19.18 Log-logistic model
loglog.mod <- survreg(Surv(dar, preg == 'pregnant') ~ herd + tx + lact +
                   thin, data = pgtrial, dist = "loglogistic")
summary(loglog.mod)


# Ex. 19.20 Individual frailty model
# Weibull model with gamma individual frailty
library(parfm)
pgtrial$preg.bin <- as.numeric(pgtrial$preg) - 1
indfr.mod <- parfm(Surv(dar, preg.bin) ~ herd + tx + lact + thin,
             cluster = "cow", data = pgtrial, dist = "weibull",
             frailty = "gamma")
indfr.mod


# Ex. 19.22 Shared frailty Cox model
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/culling.rdata"))
#unlink(temp)
load("../../data/ver2_data_R/culling.rdata")

dim(culling)
head(culling)
table(culling$culled)
summary(culling$culled)

library(frailtypack)
shfrw.mod <- frailtyPenal(Surv(dar, culled) ~ cluster(herd) +
                          as.factor(lact_c3) + johnes,
                          hazard = 'Weibull', data = culling)
summary(shfrw.mod)
shfrw.sum <- cbind(shfrw.mod$coef, sqrt(diag(shfrw.mod$varH)), 
                 shfrw.mod$coef / sqrt(diag(shfrw.mod$varH)), 
        signif(1 - pchisq((shfrw.mod$coef/sqrt(diag(shfrw.mod$varH)))^2, 1)),
                   exp(shfrw.mod$coef),
 exp(shfrw.mod$coef - abs(qnorm((1 - 0.95) / 2)) * sqrt(diag(shfrw.mod$varH))), 
 exp(shfrw.mod$coef + abs(qnorm((1 - 0.95) / 2)) * sqrt(diag(shfrw.mod$varH))))
row.names(shfrw.sum) <- c("Lactation 2", "Lactation 3+", "Johnes")
colnames(shfrw.sum) <- c("Coef.", "Std. Err.", "z", "p-value", "Hazard Ratio", 
                       "Lower CI", "Upper CI")
shfrw.sum


# Ex. 19.22 Shared frailty Cox model
shfrc.mod <- coxph(Surv(dar, culled) ~ as.factor(lact_c3) + johnes +
                   frailty(herd, distribution = "gamma"), data = culling)
summary(shfrc.mod)
