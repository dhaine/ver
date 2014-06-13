# VER Ch.13: Linear Regression

# either download data from the web:
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/daisy2.rdata"))
#unlink(temp)  # remove the temporary file

# or use the unzipped archive
load("../../data/ver2_data_R/daisy2.rdata")

dim(daisy2)
daisy2 <- daisy2[daisy2$h7 == 1, ]  # we only use a subset of the data
dim(daisy2)

summary(daisy2$milk120, na.rm = TRUE)
sd(daisy2$milk120, na.rm = TRUE)
daisy2 <- daisy2[!is.na(daisy2$milk120), ]  # get rid of missing observations for milk production

# Ex. 14.1 simple linear regression model
lm.milk <- lm(milk120 ~ parity, data = daisy2)
(lm.milk.sum <- summary(lm.milk))
anova(lm.milk)  # anova table
# table of fitted values and residuals
head(data.frame(daisy2[, c(7:8)], fitted.value = fitted(lm.milk),
     residual = resid(lm.milk)), n = 10)
# predicted values (reflect the uncertainty about the regression line)
head(predict(lm.milk, interval = "confidence"), n = 10)
### confidence bands (include also the uncertainty about future observations)
head(predict(lm.milk, interval = "prediction"), n = 10)
# plot actual values and intervals for prediction both for the mean and new observations. First create a new data frame to hold the values of parity for which we want the predictions, then compute the prediction and confidence bands and plot them.
confidence <- data.frame(parity = 1:7)
# confidence bands
confidence.band <- predict(lm.milk, int = "c" , newdata = confidence)
confidence.band <- as.data.frame(cbind(confidence.band, parity = c(1:7)))
# prediction bands
prediction.band <- predict(lm.milk, int = "p" , newdata = confidence)
prediction.band <- as.data.frame(cbind(prediction.band, parity = c(1:7)))
library(ggplot2)
ggplot(data = daisy2, aes(x = parity, y = milk120)) +
  geom_point() +
  geom_smooth(data = confidence.band, aes(x = parity, y = lwr), method = lm,
              se = FALSE, colour = "blue") + 
  geom_smooth(data = confidence.band, aes(x = parity, y = upr), method = lm,
              se = FALSE, colour = "blue") + 
  geom_smooth(data = prediction.band, aes(x = parity, y = lwr), method = lm,
              se = FALSE, colour = "green") + 
  geom_smooth(data = prediction.band, aes(x = parity, y = upr), method = lm,
              se = FALSE, colour = "green") + 
  xlab("Parity") + ylab("Milk volume in first 120 days of lactation")


# Ex. 14.2 regression with dichotomous predictor
lm.dyst <- lm(milk120 ~ as.factor(dyst), data = daisy2)
(lm.dyst.sum <- summary(lm.dyst))
anova(lm.dyst)  # anova table


# Could also treat parity as categorical
lm.milk2 <- lm(milk120 ~ as.factor(parity), data = daisy2)
(lm.milk2.sum <- summary(lm.milk2))
anova(lm.milk2)
# F-test
(src.var.total <- sum((daisy2$milk120 - mean(daisy2$milk120))^2))  # source of variation: total
(src.var.res <- sum(lm.milk2$res^2))  # source of variation: error (or residual)
(f.stat <- ((src.var.total - src.var.res) / (lm.milk2.sum$fstatistic[2])) /
  (src.var.res / lm.milk2.sum$fstatistic[3]))  # F-statistic
1 - pf(f.stat, lm.milk2.sum$fstatistic[2], lm.milk2.sum$fstatistic[3])
# R-squared
1 - (var(residuals(lm.milk2)) / var(daisy2$milk120, na.rm = TRUE))


# Ex. 14.3 Multiple variables
daisy2$milk120.sq <- daisy2$milk120^2
daisy2$parity.cat[as.numeric(daisy2$parity) >= 3] <- "3+"
daisy2$parity.cat[daisy2$parity == 2] <- "2"
daisy2$parity.cat[as.numeric(daisy2$parity) < 2] <- "1"
 
lm.milk2 <- lm(milk120 ~ parity + twin + dyst + rp + vag_disch,
               data = daisy2)
summary(lm.milk2)
anova(lm.milk, lm.milk2)


# Ex 14.4 Scaling predictor variables
daisy2$parity1 <- daisy2$parity - 1
lm.milk3 <- lm(milk120 ~ parity1,
               data = daisy2)
summary(lm.milk3)


# Ex 14.6 indicator variables
lm.milk4 <- lm(milk120 ~ parity + as.factor(herd),
               data = daisy2)
summary(lm.milk4)


# subsetting dataset
daisy2 <- daisy2[daisy2$h7 == 1, ]
dim(daisy2)
library(lubridate)
daisy2$date <- as.character(daisy2$calv_dt)
daisy2$date <- ymd(daisy2$date)
daisy2$mth <- month(daisy2$date)
daisy2$aut_calv <- with(daisy2, ifelse(mth %in% c(9:12), "fall", "other"))
daisy2$hs100 <- daisy2$herd_size / 100  # herd size scaled by dividing by 100
daisy2$hs100_sq <- daisy2$hs100^2  # scaled herd size scales
daisy2$hs100_ct <- daisy2$hs100 - mean(daisy2$hs100)  # centered
daisy2$hs100_ctsq <- daisy2$hs100_ct^2  # squared centered variable
daisy2$parity_sc <- daisy2$parity - mean(daisy2$parity)

daisy3 <- daisy2[complete.cases(daisy2), ]  # data frame with only complete cases

# ex. 14.8 centring to avoid collinearity
lm.wpc <- lm(wpc ~ hs100 + hs100_sq, data = daisy3)
summary(lm.wpc)
lm.wpc2 <- lm(wpc ~ hs100_ct + hs100_ctsq, data = daisy3)
summary(lm.wpc2)


# Ex. 14.9 interaction
lm.wpc3 <- lm(wpc ~ rp + vag_disch, data = daisy3)
summary(lm.wpc3)
lm.wpc4 <- lm(wpc ~ rp*vag_disch, data = daisy3)
summary(lm.wpc4)


# Ex. 14.10 interaction between dichotomous and continuous variable
daisy3$milk120k <- daisy3$milk120 / 1000
lm.wpc5 <- lm(wpc ~ as.factor(dyst)*milk120k, data = daisy3)
summary(lm.wpc5)


# Ex. 14.11 interaction between continuous variable
lm.wpc6 <- lm(wpc ~ parity*milk120k, data = daisy3)
summary(lm.wpc6)


# Ex. 14.12 impact of reproductive diseases on wpc
lm.wpc7 <- lm(wpc ~ aut_calv + hs100_ct + hs100_ctsq + parity_sc + twin +
             dyst + rp + vag_disch + rp*vag_disch, data = daisy3)
summary(lm.wpc7)
###Distribution of residuals should be normal with mean 0 (and thus medians should not be far from 0)
#anova(lm.wpc)

# Ex. 14.13 Homoscedasticity
# stat test: Cook-Weisberg test (or Breusch-Pagan test) (large p-values: constant variance)
library(car)
ncvTest(lm.wpc7)
library(lmtest)
bptest(lm.wpc7)
# bp.test() performs the same score test as ncvTest(), except that the default alternative hypothesis is different -- in bp.test() the error variance is a function of a linear combination of the regressors and in ncvTest() the error variance is a function of the fitted values (i.e., a *particular* linear combination of regressors). Testing against the fitted values with 1 df will have greater power if this is the real pattern of heteroscedasticity. (cf. John Fox https://stat.ethz.ch/pipermail/r-help/2009-September/212054.html)
# plot of standardised residuals vs. fitted values (which = 1)
plot(lm.wpc7, which = 1)  # diag plot in base stats package
ggplot(lm.wpc7, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_smooth(se = FALSE)

# Ex. 14.14 normality of residuals
# normal probability plot (Q-Q plot) (which = 2)
plot(lm.wpc7, which = 2)
ggplot(lm.wpc7, aes(sample = .stdresid)) +
  stat_qq() +
  geom_abline()
# p < 0.05 = non-normality
shapiro.test(resid(lm.wpc7))

# Scale-Location plot of sqrt(|residuals|) against fitted values (which = 3)
plot(lm.wpc7, which = 3)
ggplot(lm.wpc7, aes(.fitted, abs(.stdresid))) +
       geom_point() +
       geom_smooth(se = FALSE) +
       scale_y_sqrt()

# Cook's distance plot (which = 4)
plot(lm.wpc7, which = 4)
ggplot(lm.wpc7, aes(seq_along(.cooksd), .cooksd)) +
  geom_bar(stat = "identity")

# Residuals vs. leverages (which = 5)
plot(lm.wpc7, which = 5)
ggplot(lm.wpc7, aes(.hat, .stdresid)) +
  geom_vline(size = 2, colour = "white", xintercept = 0) +
  geom_hline(size = 2, colour = "white", yintercept = 0) +
  geom_point() +
  geom_smooth(se = FALSE)

# Cook's distances  vs leverages/(1-leverages) (which = 6)
plot(lm.wpc7, which = 6)
ggplot(lm.wpc7, aes(.hat, .cooksd)) +
  geom_vline(colour = NA) +
  geom_abline(slope = seq(0, 3, by = 0.5), colour = "white") +
  geom_smooth(se = FALSE) +
  geom_point()

# Ex. 14.15 linearity of predictor-outcome association
# plotting residuals against each of the continuous predictor variables
ggplot(lm.wpc7, aes(hs100_ct, .resid)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE)

# correcting distribution problems using robust standard errors
library(MASS)
rlm.wpc <- rlm(wpc ~ aut_calv + hs100_ct + hs100_ctsq + parity_sc + twin +
             dyst + rp + vag_disch + rp*vag_disch, data = daisy3)
summary(rlm.wpc)

# log transformation
daisy3$milk120.sq <- daisy3$milk120^2
lm.wpc8 <- lm(wpc ~ vag_disch + milk120 + milk120.sq, data = daisy3)
boxcox(lm.wpc8, lambda = seq(0, 0.2, by = 0.01), plotit = TRUE)
lm.lwpc <- lm(log(wpc) ~ vag_disch + milk120 + milk120.sq, data = daisy3)
summary(lm.lwpc)
(ci.lwpc <- confint(lm.lwpc))
exp(lm.lwpc$coefficients[2])
exp(ci.lwpc[2, 1])
exp(ci.lwpc[2, 2])
