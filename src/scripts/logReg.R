# VER Ch.16: Logistic Regression

# either download data from the web:
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/nocardia.rdata"))
#unlink(temp)  # remove the temporary file

# or use the unzipped archive
load("../../data/ver2_data_R/nocardia.rdata")

dim(nocardia)
str(nocardia)

# Ex. 16.1 logistic regression
mod1 <- glm(casecont ~ dcpct + dneo + dclox, family = binomial("logit"),
            data = nocardia)  # "logit" can be omitted as it is the default
(mod1.sum <- summary(mod1))
# confidence interval:
c(0.022668-1.96*0.007187, 0.022668+1.96*0.007187)
# or can also do:
confint.default(mod1, parm = "dcpct")  # CI constructed using normal approximation
                                       # for the parameter estimate
confint(mod1, parm = "dcpct")  # but better to construct a profile likelihood-based
                               # confidence interval (from MASS package)
exp(coef(mod1.sum)[ , 1])  # exponentiated coefficients
cbind(exp(coef(mod1)), exp(confint(mod1)))  # with CIs

mod1.null <- glm(casecont ~ 1, family = binomial, data = nocardia)
lr.mod1 <- -(deviance(mod1) / 2)
lr.mod1.null <- -(deviance(mod1.null) / 2)
(lr <- 2 * (lr.mod1 - lr.mod1.null))
1 - pchisq(lr, 2)  # LR test

# can also compare likelihood of the model under investigation to the likelihood
# of fully saturated model (one in which there would be one parameter fit for each
# data point)
pchisq(deviance(mod1), df.residual(mod1), lower = FALSE)

# can evaluate a single coefficient with a Wald test. If we add barn type to the
# model, we can test for an overall effect of barn type (b gives the coefficient,
# Sigma the variance covariance matrix of the error terms, Terms indicates which
# terms in the model to test (in the order they appear in the table of coefficients).
# The overall effect of barn is significant. We can also test for other barn against
# tie-stall barn.
mod1b <- glm(casecont ~ dcpct + dneo + dclox + as.factor(dbarn),
             family = binomial, data = nocardia)
library(aod)
wald.test(b = coef(mod1b), Sigma = vcov(mod1b), Terms = 4:5)
comp <- cbind(0, 0, 0, 0, 1, -1)
wald.test(b = coef(mod1b), Sigma = vcov(mod1b), L = comp)


# compare a full and nested model to test the contribution of any subset of parameters
mod1.nest <- glm(casecont ~ dcpct, family = binomial("logit"),
            data = nocardia)
lr.mod1.nest <- -(deviance(mod1.nest) / 2)
(lr <- 2 * (lr.mod1 - lr.mod1.nest))
1 - pchisq(lr, 2)
anova(mod1, mod1.nest, test = "Chisq")

# Ex. 16.2 interpreting coefficients
nocardia$dbarn <- as.factor(nocardia$dbarn)
mod2 <- glm(casecont ~ dcpct + dneo + dclox + dbarn,
              family = binomial("logit"), data = nocardia)
(mod2.sum <- summary(mod2))
cbind(exp(coef(mod2)), exp(confint(mod2)))
# note: book do not report the profile likelihood-based CI.
# Profile likelihood-based CI is preferable due to the Hauck-Donner
# effect (http://www.jstor.org/stable/2286473) (the SE are overestimated and the
# z-value is too small)
### interpretation of odds ratios: using neomycin in the herd increased the log
### odds of Nocardia mastitis by 2.7 units (or using neomycin increases the odds
### 14.7 times). Since Nocardia mastitis is a rare condition, we can interpret
### the OR as a risk ratio (RR) and say neomycin increases the risk of Nocardia
### mastitis by 15 times.
### Changing the percentage of dry cows treated from 50 to 75% increases the log
### odds of disease by (75-50)*0.022=0.55 units (or 1.022^(75-50)=1.72).
### An increase of 25% in the percentage of cows dry-treated increases the risk
### of disease by about 72% (i.e. 1.72 times).
### Tiestall barns and other barn types both have lower risks of Nocardia mastitis
### than freestall barns.

# Mutiple comparisons
library(multcomp)
mod2.mc <- glht(mod2, linfct = mcp(dbarn = "Tukey"))
(mod2.mc.sum <- summary(mod2.mc))

# Significance of the main effects can be tested with:
drop1(mod2, test = "Chisq")
# the drop1 function tests each predictor relative to the full model.

# Ex. 16.3 Effect of factors on the probability scale
nocardia$neo.pred <- predict(mod1, type = "response", se.fit = FALSE)
library(ggplot2)
ggplot(nocardia, aes(x = dcpct, y = neo.pred, group = dneo,
                     colour = factor(dneo))) + 
  stat_smooth(method = "glm", family = "binomial", se = FALSE) +
  labs(colour = "Neomycin", x = "Percent of cows dry treated",
       y = "Probability of Nocardia")


# Ex. 16.6 Covariate pattern
mod3 <- glm(casecont ~ dcpct + dneo + dclox + dneo*dclox,
            family = binomial("logit"), data = nocardia)
summary(mod3)
library(epiR)
mod3.mf <- model.frame(mod3)
(mod3.cp <- epi.cp(mod3.mf[-1]))

# Ex. 16.7 Pearson and deviance residuals
(res.dev <- residuals(mod3))  # deviance residuals
residuals(mod3, "pearson")  # pearson residuals
residuals(mod3, "response")  # response residuals
residuals(mod3, "working")  # working residuals

plot(predict(mod3), res.dev, xlab = "Fitted values", ylab = "Residuals",
     ylim = max(abs(res.dev)) * c(-1, 1))
abline(h = 0, lty = 2)
# Pearson and deviance chi-sq tests
sum(residuals(mod3, type = "pearson")^2)
deviance(mod3)
# based on individual observations:
1 - pchisq(deviance(mod3), df.residual(mod3))
# p-value is large indicating no evidence of lack of fit
1 - pchisq(sum(residuals(mod3, type = "pearson")^2), df.residual(mod3))

# Hosmer-Lemeshow test
library(Hmisc)
hosmerlem <- function (y, yhat, g = 10) {
    cutyhat <- cut2(yhat, g = g)
    obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    chisq <- sum((obs - expect)^2 / expect)
    P <- 1 - pchisq(chisq, g - 2)
    c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}
hosmerlem(y = nocardia$casecont, yhat = fitted(mod3))
# The model used has many ties in its predicted probabilities (too few covariate
# values?) resulting in an error when running the Hosmer-Lemeshow test. Using
# fewer cut-points (g = 5 or 7) does not solve the problem. This is a typical
# example when not to use this test. A better goodness-of-fit test than
# Hosmer-Lemeshow and Pearson / deviance chi^2 tests is the
# le Cessie-van Houwelingen-Copas-Hosmer unweighted sum of squares test
# for global goodness of fit (also here) implemented in the rms package (but you
# have to implement your model with the lrm function of this package):
# @ARTICLE{hos97com,
#        author = {Hosmer, D. W. and Hosmer, T. and {le Cessie}, S. and 
#     Lemeshow, S.},
#        year = 1997,
#        title = {A comparison of goodness-of-fit tests for the logistic 
#     regression model},
#        journal = Statistics in Medicine,
#        volume = 16,
#        pages = {965-980},
#        annote = {goodness-of-fit for binary logistic model;difficulty with
# Hosmer-Lemeshow statistic being dependent on how groups are defined; sum of
# squares test; cumulative sum test; invalidity of naive test based on deviance;
# goodness-of-link function;simulation setup}
library(rms)
mod3b <- lrm(casecont ~ dcpct + dneo + dclox + dneo*dclox, nocardia,
             method = "lrm.fit", model = TRUE, x = TRUE, y = TRUE,
             linear.predictors = TRUE, se.fit = FALSE)
residuals(mod3b, type = "gof")
# p-value is 0.16 so there's no evidence the model is incorrect. Even better
# than these tests would be to check for linearity of the predictors.


# Ex. 16.8 Overdispersion
library(faraway)
halfnorm(residuals(mod1))
# dispesion parameter of mod1:
(sigma2 <- sum(residuals(mod1, type = "pearson")^2) / 104)
drop1(mod1, scale = sigma2, test = "F")
# The dispersion parameter is not very different than one (no dispersion).
# If dispersion was present, you could use it in the F-tests for the predictors,
# adding scale to drop1.

# Ex. 16.10 Predictive ability
nocardia <- upData(nocardia,
                   levels = list(casecont = list('Control' = 0, 'Case' = 1)))
predicted <- predict(mod3)
library(ROCR)
prob <- prediction(predicted, nocardia$casecont, 
                   label.ordering = c('Control', 'Case'))
tprfpr <- performance(prob, "tpr", "fpr")
tpr <- unlist(slot(tprfpr, "y.values"))
fpr <- unlist(slot(tprfpr, "x.values"))
roc <- data.frame(tpr, fpr)
ggplot(roc) + 
  geom_line(aes(x = fpr, y = tpr)) + 
  geom_abline(intercept = 0, slope = 1, colour = "gray") + 
  ylab("Sensitivity") + 
  xlab("1 - Specificity")


# Ex. 16.11 Identifying important obsevrations
nocardia$casecont.num <- as.numeric(nocardia$casecont) - 1
mod1.mf <- model.frame(mod1)
mod1.cp <- epi.cp(mod1.mf[-1])
nocardia.cp <- as.data.frame(cbind(cpid = mod1.cp$id,
                                   nocardia[ , c(1, 9:11, 13)],
                                   fit = fitted(mod1)))
# Residuals and delta betas based on covariate pattern:
mod1.obs <- as.vector(by(as.numeric(nocardia.cp$casecont.num),
                         as.factor(nocardia.cp$cpid), FUN = sum))
mod1.fit <- as.vector(by(nocardia.cp$fit, as.factor(nocardia.cp$cpid),
                         FUN = min))
mod1.res <- epi.cpresids(obs = mod1.obs, fit = mod1.fit,
                         covpattern = mod1.cp)
# Plot of Pearson residuals as a function of the linear predictor (the log-odds of success). Ideally, this plot should resemble a horizontal band with most observations within the range -3 to +3:
mod1.lodds <- as.vector(by(predict(mod1), as.factor(nocardia.cp$cpid),
                           FUN = min))
plot(mod1.lodds, mod1.res$spearson,
     type = "n", ylab = "Pearson Residuals", xlab = "Logit")
text(mod1.lodds, mod1.res$spearson, labels = mod1.res$cpid, cex = 0.8)
symbols(mod1.lodds, mod1.res$pearson, circles = mod1.res$n, add = TRUE)
abline(h = 0, lty = 2)

# this one is in Dohoo (but not quite the same, guess Ian did not use the model
# he's pretending)
plot(mod1.cp$cov.pattern$id, mod1.res$spearson, type = "n",
     ylab = "Pearson Residuals", xlab = "Covariate Pattern")
text(mod1.cp$cov.pattern$id, mod1.res$pearson, labels = mod1.res$cpid,
     cex = 0.8)
symbols(mod1.cp$cov.pattern$id, mod1.res$spearson, circles = mod1.res$n,
        add = TRUE)


# Ex. 16.12 Identify the covariate patterns with the largest leverage values:
mod1.res[sort.list(mod1.res$leverage, decreasing = TRUE), ]
# The covariate patterns with the largest leverage values are 2, 14, 12, and 5.
# Identify the covariate patterns with the largest delta-betas:
mod1.res[sort.list(mod1.res$sdeltabeta, decreasing = TRUE), ]
# Covariate pattern 2 has the largest delta beta. This is not surprising because
# it contained 38 observations (around one-third of the data).


# Ex 16.13 Exact logistic regression
# robust SE
library(robust)
mod4 <- glmrob(casecont ~ dcpct + dneo + dclox + dneo*dclox,
               family = binomial, data = nocardia, method = "Mqle",
               control= glmrobMqle.control(tcc = 1.2))
summary(mod4)

# Ex. 16.15 Conditional logistic regression
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/sal_outbrk.rdata"))
#unlink(temp)
load("../../data/ver2_data_R/sal_outbrk.rdata")

library(survival)
mod5 <- clogit(casecontrol ~ slt_a + strata(match_grp), data = sal_outbrk)
summary(mod5)

