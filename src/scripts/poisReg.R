# VER Ch.18: Modelling Count and Rate Data

# either download data from the web:
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/tb_real.rdata"))
#unlink(temp)  # remove the temporary file

# or use the unzipped archive
load("../../data/ver2_data_R/tb_real.rdata")

dim(tb_real)
str(tb_real)
library(Hmisc)
tb_real<- upData(tb_real, labels = c(farm_id = 'Farm identification',
                            type = 'Type of animal',
                            sex = 'Sex', age = 'Age category',
                            reactors = 'Number of pos/reactors in the group',
                            par = 'Animal days at risk in the group'),
                 levels = list(type = list('Dairy' = 1, 'Beef' = 2,
                                 'Cervid' = 3, 'Other' = 4),
                   sex = list('Female' = 0, 'Male' = 1),
                   age = list('0-12 mo' = 0, '12-24 mo' = 1, '>24 mo' = 2)))

# histogram
hist(tb_real$reactors, breaks = 0:30 - 0.5)

# Ex. 18.1 Poisson regression
mod1 <- glm(reactors ~ type + sex + age + offset(log(par)), family = poisson,
            data = tb_real)
summary(mod1)

mod1b <- glm(reactors ~ type + sex + age, offset = log(par), family = poisson,
            data = tb_real)
# IRR and CI:
cbind(exp(coef(mod1)), exp(confint(mod1)))
# deviance gof
mod1$deviance
pchisq(deviance(mod1), df.residual(mod1), lower = FALSE)
# Pearson statistic
mod1.pearson <- residuals(mod1, type = "pearson")
sum(mod1.pearson^2)
1-pchisq(sum(mod1.pearson^2), mod1$df.residual)


# Ex. 18.2 Diagnostics
plot(mod1, which = 1)  # std residuals vs fitted values
plot(mod1, which = 3)  # scale location
plot(mod1, which = 4)  # cooks
# Anscombe residuals
library(wle)
mod1.ansc <- residualsAnscombe(tb_real$reactors, mod1$fitted.values,
                               family = poisson())
plot(predict(mod1, type = "response"), mod1.ansc)
plot(cooks.distance(mod1), mod1.ansc)

# can perform chi-square tests for each term by using drop1 and looking at the
# changes in the deviance.
drop1(mod1, test = "Chisq")

# Investigate residuals
# Pearson residuals from Poisson regression model (w/o using offset)
mod1.pearson <- residuals(mod1, type = "pearson")
temp <- lm.influence(mod1)
names(temp)
# Standardized Pearson residuals
mod1.h <- lm.influence(mod1)$h  # also can use hatvalues(model = mod.fit)
head(mod1.h)
mod1.std.pearson <- mod1.pearson / sqrt(1 - mod1.h)
head(mod1.std.pearson)

X <- model.matrix(mod1)
mu.hat <- predict(mod1, type = "response")  # also could use mu.hat <- mod.fit$fitted.values
H <- diag(sqrt(mu.hat)) %*% X %*% solve(t(X) %*% diag(mu.hat) %*% X) %*% t(X) %*% diag(sqrt(mu.hat))
diag(H)[1:5]

par(mfrow = c(2,1))  # 2x1 grid of plots
# Pearson residual vs observation number plot
plot(x = 1:length(mod1.pearson), y = mod1.pearson, xlab = "Observation number",
     ylab = "Pearson residuals",
     main = "Pearson residuals vs. observation number")
abline(h = qnorm(c(0.975, 0.995, 0.025, 0.005)), lty = 3, col = "red")
# Standardized residual vs observation number plot
plot(x = 1:length(mod1.std.pearson), y = mod1.std.pearson,
     xlab = "Observation number", ylab = "Standardized residuals",
     main = "Standardized residuals vs. observation number")
abline(h = qnorm(c(0.975, 0.995, 0.025, 0.005)), lty = 3, col = "red")


# Ex. 18.3 Negative binomial
library(MASS)
mod2 <- glm.nb(reactors ~ type + sex + age + offset(log(par)), data = tb_real)
summary(mod2)

# Ex. 18.7 zinb
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/fec.rdata"))
#unlink(temp)
load("../../data/ver2_data_R/fec.rdata")

library(pscl)
mod3 <- zeroinfl(fec ~ lact + past_lact + man_heif + man_lact, data = fec,
                 dist = "negbin")
summary(mod3)
mod3.null <- update(mod3, . ~ 1)
# fit same model with negative binomial
mod3.nb <- glm.nb(fec ~ lact + past_lact + man_heif + man_lact, data = fec)
# vuong test
vuong(mod3, mod3.nb)
# Vuong statistic is 3.3 (p < 0.001) indicating the first model, the zero-inflated
# one, is superior to the regular negative binomial. Note that R does not estimate
# alpha but its inverse, theta . alpha is 3.9, suggesting a negative binomial is
# preferable to a Poisson model.
1/mod3$theta

# Ex. 18.8 hurdle
mod4 <- hurdle(fec ~ lact + past_lact + man_heif + man_lact, data = fec,
               dist = "negbin")
summary(mod4)
# can compare zero-inflated and hurdle models by their log-likelihood.
# Hurdle model fits better:
logLik(mod4)
logLik(mod3)
logLik(mod3.null)


# Ex. 18.9 Zero-truncated
#temp <- tempfile()
#download.file(
#  "http://ic.upei.ca/ver/sites/ic.upei.ca.ver/files/ver2_data_R.zip", temp)
#load(unz(temp, "ver2_data_R/daisy2.rdata"))
#unlink(temp)
load("../../data/ver2_data_R/daisy2.rdata")

library(VGAM)
mod5 <- vglm(spc ~ parity + cf + vag_disch, family = posnegbinomial(),
             data = daisy2)
summary(mod5)


#############################################################################
# Poisson regression and risk ratio

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
# recoding number of lactation
culling$lact <- with(culling, ifelse(lact_c3 > 1, 2, lact_c3))

log.reg <- glm(culled ~ lact, family = binomial("logit"), data = culling)
summary(log.reg)
confint(log.reg, parm = "lact")
exp(coef(log.reg.sum)[ , 1])
cbind(exp(coef(log.reg)), exp(confint(log.reg)))

with(culling, table(culled, lact))
## OR is 2.19:
(364 / 158) / (102 / 97)
## or the ratio between the cross-products
(364 * 97) / (158 * 102)
# or with epicalc
library(epicalc)
with(culling, cc(culled, lact, graph = FALSE))
# RR is 1.36:
(364 / 522) / (102 / 199)

# Common outcome: log link, binomial family
log.reg2 <- glm(culled ~ lact, family = binomial(link = "log"),
              data = culling)
summary(log.reg2)
cbind(exp(coef(log.reg2)), exp(confint(log.reg2)))

# Common outcome: log link, poisson family, robust estimator
# (modified Poisson with robust estimator by Zou, https://www.ncbi.nlm.nih.gov/pubmed/15033648)
# create id for each observation
culling$id <- 1:length(culling$herd)
library(geepack)
zou.mod <- geeglm(culled ~ lact, family = poisson(link = "log"), id = id,
                             corstr = "exchangeable", data = culling)
summary(zou.mod)
# exponentiated coefficients
exp(coef(zou.mod))
# getting confidence intervals
library(doBy)
zou.mod.coefci <- esticon(zou.mod, diag(2))
zou.mod.expci <- exp(cbind(zou.mod.coefci$Estimate, zou.mod.coefci$Lower,
                           zou.mod.coefci$Upper))
rownames(zou.mod.expci) <- names(coef(zou.mod))
colnames(zou.mod.expci) <- c("RR", "Lower RR", "Upper RR")
zou.mod.expci

# Fit by glm() then test using robust SE estimator
pois.reg <- glm(culled ~ lact, family = poisson(link = "log"),
                data = culling)
# Load sandwich package for robust estimator
library(sandwich) # to get robust estimator
## Load lmtest package for coeftest
library(lmtest) # to test coefficients
# Poisson model with SE estimated via robust variance estimator
coeftest(pois.reg, vcov = sandwich)
## RR
exp(coef(pois.reg))
## CI
0.3078 + qnorm(0.05 / 2) * 0.0749 # upper 95% CI
exp(0.3078 + qnorm(0.05 / 2) * 0.0749)
0.3078 + qnorm(1 - 0.05 / 2) * 0.0749 # lower 95% CI
exp(0.3078 + qnorm(1 - 0.05 / 2) * 0.0749)



zou.mod2 <- geeglm(culled ~ lact + johnes, family = poisson(link = "log"),
                   id = id, corstr = "exchangeable", data = culling)
summary(zou.mod2)
zou.mod2.coefci <- esticon(zou.mod2, diag(3))
zou.mod2.expci <- exp(cbind(zou.mod2.coefci$Estimate, zou.mod2.coefci$Lower,
                           zou.mod2.coefci$Upper))
rownames(zou.mod2.expci) <- names(coef(zou.mod2))
colnames(zou.mod2.expci) <- c("RR", "LowerCI", "UpperCI")
zou.mod2.expci

zou.df <- as.data.frame(zou.mod2.expci)
zou.df$var <- row.names(zou.df)
library(ggplot2)
ggplot(zou.df[2:3, ], aes(y = RR, x = reorder(var, RR))) +
  geom_point() +
  geom_pointrange(aes(ymin = LowerCI, ymax = UpperCI)) + 
  scale_x_discrete(labels = c("Lactation", "Johnes")) + 
  scale_y_log10(breaks = seq(1, 2, by = 0.1)) +
  xlab("") + 
  geom_hline(yintercept = 1, linetype = "dotdash", colour = "blue") +
  coord_flip()
library(gridExtra)
grid.table(round(zou.df[ , -4], digits = 3))
