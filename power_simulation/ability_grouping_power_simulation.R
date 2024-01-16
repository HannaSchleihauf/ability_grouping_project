# title: "ability_grouping_power_simulation"
# author: "Hanna Schleihauf"
# date: "16/01/2024"

options(scipen = 9999)

## LOAD PACKAGES AND FUNCTIONS ################################################
library("tidyverse")
library("cowplot")
library("gghalves")
library("lme4")
library("kyotil") # we want to store info about convergence issues
library(glmmTMB)
source("./functions/beta_par_to_par_transf.r")

## GENERATE DATA ##############################################################
set.seed(1)
n.subj <- 80 # number of subj
n.simus <- 1000 # small number for testing
# range of age
min.age <- 6
max.age <- 11

n.per.subj <- 4 # observations per subj
n.group.group <- 2 # number of groups
n.per.group.group <- 1 # observations per subj and group (between-subj)
group <- ordered(c(".no_group", ".group"))
group <- ordered(group, levels = c(".no_group", ".group"))
n.ability.group <- 2 # number of groups
n.per.ability.group <- 2 # observations per subj and group (within-subj)
ability <- ordered(c(".low", ".high"))
ability <- ordered(ability, levels = c(".low", ".high"))
subj.id <- as.factor(paste("subj", 1:n.subj, sep = ".")) # creating subj ids
task <- c("krumen", "zarpen")

# expected performance levels
no_group.low_ability.6 <- 0.60 # <- reference level
no_group.high_ability.6 <- 0.60
group.low_ability.6 <- 0.75
group.high_ability.6 <- 0.75

no_group.low_ability.11 <- 0.50
no_group.high_ability.11 <- 0.50
group.low_ability.11 <- 0.80
group.high_ability.11 <- 0.80

start.data <- data.frame(subj.id)

# duplicate rows according to the number obs. per subj:
start.data <- start.data[rep(x = 1:nrow(start.data), times = n.per.subj), ]
start.data <- as.data.frame(start.data)
names(start.data) <- "subj.id"

# create within-subj predictor
# add ability and task type
start.data <- data.frame(expand.grid(
  subj.id = subj.id,
  ability = ability,
  task = task
))
# predictor trial number
start.data$trial <-
  ave(seq_along(start.data$subj.id), start.data$subj.id, FUN = seq_along)

# create between-subj predictors
# predictor group
start.data$group <-
  as.factor(rep(x = group, each = n.subj /
    2))[as.numeric(start.data$subj.id)]
# predictor sex
start.data$gender <-
  as.factor(rep(
    x = c(".male", ".female", ".male", ".female"),
    each = n.subj / 4
  ))[as.numeric(start.data$subj.id)]
# predictor age
start.data$age <-
  rep(x = runif(
    n = n.subj,
    min = min.age, max = max.age
  ))[as.numeric(start.data$subj.id)]
# check whether it worked
ftable(group ~ gender, start.data) / 4 # worked!

# z-transformation of covariates
# z-transform age
start.data$z.age <- as.vector(scale(start.data$age))
start.data$z.trial <- as.vector(scale(as.numeric(start.data$trial)))

# dummy code factors and center them for random slopes
start.data$group.group <-
  as.numeric(start.data$group == levels(start.data$group)[2])
start.data$group.c <-
  as.numeric(start.data$group.group) -
  mean(as.numeric(start.data$group.group)) # centering

start.data$ability.high <-
  as.numeric(start.data$ability == levels(start.data$ability)[2])
start.data$ability.c <-
  as.numeric(start.data$ability.high) -
  mean(as.numeric(start.data$ability.high)) # centering

# dummy code and center gender to make the estimates ungroupal of
# the reference category
start.data$gender.male <-
  as.numeric(start.data$gender == levels(start.data$gender)[2])
start.data$gender.male.c <-
  start.data$gender.male - mean(start.data$gender.male)

# checks whether we created the correct data drame
# does each subj have only one sex and age?
xx <- table(start.data$subj.id, start.data$gender)
range(apply(X = xx > 0, MARGIN = 1, sum)) # should be 1 and 1

xx <- table(start.data$subj.id, start.data$age)
range(apply(X = xx > 0, MARGIN = 1, sum)) # should be 1 and 1

xx <- table(start.data$subj.id, start.data$group)
range(apply(X = xx > 0, MARGIN = 1, sum))

xx <- table(start.data$subj.id, start.data$ability)
range(apply(X = xx > 0, MARGIN = 1, sum))

xx <- table(start.data$subj.id, start.data$trial)
range(apply(X = xx > 0, MARGIN = 1, sum))

## CALCULATE ESTIMATES / SLOPES BASED ########################################

# to calculate slope between two point one need to (y2-y2)/(x2-x1)
# reference levels
intercept <-
  mean(c(qlogis(no_group.low_ability.6), qlogis(no_group.low_ability.11)))

# slope of age is the effect with the other factors being on
# their reference levels (ability = low, group = no group)
s.age <-
  (qlogis(no_group.low_ability.11) - qlogis(no_group.low_ability.6)) /
    (max(start.data$z.age) - min(start.data$z.age))

# slope of group is the effect with age being at its average (0) or
# the average of the slopes against group at age
# being at its minimum and maximum, respectively
# ability is at its reference level
s.group <-
  mean(c(
    (qlogis(group.low_ability.11) - qlogis(no_group.low_ability.11)),
    (qlogis(group.low_ability.6) - qlogis(no_group.low_ability.6))
  ))

# slope of ability is the effect with age being at its average (0) or
# the average of the slopes against ability at age
# being at its minimum and maximum, respectively
# group is at its reference level
s.ability <-
  mean(c(
    (qlogis(no_group.high_ability.11) - qlogis(no_group.low_ability.11)),
    (qlogis(no_group.high_ability.6) - qlogis(no_group.low_ability.6))
  ))

# two-way interactions
# slope of age in the group group with ability at the reference level
# We need to calculate:
# ((slope of age with group being at its maximum(group)  -
# slope of age with group being at its minimum(no.group)) /
# (maximum of group(group) - min of group(no.group))
s.age.group <-
  ((qlogis(group.low_ability.11) - qlogis(group.low_ability.6)) /
    (max(start.data$z.age) - min(start.data$z.age))) - # slope of age at group being at its maximum (group)
  ((qlogis(no_group.low_ability.11) - qlogis(no_group.low_ability.6)) /
    (max(start.data$z.age) - min(start.data$z.age))) # slope of age with group being at its minimum (no.group)

# since it also needs to work the other way around
s.group.age <-
  ((qlogis(group.low_ability.11) - qlogis(no_group.low_ability.11)) -
    (qlogis(group.low_ability.6) - qlogis(no_group.low_ability.6))) /
    (max(start.data$z.age) - min(start.data$z.age))
# test whether both versions lead to the same result
round(s.age.group, 6) == round(s.group.age, 6)

# slope of age in ability with group at the reference level
# We need to calculate:
# ((slope of age with ability being at its maximum(high)  -
# slope of age with ability being at its minimum(low)) /
# (maximum of ability(high) - min of ability(low))
s.age.ability <-
  ((qlogis(no_group.high_ability.11) - qlogis(no_group.high_ability.6)) /
    (max(start.data$z.age) - min(start.data$z.age))) - # slope of age at ability being at its maximum (high)
  ((qlogis(no_group.low_ability.11) - qlogis(no_group.low_ability.6)) /
    (max(start.data$z.age) - min(start.data$z.age))) # slope of age with ability being at its minimum (low)

# since it also needs to work the other way around
s.ability.age <-
  ((qlogis(no_group.high_ability.11) - qlogis(no_group.low_ability.11)) -
    (qlogis(no_group.high_ability.6) - qlogis(no_group.low_ability.6))) /
    (max(start.data$z.age) - min(start.data$z.age))
# test whether both versions lead to the same result
round(s.age.ability, 6) == round(s.ability.age, 6)

# ((slope of group with ability being at its maximum(high)  -
# slope of group with ability being at its minimum(low)) /
# (maximum of ability(high) - min of ability(low))
s.group.ability <-
  mean(c(
    ((qlogis(group.high_ability.6) - qlogis(no_group.high_ability.6)) -
      (qlogis(group.low_ability.6) - qlogis(no_group.low_ability.6))),
    ((qlogis(group.high_ability.11) - qlogis(no_group.high_ability.11)) -
      (qlogis(group.low_ability.11) - qlogis(no_group.low_ability.11)))
  ))
# since it also needs to work the other way around
s.ability.group <-
  mean(c(
    ((qlogis(group.high_ability.6) - qlogis(group.low_ability.6)) -
      (qlogis(no_group.high_ability.6) - qlogis(no_group.low_ability.6))),
    ((qlogis(group.high_ability.11) - qlogis(group.low_ability.11)) -
      (qlogis(no_group.high_ability.11) - qlogis(no_group.low_ability.11)))
  ))
# test whether both versions lead to the same result
round(s.ability.group, 6) == round(s.group.ability, 6)

# three-way-interaction
s.group.ability.age <-
  # slope for group at max age for high_ability
  ((((qlogis(group.high_ability.11) - qlogis(no_group.high_ability.11)) -
    # slope for group at min age for high_ability
    (qlogis(group.high_ability.6) - qlogis(no_group.high_ability.6))) /
    (max(start.data$z.age) - min(start.data$z.age))) -
    # slope for group at max age for us
    (((qlogis(group.low_ability.11) - qlogis(no_group.low_ability.11)) -
      # slope for group at min age for us
      (qlogis(group.low_ability.6) - qlogis(no_group.low_ability.6))) /
      (max(start.data$z.age) - min(start.data$z.age))))

# create vector with coefficients
coefs <- c(
  "(Intercept)" = intercept,
  "z.age" = s.age,
  "group.L" = s.group,
  "ability.L" = s.ability,
  "gender.male.c" = 0,
  "z.trial" = 0,
  "z.age:group.L" = s.group.age,
  "z.age:ability.L" = s.ability.age,
  "group.L:ability.L" = s.group.ability,
  "z.age:group.L:ability.L" = s.group.ability.age
)
coefs.t <- as.data.frame(t(as.data.frame(coefs)))

str(start.data)
## DEFINE RANDOM EFFECTS AND SLOPES ########################################

# random effect
# educated guess of what the random effect could be (based on the qlogis of the reference level performance)
tiny.re <- abs(qlogis(no_group.low_ability.6) / 8)
moderate.re <- abs(qlogis(no_group.low_ability.6) / 4)
strong.re <- abs(qlogis(no_group.low_ability.6) / 2)
extrem.re <- abs(qlogis(no_group.low_ability.6) * 1)
r.effects <- c(tiny.re, moderate.re, strong.re)

# random slope for trial
r.slope.trial <- 0.1
# random slope for trial
r.slope.ability <- 0.2

subj.sd <- 1.75 # random effect of subject (standard deviation)
resid.sd <- 0.5 # residual standard deviation

## CREATE OBJECT TO STORE RESULTS #############################################

# create object to store the simulation parameters and results:
all.res <-
  data.frame(expand.grid(
    r.effect = r.effects,
    r.slope.trial = r.slope.trial,
    r.slope.ability = r.slope.ability,
    simu = 1:n.simus
  ))

# add columns for estimates
all.res$icpt <- NA
all.res$z.age <- NA
all.res$group.group <- NA
all.res$ability.high <- NA
all.res$gender.male.c <- NA
all.res$z.trial <- NA
all.res$z.age.group.group <- NA
all.res$z.age.ability.high <- NA
all.res$group.group.ability.high <- NA
all.res$z.age.group.group.ability.high <- NA

# add columns for re.sd and warnings for full model and null model
all.res$re.sd <- NA
all.res$re.ability <- NA
all.res$re.trial <- NA
all.res$warns.full <- NA
all.res$warns.null <- NA
# add columns for likelihood ratio test results (p-values)
all.res$full.null.p <- NA
all.res$lrt.p.group <- NA
all.res$lrt.p.z.age <- NA
all.res$lrt.p.ability <- NA
all.res$lrt.p.gender.male.c <- NA
all.res$lrt.p.z.trial <- NA
all.res$lrt.p.z.age.group <- NA
all.res$lrt.p.z.age.ability <- NA
all.res$lrt.p.group.ability <- NA
all.res$lrt.p.z.age.group.ability <- NA

## START SIMULATION ########################################################
xdata <- start.data # change start.data to xdata

# create model matrix
m.mat <- model.matrix(object = ~ z.age * group * ability + gender.male.c + z.trial, data = xdata) # create model matrix
# create LP wrt fixed effects
LP <- m.mat[, names(coefs)] %*% coefs

# define control structure to make convergence more likely:
# define control structure to make convergence more likely:
contr <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))

# run simulation
for (i in 1:nrow(all.res)) {
  # i=1
  set.seed(i) # allows to later replicate individual simulations

  # add random effect to linear predictor:
  xdata$response <- LP +
    # random intercept:
    rnorm(n = n.subj, sd = all.res[i, "r.effect"])[as.numeric(xdata$subj.id)] +
    # random slope of trial number:
    rnorm(
      n = n.subj,
      sd = all.res[i, "r.slope.trial"]
    )[as.numeric(xdata$subj.id)] *
      xdata$z.trial +
    # random slope of ability:
    rnorm(
      n = n.subj,
      sd = all.res[i, "r.slope.ability"]
    )[as.numeric(xdata$subj.id)] *
      xdata$ability.c  +
    # random intercept of subject
    rnorm(n = n.subj, sd = subj.sd)[as.numeric(xdata$subj.id)] +
    rnorm(n = length(xdata$subj.id), sd = resid.sd) # residual variation

  # fit full model:
  full <- keepWarnings(lmer(response ~
    (z.age + group + ability)^3 + gender.male.c + z.trial +
    (1 + ability.c + z.trial | subj.id),
  data = xdata
  ))

  # fit reduced model:
  red1 <- keepWarnings(lmer(response ~
    (z.age + group + ability)^2 + gender.male.c + z.trial +
    (1 + ability.c + z.trial | subj.id),
  data = xdata
  ))

  # fit main effects model:
  main <- keepWarnings(lmer(response ~
    (z.age + group + ability) + gender.male.c + z.trial +
    (1 + ability.c + z.trial | subj.id),
  data = xdata
  ))

  # fit null model:
  null <- keepWarnings(lmer(response ~
    1 +
    (1 + ability.c + z.trial | subj.id),
  data = xdata
  ))

  # store results:
  all.res[i, c(
    "icpt",
    "z.age", "group.group", "ability.high", "gender.male.c",
    "z.trial", "z.age.group.group", "z.age.ability.high",
    "group.group.ability.high",
    "z.age.group.group.ability.high"
  )] <- fixef(full$value)

  xx <- sapply(summary(full$value)$varcor, function(x) attr(x, "stddev"))
  all.res[i, "re.sd"] <- xx[1]
  all.res[i, "re.ability"] <- xx[2]
  all.res[i, "re.trial"] <- xx[3]

  all.res[i, "warns.full"] <- nchar(paste(full$warnings, collapse = ""))
  all.res[i, "warns.null"] <- nchar(paste(null$warnings, collapse = ""))
  all.res[i, "full.null.p"] <-
    as.data.frame(anova(null$value, full$value, test = "Chisq"))[2, "Pr(>Chisq)"]
  all.res[i, "lrt.p.z.age.group.ability"] <-
    as.data.frame(drop1(full$value, test = "Chisq"))["z.age:group:ability", "Pr(Chi)"]

  xx <- drop1(red1$value, test = "Chisq")
  all.res[i, "lrt.p.z.age.ability"] <- as.data.frame(xx)["z.age:group", "Pr(Chi)"]
  all.res[i, "lrt.p.group.ability"] <- as.data.frame(xx)["z.age:ability", "Pr(Chi)"]
  all.res[i, "lrt.p.z.age.group"] <- as.data.frame(xx)["group:ability", "Pr(Chi)"]

  xx <- drop1(main$value, test = "Chisq")
  all.res[i, "lrt.p.group"] <- as.data.frame(xx)["group", "Pr(Chi)"]
  all.res[i, "lrt.p.ability"] <- as.data.frame(xx)["ability", "Pr(Chi)"]
  all.res[i, "lrt.p.z.age"] <- as.data.frame(xx)["z.age", "Pr(Chi)"]
  all.res[i, "lrt.p.z.trial"] <- as.data.frame(xx)["z.trial", "Pr(Chi)"]
  all.res[i, "lrt.p.gender.male.c"] <- as.data.frame(xx)["gender.male.c", "Pr(Chi)"]

  print(i)
}

save.image("ability_grouping_power_simulation.RData")
load("ability_grouping_power_simulation.RData")


## EVALUATION ################################################################

# full model
all.res1 <- all.res
tapply(
  X = all.res1[, "warns.full"] > 0, INDEX = all.res1[, c("r.effect")],
  FUN = sum
)
sum(all.res1[, "warns.full"] > 0)

# null model:
tapply(
  X = all.res1[, "warns.null"] > 0, INDEX = all.res1[, c("r.effect")],
  FUN = sum
)
all.res1

lrt.data <- all.res1 %>%
  filter(full.null.p < 0.05) %>%
  group_by(r.effect) %>%
  summarise(
    proportion.lrt.group = length(lrt.p.group[lrt.p.group <= 0.05]) / n.simus,
    proportion.lrt.z.ability.group = length(lrt.p.group.ability[lrt.p.group.ability < 0.05]) / n.simus,
    proportion.lrt.z.age.group = length(lrt.p.z.age.group[lrt.p.z.age.group < 0.05]) / n.simus,
    proportion.lrt.z.age.group.ability = length(lrt.p.z.age.group.ability[lrt.p.z.age.group.ability < 0.05]) / n.simus

  )
lrt.data

