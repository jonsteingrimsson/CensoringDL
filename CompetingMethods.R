library(randomForestSRC)
library(survival)
library(glmnet)

job.number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# A function that takes in observed time, failure time indicator, and covariate
# vector and returns a cox model selected using glmnet
cov.mod.select = function(obs, delta, cov.vec){
# Creating a survival data frame
data.surv = data.frame(obs, delta, cov.vec)
  # Doing a cross validation across lambda
  cvfit = cv.glmnet(as.matrix(cov.vec), Surv(obs, delta), family = "cox")
  # Finding optimal lambda
  lambda.min = cvfit[['lambda.min']]
  # Finding the active set
  coef.min = coef(cvfit, s = "lambda.min")
  active.min = which(coef.min != 0)
  
return(active.min = active.min)
}

# A function that takes in a cox model and a test set covariate vector and outputs an estimator
# for $E[\log(T)|W] for W in the test set.

cox.calc.exp = function(cox.mod, cov.vec.test, active.min = NULL, tau1)
{
# If the cox model has at least one covariate
if(is.null(active.min) | length(active.min) > 0){
# Getting the Survival Curves. 
cox.surv <- survfit(cox.mod, newdata= cov.vec.test)
nu <- cox.mod$n
nu.test = nrow(cov.vec.test)

# Calculation of the survivor function
m1 = rep(NA, nu.test)

for(i in 1:nu.test){
# Getting the properties of the survival cox function for observation i.
time.used <- cox.surv[i][['time']]
# Calculating the jumps in the Cox model survival curve estimator
surv.diff <- c(1, cox.surv[i][['surv']][-length(cox.surv[i][['surv']])]) - cox.surv[i][['surv']]

# Calculating the conditional expectation
m1[i]= sum((log(time.used) * surv.diff)[time.used < tau1] )
}
}

# If there is no covariate in the model
if(length(active.min) == 0){
nu <- cox.mod$n
nu.test = nrow(cov.vec.test)

# Calculation of the survivor function
m1 = rep(NA, nu.test)

for(i in 1:nu.test){
# Getting the properties of the survival cox function for observation i.
time.used <- cox.mod$time

# Calculating the jumps in the Cox model survival curve estimator
surv.diff <- c(1, cox.mod$surv[-length(cox.mod[['surv']])]) - cox.mod[['surv']]
# Calculating the conditional expectation
m1[i]= sum((log(time.used) * surv.diff)[time.used < tau1])
}
}

return(m1)
}


compare.methods.brier <- function(){
pred.err.cox.cv = rep(NA, 5)
pred.err.cox.pen.cv = rep(NA, 5)
pred.err.rand.cv = rep(NA, 5)

for(j in 1:5){
# Creating data
file.name = paste("x_train", toString(job.number), "CV", toString(j), ".rda", sep = "")
load(file.name)
file.name = paste("x_test", toString(job.number),"CV", toString(j), ".rda", sep = "")
load(file.name)
file.name = paste("y_test", toString(job.number), "CV", toString(j), ".rda", sep = "")
load(file.name)
file.name = paste("data_surv", toString(job.number), "CV", toString(j), ".rda", sep = "")
load(file.name)
load("Tau1.rda")
file.name = paste("ipcw_weights", toString(job.number), "CV", toString(j), ".rda", sep = "")
load(file.name)


y.brier = (y_test > time.point)
data.test = data.frame(x_test)

if(dim(x_train)[1] > dim(x_train)[2]){
# Calculating the prediction error using Cox-model
cox.mod = coxph(Surv(obs, delta)~., data = data.surv)
data.test = data.frame(x_test)
cox.surv <- survfit(cox.mod, newdata= data.test)
pred.cox = rep(NA, nrow(data.test))
for(i in 1:nrow(data.test)){
pred.cox[i] = cox.surv[i]$surv[sum(cox.surv[i]$time < time.point)]
}
pred.err.cox.cv[j] = mean(data.ipcw * (y.brier - pred.cox)^2)
}

not.active = TRUE
if(not.active == TRUE){
# Calculate prediction error using a penalized cox model
cov.vec.dummy = model.matrix( ~ ., data=data.surv[, -c(1,2)])

active.min = cov.mod.select(data.surv$obs, data.surv$delta, cov.vec.dummy)
data.surv.new = data.frame(data.surv$obs, data.surv$delta, cov.vec.dummy)
data.surv.used = data.surv.new[, c(1:2, active.min +2)]
names(data.surv.used)[1:2] = c("obs", "delta")
test_data.dummy = model.matrix( ~ ., data=data.test)
test.data.used = data.frame(test_data.dummy[, active.min])
print(names(test.data.used))
print(names(data.surv.used)[-c(1,2)])
names(test.data.used) = names(data.surv.used)[-c(1,2)]
print(names(test.data.used))
print(names(data.surv.used)[-c(1,2)])

if(length(active.min) > 0){
cox.mod.pen <- coxph(Surv(obs, delta) ~., data =data.surv.used)
cox.surv <- survfit(cox.mod.pen, newdata= test.data.used)
pred.cox.pen = rep(NA, nrow(data.test))
for(i in 1:nrow(data.test)){
pred.cox.pen[i] = cox.surv[i]$surv[sum(cox.surv[i]$time < time.point)]
}
}
if(length(active.min) == 0){
cox.mod.pen = survfit(Surv(data.surv$obs, data.surv$delta)~1)
pred.cox.pen = rep(cox.mod.pen$surv[sum(cox.mod.pen$time < time.point)], nrow(data.test))
}
pred.err.cox.pen.cv[j] = mean(data.ipcw * (pred.cox.pen - y.brier)^2)
}

# Fitting random survival forest
rand.for = rfsrc(Surv(obs, delta)~., data = data.surv, ntree = 1000)
pred.rand.for = predict(rand.for, newdata = data.test)
used.time <- sum(pred.rand.for[['time.interest']] <= time.point) + 1
surv.est.default <- rep(0, nrow(data.test))
  for(i in 1:nrow(data.test)){
    # Calculating P(T > tau1|W_i)
    surv.est.default[i] <- c(1, pred.rand.for$survival[i, ])[used.time]
  }
pred.err.rand.cv[j] = mean(data.ipcw * (y.brier - surv.est.default)^2)
print(j)
}

if(dim(x_train)[1]> dim(x_train)[2]){object.return = list(c(pred.err.cox.pen = mean(pred.err.cox.pen.cv), pred.err.cox = mean(pred.err.cox.cv), pred.err.rand=mean(pred.err.rand.cv)))}
#object.return = list(c(pred.err.cox.pen = pred.err.cox.pen, pred.err.rand=pred.err.rand))
#if(not.active != TRUE){object.return = list(c(pred.err.cox = pred.err.cox, pred.err.rand=pred.err.rand))}
return(object.return)
}

other.methods = compare.methods.brier()


