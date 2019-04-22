# This function creates the imputations used for the deep learning 
# algorithm.

# It includes 4 key functions creating imputations for each combination of
# loss function and full data loss function
library(rpart)
library(survival)
library(party)
library(randomForestSRC)
library(MASS)
library(pec)
library(glmnet)


GfuncKM=function(obs,delta,dtype)
{ 
  n <- length(obs)
  # Changing the dataset to account for truncation. 
  aa=datach(obs,delta,dtype)
  # Observed time after truncation
  obs=aa[[1]]
  # Failure indicator after truncation.
  delta=aa[[2]]
  
  #Calculating the hazard estimator. 
  hazC=mapply(function(xx,dd){dd/sum((obs>=xx))},xx=obs,dd=1-delta)
  surC_km=mapply(function(xx){prod(1-hazC[obs<=xx])},xx=obs)
  return(list(surC_km,obs,delta))
}


# Calculates The estimator for the censoring distribution
# Input: obs = observed time, delta = failure indicator, dtype = truncation level and method.
# Output: List of bar -(tilde T|W), failure indicator changed to account for truncation
# and observed time also changed to account for truncation.
GfuncSurvivalTree=function(obs,delta,dtype,xx)
{

num = length(obs)
nu = num

# Fitting the Cox model
# Creating the data frame
data.used <- data.frame(obs, 1 - delta, xx)
names(data.used)[1:2] <- c("obs", "delta.cens")
surv.tree = rpart(Surv(obs,delta.cens)~., data = data.used, minbucket = 30)

# Getting the Survival Curves. 
pred.surv.tree <- predict(surv.tree, proximity = FALSE)

# Finding the terminal nodes
sett=unique(surv.tree[['where']])
nset=length(sett)

cens.est = matrix(0, nrow = num, ncol = num)
obs.used = rep(NA, num)
delta.used = rep(NA, num)

for (i in 1:nset){
# Finding the subset corresponding the ith node of the tree
subset=(1:nu)[surv.tree[['where']]==sett[i]]
nlen=length(subset)
# Observed time within the node
sobs=obs[subset]
# Failure indicators within each node.
sdelta=delta[subset]

# Doing truncation within that subset
# Changing the dataset to account for truncation. 
aa=datach(sobs,sdelta,dtype = dtype)
# Observed time after truncation
sobs=aa[[1]]
# Failure indicator after truncation.
sdelta = aa[[2]]

obs.used[subset] = sobs
delta.used[subset] = sdelta

# Calculating the KM estimator censoring curve within a node
# Calculating the jumps in the KM estimator
hazC=mapply(function(xx,dd){dd/sum((sobs>=xx))},xx=sobs,dd=1-sdelta)
surC_km=mapply(function(xx){prod(1-hazC[sobs<=xx])},xx=obs)
cens.est[subset, ] = matrix(surC_km,nrow=length(subset),ncol=length(surC_km),byrow=TRUE)
}

return(list(cens.est, obs.used, delta.used, surv.tree[['where']]))
}


# external calculates the parameter vector using the L_2 loss and no simulation from martingale.
# Inputs: obs = observed time, delta = failure indicator, x1-x5 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
externalRegularTreeKM = function(obs,delta,xx, mtype,dtype)
{
  n = length(obs)
  nu = length(obs)

  # Calculating the conditional expectation    
  m1 = mfunc(obs,delta,xx, mtype)

  # Calculating the conditional censoring distribution.
  tem=GfuncSurvivalTree(obs,delta,dtype,xx)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]
  # Observed event times for adjusted for truncation
  obs=tem[[2]]
  # Failure indicator adjusted for truncation
  delta=tem[[3]]
  # Finding which observations fall in which terminal node
  term.node = tem[[4]]

  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta/diag(surC_rf)
  a1=a0*log(obs)
  
  b0=(1-delta)/diag(surC_rf)
  b1=b0 * diag(m1)

c0 <- rep(NA, n)
c1 <- rep(NA, n)

# Creating the ordered data
ord.used = order(obs)
obs.order = obs[ord.used]
delta.order = delta[ord.used]

# Finding the terminal nodes
sett=unique(term.node)
nset=length(sett)

for (i in 1:nset){
# Finding the subset corresponding the ith node of the tree
subset=(1:nu)[term.node==sett[i]]
nlen=length(subset)
# Observed time within the node
sobs=obs[subset]
# Failure indicators within each node.
sdelta=delta[subset]


  kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
  c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
  c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*m1[subset,i]/kk)},tt=sobs,i=1:nlen)
}

parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),nu)
  
  return(parms)
}


# external calculates the parameter vector using the L_2 loss and no simulation from martingale.
# Inputs: obs = observed time, delta = failure indicator, x1-x5 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
externalBrierTreeKM = function(obs,delta,xx, mtype,dtype, time.point)
{
  n = length(obs)
  nu = length(obs)

  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  n = length(obs)
  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  # number of observations

  
  # Calculating the conditional expectation    
  m1=mfuncBrier(obs,delta,xx, mtype, time.point)

  # Calculating the conditional censoring distribution.
  tem = GfuncSurvivalTree(obs.t,delta.t,dtype, xx)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]
  # Finding which observations fall in which terminal node
  term.node = tem[[4]]

  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta.t/diag(surC_rf)
  a1 = a0 * (obs > time.point)
  
  b0=(1-delta.t)/diag(surC_rf)
  b1=b0 * diag(m1)

c0 <- rep(NA, n)
c1 <- rep(NA, n)

# Finding the terminal nodes
sett=unique(term.node)
nset=length(sett)

for (i in 1:nset){
# Finding the subset corresponding the ith node of the tree
subset=(1:nu)[term.node==sett[i]]
nlen=length(subset)
# Observed time within the node
sobs=obs.t[subset]
# Failure indicators within each node.
sdelta=delta.t[subset]


  kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
  c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
  c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*m1[subset,i]/kk)},tt=sobs,i=1:nlen)
}

parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),nu)
  
  return(parms)
}


# external calculates the parameter vector using the L_2 loss and no simulation from martingale.
# Inputs: obs = observed time, delta = failure indicator, x1-x5 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
externalBrierKM = function(obs,delta,xx, mtype,dtype, time.point)
{
  n = length(obs)
  nu = length(obs)

  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  n = length(obs)
  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  # number of observations

  
  # Calculating the conditional expectation    
  m1=mfuncBrier(obs,delta,xx, mtype, time.point)


  tem = GfuncKM(obs.t,delta.t,dtype)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]


  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta.t/surC_rf
  a1 = a0 * (obs > time.point)
  
  b0=(1-delta.t)/surC_rf
  b1=b0 * diag(m1)

c0 <- rep(NA, n)
c1 <- rep(NA, n)

  sobs=obs.t
  kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
  c0=mapply(function(tt){sum(b0*(sobs<=tt)/kk)},tt=obs)
  c1=mapply(function(tt,i){sum(b0*(sobs<=tt)*m1[,i]/kk)},tt=sobs,i=1:n)

parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),nu)
  
  return(parms)
}


# Calculating the model for the Conditional Expectations using random forests
random.forest = function(obs,delta,xx)
{
nu = length(obs)

# Fitting the Cox model
# Creating the data frame
data.used <- data.frame(obs, delta, xx)
names(data.used)[1:2] <- c("obs", "delta")
rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none")

# Getting the Survival Curves.
pred.rf <- predict(rand.for, proximity = FALSE)

# Calculation of the survivor function
m1=matrix(0,nu, nu)

# Finding unique event times
time.used <- pred.rf$time.interest

# Finding the jumps in the estimator \hat P(T >t|W_i)
surv.diff <- matrix(0, ncol = sum(delta), nrow = nu)

for(i in 1:nu){
# Calculating the jumps in the random forest model survival curve estimator
surv.diff[i, ] <- c(1, pred.rf$survival[i, ][-length(pred.rf$survival[i, ])]) - pred.rf$survival[i, ]
}

for(j in 1:nu)
{
if(delta[j]==FALSE){

for(i in 1:nu)
{
if(obs[j]<=obs[i]){

if(sum(surv.diff[i, ][time.used > obs[j]]) != 0){
# Calculating the conditional expectation
m1[j,i]=  sum(log(time.used[time.used > obs[j]]) * surv.diff[i, ][time.used > obs[j]])/sum(surv.diff[i, ][time.used > obs[j]])
}
}
}
if (sum(surv.diff[i, ][time.used > obs[j]]) == 0){
m1[j,]=log(obs[j])}
}
}

return(m1)
}


# Calculating the model for the Conditional Expectations using the random forest algorithm.
randomForestBrier = function(obs,delta,xx, time.point){


n = length(obs)
nu = length(obs)

# Fitting the Cox model
# Creating the data frame
data.used <- data.frame(obs, delta, xx)
names(data.used)[1:2] <- c("obs", "delta")
rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none")
  
# Getting the Survival Curves. 
pred.rf <- predict(rand.for, proximity = FALSE)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
#P(T > tau|T >u,W) = P(T>tau|W)/P(T>u|W)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvRF = function(time, cov.index){
time.point = sum(pred.rf$time.interest < time) + 1
surv = c(1, pred.rf$survival[cov.index, ])[time.point]
return(surv)
}

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(0, ncol = n, nrow = n)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
for(i in 1:n){
# Calculate the denominator and the numerator in the desired probability
prob.est.den = predsurvRF(obs[j], i)
prob.est.num = predsurvRF(time.point, i)
if(prob.est.den != 0){
m1[j, i] = prob.est.num/prob.est.den
}
if(prob.est.den == 0){
  m1[j, i] = 0.5
}

}
}
}

# Return the probabilities
return(m1)
}

# external calculates the parameter vector.
# Inputs: obs = observed time, delta = failure indicator, x1-x25 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
BJexternalBrier=function(obs,delta,xx, mtype, dtype, time.point)
{
  n = length(obs)
  nu = n
  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  # number of observations
  
  # Calculating the conditional expectation    
  m1 = mfuncBrier(obs,delta,xx, mtype, time.point)

  a1 = delta.t *  (obs > time.point) + (1 - delta.t) * diag(m1)
  
  parms = a1
  
return(parms)
}

# external calculates the parameter vector.
# Inputs: obs = observed time, delta = failure indicator, x1-x25 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
BJexternal=function(obs,delta,xx, mtype,dtype, time.point)
{
  n = length(obs)
  nu = n
  # Creating the new T(t) dataset
  
  # Calculating the conditional expectation    
  m1 =   m1=mfunc(obs,delta,xx, mtype)

  a1 = delta *  log(obs) + (1 - delta) * diag(m1)
  
  parms = a1
  
return(parms)
}

# Calculating the Tree Based conditional expectation
# Input: obs = observed time, delta = failure indicator, x1, ldots, x25 = vector of covariates
# mtype = what type of cond exp is being calculated. 
# Output: Matrix of conditional expectations m_{1i}(tilde T_j) for j, i = 1, ldots ,n
# time point is the t in the equation m1[j, i] = P(T >t|T >T_j, W_i)
mfuncBrier=function(obs,delta,xx, mtype, time.point)
{
  num=length(obs)

  if(mtype=="rand.for")
  {
    m1 = randomForestBrier(obs,delta,xx, time.point)
  }
 
  return(m1)  
}


####################################################################################################
###### The following functions calcualate the estimator for P(T >t|W)
#####################################################################################################

# Calculating the Tree Based conditional expectation
# Input: obs = observed time, delta = failure indicator, x1, ldots, x25 = vector of covariates
# xx = matrix of covariates, mtype = what type of cond exp is being calculated. 
# Output: Matrix of conditional expectations m_{1i}(tilde T_j) for j, i = 1, ldots ,n
mfunc=function(obs,delta,xx, mtype)
{

  if (mtype=="rand.for")
  {
    m1=random.forest(obs,delta,xx)
  }
  return(m1)
}

# Changes the dataset to account for truncation
# Input: obs = observed time, delta = failure indicator, dtype = how truncation is done.
# b is method 2 and a is method 1 and the number afterwards indicates what the truncation level is.
# Output: obs = observed failure time adjusted for truncation, delta = failure indicator adjusted for truncaiton.
datach=function(obs,delta,dtype, tau1 = NULL)
{
  nu = length(obs)
  if(dtype=="b0")
  {
    delta[obs==max(obs)]=TRUE
  }
  
  
  
  if(dtype=="b5")
  {
    delta[order(obs)][floor(nu-0.05*nu):nu]=TRUE
  }
  
  
  if(dtype=="b10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=TRUE
  }
  
  
  if(dtype=="b15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=TRUE
  }
  
  if(dtype=="a5")
  {    
    delta[order(obs)][floor(nu-0.05*nu):nu]=T
    obs[order(obs)][floor(nu-0.05*nu):nu]=obs[order(obs)][floor(nu-0.05*nu)]    
  } 
  if(dtype=="a10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=T
    obs[order(obs)][floor(nu-0.10*nu):nu]=obs[order(obs)][floor(nu-0.10*nu)]    
  } 
  
  
  if(dtype=="a15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=T
    obs[order(obs)][floor(nu-0.15*nu):nu]=obs[order(obs)][floor(nu-0.15*nu)]    
  }   
  
  if(dtype=="none")
  {
    obs=obs
    delta=delta
  }   
  if(dtype=="c")
  {
    obs=pmin(obs, 4.75)
    delta=delta
  }   
  return(list(obs,delta)) 
} 

# Create imputation for doubly robust and l2 loss
create.imp.dr.l2 = function(obs,delta,cov.vec,mtype,dtype){
n = length(obs)

  parms = externalRegularTreeKM(obs=obs,delta=delta,cov.vec, mtype=mtype, dtype=dtype)

  a1 <- parms[(n+1):(2*n)]
  b1 <- parms[(3*n+1):(4*n)]
  c1 <- parms[(5*n+1):(6*n)]
  y.imp.all <- a1 + b1 - c1

return(y.imp.all)
}

# Create imputation for buckley James and l2 loss
create.imp.bj.l2 = function(obs,delta,cov.vec,mtype,dtype){

  y.imp.all = BJexternal(obs=obs,delta=delta,cov.vec, mtype=mtype, dtype=dtype)

return(y.imp.all)
}

# Create imputation for dr and brier loss
create.imp.dr.brier = function(obs,delta,cov.vec,mtype,dtype, time.point){
   n = length(obs)
   parms = externalBrierTreeKM(obs,delta,cov.vec, mtype = mtype,dtype = dtype, time.point)
   a1 <- parms[(n+1):(2*n)]
   b1 <- parms[(3*n+1):(4*n)]
   c1 <- parms[(5*n+1):(6*n)]
   y.brier <- a1 + b1 - c1
return(y.brier)
}

# Create imputation for dr and brier loss
create.imp.dr.brier.km = function(obs,delta,cov.vec,mtype,dtype, time.point){
   n = length(obs)
   parms = externalBrierKM(obs,delta,cov.vec, mtype = mtype,dtype = dtype, time.point)
   a1 <- parms[(n+1):(2*n)]
   b1 <- parms[(3*n+1):(4*n)]
   c1 <- parms[(5*n+1):(6*n)]
   y.brier <- a1 + b1 - c1
return(y.brier)
}


# Create imputation for dr and brier loss
create.imp.bj.brier = function(obs,delta,cov.vec,mtype,dtype, time.point){
   y.brier = BJexternalBrier(obs,delta,cov.vec, mtype = mtype,dtype = dtype, time.point)

return(y.brier)
}

# Create imputation for dr and brier loss
create.ipcw.weights = function(obs,delta,cov.vec,dtype, time.point){

    # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  n = length(obs)
  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  # number of observations

  # Calculating the conditional censoring distribution.
  tem = GfuncSurvivalTree(obs.t,delta.t,dtype, cov.vec)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]
  # Finding which observations fall in which terminal node
  term.node = tem[[4]]

  # Calculating a0, a1, b0, b1, c0, c1
  ipcw.weights=delta.t/diag(surC_rf)
  data.ipcw = data.frame(obs.t, delta.t, ipcw.weights) 
return(data.ipcw)
}

