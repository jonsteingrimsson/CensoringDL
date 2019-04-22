library(survival)
library(randomForestSRC)
library(MASS)
library(penalized)

data(nki70)
data.nki70 = nki70
data.nki70$time = data.nki70$time + runif(nrow(data.nki70), 0, 0.001)
data.used = data.nki70

source("CreateImputations.R")
load("Tau1.rda")

create.data.metabric = function(){

# Choose level of v-fold cross validation 
v <- 5
# Creating the cross validated group. 
get.sam <- rep(c(1:v),sum(data.used$event==TRUE)/v)
if(length(get.sam)<sum(data.used$event==TRUE)){
get.sam <- c(get.sam,c(1:(sum(data.used$event==TRUE)-length(get.sam))))
}
delta.1 <- sample(get.sam,sum(data.used$event==TRUE),replace=F)
delta.0 <- sample(rep(c(1:v),length(data.used$event)-length(delta.1)),sum(data.used$event==FALSE),replace=F) 
k <- 1
l <- 1
val.sample <- NULL
for(m in 1:length(data.used$event)){
if(data.used$event[m]==0){
val.sample[m] <- delta.0[k]
k <- k+1
}
if(data.used$event[m]==1){
val.sample[m] <- delta.1[l]
l <- l+1
}
}

  for(j in 1:5){
  val.sample = sample(1:5, dim(data.used)[1], replace = TRUE)
  data.train = data.used[val.sample != 1, ]
  obs.train = data.train$time
  delta.train = data.train$event
  cov.train = data.train[, -c(1,2)]
  data.test = data.used[val.sample == 1, ]
  obs.test = data.test$time
  delta.test = data.test$event
  cov.test = data.test[, -c(1,2)]
    
  # Create imputation
  imp.val.dr = create.imp.dr.brier.km(obs = obs.train, delta = delta.train, cov.vec = cov.train, mtype = "rand.for", dtype = "b10", time.point = time.point)
  imp.val.bj = create.imp.bj.brier(obs = obs.train, delta = delta.train, cov.vec = cov.train, mtype = "rand.for", dtype = "none", time.point = time.point)
  data.ipcw = create.ipcw.weights(obs = data.used$time, delta = data.used$event, cov.vec = data.used[, -c(1,2)], dtype = "b10", time.point = time.point)
  data.ipcw = data.ipcw[val.sample == 1, 3]

  x_train = cov.train
  y_train.dr = imp.val.dr
  y_train.bj = imp.val.bj
  x_test = cov.test
  y_test = obs.test

  data.surv = data.frame(data.used$time, data.used$event, data.used[, -c(1,2)])[val.sample != 1, ]
  names(data.surv)[1:2] = c("obs", "delta")

file.name = paste("x_train", toString(job.number), "CV", toString(j),".rda", sep = "")
save(x_train, file = file.name)
file.name = paste("y_trainDR", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(y_train.dr, file = file.name)
file.name = paste("y_trainBJ", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(y_train.bj, file = file.name)
file.name = paste("x_test", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(x_test, file = file.name)
file.name = paste("y_test", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(y_test, file = file.name)
file.name = paste("data_surv", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(data.surv, file = file.name)
file.name = paste("ipcw_weights", toString(job.number), "CV", toString(j), ".rda", sep = "")
save(data.ipcw, file = file.name)
save(time.point, file = "Tau1.rda")
}
}


create.data.metabric()
