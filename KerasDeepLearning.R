library(keras)

# One simulation
one.sim.dl = function(){
  pred.err.dl = rep(NA, 5)
  pred.err.dl.cv = rep(NA, 5)
  reg.par.final = rep(NA, 5)

  for(j in 1:5){
  file.name = paste("x_train", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  file.name = paste("y_trainDR", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  file.name = paste("x_test", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  file.name = paste("y_test", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  file.name = paste("data_surv", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  file.name = paste("ipcw_weights", toString(job.number), "CV", toString(j), ".rda", sep = "")
  load(file.name)
  load("Tau1.rda")


  # Start by fitting the keras model
   train_data <- x_train
   train_targets <- y_train.dr
   test_data = x_test

# Creating dummy variables for data 
train_data.dummy = model.matrix( ~ ., data=train_data )
test_data.dummy = model.matrix( ~ ., data=test_data)

# Setting up architecture
  build_model <- function() {
    model <- keras_model_sequential() %>% 
      layer_dense(units = 15, activation = "relu",
                  input_shape = dim(train_data.dummy)[2]) %>% 
      layer_dropout(rate = 0.2) %>%
  #    layer_dense(units = 10, activation = "relu") %>% 
  #        layer_dropout(rate = 0.2) %>%
      layer_dense(units = 1, activation = "sigmoid") 
      
    model %>% compile(
      optimizer = "rmsprop", 
      loss = "mean_squared_error", 
      metrics = c("mean_squared_error")
    )
  }

    # Build the Keras model (already compiled)
    model <- build_model()
    
    batch.size = 32
    # Train the model (in silent mode, verbose=0)
    model %>% fit(train_data.dummy, train_targets,
                  epochs = 100, batch_size = batch.size, verbose = 0, optimizer='rmsprop')
                  
    predict = model %>% predict(test_data.dummy, batch_size=batch.size, verbose=0)
    y.brier = (y_test > time.point)

    # Calculating IPCW Brier loss
    pred.err.dl[j] = mean(data.ipcw * (y.brier - predict)^2)

  # Creating cross validation error for 
  build_model.cv <- function(reg.par) {
    model <- keras_model_sequential() %>% 
      layer_dense(units = 15, activation = "relu",
                  kernel_regularizer = regularizer_l2(reg.par), input_shape = dim(train_data.dummy)[[2]]) %>% 
      layer_dropout(rate = 0.2) %>%
  #    layer_dense(units = 10, activation = "relu") %>% 
  #        layer_dropout(rate = 0.2) %>%
      layer_dense(units = 1, activation = "sigmoid",
                  kernel_regularizer = regularizer_l2(reg.par)) 
      
    model %>% compile(
      optimizer = "rmsprop", 
      loss = "mean_squared_error", 
      metrics = c("mean_squared_error")
    )
  }

# Set \eta parameter
  reg.par.used = c(0, 0.001, 0.01, 0.1)
  cross.val.error = rep(NA, length(reg.par.used))

  # Loop over possible penalization parameters
  for(l in 1:length(reg.par.used)){
  # Number of cross validations
  k <- 5
  indices <- sample(1:nrow(train_data))
  folds <- cut(indices, breaks = k, labels = FALSE)
  # Number of iterations
  num_epochs <- 100
  all_scores <- c()
  for (i in 1:k) {
  val_indices <- which(folds == i, arr.ind = TRUE)
  val_data <- train_data.dummy[val_indices,]
  val_targets <- train_targets[val_indices]
  partial_train_data <- train_data.dummy[-val_indices,]
  partial_train_targets <- train_targets[-val_indices]
  model <- build_model.cv(reg.par.used[l])
  model %>% fit(partial_train_data, partial_train_targets,
  epochs = num_epochs, batch_size = batch.size, verbose = 0, optimizer='rmsprop')
  results <- model %>% evaluate(val_data, val_targets, verbose = 0)
  all_scores <- c(all_scores, results$mean_squared_error)
  }
  cross.val.error[l] = mean(all_scores)
  }

  # Fitting final model
  reg.par.final[j] = reg.par.used[which.min(cross.val.error)]

    # Build the Keras model (already compiled)
    model <- build_model.cv(reg.par.final[j])
    
    # Train the model (in silent mode, verbose=0)
    model %>% fit(train_data.dummy, train_targets,
                  epochs = 100, batch_size = batch.size, verbose = 0, optimizer='rmsprop')
                  
    predict.cv = model %>% predict(test_data.dummy, batch_size=batch.size, verbose=0)

    # Truncate to fall in interval
    pred.err.dl.cv[j] = mean(data.ipcw * (y.brier - predict.cv)^2)
}

return(list(pred.err.dl.cv = mean(pred.err.dl.cv), reg.par.final = reg.par.final, pred.err.dl = mean(pred.err.dl)))
}

dl.sim.dr = one.sim.dl()

