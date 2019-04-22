n.sim = 1
err.dl.dr = rep(NA, n.sim)
err.dl.dr.cv = rep(NA, n.sim)
err.dl.bj = rep(NA, n.sim)
err.dl.bj.cv = rep(NA, n.sim)
err.rf = rep(NA, n.sim)
err.cox = rep(NA, n.sim)
err.cox.pen = rep(NA, n.sim)
time.used = rep(NA, n.sim)

time.point = 10
save(time.point, file = "Tau1.rda")
# Creating data
source("CreateData.R")
# Implementing competing methods
source("CompetingMethods.R")
# Implementing the CUDL algorithms
source("KerasDeepLearning.R")
source("KerasDeepLearningBJ.R")

performance = list(deep.learning.dr = dl.sim.dr, deep.learning.bj = dl.sim.bj, competing.methods = other.methods)
err.dl.dr.cv = performance$deep.learning.dr$pred.err.dl
err.dl.dr.cv = performance$deep.learning.dr$pred.err.dl.cv
err.dl.bj = performance$deep.learning.bj$pred.err.dl
err.dl.bj.cv = performance$deep.learning.bj$pred.err.dl.cv
err.rf = performance$competing.methods[[1]][3]
err.cox = performance$competing.methods[[1]][2]
err.cox.pen = performance$competing.methods[[1]][1]
