library(foreach)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(devtools)

##### Functions #####

traceplot <- function(modeloutput, par.mfrow = c(1,1)){
  par(mfrow = par.mfrow, mar = c(1,1,1,1), oma = c(2,2,2,2))
  mean_fun <- function(x) apply(x, 2, mean)
  SE_fun <- function(x) sqrt(apply(x, 2, var)/Paths)
  results <- modeloutput$results
  mean = lapply(results, mean_fun)
  se = lapply(results, SE_fun)
  pars <- modeloutput$pars
  plot(1:pars$n_gener[1], 1-results[[1]][1,], type = "l", ylim = c(0, 1), col = "grey", xaxt = "n",
       main = pars[1,], cex.main = 1, las = 2)
  for(i in 1:nrow(results[[1]])) lines(1:pars$n_gener[1], 1-results[[1]][i,], col = "grey")
  lines(1:pars$n_gener[1], 1-mean[[1]], col = "red")
  lines(1:pars$n_gener[1], 1-mean[[1]] + 1.96 * se[[1]], lty = 2, col = "red")
  lines(1:pars$n_gener[1], 1-mean[[1]] - 1.96 * se[[1]], lty = 2, col = "red")
  text(y = 0.8, x = pars$n_gener[1]/2, labels = pars[1,])

  if(length(results) > 1){
    for(i in 2:length(results)){
      plot(1:pars$n_gener[i], 1-results[[i]][1,], type = "l", ylim = c(0, 1), col = "grey", xaxt = "n",
           main = pars[i,], cex.main = 1, las = 2)
      for(j in 1:nrow(results[[1]])) lines(1:pars$n_gener[i], 1-results[[i]][j,], col = "grey")
      lines(1:pars$n_gener[i], 1-mean[[i]], col = "red")
      lines(1:pars$n_gener[i], 1-mean[[i]] + 1.96 * se[[i]], lty = 2, col = "red")
      lines(1:pars$n_gener[i], 1-mean[[i]] - 1.96 * se[[i]], lty = 2, col = "red")
      text(y = 0.8, x = pars$n_gener[i]/2, labels = pars[i,])
    }
  }
}

plotLeggett <- function(equi, equi_se, par_c, par_k, n_mates){
  par(mfrow = c(2,3))
  par(mar = c(1,1,1,1), oma = c(4,4,4,4))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  for(i in 1:length(par_k)){
    for(ii in 1:length(par_c)){
      plot(n_mates, equi[which(param_combi[,8] == par_k[i] & param_combi[,7] == par_c[ii])], ylim = c(0,1),
           ylab = "", xlab = "", pch = "x", type = "b",
           xaxt = "n", yaxt = "n", cex = 1.5)
      lines(n_mates, equi[which(param_combi[,8] == par_k[i] & param_combi[,7] == par_c[ii])]
            + 1.96 * equi_se[which(param_combi[,8] == par_k[i] & param_combi[,7] == par_c[ii])],
            ylim = c(0,1),
            lty = 2, cex = 1.5)
      lines(n_mates, equi[which(param_combi[,8] == par_k[i] & param_combi[,7] == par_c[ii])]
            - 1.96 * equi_se[which(param_combi[,8] == par_k[i] & param_combi[,7] == par_c[ii])],
            ylim = c(0,1),
            lty = 2, cex = 1.5)
      if(i == 1){
        mtext(paste("c = ", par_c[ii]), outer = F, line = 1)
      }
      if(i == 2){
        axis(1)
        mtext("M", side = 1, outer = F, line = 2, cex = 0.7)
      }
      if(ii == 1){
        axis(2, las = 2)
        mtext("helping, 1 - d*", side = 2, outer = F, line = 2, cex = 0.7)
      }
      if(ii == 3){
        mtext(paste("k = ", par_k[i]), side = 4, outer = F, cex = 1, line = 1)
      }
    }
  }

}



#####################
#### Parameters ################################################################
#####################

Paths = 20 # number of Paths to simulate
T = 10 # number of generations to simulate
P = 10 # number of patches
MutStep = 0.001 # mutation step
M = 2 #seq(1,9,2) # number of mates
K = 10 # number of offspring
c = 1.5 #c(0.75, 1, 1.25) # parameter
k = 1 #c(0.5, 1) # parameter

param_combi <- as.data.frame(matrix(0,
                                    nrow = length(M) * length(K) * length(c) * length(k), ncol = 8))
colnames(param_combi) <- c("paths", "n_gener","n_patches", "MutStep", "n_mates", "n_off", "par_c", "par_k")
count = 0
for(i in 1:length(M)){
  for(ii in 1:length(K)){
    for(iii in 1:length(c)){
      for(iiii in 1:length(k)){
        count <- count + 1
        param_combi[count,1] <- Paths
        param_combi[count,2] <- T
        param_combi[count,3] <- P
        param_combi[count,4] <- MutStep
        param_combi[count,5] <- M[i]
        param_combi[count,6] <- K[ii]
        param_combi[count,7] <- c[iii]
        param_combi[count,8] <- k[iiii]
      }
    }
  }
}



#####################
#### Run the model ################################################################
#####################

# start timing
ptm <- proc.time()

# use all cores available, except 1
cores <- parallel::detectCores() - 1

# set up cluster
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

# start parallelised for-loop
test1 <- foreach::foreach(i = 1:nrow(param_combi)) %dopar%{
  coopbreed::coopbreed(paths = param_combi[i,1], n_gener =  param_combi[i,2], n_patches =  param_combi[i,3],
                       MutStep = param_combi[i,4], n_mates = param_combi[i,5],
                       n_off = param_combi[i,6], par_c= param_combi[i,7], par_k= param_combi[i,8])
}
i = 1
coopbreed(paths = param_combi[i,1], n_gener =  param_combi[i,2], n_patches =  param_combi[i,3],
          MutStep = param_combi[i,4], n_mates = param_combi[i,5],
          n_off = param_combi[i,6], par_c= param_combi[i,7], par_k= param_combi[i,8])


# stop cluster, free memory form workers
parallel::stopCluster(cl = cl)

# stop timing
proc.time() - ptm

#26.85 sec















system.time({
  results <- coopbreed2(paths = Paths, n_gener =  T, n_patches =  P, MutStep = MutStep, n_mates = M,
                        n_off = K, par_c= c, par_k= k)
  save(results, file = "/home/simon/Dokumente/ST_cooperative_breeding/Results1.RData")
}) # 54224.28 sec
results <- test




# calculate mean and standart error
mean_fun <- function(x) apply(x, 2, mean)
SE_fun <- function(x) sqrt(apply(x, 2, var)/Paths)
MeanPath_list1 = lapply(test1, mean_fun)
SE_list1 = lapply(test1, SE_fun)

MeanPath_list2 = lapply(test2, mean_fun)
SE_list2 = lapply(test2, SE_fun)


#############################
#### Plot results ################################################################
#############################

# plot evolutionary trajectories of each model run


traceplot(test1, MeanPath_list1, SE_list1)
traceplot(test2, MeanPath_list2, SE_list2)




# plot from leggett et al.:
# equi = 1 - d* (mean value for last generation simulated)
equi1 <- 1-sapply(MeanPath_list1, FUN = function(x) return(x[T]))
equi_se1 <- sapply(SE_list1, FUN = function(x) return(x[T]))
equi2 <- 1-sapply(MeanPath_list2, FUN = function(x) return(x[T]))
equi_se2 <- sapply(SE_list2, FUN = function(x) return(x[T]))




plotLeggett(equi1, equi_se1, c, k, M)
plotLeggett(equi2, equi_se2, c, k, M)




install_github('simonpedrobolivar/coopbreed3','simonpedrobolivar')


#########################################
##### Schmierblatt #############################################
#########################################

######### create package
setwd("/home/simon/Dokumente/ST_cooperative_breeding/R_package")

create("my1stRpackage")
setwd("./my1stRpackage")
document()
# back to upper level
setwd("..")
# build the .tar.gz file:
build("my1stRpackage")
# install on the computer:
install("my1stRpackage")
my1stRpackage::afun


setwd("/home/simon/Dokumente/ST_cooperative_breeding/R_stuff")
Rcpp.package.skeleton("coopbreed")
setwd("./coopbreed")
Rcpp::compileAttributes()
setwd("..")
build("coopbreed")

install("coopbreed")
coopbreed::test(2,3)
coopbreed::rcpp_hello_world(2,3)
coopbreed::Sb(2,3)



setwd("/home/simon/Dokumente/ST_cooperative_breeding/test")
Rcpp.package.skeleton("testpackage")
Rcpp::compileAttributes("testpackage")
install("testpackage")
testpackage::RcppDataFrame()
testpackage::RcppDateExample()
testpackage::Sb(2,3)

#####################################
# plot ONE model run:
system.time({
  res <- coopbreed(paths = Paths, n_gener =  T, n_patches =  P, MutStep = MutStep, n_mates = M,
                   n_off = K, par_c= c, par_k= k)
})
save(res, file = "/home/simon/Dokumente/ST_cooperative_breeding/Results_m2_p15_k1.RData")

# calculate mean and standart error
MeanPath = apply(res, 2, mean)
SE=sqrt(apply(res, 2, var)/Paths)
dim(res)








unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister()












X <- runif(P)
Y <- runif(P)
PopnNext <- ifelse(X > BreederSurvival,
                   # Breeder dies
                   ifelse(Y < ProbLocal,
                          #Winner Comes from Local Patch
                          sample(1:K,1,replace = T, prob = 1-t(OffspringPhenotype[,PatchWinner])),
                          #Winner Comes from Other Patch
                          sample(1:K,1,replace = T, prob = t(OffspringPhenotype[,PatchWinner]))),
                   # Breeder survives
                   ....Todo)


myfun <- function(x){
  X = runif(1)
  if(X > BreederSurvival[x]){
    # Breeder Dies
    Y=runif(1)
    #########################################
    #TODO
    #ifelse(Y < ProbLocal[i], )

    if(Y < ProbLocal[x]){
      #Winner Comes from Local Patch
      PatchWinner=x
      # Local Offspring Compete
      Winner=sample(1:K,1,replace = T, prob = 1-t(OffspringPhenotype[,PatchWinner]))
      ################################
    }else{
      # Patches Compete
      PatchWinner=sample(1:P, 1, replace = T, prob = AvgPhenotype);
      # Offspring On Winning Patch Compete
      Winner=sample(1:K,1,replace = T, prob = t(OffspringPhenotype[,PatchWinner]))
    }
    WinnersGenes= array(Offspring[Winner, ,PatchWinner], c(1, 2))

    # Mutate Winner's Genes
    WinnersGenes = WinnersGenes + MutStep * rnorm(2)
    WinnersGenes[which(WinnersGenes <= 0)] = MutStep
    WinnersGenes[which(WinnersGenes >= 1)] = 1 - MutStep
    # Update Next Gen Population Array
    PopnNext[,x] = t(WinnersGenes)
  }else{
    # Breeder Survives
    PopnNext[,x] = Popn[,x]
  }
  return(PopnNext)

}
sapply(1:P,FUN = myfun)


