# load required packages
library(foreach)
library(doParallel)



#####################
#### Parameters ################################################################
#####################

Paths = 5 # number of Paths to simulate
T = 3000 # number of generations to simulate
P = 1000 # number of patches
MutStep = 0.001 # mutation step
M = 2 # number of mates
K = 100 # number of offspring
c = 1.5 # parameter
k = 1 # parameter


set.seed(632)

#############################
#### The model ####################################################################
#############################

coopbreed_fun <- function(Paths, T, P, MutStep, M, K, c, k){
  # breeder survival
  Sb <- function(d){
    1 - exp(-k * (1 - d))
  }
  Data = matrix(nrow = Paths, ncol = T)
  
  for(path in 1:Paths){ # for all paths
    #print(paste("Path: ", path))
    # Initialize Population
    Popn = matrix(runif(P*2), nrow = 2, ncol = P)
    #Popn[,1:10]
    # one individual per patch 
    # two chromosomes per individual
    
    for(t in 1:T){ # for all generations
      #print(paste("Generation: ", t))
      mates = matrix(sample(1:P, M * P, replace = T), nrow = M, ncol = P)
      #dim(mates)
      #mates[,1:10]
      # stores index of each mate of each individual
      #Offspring <- list()
      #Offspring[[1]] <- Popn[sample(2, K, replace = T), mates[sample(M,K,replace = T), 1:P]]
      #Offspring[[2]] <- Popn[sample(2, K, replace = T),1:P]
      
      
      Offspring = array(dim = c(K,2, P)) 
      #mates[Fathers,][1:10,1:10]
      #mates[,1:10]
      fun1 <- function(x) return(x[sample(M,1,replace = T)])
      # Genotype 1: Choose one chromosome from Father. apply: choose "Who is the father from the mates in each patch?"
      Offspring[,1,] <- Popn[sample(2, K, replace = T), apply(mates, 2, FUN = fun1)]
      # choose chromose from mother:
      Offspring[,2,] <- Popn[sample(2, K, replace = T),]
      
      # stores the genotypes of each off Offspring
      
      
      OffspringPhenotype <- array(apply(Offspring,c(1,3), mean), c(K, P)) # mean jedes einzelnen nachkommkens (in each patch)
      #length(apply(Offspring,c(1,3), mean))
      #mean(c(Offspring[1,1,1], Offspring[1,2,1]))
      #apply(Offspring,c(1,3), mean)[1]
      
      #dim(OffspringPhenotype)
      AvgPhenotype = apply(OffspringPhenotype, 2, mean)  # average phenotype per patch
      #length(AvgPhenotype)
      GlobalAvgPhenotype = mean(t(AvgPhenotype))
      
      ProbLocal=(1 - AvgPhenotype)/(c * GlobalAvgPhenotype + (1 - AvgPhenotype))
      BreederSurvival <- Sb(AvgPhenotype)
      PopnNext = matrix(nrow = 2, ncol = P)
      
      
      # TODO
      
      # start timing
      #ptm <- proc.time() 
      # use all cores available, except 1
      #cores <- parallel::detectCores() - 1
      # set up cluster
      #cl <- parallel::makeCluster(cores)
      #doParallel::registerDoParallel(cl)
      
      for(i in 1:P){
        #foreach::foreach(i = 1:P, .combine = c)%dopar%{ # for all patches
        X = runif(1)
        if(X > BreederSurvival[i]){
          # Breeder Dies 
          Y=runif(1)
          #########################################
          #TODO
          #ifelse(Y < ProbLocal[i], )
          
          if(Y < ProbLocal[i]){
            #Winner Comes from Local Patch
            PatchWinner = i
            # Local Offspring Compete
            Winner=sample(1:K,1,replace = T, prob = 1-t(OffspringPhenotype[,i]))
            ################################
          }else{
            # Patches Compete
            PatchWinner=sample(1:P, 1, replace = T, prob = AvgPhenotype)
            # Offspring On Winning Patch Compete
            Winner=sample(1:K,1,replace = T, prob = t(OffspringPhenotype[,PatchWinner]))
          }
          WinnersGenes= array(Offspring[Winner,,PatchWinner], c(1, 2))
          
          # Mutate Winner's Genes
          WinnersGenes = WinnersGenes + MutStep * rnorm(2) 
          WinnersGenes[which(WinnersGenes <= 0)] = MutStep
          WinnersGenes[which(WinnersGenes >= 1)] = 1 - MutStep
          # Update Next Gen Population Array
          PopnNext[,i] = t(WinnersGenes)
        }else{
          # Breeder Survives
          PopnNext[,i] = Popn[,i]
        }
      }
      
      Popn = PopnNext
      # Record Data
      Data[path,t] = mean(apply(Popn, 1, mean))
    }
  }  
  MeanPath = apply(Data, 2, mean)
  SE=sqrt(apply(Data, 2, var)/Paths)
  
  return(list("Data" = Data, "Mean" = MeanPath, "SE" = SE))
}




system.time({
  res <- coopbreed_fun(Paths = 1, T = 0.1 *T, P = 0.1 *P, MutStep = MutStep, M = M, K = K, c= c, k= k)
})
coopbreed_fun2 <- function(x) return(coopbreed_fun(Paths = 1, T = T, P = P, MutStep = MutStep, M = M, K = K, c= c, k= k))
Data2 = matrix(nrow = Paths, ncol = T)

apply(Data2, 1, coopbreed_fun2)

#############################
#### Plot results ################################################################
#############################
plot(1:(0.1*T), 1-res$Data[1,], type = "l", ylim = c(0, 1), col = "grey")
for(i in 1:nrow(res$Data)) lines(1:(0.1*T), 1-res$Data[i,], col = "grey")

lines(1:(0.1*T), 1-res$Mean, col = "red")

dim(Data)



