library(foreach)
library(doParallel)

coopbreedR <- function(Paths, T, P, MutStep, M, K, c, k){
  # breeder survival
  Sb <- function(d){
    1 - exp(-k * (1 - d))
  }

  Data_vec = numeric(length = T)
  # use all cores available, except 1
  cores <- parallel::detectCores() - 1

  # set up cluster
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  Data <- foreach::foreach(path = 1:Paths, .combine = cbind) %dopar%{
  # for all paths

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
      Data_vec[t] = mean(apply(Popn, 1, mean))
    }
    Data_vec
  }
  #MeanPath = apply(Data, 2, mean)
  #SE=sqrt(apply(Data, 2, var)/Paths)
  return(Data)
  #return(list("Data" = Data, "Mean" = MeanPath, "SE" = SE))
}

