library(rbenchmark)

Paths = 4 # number of Paths to simulate
T = 300 # number of generations to simulate
P = 100 # number of patches
MutStep = 0.001 # mutation step
M = c(1,3,5) # number of mates
K = 10 # number of offspring
c = 1.25 # parameter
k = 0.5 # parameter


system.time({
  test_results <- coopbreed2(paths = Paths, n_gener =  T, n_patches =  P, MutStep = MutStep, n_mates = M, 
                             n_off = K, par_c= c, par_k= k)
  
}) # 54224.28 sec



pars_mat <- matrix(nrow = 3, ncol = 8)
for(i in 1:nrow(pars_mat)){
  pars_mat[i,1] <- M[i]
  pars_mat[i,2] <- K
  pars_mat[i,3] <- c
  pars_mat[i,4] <- k
  pars_mat[i,5] <- Paths
  pars_mat[i,6] <- T
  pars_mat[i,7] <- P
  pars_mat[i,8] <- MutStep
}

system.time({
  test_results_parallel <- coopbreed_parallel(pars_mat)
})







a <- benchmark(coopbreed_fun(Paths = Paths, T = T, P = P, MutStep = MutStep, M = M, K = K, c = c, k = k), 
               coopbreed2(paths = Paths, n_gener =  T, n_patches =  P, MutStep = MutStep, n_mates = M, 
                          n_off = K, par_c= c, par_k = k))

res <- benchmark(matrixSqrt(m),
                 parallelMatrixSqrt(m),
                 order="relative")
res[,1:4]
