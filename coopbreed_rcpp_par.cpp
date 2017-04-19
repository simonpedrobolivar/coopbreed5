#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <iostream>
#include <progress.hpp>
#include <list>
#include <vector>
#include <string> 
#include <sstream>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

// function for casting integer to string
#define SSTR( x ) static_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()

// [[Rcpp::export]]
NumericVector Sb(NumericVector d, double k){ // breeder survival function
  int n;
  n = d.size();
  NumericVector res(n);
  for(int i = 0; i < n; i++){
    res(i) = 1 - exp(-k * (1 - d(i)));
  }
  return res;
}

// [[Rcpp::export]]
double mean_rcpp(NumericVector x){ // mean function
  int n = x.size();
  //Size of vector
  double sum = 0;
  //Sum value
  //For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    //Shorthand for sum=sum+
    sum += x[i];
  }
  return sum/n;
  //Obtain and return the Mean
}




struct coopbreed3 : public Worker{
  // input matrix to read from
  const RMatrix<double> pars;
  
  // output matrix to write to
  RMatrix<double> results;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  coopbreed3(const NumericMatrix pars, NumericMatrix results)
    : pars(pars), results(results) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    
    
    //int n_simulations = n_mates_vec.size() * n_off_vec.size() * par_c_vec.size() * par_k_vec.size();
    //vector<string> names;
    
    for(std::size_t par = begin; par < end; par++){
      int n_mates = pars(par, 0);
      int n_off = pars(par, 1);
      int par_c = pars(par, 2);
      int par_k = pars(par, 3);
      int paths =  pars(par, 4);
      int n_gener=  pars(par, 5);
      int n_patches =  pars(par, 6);
      double MutStep =  pars(par, 7);
      
      // create string for naming the model outcome
      //  names.push_back(SSTR("n_mates = " << n_mates_vec(mates) << "/ n_off = " << n_off_vec(off) 
      //                                  << "/ par_c = " << par_c_vec(pc) << "/ par_k = " << par_k_vec(pk)));
      
      
      NumericVector n_mates_v(n_mates);
      for(std::size_t n = 0; n < n_mates; n++){
        n_mates_v(n) = n;
      }
      
      NumericVector n_patches_v(n_patches);
      for(std::size_t n = 0; n < n_patches;n++){
        n_patches_v(n) = n;
      }
      NumericVector n_off_v(n_off);
      for(std::size_t n = 0; n < n_off;n++){
        n_off_v(n) = n;
      }
      
      NumericMatrix Data(paths,n_gener); // matrix for storing the results of one model run
      
      for (std::size_t path = 0 ; path < paths; path++) {
        // Initialize Population
        NumericMatrix Popn(2,n_patches);
        Popn(0,_) = runif(n_patches);
        Popn(1,_) = runif(n_patches);
        //cout << "popn: "; cout << Popn; cout << " "; cout << ", ";
        //one individual per patch 
        //two chromosomes per individual 
        
        //for (int i = 0; i < n_patches; i++){
        //  Popn(0,i) = R::runif(0,1);
        //  Popn(1,i) = R::runif(0,1);
        //
        //}
        
        for(std::size_t t = 0; t < n_gener; t++){ // for all generations
          NumericMatrix mates(n_mates, n_patches);// stores index of each mate of each individual
          NumericVector sample1(n_mates);
          for(std::size_t j = 0; j < n_patches; j++){
            sample1 = RcppArmadillo::sample(n_patches_v, n_mates, TRUE);  
            for(std::size_t jj = 0; jj < n_mates; jj++){
              mates(jj,j) = sample1(jj);
            }
          }
          
          double Offspring[n_off][2][n_patches];// stores the genotypes of each off Offspring
          NumericVector vec = NumericVector::create(0,1);
          
          for(std::size_t i = 0; i < n_patches; i++){
            NumericVector Fathers(n_off);
            Fathers = RcppArmadillo::sample(n_mates_v, n_off, TRUE);
            NumericVector PatContribs(n_off);
            PatContribs = RcppArmadillo::sample(vec, n_off, TRUE);
            NumericVector MatContribs(n_off);
            MatContribs = RcppArmadillo::sample(vec, n_off, TRUE);
            
            for(std::size_t j = 0; j < n_off; j++){ // for all offspring
              Offspring[j][0][i] = Popn(PatContribs(j), mates(Fathers(j), i));
              Offspring[j][1][i] = Popn(MatContribs(j),i);
            }
          }
          
          NumericMatrix OffspringPhenotype(n_off, n_patches);
          NumericVector AvgPhenotype(n_patches); // average phenotype per patch
          NumericVector vec2(2);
          
          for(std::size_t i = 0; i < n_off; i++){
            for(std::size_t j = 0; j < n_patches; j++){
              vec2(0) = Offspring[i][0][j];
              vec2(1) = Offspring[i][1][j];
              OffspringPhenotype(i,j) = mean_rcpp(vec2); // mean of maternal and paternal genes
            }
          }
          for(std::size_t i = 0; i < n_patches; i++){
            AvgPhenotype(i) = mean_rcpp(OffspringPhenotype(_,i)); 
          }
          double GlobalAvgPhenotype = mean_rcpp(AvgPhenotype);
          NumericVector ProbLocal(n_patches); 
          for(std::size_t n = 0; n < n_patches; n++){
            ProbLocal(n) = (1 - AvgPhenotype(n))/(par_c * GlobalAvgPhenotype + (1 - AvgPhenotype(n)));//check if results is rigth
          }
          NumericVector BreederSurvival = Sb(AvgPhenotype, par_k);
          NumericMatrix PopnNext(2, n_patches);
          
          
          for(std::size_t i = 0; i < n_patches; i++){ // for all patches
            double X;
            X = R::runif(0,1);
            if(X > BreederSurvival(i)){ // Breeder Dies
              double Y;
              Y = R::runif(0,1);
              int PatchWinner;
              int Winner;
              if(Y < ProbLocal(i)){ 
                PatchWinner = i;//Winner Comes from Local Patch
                NumericVector prob1(n_off);// Local Offspring Compete
                for(std::size_t n=0; n < n_off; n++){
                  prob1(n) = 1.0 - OffspringPhenotype(n,PatchWinner);
                }
                Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob1)(0);
                //cout << "Winner: "; cout << Winner; cout << ", ";
              }else{
                // Patches Compete
                PatchWinner = RcppArmadillo::sample(n_patches_v, 1, TRUE, AvgPhenotype)(0); 
                NumericVector prob2(n_off);
                for(std::size_t n=0; n < n_off; n++){// Offspring On Winning Patch Compete
                  prob2(n) = OffspringPhenotype(n,PatchWinner);
                }
                Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob2)(0);
              }
              
              NumericVector WinnersGenes(2);
              for(std::size_t n = 0; n < 2; n++){// Mutate Winner's Genes
                WinnersGenes(n) = Offspring[Winner][n][PatchWinner] + (MutStep * R::rnorm(0,1));
                if(WinnersGenes(n) <= 0){
                  WinnersGenes(n) = MutStep;
                  //cout << "winnersgenes: "; cout << WinnersGenes(n);
                }
                if(WinnersGenes(n) >= 1.0){
                  WinnersGenes(n) = 1.0 - MutStep;
                  
                }
              }
              PopnNext(_,i) = WinnersGenes;// Update Next Gen Population Array
            }else{// Breeder Survives
              PopnNext(_,i) = Popn(_,i);
            }
          }
          Popn = PopnNext;
          // Record Data
          NumericVector mean_perind(n_patches);
          for(std::size_t n = 0; n < n_patches; n++){
            mean_perind(n) = mean_rcpp(Popn(_,n));
          }
          Data(path,t) = mean_rcpp(mean_perind);
        }
      }  
      //MeanPath = apply(Data, 2, mean)
      //  SE=sqrt(apply(Data, 2, var)/paths)
      
      //return(list("Data" = Data, "Mean" = MeanPath, "SE" = SE))
      //return Data;
      for(std::size_t i = 0; i < paths; i++){
        for(std::size_t ii = 0; ii < n_gener; ii++){
          results(par * paths + i, ii) = Data(i, ii);
        }
      }
      
      //results.push_back(Data);
      //results.attr("names") = names;
    }
  }
};


// [[Rcpp::export]]
NumericMatrix coopbreed_parallel(NumericMatrix pars){
  // allocate the matrix we will return
  NumericMatrix results(pars.nrow() * pars(1, 4), pars(1, 5));
  //cout << pars.nrow()* pars(1, 4); cout << ", "; cout << pars(1, 5);
  // create the worker
  coopbreed3 coopbreed3(pars, results);
  cout << "here i am first";
  // call it with parallelFor
  parallelFor(0, pars.nrow(), coopbreed3);
  cout << "here i am";
  return results;
}







