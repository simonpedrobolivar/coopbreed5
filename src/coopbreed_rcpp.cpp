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

// [[Rcpp::export]]
Rcpp::List coopbreed2(int paths, int n_gener, int n_patches, double MutStep,
                      NumericVector n_mates_vec, NumericVector n_off_vec,
                      NumericVector par_c_vec, NumericVector par_k_vec){


  int n_simulations = n_mates_vec.size() * n_off_vec.size() * par_c_vec.size() * par_k_vec.size();
  vector<string> names;
  Rcpp::List Datalist2;
  Progress p(n_simulations, TRUE);
  for(int mates = 0; mates < n_mates_vec.size(); mates++){
    int n_mates = n_mates_vec(mates);
    for(int off = 0; off < n_off_vec.size(); off++){
      int n_off = n_off_vec(off);
      for(int pc = 0; pc < par_c_vec.size(); pc++){
        int par_c = par_c_vec(pc);
        for(int pk = 0; pk < par_k_vec.size(); pk++){
          int par_k = par_k_vec(pk);


          if (Progress::check_abort()) // check if interrupted by user
            return -1.0;

          // create string for naming the model outcome
          names.push_back(SSTR("n_mates = " << n_mates_vec(mates) << "/ n_off = " << n_off_vec(off)
                                            << "/ par_c = " << par_c_vec(pc) << "/ par_k = " << par_k_vec(pk)));


          NumericVector n_mates_v(n_mates);
          for(int n = 0; n < n_mates; n++){
            n_mates_v(n) = n;
          }

          NumericVector n_patches_v(n_patches);
          for(int n = 0; n < n_patches;n++){
            n_patches_v(n) = n;
          }
          NumericVector n_off_v(n_off);
          for(int n = 0; n < n_off;n++){
            n_off_v(n) = n;
          }

          NumericMatrix Data(paths,n_gener); // matrix for storing the results of one model run

          for (int path = 0 ; path < paths; path++) {
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

            for(int t = 0; t < n_gener; t++){ // for all generations
              NumericMatrix mates(n_mates, n_patches);// stores index of each mate of each individual
              NumericVector sample1(n_mates);
              for(int j = 0; j < n_patches; j++){
                sample1 = RcppArmadillo::sample(n_patches_v, n_mates, TRUE);
                for(int jj = 0; jj < n_mates; jj++){
                  mates(jj,j) = sample1(jj);
                }
              }

              double Offspring[n_off][2][n_patches];// stores the genotypes of each off Offspring
              NumericVector vec = NumericVector::create(0,1);

              for(int i = 0; i < n_patches; i++){
                NumericVector Fathers(n_off);
                Fathers = RcppArmadillo::sample(n_mates_v, n_off, TRUE);
                NumericVector PatContribs(n_off);
                PatContribs = RcppArmadillo::sample(vec, n_off, TRUE);
                NumericVector MatContribs(n_off);
                MatContribs = RcppArmadillo::sample(vec, n_off, TRUE);

                for(int j = 0; j < n_off; j++){ // for all offspring
                  Offspring[j][0][i] = Popn(PatContribs(j), mates(Fathers(j), i));
                  Offspring[j][1][i] = Popn(MatContribs(j),i);
                }
              }

              NumericMatrix OffspringPhenotype(n_off, n_patches);
              NumericVector AvgPhenotype(n_patches); // average phenotype per patch
              NumericVector vec2(2);

              for(int i = 0; i < n_off; i++){
                for(int j = 0; j < n_patches; j++){
                  vec2(0) = Offspring[i][0][j];
                  vec2(1) = Offspring[i][1][j];
                  OffspringPhenotype(i,j) = mean_rcpp(vec2); // mean of maternal and paternal genes
                }
              }
              for(int i = 0; i < n_patches; i++){
                AvgPhenotype(i) = mean_rcpp(OffspringPhenotype(_,i));
              }
              double GlobalAvgPhenotype = mean_rcpp(AvgPhenotype);
              NumericVector ProbLocal(n_patches);
              for(int n = 0; n < n_patches; n++){
                ProbLocal(n) = (1.0 - AvgPhenotype(n)) / ((par_c * GlobalAvgPhenotype)
                                                            + (1.0 - AvgPhenotype(n)));//check if results is rigth
                //cout << "problocal: "; cout << ProbLocal(n); cout << ", ";
              }

              //ProbLocal=(1 - AvgPhenotype)/(c * GlobalAvgPhenotype + (1 - AvgPhenotype))

              NumericVector BreederSurvival = Sb(AvgPhenotype, par_k);
              NumericMatrix PopnNext(2, n_patches);


              for(int i = 0; i < n_patches; i++){ // for all patches
                double X = R::runif(0,1);
                if(X > BreederSurvival(i)){ // Breeder Dies
                  double Y = R::runif(0,1);
                  int PatchWinner;
                  int Winner;
                  if(Y < ProbLocal(i)){
                    PatchWinner = i; //Winner Comes from Local Patch
                    NumericVector prob1(n_off);// Local Offspring Compete
                    for(int n=0; n < n_off; n++){
                      prob1(n) = 1.0 - OffspringPhenotype(n,PatchWinner);
                    }
                    Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob1)(0);
                    //cout << "Winner: "; cout << Winner; cout << ", ";
                  }else{
                    // Patches Compete
                    PatchWinner = RcppArmadillo::sample(n_patches_v, 1, TRUE, AvgPhenotype)(0);
                    NumericVector prob2(n_off);
                    for(int n=0; n < n_off; n++){// Offspring On Winning Patch Compete
                      prob2(n) = OffspringPhenotype(n,PatchWinner);
                    }
                    Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob2)(0);
                  }

                  NumericVector WinnersGenes(2);
                  for(int n = 0; n < 2; n++){// Mutate Winner's Genes
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
              for(int n = 0; n < n_patches; n++){
                mean_perind(n) = mean_rcpp(Popn(_,n));
              }
              Data(path,t) = mean_rcpp(mean_perind);
            }
          }
          //MeanPath = apply(Data, 2, mean)
          //  SE=sqrt(apply(Data, 2, var)/paths)

          //return(list("Data" = Data, "Mean" = MeanPath, "SE" = SE))
          //return Data;
          Datalist2.push_back(Data);
          Datalist2.attr("names") = names;
          p.increment();
        }

      }
    }
  }

  return Datalist2;
}









// [[Rcpp::export]]
NumericMatrix coopbreed(int paths, int n_gener, int n_patches, double MutStep, int n_mates, int n_off, double par_c, double par_k){


  Progress p(paths * n_gener * n_patches, TRUE);

  NumericVector n_mates_v(n_mates);
  for(int n = 0; n < n_mates; n++){
    n_mates_v(n) = n;
  }

  NumericVector n_patches_v(n_patches);
  for(int n = 0; n < n_patches;n++){
    n_patches_v(n) = n;
  }
  NumericVector n_off_v(n_off);
  for(int n = 0; n < n_off;n++){
    n_off_v(n) = n;
  }
  NumericMatrix Data(paths,n_gener);

  for (int path = 0 ; path < paths; path++) {
    // cout << "Path: "; cout << path;
    // Initialize Population
    NumericMatrix Popn(2,n_patches);
    for (int i = 0; i < n_patches; i++){
      Popn(0,i) = R::runif(0,1);
      Popn(1,i) = R::runif(0,1);
      //cout << Popn(1,i); cout << ", ";
      //one individual per patch
      // two chromosomes per individual
    }

    //Popn(1,_) = runif(n_patches);
    //Popn(2,_) = runif(n_patches);
    for(int t = 0; t < n_gener; t++){ // for all generations
      //cout << "Generation: "; cout << t;
      NumericMatrix  mates(n_mates, n_patches);// stores index of each mate of each individual
      //for(int i = 0; i < n_mates; i++){
      //  for(int j = 0; j < n_patches; j++){
      //    mates(i,j) = 1; //RcppArmadillo::sample(n_patches_v, 1, TRUE, NumericVector::create());
      //  }

      //}
      NumericVector sample1(n_mates);
      for(int j = 0; j < n_patches; j++){
        sample1 = RcppArmadillo::sample(n_patches_v, n_mates, TRUE);
        for(int jj = 0; jj < n_mates; jj++){
          mates(jj,j) = sample1(jj);
          //cout << mates(jj,j); cout << ", ";
        }
      }


      //(n_off,2,n_patches);
      double Offspring[n_off][2][n_patches];// stores the genotypes of each off Offspring
      NumericVector vec = NumericVector::create(0,1);

      for(int i = 0; i < n_patches; i++){
        p.increment();
        //cout << "Patch nr. "; cout << i;
        NumericVector Fathers(n_off);
        Fathers = RcppArmadillo::sample(n_mates_v, n_off, TRUE);
        //cout << Fathers; cout << ", ";
        NumericVector PatContribs(n_off);
        PatContribs = RcppArmadillo::sample(vec, n_off, TRUE);
        NumericVector MatContribs(n_off);
        MatContribs = RcppArmadillo::sample(vec, n_off, TRUE);

        for(int j = 0; j < n_off; j++){ // for all offspring
          Offspring[j][0][i] = Popn(PatContribs(j), mates(Fathers(j), i));
          Offspring[j][1][i] = Popn(MatContribs(j),i);
          //cout << Offspring[j][0][i]; cout << ", ";
        }
      }

      NumericMatrix OffspringPhenotype(n_off, n_patches);
      NumericVector AvgPhenotype(n_patches);
      NumericVector vec2(2);

      for(int i = 0; i < n_off; i++){
        for(int j = 0; j < n_patches; j++){
          // TODO
          vec2(0) = Offspring[i][0][j];
          vec2(1) = Offspring[i][1][j];
          OffspringPhenotype(i,j) = mean_rcpp(vec2);
          //cout << OffspringPhenotype(i,j); cout << ", ";
        }
      }

      //for(int i = 0; i < n_patches; i++){
      //  NumericVector patch_vec(n_off);
      //  for(int ii = 0; ii < n_off; ii++){
      //    patch_vec(ii) = OffspringPhenotype(ii, i);
      //  }
      //  AvgPhenotype(i) = mean_rcpp(patch_vec);
      //}

      for(int i = 0; i < n_patches; i++){
        AvgPhenotype(i) = mean_rcpp(OffspringPhenotype(_,i));
      }
      double GlobalAvgPhenotype = mean_rcpp(AvgPhenotype);
      //cout << "avg = "; cout << AvgPhenotype;
      NumericVector ProbLocal(n_patches);
      for(int n = 0; n < n_patches; n++){
        ProbLocal(n) = (1 - AvgPhenotype(n))/(par_c * GlobalAvgPhenotype + (1 - AvgPhenotype(n))); //check if results is rigth
        //cout << "probloca: "; cout << ProbLocal(n); cout << ", ";
      }
      NumericVector BreederSurvival = Sb(AvgPhenotype, par_k);
      //cout << "breddersurv: "; cout << BreederSurvival; cout << ", ";
      NumericMatrix PopnNext(2, n_patches);


      for(int i = 0; i < n_patches; i++){ // for all patches
        //cout << "Patch: ";cout << i;cout << "  ";
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
            for(int n=0; n < n_off; n++){
              prob1(n) = 1.0 - OffspringPhenotype(n,PatchWinner);
            }
            Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob1)(0);
          }else{
            PatchWinner = RcppArmadillo::sample(n_patches_v, 1, TRUE, AvgPhenotype)(0); // Patches Compete
            NumericVector prob2(n_off);
            for(int n=0; n < n_off; n++){// Offspring On Winning Patch Compete
              prob2(n) = OffspringPhenotype(n,PatchWinner);
            }
            Winner = RcppArmadillo::sample(n_off_v,1,TRUE, prob2)(0);
          }
          NumericVector WinnersGenes(2);
          for(int n = 0; n < 2; n++){// Mutate Winner's Genes
            WinnersGenes(n) = Offspring[Winner][n][PatchWinner] + (MutStep * R::rnorm(0,1));
            if(WinnersGenes(n) <= 0){
              WinnersGenes(n) = MutStep;
            }
            if(WinnersGenes(n) >= 1){
              WinnersGenes(n) = 1 - MutStep;
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
      for(int n = 0; n < n_patches; n++){
        mean_perind(n) = mean_rcpp(Popn(_,n));
      }
      Data(path,t) = mean_rcpp(mean_perind);

    }
  }
  //MeanPath = apply(Data, 2, mean)
  //  SE=sqrt(apply(Data, 2, var)/paths)

  //return(list("Data" = Data, "Mean" = MeanPath, "SE" = SE))
  return Data;
}












