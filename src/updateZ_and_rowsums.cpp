#include <Rcpp.h>
using namespace Rcpp;

//'  @keywords internal
 // [[Rcpp::export]]
 Rcpp::List updateZ_and_rowsums(const Rcpp::IntegerMatrix& n_ij,
                                const Rcpp::IntegerVector& x,
                                const Rcpp::NumericVector& lambda) {
   const int K = n_ij.nrow();
   Rcpp::NumericMatrix Z(K, K);
   Rcpp::NumericVector rowSumsZ(K);

   Rcpp::NumericVector lam_i(K);
   for (int i = 0; i < K; ++i) {
     int lab = x[i];
     if (lab <= 0 || lab > lambda.size() || Rcpp::NumericVector::is_na(lambda[lab - 1])) {
       lam_i[i] = NA_REAL;
     } else {
       lam_i[i] = lambda[lab - 1];
     }
   }

   for (int i = 0; i < K; ++i) {
     for (int j = i + 1; j < K; ++j) {
       int nij = n_ij(i, j);
       if (nij > 0) {
         double rate = lam_i[i] + lam_i[j];
         if (!R_finite(rate) || rate <= 0.0) {
           Z(i, j) = 0.0; Z(j, i) = 0.0;
         } else {
           double draw = R::rgamma((double)nij, 1.0 / rate); // shape, scale
           Z(i, j) = draw; Z(j, i) = draw;
           rowSumsZ[i] += draw; rowSumsZ[j] += draw;
         }
       } else {
         Z(i, j) = 0.0; Z(j, i) = 0.0;
       }
     }
   }

   return Rcpp::List::create(_["Z"] = Z, _["rowSumsZ"] = rowSumsZ);
 }
