#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Draw exponential-race auxiliary variables and calculate each ranking's
// item-level exposure.  `ranking_cluster` is one-based and selects a row of
// `latent_strength`; the same kernel covers PL mixtures and PL--LBM after an
// LBM cell-strength matrix has been expanded to individual items.
// [[Rcpp::export]]
List pl_grouped_augmentation_cpp(const IntegerMatrix& rankings,
                                 const IntegerVector& ranking_cluster,
                                 const NumericMatrix& latent_strength) {
  const int n_rankings = rankings.nrow();
  const int rank_length = rankings.ncol();
  const int n_groups = latent_strength.nrow();
  const int n_items = latent_strength.ncol();

  if (ranking_cluster.size() != n_rankings) {
    stop("`ranking_cluster` must have one label per ranking.");
  }
  if (n_groups < 1 || n_items < 1) {
    stop("`latent_strength` must have at least one row and one column.");
  }

  NumericMatrix waiting_time(n_rankings, rank_length);
  NumericMatrix exposure(n_rankings, n_items);

  for (int ranking_index = 0; ranking_index < n_rankings; ++ranking_index) {
    const int group_index = ranking_cluster[ranking_index] - 1;
    if (group_index < 0 || group_index >= n_groups) {
      stop("`ranking_cluster` contains an out-of-range label.");
    }

    double total_strength = 0.0;
    for (int item_index = 0; item_index < n_items; ++item_index) {
      const double value = latent_strength(group_index, item_index);
      if (!R_finite(value) || value <= 0.0) {
        stop("`latent_strength` must be finite and strictly positive.");
      }
      total_strength += value;
    }

    double selected_strength = 0.0;
    double elapsed_time = 0.0;
    for (int position = 0; position < rank_length; ++position) {
      const int item_index = rankings(ranking_index, position) - 1;
      if (item_index < 0 || item_index >= n_items) {
        stop("`rankings` contains an out-of-range item id.");
      }
      const double risk_set_strength = total_strength - selected_strength;
      if (!R_finite(risk_set_strength) || risk_set_strength <= 0.0) {
        stop("PL risk-set strength became non-positive.");
      }
      const double draw = R::rexp(1.0 / risk_set_strength);
      waiting_time(ranking_index, position) = draw;
      elapsed_time += draw;
      exposure(ranking_index, item_index) += elapsed_time;
      selected_strength += latent_strength(group_index, item_index);
    }

    for (int item_index = 0; item_index < n_items; ++item_index) {
      bool ranked = false;
      for (int position = 0; position < rank_length; ++position) {
        if (rankings(ranking_index, position) - 1 == item_index) {
          ranked = true;
          break;
        }
      }
      if (!ranked) exposure(ranking_index, item_index) += elapsed_time;
    }
  }

  return List::create(
    Named("augmentation") = waiting_time,
    Named("exposure") = exposure
  );
}
