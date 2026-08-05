# Toy datasets and model workflows

~{r setup, include=FALSE} knitr::opts_chunk\$set(collapse = TRUE,
comment = “#\>”) ~

BTSBM includes small, reproducible synthetic datasets for learning the
interfaces and writing regression tests. They are deliberately not
copies of respondent-level or match-level data from real applications.

## Pairwise tennis league

tennis_toy() creates a round-robin paired-comparison matrix with three
generating player tiers. It is useful for the existing BT and BT–SBM
models.

~{r tennis} library(BTSBM)

tennis \<- tennis_toy(matches_per_pair = 8, seed = 1) tennis

tennis_fit \<- fit_btsbm( tennis, bt_sbm_model( item_clustering =
gnedin_prior(gnedin_hyperparameter = 0.5) ), mcmc_control(iter = 100,
warmup = 50, seed = 2) )

posterior_strength(tennis_fit, summary = “mean”)
posterior_similarity(tennis_fit, target = “item”) ~

## Eurovision-inspired ranking mixture

eurovision_toy() represents juror/viewer ranking profiles over fictional
acts. The useful inferential target is the partition of ranking rows, so
a PL mixture is natural.

~{r eurovision} votes \<- eurovision_toy(n_rankings = 30, rank_length =
5, seed = 3) votes

mixture_fit \<- fit_btsbm( votes, pl_mixture_model( ranking_clustering =
dirichlet_process_prior(concentration = 1) ), mcmc_control(iter = 120,
warmup = 60, seed = 4) )

posterior_similarity(mixture_fit, target = “ranking”)
mcmc_diagnostics(mixture_fit) ~

The matrix returned by posterior_similarity() is preferable to raw
cluster labels because the labels themselves are exchangeable across
MCMC draws.

## Sushi-inspired Plackett–Luce model

sushi_toy() is a synthetic dataset using Sushi Preference Dataset-style
items. It supports a simple PL model, PL–SBM, or a ranking mixture. The
original Sushi data are not bundled; see the PL ranking-model vignette
for its source and licence.

~{r sushi} sushi \<- sushi_toy(n_rankings = 30, rank_length = 5, seed =
5)

pl_fit \<- fit_btsbm( sushi, pl_model(latent_strength =
latent_strength(shape = 1, rate = 1)), mcmc_control(iter = 100, warmup =
50, seed = 6) )

posterior_strength(pl_fit, summary = “mean”) mcmc_diagnostics(pl_fit) ~

## Joint ranking/item structure

pl_lbm_toy() contains known ranking and item groups. It is a compact
fixture for a PL–LBM workflow and for testing both posterior-similarity
matrices.

~{r lbm} structured_rankings \<- pl_lbm_toy(n_rankings = 20, rank_length
= 4, seed = 7) lbm_fit \<- fit_btsbm( structured_rankings, pl_lbm_model(
ranking_clustering = finite_partition(max_clusters = 3), item_clustering
= finite_partition(max_clusters = 3) ), mcmc_control(iter = 120, warmup
= 60, seed = 8) )

mcmc_diagnostics(lbm_fit) posterior_similarity(lbm_fit, target =
“ranking”) posterior_similarity(lbm_fit, target = “item”) ~

For real analyses, run multiple chains, inspect mcmc_diagnostics(), and
perform posterior predictive checks before interpreting a latent
partition.
