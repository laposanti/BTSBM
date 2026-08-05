# Set partition-move controls

Set partition-move controls

## Usage

``` r
partition_moves(
  split_merge = 0L,
  ranking_split_merge = 0L,
  item_split_merge = 0L,
  pair_proposal = c("uniform", "balanced", "informed"),
  restricted_scans = 1L
)
```

## Arguments

- split_merge:

  Number of generic split–merge moves per iteration.

- ranking_split_merge, item_split_merge:

  PL–LBM moves for the ranking and item partitions respectively.

- pair_proposal:

  Pair-selection rule for a future split–merge kernel.

- restricted_scans:

  Number of restricted Gibbs scans.

## Value

An object of class `btsbm_partition_moves`.
