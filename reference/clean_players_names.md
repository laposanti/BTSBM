# Format player/item names as "Surname F."

Format player/item names as "Surname F."

## Usage

``` r
clean_players_names(name)
```

## Arguments

- name:

  character scalar or vector (e.g., "Roger Federer").

## Value

Character vector with formatted names (e.g., "Federer R.").

## Examples

``` r
clean_players_names(c("Roger Federer", "Nadal"))
#> [1] "Federer R ." "Nadal"      
```
