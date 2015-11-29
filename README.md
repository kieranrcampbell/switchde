
# switchde
Switchde provides distribution fitting for sigmoidal differential gene expression across pseudotime. It supports zero-inflation similar to that found in Pearson & Yau (Genome Biology, 2015), implemented through the EM algorithm.

## Usage

Access to the fitting (unless you want to go into the internals) is through two functions:

```r
fitModel(object, pseudotime = NULL, zero_inflated = FALSE)
testDE(object, pseudotime = NULL, zero_inflated = FALSE)
```

`fitModel` returns the sigmoidal model parameters while testDE returns the parameters as well as a p-value compared to the null model. The parameters are:
- `object`: Gene expression across cells. Either a vector, a matrix (columns cells, rows genes) or an `SCESet` object
- `pseudotime`: A pseudotime vector of length number of cells. Can be `NULL` if `object` is an `SCESet` and `pData(object)$pseudotime` is defined
- `zero_inflated`: Should a zero-inflated component be modelled? While this is marginally more accurate, it is far slower.