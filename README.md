[![DOI](https://zenodo.org/badge/168929183.svg)](https://zenodo.org/badge/latestdoi/168929183)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ttTensor)](
https://cran.r-project.org/package=ttTensor)
[![Downloads](https://cranlogs.r-pkg.org/badges/ttTensor)](https://CRAN.R-project.org/package=ttTensor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ttTensor?color=orange)](https://CRAN.R-project.org/package=ttTensor)
[![:name status badge](https://rikenbit.r-universe.dev/badges/:name)](https://rikenbit.r-universe.dev)
[![:registry status badge](https://rikenbit.r-universe.dev/badges/:registry)](https://rikenbit.r-universe.dev)
[![:total status badge](https://rikenbit.r-universe.dev/badges/:total)](https://rikenbit.r-universe.dev)
[![ttTensor status badge](https://rikenbit.r-universe.dev/badges/ttTensor)](https://rikenbit.r-universe.dev)
![GitHub Actions](https://github.com/rikenbit/ttTensor/actions/workflows/build_test_push.yml/badge.svg)

# ttTensor
R package for Tensor-Train Decomposition

Installation
======
~~~~
git clone https://github.com/rikenbit/ttTensor/
R CMD INSTALL ttTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/ttTensor")
~~~~

References
======
- TTSVD : [Tensor-Train Decomposition, I. V. Oseledets, 2011](https://epubs.siam.org/doi/10.1137/090752286)
- TTWOPT : [Completion of high order tensor data with missing entries via tensor-train decomposition, Yuan, Longhao, et. al., 2017](https://arxiv.org/abs/1709.02641)
- TTCross : [TT-cross approximation for multidimensional arrays, I. V. Oseledets, et. al., 2010](https://www.sciencedirect.com/science/article/pii/S0024379509003747), [On selecting a maximum volume sub-matrix of a matrix and related problems, Ali Civril, et. al., 2009](https://www.sciencedirect.com/science/article/pii/S0304397509004101)

## License
Copyright (c) 2019 Koki Tsuyuzaki and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [Artistic License 2.0](https://www.perlfoundation.org/artistic-license-20.html).

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido
