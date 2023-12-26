# LDSC Notes

<!-- badged -->
[![GitHub Pages](https://github.com/LucaJiang/ldsc_notes/actions/workflows/jekyll-gh-pages.yml/badge.svg?branch=gh-pages)](https://github.com/LucaJiang/ldsc_notes/actions/workflows/jekyll-gh-pages.yml)

Notes for LD Score Regression.

- [LDSC Notes](#ldsc-notes)
  - [Structure of the Repository](#structure-of-the-repository)
  - [How to Run](#how-to-run)
  - [Implementation Details](#implementation-details)
  - [Presentation Page](#presentation-page)
  - [Reference](#reference)

## Structure of the Repository

- 'master' branch: the source code of my implementation of LDSC.
  - 'ldscore' folder: the source code of LDSC.
  - 'ldsc.py' file: the entry point of LDSC.
  - 'data' folder: the data used by LDSC, i.e. $\ell^2$ and reference panel.
  - 'results' folder: the output log of LDSC.
- 'gh-pages' branch: the code of the presentation page.

## How to Run

- Calculate LD scores and estimate heritability:

```bash
python3 ldsc.py \
-M all \
-s ./data/full.sumstats \
-r ./data/eur_w_ld_chr/ \
-m CM -w 1e-4 \
-o test
```

- Calculate LD scores only:

```bash
python3 ldsc.py \
-M ldsc \
-s ./data/full.sumstats \
-r ./data/eur_w_ld_chr/ \
-m CM -w 1e-4 \
-o test
```

- Estimate heritability only:

```bash
python3 ldsc.py \
-M h2 \
-s ./results/test.out \
-o test
```

## Implementation Details

- [x] LD Score Calculation with MPI
  - [ ] Delete MAF < 0.01 SNPs
  - [x] Mapping from SNP in sumsstats to SNP in reference panel
  - [x] Calculate LD Score in each chromosome with sliding window
- [ ] Iterative ReWeighted Least Squares(IRWLS):
  $$\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$$
$$\text{argmin}\_{\boldsymbol{\beta}} \sum_{i=1}^n ||y_i - \mathbf{x}_i\boldsymbol{\beta}||_p$$
  where $\mathbf{X}=[\mathbf{1},\mathbf{\ell}]$ is a $n\times 2$ matrix.
  - [x] Initial with pesudo-inverse:
    $$\hat{\boldsymbol{\beta}}_0 = \text{pinv}(\mathbf{X})\mathbf{Y}$$
  - [x] Iteratively update the weight matrix, for $k=1,\ldots$, until convergence or max iteration:
      1. Calculate the residuals:
        $$\mathbf{r}_k = \mathbf{Y} - \mathbf{X}\hat{\boldsymbol{\beta}}_k$$
      2. Calculate the weight matrix:
        $$\mathbf{w}_k = | \mathbf{r}_k |^{p-2}$$
        $$\mathbf{W}\_k = \text{diag}(\mathbf{w}\_k/\sum \mathbf{w}\_k)$$
      3. Update the regression coefficient:
        $$\hat{\boldsymbol{\beta}}_{k+1} = (\mathbf{X}^T\mathbf{W}_k\mathbf{X})^{-1}\mathbf{X}^T\mathbf{W}_k\mathbf{Y}$$
      4. Check convergence:
        $$||\mathbf{r}\_k - \mathbf{r}\_{k+1}||_2 < \epsilon$$

  - [ ] jackknife estimate of the variance:
      <!-- $$\hat{\sigma}^2 = \frac{1}{n}\sum_{i=1}^n \frac{r_i^2}{w_i}$$ -->

## Presentation Page

The presentation on LDSC is [here](https://lucajiang.github.io/ldsc_notes/#/), which is generated base on the 'gh-pages' branch and deployed by github action.

## Reference

1. Bulik-Sullivan, Brendan K., et al. "LD score regression distinguishes confounding from polygenicity in genome-wide association studies." Nature genetics 47.3 (2015): 291-295. doi: [10.1038/ng.3211](https://doi.org/10.1038/ng.3211)
2. Burrus C S. Iterative reweighted least squares[J]. OpenStax CNX. Available online: [web.archive](https://web.archive.org/web/20221017041048/https://cnx.org/exports/92b90377-2b34-49e4-b26f-7fe572db78a1@12.pdf/iterative-reweighted-least-squares-12.pdf), 2012, 12.
