<!-- title page -->
# Notes on LDSC

Jiang Wenxin

2023.12.22

--------------------
<!-- toc -->
- [Notes on LDSC](#notes-on-ldsc)
  - [Individual-Level Model](#individual-level-model)
  - [From Individual-level to LD Score and Heritability](#from-individual-level-to-ld-score-and-heritability)
  - [Summary-Level Model](#summary-level-model)
    - [Proof of Proposition 1](#proof-of-proposition-1)
  - [Reference](#reference)

--------------------

<!-- ## Introduction -->

<!-- $$\begin{equation}\tag{eq:1}
\begin{aligned}
\sum_{i=1}^n b&=0\\
\sum_{i=1}^n c&=0\\
\end{aligned}
\end{equation}$$ -->

<!-- -------------------- -->

## Individual-Level Model

**Model:**
$$\begin{equation}\tag{1.1}
\phi=X\beta+\epsilon
\end{equation}$$

\(\phi\): N by 1 vector of phenotypes
\(X\): N by M matrix of normalized genotypes
\(\beta\): M by 1 vector of per-normalized-genotype effect sizes
\(\epsilon\): N by 1 vector of environmental effects

**Expectation and Variance:**

<!-- table center-->

<div align="center">

| | \(\mathbb{E}\) | \(\text{Var}\) |
| --- | --- | --- |
| \(\beta\) | 0 | \(\frac{h_g^2}{M}I\) |
| \(\epsilon\) | 0 | \(\left(1-h_g^2\right)I\) |
| \(\phi\) | 0 | \(\frac{h_g^2}{M}XX^\text{T}+\left(1-h_g^2\right)I\) |

</div>

**Assumptions:**

1. The genotype at variant \(j\) for individual \(i\) is independent of other individuals' genotypes.
2. \(X,\beta\) and \(\epsilon\) are mutually independent.

--------------------

## From Individual-level to LD Score and Heritability

**Linkage disequilibrium (LD):**
Define the correlation between two variants \(j\) and \(k\) as
$$\begin{equation*}r_{jk}:=\mathbb{E}\left[X_{ij}X_{ik}\right]\end{equation*}$$

**LD Score:**
Define the LD score of variant \(j\) as
$$\begin{equation}\tag{1.2}\ell_j:=\sum_{k=1}^M r_{jk}^2\end{equation}$$

**Heritability:**
Since \(\text{Var}\left[\phi|X\right]=\frac{h_g^2}{M}XX^\text{T}+\left(1-h_g^2\right)I\), we have
$$\begin{equation}\tag{1.3}h_g^2=\frac{\text{Var}\left[\phi|X\right]-I}{XX^\text{T}/M-I}\end{equation}$$

> \(h_g^2\) is the proportion of phenotypic variance explained by the genotypes. It is also called the SNP heritability. It is the ratio of the variance of the phenotype explained by the genotypes to the total variance of the phenotype.

--------------------

## Summary-Level Model

**\(\chi^2\) statistic:**
Define the \(\chi^2\) statistic of variant \(j\) as
$$\begin{equation*}\chi^2_j:=N\hat{\beta}_j^2\end{equation*}$$
where \(\hat{\beta}_j:=X_j^\text{T}\phi/N\) is the estimated effect size of variant \(j\).

**Proposition 1:**
The expected \(\chi^2\)-value of variant \(j\) is
$$\begin{equation}\tag{1.3}\mathbb{E}\left[\chi^2_j\right]=\frac{N h_g^2}{M}\ell_j+1\end{equation}$$

--------------------

### Proof of Proposition 1

Since \(\chi^2_j=N\hat{\beta}_j^2\), we have
$$\begin{equation*}\begin{aligned}
\mathbb{E}\left[\chi^2_j\right]&=\mathbb{E}\left[N\hat{\beta}_j^2\right]\\
&=N\left(\text{Var}\left[\hat{\beta}_j\right]+\mathbb{E}^2\left[\hat{\beta}_j\right]\right)\qquad [\text{MoM}]\\
&=N\text{Var}\left[\hat{\beta}_j\right] \qquad \left[\mathbb{E}\left[\hat{\beta}_j\right]=0 \right]\\
\end{aligned}\end{equation*}$$

Obtain \(\text{Var}\left[\hat{\beta}_j\right]\) by the law of total variance:
$$\begin{equation}\tag{1.4}\begin{aligned}
\text{Var}\left[\hat{\beta}_j\right]&=\mathbb{E}\left[\text{Var}\left[\hat{\beta}_j|X\right]\right]+\text{Var}\left[\mathbb{E}\left[\hat{\beta}_j|X\right]\right]\\
&=\mathbb{E}\left[\text{Var}\left[\hat{\beta}_j|X\right]\right]+0 \qquad \left[\mathbb{E}\left[\hat{\beta}_j|X\right]=0 \right]\\
&=\mathbb{E}\left[\text{Var}\left[X_j^\text{T}\phi/N|X\right]\right] \qquad \left[\text{LSE}: \hat{\beta}_j=X_j^\text{T}\phi/N \right]\\
&=\mathbb{E}\left[X_j\text{Var}\left[\phi|X\right]X_j^\text{T}/N^2\right]\\
&=\frac{1}{N^2}\mathbb{E}\left[\left(\frac{h_g^2}{M}X_j^\text{T}XX^\text{T}X_j+N\left(1-h_g^2\right)\right)\right] \qquad \left[X_j^\text{T}X_j=N\right]
\end{aligned}\end{equation}$$

Let \(\tilde{r}_{jk}:=\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}\) be the sample correlation between additively-coded genotypes at variants \(j\) and \(k\). Then, rewrite the term on the left as
$$\begin{equation}\tag{1.6}\frac{1}{N^2}X_j^\text{T}XX^\text{T}X_j=\sum_{k=1}^M \tilde{r}_{jk}^2\end{equation}$$

Obtain the expectation of \(\tilde{r}_{jk}^2\) by delta method:
$$\begin{equation}\tag{1.7}\begin{aligned}
\mathbb{E}\left[\tilde{r}_{jk}^2\right]&=\text{Var}\left[\tilde{r}_{jk}\right]+\mathbb{E}^2\left[\tilde{r}_{jk}\right] \qquad [\text{MoM}]\\
&=\text{Var}\left[\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}\right]+\mathbb{E}^2\left[\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}\right]\\
&=\frac{1}{N^2}\sum_{i=1}^N\text{Var}\left[X_{ij}X_{ik}\right]+\frac{1}{N^2}\sum_{i=1}^N\mathbb{E}^2\left[X_{ij}X_{ik}\right] \qquad \left[\text{independence of} i\right]\\
&=\frac{1}{N^2}\sum_{i=1}^N\text{Var}\left[X_{ij}X_{ik}\right]+r_{jk}^2 \qquad \left[\mathbb{E}\left[X_{ij}X_{ik}\right]=r_{jk}\right]\\
&\approx r_{jk}^2+(1-r_{jk}^2)/N \qquad \left[\text{delta method}\right] ??
\end{aligned}\end{equation}$$

We have
$$\begin{equation}\tag{1.8}\begin{aligned}
\mathbb{E}\left[\sum_{k=1}^M \tilde{r}_{jk}^2\right] &\approx \sum_{k=1}^M r_{jk}^2+\frac{\sum_{k=1}^M\left(1-r_{jk}^2\right)}{N}\\
&\approx \ell_j+\frac{M-\ell_j}{N} \qquad \left[\text{definition of} \ell_j\right]
\end{aligned}\end{equation}$$

Therefore, we can continue the derivation of \(\text{Var}\left[\hat{\beta}_j\right]\) as
$$\begin{equation*}\begin{aligned}
\text{Var}\left[\hat{\beta}_j\right]&=\frac{1}{N^2}\mathbb{E}\left[\left(\frac{h_g^2}{M}X_j^\text{T}XX^\text{T}X_j+N\left(1-h_g^2\right)\right)\right]\\
&=\left(\frac{h_g^2}{M}\mathbb{E}\left[\sum_{k=1}^M \tilde{r}_{jk}^2\right]+\left(1-h_g^2\right)/N\right) \qquad \left[\text{Eq. (1.6)}\right]\\
&=\left(\frac{h_g^2}{M}\left(\ell_j+\frac{M-\ell_j}{N}\right)+\left(1-h_g^2\right)/N\right) \qquad \left[\text{Eq. (1.8)}\right]\\
&=\frac{(1-1/N)h_g^2}{M}\ell_j+\frac{1}{N}\\
&\approx\frac{h_g^2}{M}\ell_j+\frac{1}{N} \qquad \left[1-1/N\approx 1\right]
\end{aligned}\end{equation*}$$

Finally, we can obtain the expectation of \(\chi^2_j\) as
$$\begin{equation}\tag{1.9}\mathbb{E}\left[\chi^2_j\right]=N\text{Var}\left[\hat{\beta}_j\right]\approx\frac{N h_g^2}{M}\ell_j+1\end{equation}$$

--------------------

## Reference

1. B. K. Bulik-Sullivan et al., LD Score regression distinguishes confounding from polygenicity in genome-wide association studies, Nat Genet, vol. 47, no. 3, Art. no. 3, Mar. 2015, doi: 10.1038/ng.3211.

<!-- 2. ... -->
