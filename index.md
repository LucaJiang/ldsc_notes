<!-- title page -->
# Notes on LDSC

Jiang Wenxin

2023.12.22

--------------------
<!-- toc -->
- [Notes on LDSC](#notes-on-ldsc)
  - [Motivation, Aims and Formula](#motivation-aims-and-formula)
  - [Individual-Level Model](#individual-level-model)
  - [From Individual-level to LD Score and Heritability](#from-individual-level-to-ld-score-and-heritability)
  - [Summary-Level Model](#summary-level-model)
    - [Proof of Proposition 1](#proof-of-proposition-1)
  - [Algorithm for Estimating chi-squared and Heritability](#algorithm-for-estimating-chi-squared-and-heritability)
  - [Limitations](#limitations)
  - [Reference](#reference)
  - [Cross-trait LDSC](#cross-trait-ldsc)
    - [Shared Samples](#shared-samples)

--------------------

<!-- ## Introduction -->

<!-- $$\begin{equation}\tag{eq:1}
\begin{aligned}
\sum_{i=1}^n b&=0\\
\sum_{i=1}^n c&=0\\
\end{aligned}
\end{equation}$$ -->

<!-- -------------------- -->

## Motivation, Aims and Formula

- **Individual-level model:** $\phi=X\beta+\epsilon$. If individual data is available, we can estimate the effect size $\beta$ directly.
- **Summary-level data:** If only the z-score of variant $j$
$$z_j=\left(\hat{\beta}_j-\bar{\beta}_j\right)/\text{SE}_j=\hat{\beta}_j$$
is available, where $\hat{\beta}_j$ is the estimated effect size and $\text{SE}_j$ is the standard error of $\hat{\beta}_j$, then calculate the $\chi^2$ statistic of variant $j$ with $$\chi^2_j=N\hat{\beta}_j^2=Nz_j^2.$$

--------------------

- **Aim:** Estimate and quantify the contribution of the polygenic signal and the confounding factors to the $\chi^2$ statistics of a variant in GWAS summary statistics. Then we can estimate the SNP heritability $h_g^2$ by regressing the $\chi^2$ statistics on the LD score of each variant. And we can calculate the genetic correlation between two traits.
- **Formula**: Examining the relationship between the $\chi^2$ statistics and the LD score of a variant, $$\begin{equation*}\mathbb{E}\left[\chi^2_j|\ell_j\right]=\frac{N h_g^2}{M}\ell_j+Na+1\end{equation*}$$
where $\ell_j$ is the LD score of variant $j$, $h_g^2$ is the SNP heritability, $N$ is the sample size, $M$ is the number of variants, and $a$ measures the confounding factors.

--------------------

## Individual-Level Model

**Model:**
$$\begin{equation}
\phi=X\beta+\epsilon\tag{1.1}
\end{equation}$$

<!-- $\phi$: N by 1 vector of phenotypes

$X$: N by M matrix of normalized genotypes

$\beta$: M by 1 vector of per-normalized-genotype effect sizes

$\epsilon$: N by 1 vector of environmental effects

**Expectation and Variance:** -->

<!-- table -->

<!-- | | $\mathbb{E}$ | $\text{Var}$ |
| --- | --- | --- |
| $\beta$ | 0 | $\frac{h_g^2}{M}I$ |
| $\epsilon$ | 0 | $\left(1-h_g^2\right)I$ |
| $\phi$ | 0 | $\frac{h_g^2}{M}XX^\text{T}+\left(1-h_g^2\right)I$ | -->

**Notations:**

| Variable | Size | Description | $\mathbb{E}$ | $\text{Var}$ |
| --- | --- | --- | --- | --- |
| $\phi$ | N*1 | Phenotype vector | $0$  | $\frac{h_g^2}{M}XX^\text{T}+\left(1-h_g^2\right)I$ |
| $X$ | N*M | Normalized Genotype matrix | $0$ | $I$ |
| $\beta$ | M*1 | Per-normalized Effect size vector | $0$ | $\frac{h_g^2}{M}I$ |
| $\epsilon$ | N*1 | Environmental effect vector | $0$ | $\left(1-h_g^2\right)I$ |

**Assumptions:**

1. The genotype at variant $j$ for individual $i$ is independent of other individuals' genotypes.
2. $X,\beta$ and $\epsilon$ are mutually independent.

--------------------

## From Individual-level to LD Score and Heritability

**Linkage disequilibrium (LD):**
Define the correlation between two variants $j$ and $k$ as
$$\begin{equation*}r_{jk}:=\mathbb{E}\left[X_{ij}X_{ik}\right]\end{equation*}$$

*Note:* $r_{jk}=\mathbb{E}\left[\frac{\text{Cov}\left[X_{ij},X_{ik}\right]}{\sqrt{\text{Var}\left[X_{ij}\right]\text{Var}\left[X_{ik}\right]}}\right]=\mathbb{E}\left[X_{ij}X_{ik}\right] - \mathbb{E}\left[X_{ij}\right]\mathbb{E}\left[X_{ik}\right]=\mathbb{E}\left[X_{ij}X_{ik}\right]$

**LD Score:**
Define the LD score of variant $j$ as
$$\begin{equation}\ell_j:=\sum_{k=1}^M r_{jk}^2\tag{1.2}\end{equation}$$

**Heritability:**
Since $\text{Var}\left[\phi|X\right]=\frac{h_g^2}{M}XX^\text{T}+\left(1-h_g^2\right)I$, we have
$$\begin{equation*}h_g^2=\frac{\text{Var}\left[\phi|X\right]-I}{XX^\text{T}/M-I}\end{equation*}$$

*Note:* $h_g^2$ is the proportion of phenotypic variance explained by the genotypes. It is also called the SNP heritability. It is the ratio of the variance of the phenotype explained by the genotypes to the total variance of the phenotype.

--------------------

## Summary-Level Model

**$\chi^2$ statistic:**
Define the $\chi^2$ statistic of variant $j$ as
$$\begin{equation*}\chi^2_j:=N\hat{\beta}_j^2\end{equation*}$$

where $\hat{\beta}_j:=X_j^\text{T}\phi/N$ is the estimated effect size of variant $j$. Here $H_0: \hat{\beta}_j\sim\mathcal{N}\left(0,\text{Var}\left[\hat{\beta}_j\right]\right)$, when $N$ is large enough.

**Proposition 1:**
The expected $\chi^2$-value of variant $j$ is
$$\begin{equation}\mathbb{E}\left[\chi^2_j\right]=\frac{N h_g^2}{M}\ell_j+1\tag{1.3}\end{equation}$$

--------------------

### Proof of Proposition 1

Since $\chi^2_j=N\hat{\beta}_j^2$, we have
$$\begin{equation*}\begin{aligned}
\mathbb{E}\left[\chi^2_j\right]&=\mathbb{E}\left[N\hat{\beta}_j^2\right]\\
&=N\left(\text{Var}\left[\hat{\beta}_j\right]+\mathbb{E}^2\left[\hat{\beta}_j\right]\right)\qquad [\text{MoM}]\\
&=N\text{Var}\left[\hat{\beta}_j\right] \qquad \left[\mathbb{E}\left[\hat{\beta}_j\right]=0 \right]\\
\end{aligned}\end{equation*}$$

Obtain $\text{Var}\left[\hat{\beta}_j\right]$ by the law of total variance:
$$\begin{equation}\begin{aligned}
\text{Var}\left[\hat{\beta}_j\right]&=\mathbb{E}\left[\text{Var}\left[\hat{\beta}_j|X\right]\right]+\text{Var}\left[\mathbb{E}\left[\hat{\beta}_j|X\right]\right]\\
&=\mathbb{E}\left[\text{Var}\left[\hat{\beta}_j|X\right]\right]+0 \qquad \left[\mathbb{E}\left[\hat{\beta}_j|X\right]=0 \right]\\
&=\mathbb{E}\left[\text{Var}\left[X_j^\text{T}\phi/N|X\right]\right] \qquad \left[\text{LSE}: \hat{\beta}_j=X_j^\text{T}\phi/N \right]\\
&=\mathbb{E}\left[X_j\text{Var}\left[\phi|X\right]X_j^\text{T}/N^2\right]\\
&=\mathbb{E}\left[\frac{h_g^2}{MN^2}X_j^\text{T}XX^\text{T}X_j+\frac{1-h_g^2}{N}\right] \quad \left[X_j^\text{T}X_j=N\right]
\end{aligned}\tag{1.4}\end{equation}$$

--------------------

Let $\tilde{r}_{jk}:=\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}$ be the sample correlation between additively-coded genotypes at variants $j$ and $k$. Then, rewrite the term on the left as
$$\begin{equation}\tag{1.6}\frac{1}{N^2}X_j^\text{T}XX^\text{T}X_j=\sum_{k=1}^M \tilde{r}_{jk}^2\end{equation}$$

Obtain the expectation of $\tilde{r}_{jk}^2$ by delta-method:
$$\begin{equation}\begin{aligned}
\mathbb{E}\left[\tilde{r}_{jk}^2\right]&=\text{Var}\left[\tilde{r}_{jk}\right]+\mathbb{E}^2\left[\tilde{r}_{jk}\right] \qquad [\text{MoM}]\\
% &=\text{Var}\left[\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}\right]+\mathbb{E}^2\left[\frac{1}{N}\sum_{i=1}^N X_{ij}X_{ik}\right]\\
% &=\frac{1}{N^2}\sum_{i=1}^N\text{Var}\left[X_{ij}X_{ik}\right]+\frac{1}{N^2}\left(\sum_{i=1}^N\mathbb{E}\left[X_{ij}X_{ik}\right]\right)^2 \\
% &=\frac{1}{N^2}\sum_{i=1}^N\text{Var}\left[X_{ij}X_{ik}\right]+r_{jk}^2 \qquad \left[\mathbb{E}\left[X_{ij}X_{ik}\right]=r_{jk}\right]\\
&\approx r_{jk}^2+(1-r_{jk}^2)/N \qquad \left[\text{delta method}\right] \textbf{??}
\end{aligned}\tag{1.7}\end{equation}$$

> [Discussions on google group](https://groups.google.com/g/ldsc_users/c/mxbnbodkGj0):
>
> $$\begin{equation*}\begin{aligned}\mu_2^{\prime}&=1-\frac{(n-2)\left(1-\rho^2\right)}{n-1} F\left(1,1, \frac{1}{2}(n+1) ; \rho^2\right)\\&\rightarrow 1-\left(1-\frac{1}{n-1}\right)\left(1-\rho^2\right), n\rightarrow \infty\end{aligned}\end{equation*}$$
>
> Here $\mu_2^{\prime}$ is the 2nd moment of sample correlation, $\rho$ is the population correlation coefficient, $\mathrm{n}$ is # observations.

<!-- $$
F\left(\alpha, \beta, \delta ; \rho^2\right)=\frac{\Gamma(\delta)}{\Gamma(\alpha) \Gamma(\beta)} \sum_{j=0}^{\infty} \frac{\Gamma(\alpha+j) \Gamma(\beta+j)}{\Gamma(\delta+j)} \frac{\left(\rho^2\right)^j}{j !}
$$ -->

--------------------

We have
$$\begin{equation}\begin{aligned}
\mathbb{E}\left[\sum_{k=1}^M \tilde{r}_{jk}^2\right] &\approx \sum_{k=1}^M r_{jk}^2+\frac{\sum_{k=1}^M\left(1-r_{jk}^2\right)}{N}\\
&\approx \ell_j+\frac{M-\ell_j}{N} \qquad \left[\text{def. of } \ell_j\right]
\end{aligned}\tag{1.8}\end{equation}$$

Therefore, we can continue the derivation of $\text{Var}\left[\hat{\beta}_j\right]$ as
$$\begin{equation*}\begin{aligned}
\text{Var}\left[\hat{\beta}_j\right]&= \mathbb{E}\left[\frac{h_g^2}{MN^2}X_j^\text{T}XX^\text{T}X_j+\frac{1-h_g^2}{N}\right]\\
&=\frac{h_g^2}{MN^2}\mathbb{E}\left[X_j^\text{T}XX^\text{T}X_j\right]+\frac{1-h_g^2}{N}\\
&=\frac{h_g^2}{MN^2}\mathbb{E}\left[\sum_{k=1}^M \tilde{r}_{jk}^2\right]+\frac{1-h_g^2}{N} \qquad \left[\text{Eq. (1.6)}\right]\\
&=\frac{h_g^2}{MN^2}\left(\ell_j+\frac{M-\ell_j}{N}\right)+\frac{1-h_g^2}{N} \qquad \left[\text{Eq. (1.8)}\right]\\
&=\frac{(1-1/N)h_g^2}{M}\ell_j+\frac{1}{N}\\
&\approx\frac{h_g^2}{M}\ell_j+\frac{1}{N} \qquad \left[1-1/N\rightarrow 1\right]
\end{aligned}\end{equation*}$$

Finally, we can obtain the expectation of $\chi^2_j$ as
$$\begin{equation}\tag{1.9}\mathbb{E}\left[\chi^2_j\right]=N\text{Var}\left[\hat{\beta}_j\right]\approx\frac{N h_g^2}{M}\ell_j+1\end{equation}$$

--------------------

## Algorithm for Estimating chi-squared and Heritability

1. Calculate the LD score $\ell_j=\sum_{k=1}^M r_{jk}^2$ for each variant $j$ from the sumstat of a reference panel.
2. Calculate the expected $\chi^2$-value of each variant $j$ with $\mathbb{E}\left[\chi^2_j\right]=Nz_j^2$.
3. Estimate the $h_g^2$ by regressing the expected $\chi^2$-value on the LD score $\ell_j$ according to $\mathbb{E}\left[\chi^2_j\right]=\frac{N h_g^2}{M}\ell_j+1$.

--------------------

## Limitations

- The same LD score is used for all pairs of traits.
- The LD score is a weak instrument that explains only a small proportion of variance in the dependent variable.
- The method relies on numerous linearity and independence assumptions.
- The method does not provide a visual representation of an LD score regression analysis.
- The method does not perform standard regression diagnostics.
- The method does not always provide consistent results with known causal relationships.
- The method does not always provide a validation against known causal relationships.
- The method does not provide a causal relationship between two traits.

<!-- A criticism of LD score regression is that every analysis for each pair of traits uses the same LD scores as the dependent variable in the regression model (and as LD scores have been pre- computed by its proponents, literally the same LD scores are used in the majority of applied analyses). This means that any influential points in the regression will affect not only one LD score regression analysis, but all such analyses. LD scores are also likely to be a ‘weak instru- ment’ in the language of Mendelian randomization, as they will only explain a small propor- tion of variance in the dependent variable. Additionally, due to the scale of the data, it is not possible to provide a visual representation of an LD score regression analysis. Standard regres- sion diagnostics are rarely, if ever, performed. Finally, results from LD score regression are not always consistent with known causal relationships; for example, the method did not find evidence for a genetic correlation between LDL cholesterol and CHD risk that survived a multiple- testing correction (Bulik-Sullivan et al., 2015). The method has utility in mapping the genetic distance between related phenotypes, such as determining how closely related different psychiatric disorders are in terms of their genetic predictors (Cross-Disorder Group of the Psychiatric Genomics Consortium, 2013). However, the reliance of the method on numerous linearity and independence assumptions, incorrect weighting in the linear regression model (correct weights would require computation of the Cholesky decomposition of a matrix with dimension equal to the number of genetic variants in the model – misspecified weights are recommended for use in practice), and lack of validation against known causal relationships mean that results from the method should not be treated too seriously as an assessment of causality. -->
--------------------
- No constraint is applied to heritability estimates, therefore it is possible to get non-interpretable estimates that are below 0 or above 1.
- Variance of estimators is calculated using a resampling method, the blockwise jackknife. For moderate GWAS sample size, this approach often leads to very wide confidence intervals for enrichment, thereby reducing statistical power.
- z-scores for SNPs in LD may be strongly correlated. LD score regression reduces the impact of dependence among data points by using specially designed weights in weighted least squares estimation and only including HapMap SNPs in the model. Still, it remains unclear if such empirical approaches are sufficient to remove bias.
<!-- Despite its success, LD score regression has lim-itations. First, no constraint is applied to heritability estimates, therefore it is possible to get non-interpretable estimates that are below 0 or above 1. Second, variance of estimators is cal- culated using a resampling method, the blockwise jackknife. For moderate GWAS sample size, this approach often leads to very wide confidence intervals for enrichment, thereby reduc- ing statistical power. Finally, z-scores for SNPs in LD may be strongly correlated. LD score regression reduces the impact of dependence among data points by using specially designed weights in weighted least squares estimation and only including HapMap SNPs in the model. Still, it remains unclear if such empirical approaches are sufficient to remove bias. -->

--------------------

## Reference

1. B. K. Bulik-Sullivan et al., LD Score regression distinguishes confounding from polygenicity in genome-wide association studies, Nat Genet, vol. 47, no. 3, Art. no. 3, Mar. 2015, doi: 10.1038/ng.3211.

2. https://github.com/YangLabHKUST/XMAP/blob/main/R/ldsc.R

3. Handbook of statistical genomics[M]. John Wiley & Sons, 2019.

--------------------

## Cross-trait LDSC

Cross-trait LD score regression is based on the following model:
$$
\begin{array}{ll}
Y_1=X \beta+\varepsilon, & \beta \sim N\left(0, \frac{h_1^2}{M} I\right), \quad \varepsilon \sim N\left(0,\left(1-h_1^2\right) I\right), \\
Y_2=Z \gamma+\delta, & \gamma \sim N\left(0, \frac{h_2^2}{M} I\right), \quad \delta \sim N\left(0,\left(1-h_2^2\right) I\right) .
\end{array}
$$

SNPs' effects on two different traits are assumed to be correlated:
$$
E\left(\beta \gamma^T\right)=\frac{\rho_g}{M} I,
$$
where $\rho_g$ is the genetic covariance parameter that quantifies the genetic sharing between traits $Y_1$ and $Y_2$.

--------------------

### Shared Samples

Importantly, LD score regression allows two GWASs to share a fraction of samples. Without loss of generality, assume the first $n_s$ samples (i.e. rows) in $X$ and $Z$ are identical and the non-genetic effects for shared samples are correlated:
$$
E\left(\varepsilon_i \delta_i\right)=\rho_e, \quad i=1, \ldots, n_s
$$

Then, the cross-trait LD score regression formula is
$$
E\left(\left(z_1\right)_j\left(z_2\right)_j\right)=\frac{\sqrt{n_1 n_2} \rho_g}{M} l_j+\frac{\left(\rho_g+\rho_e\right) n_s}{\sqrt{n_1 n_2}}
$$

Similar to single-trait analysis, regression coefficients can be used to estimate genetic covariance, or a close-related but more interpretable metric - genetic correlation:
$$
\text { corr }=\frac{\rho_g}{\sqrt{h_1^2 h_2^2}}
$$

--------------------

A recent method, GNOVA (Lu et al., 2017a), has extended cross-tissue LD score regression to model annotation-stratified genetic covariance. When $K$ functional annotations are present in the model,
$$
\begin{aligned}
Y_1  =\sum_{i=1}^K X_i \beta_i+\varepsilon,\quad Y_2  =\sum_{i=1}^K Z_i \gamma_i+\chi, \quad
E\left(\beta_i \gamma_i^T\right) & =\frac{\rho_i}{m_i} I .
\end{aligned}
$$

Parameters $\rho_i(i=1, \ldots, K)$ quantify the genetic covariance components for each functional annotation. GNOVA used an estimator based on the method of moments to estimate genetic covariance:
$$
\left(\begin{array}{c}
\hat{\rho}_1 \\
\vdots \\
\hat{\rho}_K
\end{array}\right)=\frac{1}{\sqrt{n_1 n_2}}\left(\begin{array}{ccc}
\frac{1}{m_1 m_1} l_{11} & \cdots & \frac{1}{m_K m_1} l_{K 1} \\
\vdots & \ddots & \vdots \\
\frac{1}{m_1 m_K} l_{1 K} & \cdots & \frac{1}{m_K m_K} l_{K K}
\end{array}\right)^{-1}\left(\begin{array}{c}
\frac{1}{m_1}\left(z_1\right)_1^T\left(z_2\right)_1 \\
\vdots \\
\frac{1}{m_K}\left(z_1\right)_K^T\left(z_2\right)_K
\end{array}\right),
$$
where $l_{i j}$ denotes the sum of all pairwise LD between SNPs in the $i$ th and $j$ th functional annotations; $\left(z_j\right)_i$ denotes the vector of $z$-scores from the $j$ th study for SNPs in the $i$ th functional annotation.

--------------------

To date, LD score regression is not able to estimate annotation-stratified genetic covariance. GNOVA also shows higher statistical power than LD score regression in both simulations and real data (Lu et al., 2017a).

--------------------

END
