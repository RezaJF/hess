# Frequently Asked Questions

This page addresses some frequently asked questions from users.
Most of the question that applies to HESS also applies to \\(\rho\\)-HESS.

## What is the sample size requirement for HESS?

We recommend to apply HESS on GWAS with sample size greater
than 50,000. HESS tends to be downwardly biased at smaller
GWAS since fewer eigenvectors are necessary in the truncated-SVD
regularization to obtain a stable estimate.

On the other hand, if genome-wide SNP-heritability is available for the GWAS
(e.g. estimated from individual-level data), then it may still be possible to
obtain the local estimates
(see [here](http://huwenboshi.github.io/hess/local_hsqg/#running-the-tool-using-total-snp-heritability)).

## Why is it necessary to re-inflate genomic control factor?

Most GWAS studies apply genomic control factor (\\(\lambda_{GC}\\)) on the
association statistics. This could result in a downward bias in the estimated
heritability since the second step of HESS involves an estimation of the
environmental effect variance (\\(\sigma_e^2\\)) -- with \\(\lambda_{GC}\\)
corrected GWAS summary stats, HESS tends to overestimate \\(\sigma_e^2\\),
resulting in downward bias in local SNP-heritability.

## Why do I get negative variance estimates for local SNP-heritability?

The variance estimates for local SNP-heritbility (genetic covariance) is a
random variable that is not constrained to be positive. For some loci, the
local SNP-heritability variance estimates may be negative.

Usually, negative variance estimates are caused by relatively small sample
size of the GWAS.
