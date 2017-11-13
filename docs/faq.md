# Frequently Asked Questions

This page addresses some frequently asked questions from users.

## What is the sample size requirement for HESS and \\(\rho\\)-HESS?

We recommend to apply HESS and \\(\rho\\)-HESS on GWAS with sample size greater
than 50,000. HESS and \\(\rho\\)-HESS tend to be downwardly biased at smaller
GWAS since fewer eigenvectors are necessary in the truncated-SVD
regularization to obtain a stable estimate.

On the other hand, if genome-wide SNP-heritability is available for the GWAS
(e.g. estimated from individual-level data), then it may still be possible to
obtain the local estimates ( see [here](http://huwenboshi.github.io/hess/local_hsqg/#running-the-tool-using-total-snp-heritability)).
