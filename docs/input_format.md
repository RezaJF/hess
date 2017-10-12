# Input Format

This page describes the format of the GWAS summary association data (supplied
through the `--local-hsqg` or `--local-rhog` flag), genome partition file
(supplied through the `--partition` flag) and the reference panel required
by HESS (supplied through the `--bfile` format).

## GWAS summary association data

HESS requires a single file in plain text or gzipped text containing the
following columns to be present in the GWAS summary association data:

* SNP - rs ID of the SNP (e.g. rs62442).
* CHR - Chromosome number of the SNP. This should be a number between 1 and 22.
* BP - Base pair position of the SNP.
* A1 - Effect allele of the SNP. The sign of the Z-score is with respect to this allele.
* A2 - The other allele of the SNP.
* Z - The Z-score of the SNP.
* N - Sample size of the SNP.

For estimating local SNP-heritability, one summary statistics file should be
supplied through the `--local-hsqg` flag. For estimating local genetic
covariance, two summary statistics files should be supplied through
the `--local-rhog` flag.

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>:  This file format is compatible with <a href="https://github.com/bulik/ldsc">LDSC</a>.
However, HESS requires two additional columns than LDSC [CHR and BP].
Therefore, one may not directly run HESS on summary stats in LDSC format
unless these two columns are provided. Other columns [e.g. MAF, INFO, etc.]
may be included in the file, but will not be used. It is also recommended
[although not required] that the summary data files are sorted by their
chromosomes and base pairs. All SNPs with either duplicate ID or position
will be removed before any analysis. )
</div>

<div style="background-color:rgba(240, 128, 128, 0.2);">
( <b>Note</b>: We recommend to run HESS on imputed summary statistics data.
HESS may yield downwardly biased estimates when only genotyped SNPs are used. )
</div>

## Genome partition

The genome partition file should be in bed format, one for each chromosome. We
recommend to use the approximately LD-independent loci provided by Berisa
et al., which can be downloaded [here](https://bitbucket.org/nygcresearch/ldetect-data/src/ac125e47bf7f?at=master).
Please choose the genome partition that corresponds to the GWAS population.
The genome partition file is supplied through the `--partition` flag.

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: In order to obtain unbiased estimates, the partitions should be as
LD-independent as possible. Using a small window may result in LD leakage and
biased estimates. )
</div>

## Reference panel

Reference panels should be in [PLINK format](https://www.cog-genomics.org/plink/2.0/input#bed)
(supplied through the `--bfile` flag).

The following is a list of popular publicly available reference panels.

* [1000 Genomes Project](http://www.internationalgenome.org/data/)
* [UK10K](https://www.uk10k.org/data_access.html)

We provide 1000 Genomes reference panel for Europeans [here](https://ucla.box.com/s/4ya4lxvwjsujt6rrotn3keaiw9dkolld).
All SNPs in this reference panel have minor allele frequency greater than 1%.

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: If a SNP has missing genotypes at some individuals in the
reference panel, HESS imputes the missing value with the population average.
HESS also removes all SNPs with duplicate ID or position. )
</div>
