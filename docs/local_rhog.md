# Estimating local genetic covariance and correlation (\\(\rho\\)-HESS)

This page describes the steps to estimate local genetic covariance and
correlation from GWAS summary association data. In the first step,
HESS computes the eigenvalues of LD matrices, the squared projections of
GWAS effect size vector onto the eigenvectors of LD matrices for each trait,
and the product of projections of GWAS effect size vectors onto the
eigenvectors. In the second step, HESS uses the output from step 1 to
obtain estimates of local SNP-heritability of each trait. In the third step,
HESS, uses output from step 2 to obtain local genetic covariance estimates
and their standard errors.

## Step 1 - compute eigenvalues, squared projections, and product of projections

### Running the tool

In this step, HESS computes eigenvalues, squared and product of projections of
GWAS effect size vector of each trait onto the eigenvectors of LD matrices.
The following code provides an example of how to perform this step.

```
for chrom in $(seq 22)
do
    python hess.py \
        --local-rhog <summary stats for trait 1> <summary stats for trait 2> \
        --chrom $chrom \
        --bfile <reference panel in PLINK format for the specific chromosome> \
        --partition <genome partition file for the specific chromosome> \
        --out step1
done
```

In the command above, `--local-rhog` tells HESS to estimate local
genetic covariance, and is used to specify GWAS summary statistics data for
two traits; `--chrom` is used to specify the chromosome number; `--bfile` is
used to specify the reference panel for the corresponding chromsome;
`--partition` is used to specify the genome partition file; `--out` is used
to specify the prefix of the output file. For input file format, please refer
to [Input Format](https://huwenboshi.github.io/hess-0.5/input_format/).

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: The above for loop can be parallelized by chromosome in case a
computer cluster is available. In addition, users can specify the minimum
minor allele frequency (MAF) of the SNPs used for estimation through the
"--min-maf" flag. The default MAF threshold is 0.05. )
</div>

### Interpreting the output

After executing the command above, 9 files will be created for each
chromosome (i.e. 198 files for all 22 chromosomes in total). The following is
an example output obtained for chromosome 22. The format of the output files
are identical to the output of
[local SNP-heritability analysis](https://huwenboshi.github.io/hess-0.5/local_hsqg/#interpreting-the-output).

The following two files will be used in the 3rd step to estimate local
genetic covariance.

* **step1_chr22.eig.gz** - contains the eigenvalues of LD matrix at each locus
* **step1_chr22.prjprod.gz** - contains the product of projections

The following three files will be used in the 2nd step to estimate local
SNP-heritability of trait 1

* **step1_trait1_chr22.eig.gz** - contains the eigenvalues of LD matrix at each locus
* **step1_trait1_chr22.info.gz** - contains information of each locus
* **step1_trait1_chr22.prjsq.gz** - contains the squared projections

The following three files are the same type of files as the previous three
files, and will be used in the 2nd step to estimate local SNP-heritability
of trait 2

* **step1_trait2_chr22.eig.gz** - contains the eigenvalues of LD matrix of each locus
* **step1_trait2_chr22.info.gz** - contains information of each locus
* **step1_trait2_chr22.prjsq.gz** - contains the squared projections

In addition, a log file **step1_trait2_chr22.log** will be created to document
the details of each step.

## Step 2 - estimate local SNP-heritability of each trait

### Running the tool

In this step, we estimate local SNP-heritability for trait 1 and trait 2 using
output from step 1. The following script provide an example of how to perform
this step. Please see [local SNP-heritability analysis](https://huwenboshi.github.io/hess-0.5/local_hsqg/#step-2-estimate-local-snp-heritability-and-standard-error)
for more detail.

```
# estimate local SNP-heritability for trait 1
python hess.py --prefix step1_trait1 --out step2_trait1

# estimate local SNP-heritability for trait 2
python hess.py --prefix step1_trait2 --out step2_trait2
```

### Interpreting the output

The above command will result in 4 files, 2 for each trait, containing local
SNP-heritability estimates at each locus.

* **step2_trait1.txt** - contains local SNP-heritability estimates at each locus for trait 1
* **step2_trait2.txt** - contains local SNP-heritability estimates at each locus for trait 2

In addition, 2 log files will also be created.

* **step2_trait1.log** - contains debugging information for trait 1
* **step2_trait2.log** - contains debugging information for trait 2

Please see [local SNP-heritability analysis](https://huwenboshi.github.io/hess-0.5/local_hsqg/#interpreting-the-output_1)
for more detail.

## Step 3 - estimate local genetic covariance and standard error

### Estimate phenotypic correlation

\\(\rho\\)-HESS requires phenotypic correlation between a pair of traits to
obtain an unbiased estimates of local genetic covariance. If phenotypic data
of the GWAS is available, we recommend to obtain phenotypic correlation of
the pair of traits by taking the Pearson correlation between the phenotype
values of the two traits.

If individual-level phenotype data are not available, one can obtain an
estimate through
[cross-trait LDSC](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).
The intercept term corresponding to the genetic covariance estimates provides
an approximation of phenotypic correlation. More precisely, the estimated 
phenotypic correlation \\(r\_{pheno}\\) is
\\[
    r\_{pheno} =  \delta \times \sqrt{N\_1 N\_2} / N\_s,
\\]
where \\(\delta\\) is the intercept term, \\(N\_1\\) and \\(N\_2\\) sample
size for the two GWAS, and \\(N_s\\) number of shared samples between the
two GWASs.

We provide a simple script (`misc/estimate_phenocor.py`) for obtaining 
phenotypic correlation from cross-trait LDSC log files.
```
python misc/estimate_phenocor.py \
    --ldsc-log <cross-trait LDSC log file> \
    --n1 <sample size for GWAS 1> --n2 <sample size for GWAS 2> \
    --ns <number of shared samples>
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: If there is no sample overlap between the two GWASs, then
one does not need to estimate phenotypic correlation. This is because bias
in local genetic covariance estimate are caused by environmental covariance
coming from overlapping GWAS samples. And one needs to know the phenotypic
correlation to infer environmental covariance. When there is no sample
overlap, there is no need to correct for bias caused by environmental
covariance. )
</div>
<div style="background-color:rgba(240, 128, 128, 0.2);">
( <b>Note</b>: The cross-trait LDSC intercept here should correspond to
genetic covariance and not SNP-heritability. )
</div>

### Running the tool

The following script combines output from step 1 and step 2 to obtain local
genetic covariance estimates.

```
python hess.py \
    --prefix step1 \
    --local-hsqg-est step2_trait1.txt step2_trait2.txt \
    --num-shared <number of overlapping samples in the two GWASs> \
    --pheno-cor <phenotypic correlation between the two traits> \
    --out step3
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: When "--num-shared" is set to zero, "--pheno-cor" can be
set to any value (e.g. 0.0) and the result will not be affected. Also, note
that no for loop is required here. &rho;-HESS will automatically look
for output from all chromosomes.)
</div>

#### Other available flags

* `--max-num-eig`: Specifies the number of eigenvectors to use in the
truncated-SVD for inverting the LD matrix, default 50
* `--min-eigval`: Specifies the minimum eigenvalue cut off in
truncated-SVD, default 1.0
* `--reinflate-lambda-gc`: Genomic control factors that have been applied to
the summary statistics, default 1.0 for both
* `--gwse-thres`: A threshold to cap the standard error of total
SNP-heritability estimate

### Interpreting the output

After step3, 2 files will be created. These include
   
* **step3.txt** - contains the local genetic covariance estimates across
the entire genome. The columns are chromosome, start and end positions of
the locus, number of SNPs in the locus, number of eigenvectors used in the
truncated-SVD, local genetic covariance estimates, and variance estimates.
We also report the standard error and test statistics.

```
chr  start     end         num_snp k    local_rhog    var          se           z            p
1    10583     1892607     1286    50   1.3906e-06    2.1507e-09   4.6375e-05   0.029985     0.97608
1    1892607   3582736     3045    50   -3.2351e-05   4.0829e-09   6.3898e-05   -0.50629     0.61265
1    3582736   4380811     1622    50   0.00011446    2.9594e-09   5.44e-05     2.104        0.035379
1    4380811   5913893     3790    50   -1.898e-06    3.3276e-09   5.7685e-05   -0.032903    0.97375
...  ...       ...         ...     ...  ...           ...          ...          ...          ...
22   46470495  47596318    2444    50   -1.8303e-05   2.7353e-09   5.23e-05     -0.34997     0.72636
22   47596318  48903703    2997    50   -4.6613e-06   3.1558e-09   5.6176e-05   -0.082977    0.93387
22   48903703  49824534    3773    50   -1.261e-06    3.2769e-09   5.7245e-05   -0.022028    0.98243
22   49824534  51243298    2789    50   6.7939e-05    3.8538e-09   6.2079e-05   1.0944       0.27378
```

* **step3.log** - contains estimates of genome-wide genetic covariance and
correlation and information useful for debugging

```
[INFO] Command started at: Wed, 04 Oct 2017 00:29:13
[INFO] Command issued:
[INFO] Total SNP-heritability of trait 1: ...
[INFO] Total SNP-heritability of trait 2: ...
                    ...
[INFO] Genome-wide genetic covariance estimate: ...
[INFO] Genome-wide genetic correlation estimate: ...
[INFO] Command finished at: ...
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: We estimate the standard error of genome-wide genetic correlation
through jackknife.)
</div>
