# Estimating local SNP-heritability

This page describes the steps to estimate local SNP-heritability from GWAS
summary association data. In the first step, HESS computes the eigenvalues of
LD matrices, and the squared projections of GWAS effect size vector onto the
eigenvectors of LD matrices. In the second step, HESS uses the output from
step 1 to obtain estimates of local SNP-heritability and their standard errors.

## Step 1 - compute eigenvalues and squared projections

### Running the tool

In this step, HESS computes the eigenvalues, and the squared projections of
GWAS effect size vector onto the eigenvectors of LD matrices. The following
code provides an example of how to perform this tep.

```
for chrom in $(seq 22)
do
    python hess.py \
        --local-hsqg <summary statistics> \
        --chrom $chrom \
        --bfile <reference panel in PLINK format for the specific chromosome> \
        --partition <genome partition file for the specific chromosome> \
        --out step1
done
```

In the command above, `--local-hsqg` tells HESS to estimate local
SNP-heritability, and is used to specify the GWAS summary statistics data;
`--chrom` is used to specify the chromosome number; `--bfile` is used to
specify the reference panel for the corresponding chromsome; `--partition` is
used to specify the genome partition file; `--out` is used to specify the
prefix of the output file. For input file format, please refer to
[Input Format](https://huwenboshi.github.io/hess-0.5/input_format/).

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: The above for loop can be parallelized by chromosome in case a
computer cluster is available. In addition, users can specify the minimum
minor allele frequency (MAF) of the SNPs used for estimation through the
"--min-maf" flag. The default MAF threshold is 0.05. )
</div>

### Interpreting the output

After executing the command above, 4 files will be created for each
chromosome (i.e. 88 files for all 22 chromosomes in total), taking up ~10MB
of space for the entire genome. The following is an example output obtained
for chromosome 22.

* **step1_chr22.info.gz** - contains the locus information, including start and end
positions, number of SNPs, rank of LD matrices, and sample size

```
16050408        17674294        371     274     91273
17674295        18296087        419     306     89182
18296088        19912357        947     502     90231
...             ...             ...     ...     ...
```

* **step1_chr22.eig.gz** - contains the eigenvalues of LD matrix (up to the rank
of the LD matrix) at each locus, one line per locus

```
39.31792281  31.23990243  23.81549256  23.47296559  20.45343550  ...
48.73186142  26.95692375  25.32769526  22.11750791  20.55766423  ...
82.58157342  67.42588424  59.52766188  43.10471854  32.15181631  ...
...          ...          ...          ...          ...          ...
```

* **step1_chr22.prjsq.gz** - contains the squared projections of GWAS effect size
vector onto the eigenvectors of LD matrix (up to the rank of the LD matrix) at
each locus, one line per locus

```
0.00008940  0.00001401  0.00013805  0.00009906  0.00007841  ...
0.00054948  0.00001756  0.00008532  0.00002303  0.00004706  ...
0.00008693  0.00005737  0.00070234  0.00008411  0.00004001  ...
...         ...         ...         ...         ...         ...
```

* **step1_chr22.log** - contains helpful information for debugging, including
number of SNPs, number of SNPs filtered, etc.

```
[INFO] Command started at: ...
[INFO] Command issued: ...
[INFO] Loaded ... partitins in ...
[INFO] Average window size is ...
[INFO] ... SNPs read from reference panel
                     ...
[INFO] Command finished at: ...
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: Log files are very useful in solving issues with running the
software. Please include the log file in the email when you report any issue
with HESS. )
</div>

## Step 2 - estimate local SNP-heritability and standard error

### Running the tool (without using total SNP-heritability)

In the second step, HESS uses output from step 1 across all chromosomes
(step1_chr{1..22}.info.gz, step1_chr{1..22}.eig.gz, step1_chr{1..22}.prjsq.gz)
to obtain local SNP-heritability estimates and their standard errors. The
following command automatically looks for output from step 1 across all
chromosomes with the prefix "step1" to obtain local SNP-heritability estimates.

```
python hess.py --prefix step1 --out step2
```

<div style="background-color:rgba(240, 128, 128, 0.2);">
( <b>Note</b>: This is the recommended way of running HESS. However, if sample
size of the GWAS is relatively small and total SNP-heritability is known, e.g.
estimated from individual-level data, then one may use the script in the
next section to obtain local SNP-heritability estimates. Also note that no
for loop is needed here. )
</div>

### Note on re-inflating \\(\lambda_{GC}\\)

Most GWAS summary stats are corrected for genomic control factor
\\(\lambda_{GC}\\). This could result in a downward bias in the estimated
SNP-heritability. If the GWAS summary stats has been corrected for
\\(\lambda_{GC}\\), it is recommended to use the following code to perform
step 2.

```
python hess.py --prefix step1 --reinflate-lambda-gc <lambda gc to reinflate> --out step2
```

### Running the tool (using total SNP-heritability)

If sample size of the GWAS is relatively small, the standard error of total
SNP-heritability estimates of HESS can be large, which could then lead to
large standard error in the local estimates. However, if an estimate of total
SNP-heritability is available, then one may use the following script to
decompose the total SNP-heritability into each locus.

```
python hess.py --prefix step1 --tot-hsqg <total SNP-heritability> <SE> --out step2
```

#### Other available flags

* `--max-num-eig`: Specifies the number of eigenvectors to use in the
truncated-SVD for inverting the LD matrix, default 50
* `--min-eigval`: Specifies the minimum eigenvalue cut off in
truncated-SVD, default 1.0
* `--reinflate-lambda-gc`: Specifies the genomic control factor that has been
applied to the summary statistics, default 1.0
* `--gwse-thres`: A threshold to cap the standard error of total
SNP-heritability estimate

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: We recommend to run HESS on GWAS summary statistics that has
not been applied the genomic control factor, which could results in downward
bias. )
</div>

We provide a simple script (`misc/estimate_lambdagc.py`) to infer the genomic
control factor to reinflate the summary statistics if it is unknown. The
following script should be executed after step 1 completes.
```
python misc/estimate_lambdagc.py --prefix step1
```

### Interpreting the output

After step 2, 2 files will be created. These include

* **step2.txt**: contains local SNP-heritability estimates for all loci across all
chromosomes. The columns are chromosome number, locus start position, locus
end position, number of SNPs in locus, number of eigenvectors used
in truncated-SVD, estimates of local SNP-heritability, estimated variance
of the local SNP-heritability. We also report standard error and test
statistics.

```
chr start       end         num_snp k   local_h2g   var         se          z           p
1   10583       1892607     1286    50  0.0001293   1.2819e-08  0.00011322  1.142       0.12673
1   1892607     3582736     3045    50  0.0003232   2.7102e-08  0.00016463  1.9632      0.024809
1   3582736     4380811     1622    50  6.7587e-05  1.282e-08   0.00011323  0.59692     0.27528
1   4380811     5913893     3790    50  3.1445e-05  2.6111e-08  0.00016159  0.1946      0.42285
... ...         ...         ...     ... ...         ...         ...         ...         ...
22  46470495    47596318    2444    50  0.00018004  2.0125e-08  0.00014186  1.2691      0.1022
22  47596318    48903703    2997    50  8.389e-05   2.4784e-08  0.00015743  0.53287     0.29706
22  48903703    49824534    3773    50  6.9487e-06  2.5931e-08  0.00016103  0.043151    0.48279
22  49824534    51243298    2789    50  0.00015227  2.4024e-08  0.000155    0.98241     0.16295
```

* **step2.log**: contains useful information for debugging

```
[INFO] Command started at: Tue, 03 Oct 2017 21:01:34
[INFO] Command issued:
    /u/project/pasaniuc/shihuwen/software/hess/hess.py \
        --prefix ... \
        --out ...
[INFO] Loaded results for ... loci from step 1
[INFO] Using 4381458 SNPs with average sample size ...
[INFO] Re-inflate the summary statistics with lambda_gc: ...
[INFO] Total SNP-heritability estimate: ...
[INFO] Command finished at: ...
```
