## HESS (Heritability Estimation from Summary Statistics v0.1-beta)

---

#### Overview

HESS estimates the amount of variance in trait explained by typed SNPs at
each single locus on the genome (local SNP-heritability) from GWAS summary
statistics, while accounting for linkage disequilibrium (LD).

---

#### <a name="input_file_format"></a> Input file format

HESS requires as input
(1) GWAS summary statistics
(2) reference panel matching the GWAS population
(3) bed files specifying start and end positions of each locus.

###### Summary statistics

To improve computational efficiency and parallelizability, HESS requires that
users split summary statistics into chromosomes and that SNPs are sorted by
position. For each SNP, HESS requires 6 information (in the listed order):
(1) rs ID (2) position (3) reference allele (4) alternative allele
(5) Z-score (6) sample size. HESS internally filters out strand-ambiguous
SNPs and flips signs of Z-scores based on alleles in the reference panel.
However, user awareness of these details are highly recommended. The
following is an example of summary statistics file.

```
rsID pos A0 A1 Z-score N
rs1000 29321 G A -1.6434 91021
rs1001 29478 T C -0.0152 89834
rs1002 30500 G A 0.7238 95831
```

###### Reference panel

1000 Genomes Project (phase 3) reference panel for SNPs with MAF > 5% in the
EUR population can be downloaded
[here](https://drive.google.com/open?id=0B0OmLzMQAvWqc3FPcVRDWkdvc2c). 

###### Partition file (bed format)

Can be downloaded [here](https://bitbucket.org/nygcresearch/ldetect-data/src).

---

#### Pipeline

HESS estimates local heritability in 2 steps. In step 1, HESS computes
the eigenvalues of LD matrices, and the squared projections of GWAS effect
size vector onto the eigenvectors of LD matrices. In step 2, HESS computes
local SNP heritability estimates and their standard errors, using results
from step 1.

###### Step 1 - compute eigenvalues and projections

In this step, HESS computes the eigenvalues of LD matrices, and the squared
projections of GWAS effect size vector onto the eigenvectors of LD matrices.
The following code snippet illustrates the 1st step of HESS.

```{r, engine='sh', count_lines}
# this for loop can be parallelized, i.e. one CPU for each chromosome
for i in $(seq 22)
do
    python hess.py \
        --chrom $i \
        --zscore-file zscore.chr"$i" \
        --reference-panel refpanel_genotype_chr"$i".gz \
        --legend-file refpanel_legend_chr"$i".gz \
        --partition-file partition_chr"$i".bed
        --output-file step1
done
```

In the command above, `--chrom` specifies the chromosome number;
`--zscore-file` specifies the summary statistics for SNPs in the
corresponding chromosome; `--reference-panel` specifies the genotype file
for the reference panel; `--legend-file` specifies the legend file for the
reference panel; `--partition-file` specifies start and end positions
of the loci; `--output-file` specifies the prefix of the output for step 1.
For input file format, please refer to
[Input file format](#input_file_format).

After executing the command above, 4 files will be created for each
chromosome, taking up ~10MB of space for the entire genome. Here's an example
obtained for chromosome 22.

* step1\_chr22.info.gz - contains the information of each locus (start and
   end positions,  number of SNPs, rank of LD matrices, sample size)
```
16050408        17674294        371     274     91273
17674295        18296087        419     306     89182
18296088        19912357        947     502     90231
  ...             ...           ...     ...      ...
```
* step1\_chr22.eig.gz - contains the positive eigenvalues of LD matrix at
each locus, one line per locus
```
39.31792281  31.23990243  23.81549256  23.47296559  20.45343550  ...
48.73186142  26.95692375  25.32769526  22.11750791  20.55766423  ...
82.58157342  67.42588424  59.52766188  43.10471854  32.15181631  ...
    ...          ...          ...          ...          ...
```
* step1\_chr22.prjsq.gz - contains the squared projections of effect
size vector onto the eigenvectors of LD matrix at each locus, one
line per locus
```
0.00008940  0.00001401  0.00013805  0.00009906  0.00007841  ...
0.00054948  0.00001756  0.00008532  0.00002303  0.00004706  ...
0.00008693  0.00005737  0.00070234  0.00008411  0.00004001  ...
   ...          ...        ...         ...         ...
```
* step1\_chr22.log - contains logging information (e.g. number of SNPs,
number of SNPs filtered, etc.)
```
Command started at: ...
Command issued: hess.py ...
Number of SNPs in reference panel: ...
Number of SNPs in Z-score file: ...
Number of SNPs in Z-score file after filtering: ...
Number of loci in partition file: ...
Command finished at: ...
```

###### Step 2 - compute local SNP heritability
In this step, HESS uses results from step 1 (step1\_chr22.info.gz,
step1\_chr22.eig.gz, step1\_chr22.prjsq.gz) to compute local SNP
heritability estimates and their standard error.

```{r, engine='sh', count_lines}
python hess.py \
    --prefix step1 \
    --k 50 \
    --out step2.txt
```

In the command above, `--prefix` specifies prefix of the files generated
during step 1, "step1", in this case; `--k`, default at 50, specifies the
maximum number of eigenvectors to use in estimating local SNP heritability;
`--out` specifies the name of the output file. 

After executing the command above, 2 files will be created.

* step2.txt - contains local SNP heritability estimates (including
chromosome number, locus start position, locus end position, number of SNPs
in locus, number of eigenvectors used, local SNP heritability, variance)
```
chr     start   end     num_snp k       local_h2g       var
1       10583   1892606 158     24      0.0001786340    0.000000011374
1       1892607 3582735 814     40      0.0004164805    0.000000039661
1       3582736 4380810 558     40      0.0001844619    0.000000027595
1       4380811 5913892 1879    40      0.0000738749    0.000000032164
...       ...     ...    ...    ...       ...             ...
```
* step2.txt.log - contains logging information (e.g. estimated genomic
control factor, total SNP heritability, etc.)
```
Command started at: ...
Command issued: ...
Number of loci from step 1: ...
Total number of SNPs: ...
Using lambda gc: ...
Estimated total h2g: ...
Command finished at: ...
```

###### Additional flags for step 2

For step 2, HESS has 4 additional flags:
* `--lambda_gc` allows users to specify their own genomic control factor to
re-inflate the summary statistics, if not specified, HESS will estimates
the genomic control factor from data
* `--tot-h2g <h2g>,<s.e.>` allows users to specify total SNP heritability
of the trait
* `--sense-threshold-joint` default at 2.0, allows users to control standard
error of the estimates when total SNP heritability is not known, the smaller
the threshold, the smaller the standard error (at the cost of downward bias)
* `--sense-threshold-indep` default at 0.5, allows users to control standard
error of the estimates when total SNP heritability is available, the smaller
the threshold, the smaller the standard error (at the cost of downward bias)

---

#### Contact

Please contact Huwenbo Shi (shihuwenbo\_AT\_ucla.edu) for questions
related to HESS.

---

#### Reference

Manuscript describing HESS can be found at the
[preprint](http://biorxiv.org/content/early/2016/01/04/035907).
