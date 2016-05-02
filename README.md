## HESS (Heritability Estimation from Summary Statistics)

---

#### Overview

HESS estimates the amount of variance in trait explained by typed SNPs at
each single locus on the genome (local SNP-heritability) from GWAS summary
statistics, while accounting for linkage disequilibrium (LD).

---

#### <a name="input_file_format"></a> Input file format

HESS requires as input
(1) GWAS summary statistics
(2) Reference panel matching the GWAS population
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
[here](https://github.com/huwenboshi/hess). 

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
        --output-file step1 \
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
chromosome:

1. step1\_chr22.info.gz - containing the information of each locus (start and
   end positions,  number of SNPs, rank of LD matrices, sample size)

```
16050408        17674294        371     274     91273
17674295        18296087        419     306     89182
18296088        19912357        947     502     90231
  ...             ...           ...     ...      ...
```

2. step1\_chr22.eig.gz - containing the positive eigenvalues of LD matrix at
each locus, one line per locus

```
39.31792281  31.23990243  23.81549256  23.47296559  20.45343550  ...
48.73186142  26.95692375  25.32769526  22.11750791  20.55766423  ...
82.58157342  67.42588424  59.52766188  43.10471854  32.15181631  ...
    ...          ...          ...          ...          ...
```

3. step1\_chr22.prjsq.gz - containing the squared projections of effect
size vector onto the eigenvectors of LD matrix at each locus, one
line per locus

```
0.00008940  0.00001401  0.00013805  0.00009906  0.00007841  ...
0.00054948  0.00001756  0.00008532  0.00002303  0.00004706  ...
0.00008693  0.00005737  0.00070234  0.00008411  0.00004001  ...
   ...          ...        ...         ...         ...
```

4. step1\_chr22.log

###### Step 2 - compute local SNP heritability

```{r, engine='sh', count_lines}
python hess.py \
    --prefix step1 \
    --k 50 \
    --out step2.txt
```

---

#### Contact

Please contact Huwenbo Shi (shihuwenbo\_AT\_ucla.edu) for questions
related to HESS.

---

#### Reference

Manuscript describing HESS can be found at the
[preprint](http://biorxiv.org/content/early/2016/01/04/035907).
