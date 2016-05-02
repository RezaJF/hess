## HESS (Heritability Estimation from Summary Statistics)

---

#### Overview

HESS estimates the amount of variance in trait explained by typed SNPs at
each single locus on the genome (local SNP-heritability) from GWAS summary
statistics, while accounting for linkage disequilibrium (LD).

---

#### File format

HESS requires as input: <br/>
1. GWAS summary statistics <br/>
2. Reference panel matching the GWAS population <br/>
3. bed files specifying start and end positions of each locus

###### Summary statistics

rsID pos A0 A1 Z-score N <br/>
rs1000 29321 G A -1.6434 91021 <br/>
rs1001 29478 T C -0.0152 89834 <br/>
rs1002 30500 G A 0.7238 95831 <br/>

###### Reference panel

Can be downloaded
[here](https://drive.google.com/open?id=0B0OmLzMQAvWqT3pnTUhtaTBKbDA).

###### Partition file (bed format)

Can be downloaded [here](https://bitbucket.org/nygcresearch/ldetect-data/src).

---

#### Pipeline

HESS estimates local heritability in 2 steps. In step 1, HESS computes
the eigenvalues of LD matrices, and the squared projections of GWAS effect
size vector onto the eigenvectors of LD matrices. In step 2, HESS computes
local heritability estimates and their standard errors, using results
from step 1.

###### Step 1
```{r, engine='sh', count_lines}
# can be parallelized
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

###### Step 2
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
