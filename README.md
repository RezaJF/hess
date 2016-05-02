## HESS (Heritability Estimation from Summary Statistics)

#### Overview

HESS estimates the amount of variance in trait explained by typed SNPs at
each single loci on the genome (local SNP-heritability) from GWAS summary
statics, while accounting for linkage disequilibrium (LD).

#### File format

HESS requires as input: <br/>
1. GWAS summary statistics <br/>
2. Reference panel matching the GWAS population <br/>
3. bed files specifying how the genome is partitioned

###### Summary statistics

rsID pos A0 A1 Z-score N <br/>
rs1000 29321 G A -1.6434 91021 <br/>
rs1001 29478 T C -0.0152 89834 <br/>
rs1002 30500 G A 0.7238 95831 <br/>

###### Reference panel

Can be downloaded [here](https://drive.google.com/open?id=0B0OmLzMQAvWqT3pnTUhtaTBKbDA).

#### Pipeline

HESS estimates local heritability in 2 step2. In step 1, HESS computes
the eigenvalues of LD matrices, and the squared projections of GWAS effect
size vector onto the corresponding eigenvectors of LD matrices. In step 2,
HESS obtain local heritability estimates and their standard errors, using
results from step 1.

###### Step 1

df

###### Step 2
