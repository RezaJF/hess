## HESS (Heritability Estimation from Summary Statistics)

#### Overview

HESS estimates the amount of variance in trait explained by typed SNPs at
each single loci on the genome (local SNP-heritability) from GWAS summary
statics, while accounting for linkage disequilibrium (LD).

#### File format

HESS requires GWAS summary statistics and reference panel matching the GWAS
population as input.

###### Summary statistics

SNP\_rsID SNP\_pos A0 A1 Z-score N <br/>
rs1000 29321 G A -1.6434 91021
rs1001 29478 T C -0.0152 89834
rs1002 30500 G A 0.7238 95831

###### Reference panel

Can be downloaded [here](https://drive.google.com/open?id=0B0OmLzMQAvWqT3pnTUhtaTBKbDA).

#### Pipeline

###### Step 1

###### Step 2
