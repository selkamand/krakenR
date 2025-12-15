## Simulated Data

Paired reads fastqs were simulated from reference *E. coli* (GCF_000005845.2_ASM584v2) and *H. pylori* (GCF_000307795.1_ASM30779v1) 
using the micritebenchmark tool.

They were classified using Kraken2 against the k2_standard_08gb_20240605 database.

```bash
db=$1
kraken2 --db ${db} --output --report ecoli_0_hpylori_1000.kreport --paired ecoli_0_hpylori_1000.R1.fq ecoli_0_hpylori_1000.R2.fq
kraken2 --db ${db} --output --report ecoli_1000_hpylori_0.kreport --paired ecoli_1000_hpylori_0.R1.fq ecoli_1000_hpylori_0.R2.fq
kraken2 --db ${db} --output --report ecoli_100_hpylori_900.kreport --paired ecoli_100_hpylori_900.R1.fq ecoli_100_hpylori_900.R2.fq
kraken2 --db ${db} --output --report ecoli_500_hpylori_500.kreport --paired ecoli_500_hpylori_500.R1.fq ecoli_500_hpylori_500.R2.fq
kraken2 --db ${db} --output --report ecoli_900_hpylori_100.kreport --paired ecoli_900_hpylori_100.R1.fq ecoli_900_hpylori_100.R2.fq
```
