# Overview
This repository contains a script for formatting a .sync file (Kofler et al., 2011) to 
an input file that will work for the software BayeScan (Foll and Gaggiotti 2008) or 
BayeScEnv (de Villemereuil et al., 2015).

# Example:

- sample.sync: A .sync file with two SNPs (Chr1 position 3 and Chr1 position 9).

The command: 
```
sync2Bscan.py -filename sample.sync -n 50 -N 6 -min_cov 10 -max_cov 500 -min_count 4 -rescale 0 -prefix sample
```
Creates the files:

- sample_all_bscanin.txt: Is the input file formated for BayeScan or BayeScEnv
```
[loci]=2
[populations]=6
[pop]=1
1	68	2	28	40
2	68	2	28	40
[pop]=2
1	74	2	74	0
2	74	2	74	0
[pop]=3
1	61	2	26	35
2	61	2	26	35
...
```

- sample_all_snps.txt: Is a file which maps the SNP numbers in the input file (above)
to their genomic locations as given by the .sync file
```
Chr1	3	1
Chr1	9	2
...
```

# References
* PoPoolation2 (and .sync files): Kofler et al., 2011 Bioinformatics 27:3435-3436
* BayeScan: Foll and Gaggiotti 2008 Genetics 180: 977-993
* BayeScEnv: de Villemereuil et al., 2015 Methods in Ecology and Evolution 6:1248-1258
