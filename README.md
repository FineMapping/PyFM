# PyFM

**PyFM** uses Bayes Factor-based modeling of different configurations of given variants for fine mapping. The framework is essentialy based on CaviarBF [[1]](#1). But like FINEMAP [[2]](#2) our method further uses Shotgun Stochastic Search to explore the configurations-space more intelligently.

We copied the example data from CaviarBF into the `example` subdirectory.

Example run from command-line:
```shell
python src/fm.py \
  -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
  -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
  -o pyfm_results \
  -n 471 -c 2 -t 0 -e 0.1 -a 1.6
```

#### Arguments:
For furtherr explanation of the arguments, check CaviarBF arguments since they are very similar.
| Argument | Type    | Description                                                    |
|----------|---------|----------------------------------------------------------------|
| `-z`     | FILE    | zfile                                                          |
| `-r`     | FILE    | rfile containing pairwise LD correlation matrix                |
| `-o`     | DIR     | output directory                                               |
| `-n`     | int     | (default: 0) sample number                                     |
| `-c`     | int     | (default: 3) maximum number of causal variants considered      |
| `-t`     | int     | (default: 0) prior type                                        |
| `-e`     | float   | (default: 0) epsilon, noise factor added to correlation matrix |
| `-a`     | [float] | (default: 0.1,0.2,0.4) priors to run on; For GWAS, (0.1, 0.2, 0.4) is recommeneded [[3]](#3) For eQTL, (0.1, 0.2, 0.4, 0.8, 1.6) is recommended |
| `-p`     | float   | (default: 1) rho cutoff to be used                             |


## Compare results with CaviarBF

Download CaviarBF from by 
```
git clone https://bitbucket.org/Wenan/caviarbf.git
```

Install by running the Makefile using
```
make
```

Run CaviarBF on the example files
```shell
../caviarbf/caviarbf \
  -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
  -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
  -o caviarbf_results/test \
  -n 471 -c 2 -t 0 -e 0.1 -a 0.1,0.2,0.4,0.8,1.6
```


## References
* <a id="1">[1]</a> Chen W, Larrabee BR, Ovsyannikova IG, Kennedy RB, Haralambieva IH, Poland GA, Schaid DJ. Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics. Genetics. 2015 Jul;200(3):719-36. doi: 10.1534/genetics.115.176107. Epub 2015 May 6. PMID: 25948564; PMCID: PMC4512539.
* <a id="2">[2]</a> Christian Benner, Chris C.A. Spencer, Aki S. Havulinna, Veikko Salomaa, Samuli Ripatti, Matti Pirinen, FINEMAP: efficient variable selection using summary data from genome-wide association studies, Bioinformatics, Volume 32, Issue 10, May 2016, Pages 1493–1501, https://doi.org/10.1093/bioinformatics/btw018
