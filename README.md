# PyFM

**PyFM** uses Bayes Factor-based modeling of different configurations of given variants for fine mapping. The framework is essentialy based on CaviarBF [[1]](#1). But like FINEMAP [[2]](#2) our method further uses Shotgun Stochastic Search to explore the configurations-space more intelligently.

We copied the example data from CaviarBF and FineMap into the `example` subdirectory.

*We prepared 2 example runs from command-line. Example 1's PIP and Rho match up
nearly exactly with the CaviarBF results (difference arise from the Cholesky
Decomposition round-offs). Example 2's Rho values don't match up quite
well, though the SNP selections are the same (under investigation).*

**Example 1**
```shell
python src/fm.py \
  -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
  -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
  -o pyfm_results \
  -n 471 -c 2 -t 0 -e 0.1 -a 1.6
```

**Example 2**
```shell
python src/fm.py \
  -z example/finemap_examples/finemap_data.z \
  -r example/finemap_examples/finemap_data.ld \
  -o pyfm_results \
  -n 5363 -c 2 -t 0 -e 0.1 -a 1.6
```

**Example 3 (Shotgun Stochastic Search)**
```shell
python src/fm.py \
  -z example/finemap_examples/finemap_data.z \
  -r example/finemap_examples/finemap_data.ld \
  -o pyfm_results \
  --configs-method SSSConfigurations \
  -n 5363 -c 2 -t 0 -e 0.1 -a 1.6
```


#### Arguments:
For further explanation of the arguments, check CaviarBF arguments since they are very similar:
 [CaviarBF Manual](CaviarBF_Manual.pdf) 
(if does not open on webpage, please download the PDF file from above and open locally)

| Argument | Type    | Description                                                    |
|----------|---------|----------------------------------------------------------------|
| `-z`     | FILE    | zfile                                                          |
| `-r`     | FILE    | rfile containing pairwise LD correlation matrix                |
| `-o`     | DIR     | output directory                                               |
| `-n`     | int     | (default: 0) sample number                                     |
| `-c`     | int     | (default: 3) maximum number of causal variants considered      |
| `-t`     | int     | (default: 0) prior type                                        |
| `-e`     | float   | (default: 0) epsilon, noise factor added to correlation matrix |
| `-a`     | [float] | (default: 0.1 0.2 0.4 0.8 1.6) priors to run on                |
| `-p`     | float   | (default: 1; i.e. reports all) rho cutoff to be used           |


## Compare results with CaviarBF

Download CaviarBF from by 
```
git clone https://bitbucket.org/Wenan/caviarbf.git
```

Install by running the Makefile using 
```
make
```

Run CaviarBF on the example files. Note: CaviarBF has two modules, 
`caviarbf` and `model_search`

`caviarbf` builds the Bayes factors for each SNP, and `model_search` find the best
model of SNP combinations based on exhaustive/greedy search. Again, please refer
to [CaviarBF Manual](CaviarBF_Manual.pdf).

We had prepared the example code for the following file structures
```
caviarbf/
    src/
    Makefile
    caviarbf (executable)
    ...
PyFM/
    src/
        fm.py (executable)
        ...
    example/
    pyfm_results/
    caviarbf_results/
    ...
```

### CaviarBF: caviarbf module
*Don't know why sometimes it returns killed: 9, but just run it again.*

Similar to PyFM, but `-o` is a PATH to the FILE, instead of PATH to DIR

**Example 1a**
```shell
../caviarbf/caviarbf \
  -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
  -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
  -o caviarbf_results/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.bf \
  -n 471 -c 2 -t 0 -e 0.1 -a 1.6
```

**Example 2a**
```shell
../caviarbf/caviarbf \
  -z example/finemap_examples/finemap_data.z \
  -r example/finemap_examples/finemap_data.ld \
  -o caviarbf_results/finemap_example.bf \
  -n 5363 -c 2 -t 0 -e 0.1 -a 1.6
```

### CaviarBF: build_model module

`-e` for exhaustive search, and `-s` for greedy stepwise search

**Example 1b**
```shell
../caviarbf/model_search \
  -i caviarbf_results/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.bf \
  -o caviarbf_results/test_stepwise \
  -s -m 237 -p 0 > caviarbf_results/log.txt
```

**Example 2b**
```shell
../caviarbf/model_search \
  -i caviarbf_results/finemap_example.bf \
  -o caviarbf_results/test_finemap_stepwise \
  -s -m 55 -p 0 > caviarbf_results/log.txt
```


## References
* <a id="1">[1]</a> Chen W, Larrabee BR, Ovsyannikova IG, Kennedy RB, Haralambieva IH, Poland GA, Schaid DJ. Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics. Genetics. 2015 Jul;200(3):719-36. doi: 10.1534/genetics.115.176107. Epub 2015 May 6. PMID: 25948564; PMCID: PMC4512539.
* <a id="2">[2]</a> Christian Benner, Chris C.A. Spencer, Aki S. Havulinna, Veikko Salomaa, Samuli Ripatti, Matti Pirinen, FINEMAP: efficient variable selection using summary data from genome-wide association studies, Bioinformatics, Volume 32, Issue 10, May 2016, Pages 1493–1501, https://doi.org/10.1093/bioinformatics/btw018
