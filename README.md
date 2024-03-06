# PyFM

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
  -n 471 -c 2 -t 0 -e 0.1 -a 1.6
```


