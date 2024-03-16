Simulation data was generated using Elison, Westone's tool (https://github.com/westonelison/CSE_284_Finemapping),
which is a wrapper of the simGWAS tool (https://github.com/chr1swallace/simGWAS).

In the metadata folder, we have info on the selected causal SNPs for each simulation. LD files need to be generated
by following instruction on Elison's tool (using plink).

Batch code are provided to easily generate many simulation, and run them using either PyFM or 
FINEMAP.
- Please change the file paths when using the batch codes

Basically, encode your PyFM run parameters in `batch_run.sh`, and in the `job_assigner.sh`
specify input and output DIR. Call on `job_assigner.sh` to run batch.