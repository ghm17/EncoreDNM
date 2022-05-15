# EncoreDNM
EncoreDNM is a powerful method to estimate concordant de novo mutation (DNM) associations between disorders. 

## Before starting
EncoreDNM is built upon `R`. The R-package `mvtnorm` is required to generate multivariate normal random vectors, and the R-package `snowfall` is required for parallel computing. 

Suppose we would like to estimate enrichment correlation between two disorders, e.g. developmental disorder (DD) and autism spectrum disorder (ASD). We need to prepare the following files:
* Two DNM summary data files: We have prepared the example data for you in `./example/DD.txt` and `./example/ASD.txt`. The DNM summary data files need to be transformed into the standard format by yourself. The first few lines should look like this: 

      head DD.txt
      
      chr gene lof Dmis Tmis syn
      10 A1CF 0 0 0 0
      10 ABCC2 0 0 5 1
      10 ABI1 0 0 0 0

The first two columns denote chromosome and gene name, the following columns denote different variant classes (must be consistent with variant classes in the mutability table). 
* Mutability table: We have prepared the example data for you in `./example/mutability_table.txt`. The mutability table also needs to be transformed into the standard format by yourself. The first few lines should look like this: 

      head mutability_table.txt

      chr gene lof Dmis Tmis syn
      10 A1CF 3.89954925e-06 1.099817e-06 1.5858302e-05 7.031111e-06
      10 ABCC2 6.98155025e-06 2.1032377e-05 2.0224272e-05 1.725152e-05
      10 ABI1 2.69662075e-06 6.712038e-06 7.31911e-06 5.48011e-06

## Single disorder analysis
To perform single-disorder EncoreDNM analysis, run the following commands:
```bash
Rscript EncoreDNM_single_disorders.R \
--dat DNM_FILE \
--N N \
--mut MUTABILITY TABLE \
--out OUTFILE \
--n_cores N_CORE
```

Here the inputs are
* `DNM_FILE` specifies the input DNM summary data file.
* `N` specifies the trio sample size.
* `MUTABILITY_TABLE` specifies the mutability table file.
* `OUTFILE` specifies the full path of the output file.
* `n_cores` specifies the number of cores used in parallel computing.

A concrete example is shown below
```bash
cd ./EncoreDNM

Rscript ./EncoreDNM_single_disorders.R \
--dat ./example/DD.txt \
--N 31058 \
--mut ./example/mutability_table.txt \
--out ./example/DD_single_disorder_output.txt \
--n_cores 20
```

### Output for single disorder analysis
After running the above step, EncoreDNM outputs a text file `DD_single_disorder_output.txt`, with the following columns:

* `variant_class`: e.g. lof, Dmis, Tmis, syn. 
* `beta_est`: Estimated value for elevation parameter beta. 
* `beta_se`: Standard error for elevation parameter beta. 
* `beta_p`: Wald test p-value for elevation parameter beta. 
* `sigma_est`: Estimated value for dispersion parameter sigma. 
* `sigma_se`: Standard error for dispersion parameter sigma. 
* `sigma_p`: Wald test p-value for dispersion parameter sigma. 
* `deviance_null`: The deviance between the mixed-effects Poisson model and the fixed-effects Poisson model.
* `p_null`: The likelihood ratio test p-value between the mixed-effects Poisson model and the fixed-effects Poisson model. 

We have prepared the example output files for DD and ASD in `./example/DD_single_disorder_output.txt` and `./example/ASD_single_disorder_output.txt`. 


## Cross disorders analysis
To perform cross-disorders EncoreDNM analysis, run the following commands:
```bash
Rscript EncoreDNM_cross_disorders.R \
--dat1 DNM_FILE_1 \
--dat2 DNM_FILE_2 \
--N1 N1 \
--N2 N2 \
--dat1_encorednm_single SINGLE_OUTFILE_1 \
--dat2_encorednm_single SINGLE_OUTFILE_2 \
--mut MUTABILITY TABLE \
--out OUTFILE \
--n_cores N_CORE
```

Here the inputs are
* `DNM_FILE_1` and `DNM_FILE_2` specify the input DNM summary data files for two disorders.
* `N1` and `N2` specify the trio sample sizes for two disorders.
* `SINGLE_OUTFILE_1` and `SINGLE_OUTFILE_2` specify the full path of output files in single disorder analysis for two disorders.
* `MUTABILITY_TABLE` specifies the mutability table file.
* `OUTFILE` specifies the full path of the output file.
* `n_cores` specifies the number of cores used in parallel computing.

A concrete example is shown below
```bash
cd ./EncoreDNM

Rscript ./EncoreDNM_cross_disorders.R \
--dat1 ./example/DD.txt \
--dat2 ./example/ASD.txt \
--N1 31058 \
--N2 6430 \
--dat1_encorednm_single ./example/DD_single_disorder_output.txt \
--dat2_encorednm_single ./example/ASD_single_disorder_output.txt \
--mut ./example/mutability_table.txt \
--out ./example/DD_ASD_cross_disorders_output.txt \
--n_cores 20
```

### Output for cross disorders analysis
After running the above step, EncoreDNM outputs a text file `DD_ASD_cross_disorders_output.txt`, with the following columns:

* `variant_class`: e.g. lof, Dmis, Tmis, syn. 
* `beta1_est`: Estimated value for elevation parameter beta in the first disorder. 
* `beta2_est`: Estimated value for elevation parameter beta in the second disorder. 
* `sigma1_est`: Estimated value for dispersion parameter sigma in the first disorder. 
* `sigma2_est`: Estimated value for dispersion parameter sigma in the second disorder. 
* `rho_est`: Estimated value for enrichment correlation rho. 
* `beta1_se`: Standard error for elevation parameter beta in the first disorder. 
* `beta2_se`: Standard error for elevation parameter beta in the second disorder. 
* `sigma1_se`: Standard error for dispersion parameter sigma in the first disorder. 
* `sigma2_se`: Standard error for dispersion parameter sigma in the second disorder. 
* `rho_se`: Standard error for enrichment correlation rho. 
* `beta1_p`: Wald test p-value for elevation parameter beta in the first disorder. 
* `beta2_p`: Wald test p-value for elevation parameter beta in the second disorder. 
* `sigma1_p`: Wald test p-value for dispersion parameter sigma in the first disorder. 
* `sigma2_p`: Wald test p-value for dispersion parameter sigma in the second disorder. 
* `rho_p`: Wald test p-value for enrichment correlation rho. 

We have prepared the example output file in `./example/DD_ASD_cross_disorders_output.txt`. 
