# EncoreDNM
EncoreDNM is a powerful method to estimate concordant de novo mutation associations between disorders. 

## Before starting
EncoreDNM is built upon `R`, make sure that your R-version is no less than 3.5.0, and the R-package `mvtnorm` is required to generate multivariate normal random vectors. 

## Tutorial
First download EncoreDNM and the corresponding data.
        
    git clone git@github.com:ghm17/EncoreDNM.git

**All the following steps should be carried out under the `./EncoreDNM` directory!** Suppose we would like to estimate enrichment correlation between two disorders, e.g. developmental disorder (DD) and autism spectrum disorder (ASD). We need to prepare the following files:
* Two DNM summary data files: We have prepared the example data for you in the directory `./example/DD.txt` and `./example/ASD.txt`. The DNM summary data files need to be transformed into the standard format by yourself. The first few lines should look like this: 

      head DD.txt
      
      chr gene lof Dmis Tmis syn
      10 A1CF 0 0 0 0
      10 ABCC2 0 0 5 1
      10 ABI1 0 0 0 0
      10 ABLIM1 1 1 1 3
      10 ABRAXAS2 0 0 2 0
      10 ACADSB 0 2 1 1

The first two columns denote chromosome and gene name, the following columns denote different variant classes which can be user-defined (as long as to be consistent with variant classes defined in the mutability table). 
* Mutability table: We have prepared the example data for you in the directory `./example/mutability_table.txt`. The mutability table also needs to be transformed into the standard format by yourself. The first few lines should look like this: 

      head mutability_table.txt

      chr gene lof Dmis Tmis syn
      10 A1CF 3.89954925e-06 1.099817e-06 1.5858302e-05 7.031111e-06
      10 ABCC2 6.98155025e-06 2.1032377e-05 2.0224272e-05 1.725152e-05
      10 ABI1 2.69662075e-06 6.712038e-06 7.31911e-06 5.48011e-06
      10 ABLIM1 4.7877255e-06 6.321814e-06 2.3908341e-05 1.2047648e-05
      10 ABRAXAS2 2.0691095e-06 7.40892857142857e-12 1.222619e-05 5.828581e-06
      10 ACADSB 1.81531575e-06 8.744476e-06 2.301889e-06 4.834414e-06

### Single disorder analysis
To perform single-disorder EncoreDNM analysis, run the following commands:
    
    cd ./EncoreDNM
    Rscript ./EncoreDNM_single_disorder.R DD.txt N mutability_table.txt DD_single_disorder_output.txt

The argument `DD.txt` specifies the input DNM summary data file, the argument `N` specifies the trio sample size which is 31058 for DD here, the argument `mutability_table.txt` specifies the mutability table file, and the last argument `DD_single_disorder_output.txt` specifies the output file in single disorder analysis. 

### Output for single disorder analysis
After running the above step, EncoreDNM outputs a text file `DD_single_disorder_output.txt`, with the following columns:

* `variant_class`: e.g. lof, Dmis, Tmis, syn. 
* `beta_est`: Estimated value for elevation parameter beta. 
* `beta_sd`: Standard error for elevation parameter beta. 
* `beta_p`: Wald test p-value for elevation parameter beta. 
* `sigma_est`: Estimated value for dispersion parameter sigma. 
* `sigma_sd`: Standard error for dispersion parameter sigma. 
* `sigma_p`: Wald test p-value for dispersion parameter sigma. 
* `deviance_null`: The deviance between the mixed-effects Poisson model and the fixed-effects Poisson model.
* `p_null`: The likelihood ratio test p-value between the mixed-effects Poisson model and the fixed-effects Poisson model. 

We have prepared the example output files for DD and ASD in `./example/DD_single_disorder_output.txt` and `./example/ASD_single_disorder_output.txt`. 


### Cross disorders analysis
To perform cross-disorders EncoreDNM analysis, run the following commands:
    
    cd ./EncoreDNM
    Rscript ./EncoreDNM_cross_disorders.R DD.txt ASD.txt N1 N2 DD_single_disorder_output.txt ASD_single_disorder_output.txt mutability_table.txt DD_ASD_cross_disorders_output.txt

The argument `DD.txt` and `ASD.txt` specify two input DNM summary data files, the argument `N1` and `N2` specify the trio sample sizes for two disorders which are 31058 for DD and 6430 for ASD here, the argument `DD_single_disorder_output.txt` and `ASD_single_disorder_output.txt` specify two output files in single disorder analysis, the argument `mutability_table.txt` specifies the mutability table file, and the last argument `DD_ASD_cross_disorders_output.txt` specifies the output file in cross disorders analysis. 

### Output for cross disorders analysis
After running the above step, EncoreDNM outputs a text file `DD_single_disorder_output.txt`, with the following columns:

* `variant_class`: e.g. lof, Dmis, Tmis, syn. 
* `beta1_est`: Estimated value for elevation parameter beta in the first disorder. 
* `beta2_est`: Estimated value for elevation parameter beta in the second disorder. 
* `sigma1_est`: Estimated value for dispersion parameter sigma in the first disorder. 
* `sigma2_est`: Estimated value for dispersion parameter sigma in the second disorder. 
* `rho_est`: Estimated value for enrichment correlation rho. 
* `beta1_sd`: Standard error for elevation parameter beta in the first disorder. 
* `beta2_sd`: Standard error for elevation parameter beta in the second disorder. 
* `sigma1_sd`: Standard error for dispersion parameter sigma in the first disorder. 
* `sigma2_sd`: Standard error for dispersion parameter sigma in the second disorder. 
* `rho_sd`: Standard error for enrichment correlation rho. 
* `beta1_p`: Wald test p-value for elevation parameter beta in the first disorder. 
* `beta2_p`: Wald test p-value for elevation parameter beta in the second disorder. 
* `sigma1_p`: Wald test p-value for dispersion parameter sigma in the first disorder. 
* `sigma2_p`: Wald test p-value for dispersion parameter sigma in the second disorder. 
* `rho_p`: Wald test p-value for enrichment correlation rho. 

We have prepared the example output file in `./example/DD_ASD_cross_disorders_output.txt`. 
