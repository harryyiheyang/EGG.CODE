# EGG Simulation

This repository provides the codes used in our paper [EGG](https://arxiv.org/abs/2401.07767).

## Files Description

### EASGraph.R and EURGraph.R
Scripts used for estimating the genetic precision matrices in the real data analysis.

### EASZZ 
The Z-score matrix of the EAS GWAS summary data for the 20 traits, including only the instrumental variables (IVs) used in the real data analysis.

### EASZR 
The estimation error correlation matrix of the EAS GWAS summary data.

### EURZZ 
The Z-score matrix of the EUR GWAS summary data for the 20 traits, including only the instrumental variables (IVs) used in the real data analysis.

### EURZR 
The estimation error correlation matrix of the EUR GWAS summary data.

### PE_Simulation.R and PE_program.R 
Scripts for simulations comparing prediction errors.

### SE_Simulation.R and SE_program.R 
Scripts for simulations comparing standard error estimations.

### basicfunction.R and functions.R
During the simulation phase, the EGG package was not yet developed for public use. These scripts contain the original functions used in our simulations, as utilized in `PE_Simulation.R` and `PE_program.R`.

## Contact

Yihe Yang  
Email: yxy1234@case.edu  
ORCID: [0000-0001-6563-3579](https://orcid.org/0000-0001-6563-3579)

## License

This project is licensed under the MIT License.
