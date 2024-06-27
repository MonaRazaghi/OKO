
# The optimization of enzyme catalytic rates in engineering of metabolic phenotypes

## OS
* Code was tested on Windows 10

## Dependencies
* Matlab (tested with 2023a)
* [Gurobi solver](https://www.gurobi.com/) (tested with version 9.5.2)

## Setup
* Set up the Gurobi solver and connect it with the COBRA toolbox using install instructions that can be found [here](https://opencobra.github.io/cobratoolbox/stable/installation.html#gurobi)

## Run k<sub>cat</sub> optimization with OKO

### To reproduce the results presented in the paper from OKO run the following scripts
* _S. cerevisiae_: `Yeast_OKO`
* _E. coli_: `Ecoli_OKO`

### To reproduce the results presented in the paper from OKO run the following scripts
* _S. cerevisiae_: `Yeast_OKO_plus`
* _E. coli_: `Ecoli_OKO_plus`

### In the case of infeasibility, biomass/product fraction must be reduced by 0.1
* for instance lower(find(ecModel.c)) = 0.99 * opt_bio; can be changed to lower(find(ecModel.c)) = 0.98 * opt_bio; 

### To apply OKO to another model
create input data file similar to Metabolites.xlsx, where the B column contains the metabolite index in the ecGEM
