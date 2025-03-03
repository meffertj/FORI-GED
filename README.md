# FORI-GED

Code for our paper "Enhancing Graph Edit Distance Computation: Stronger and Orientation-based ILP Formulations"

## Installing

### Following software is needed:

- GUROBI
- GEDLIB https://github.com/dbblumenthal/gedlib
- LIBLSAP: https://forge.greyc.fr/projects/liblsap/ (header only, add LIBLSAP_ROOT=~/liblsap to cmake options)
need to change ```static const DataType _zero = 0;``` in liblsap/cpp/include/lsap.hh, line 50 to ```static constexpr DataType _zero = 0;```

### Variables to add to path/cmake:

#### add the following to your cmake options:

- -DGUROBI_DIR=/path/to/gurobi1200/linux64/

#### add the following to your environment:

- LIBLSAP_ROOT=/path/to/liblsap;GEDLIB_ROOT=/path/to/gedlib;GUROBI_HOME="/path/to/gurobi1200/linux64/"

## How to replicate our results

## Without installation
We provide .mps formulation files in the ```Formulations``` folder, grouped by model. In order to run formulations by bin, the scripts folder contains the files ```q1_script_execute_mps.sh``` and ```q2_script_execute_mps.sh```, which execute all the .mps file in the folder they are placed in, with the Gurobi settings used to answer research questions n.1 and n.2, respectively (cf. our paper for further details).   
By default, both scripts invoke the Gurobi command line tool. Feel free to edit the scripts in order to read and solve the formulations with your preferred ILP solver.

### Example of use:
Suppose you want to reproduce the experiments for FORI formulation on the AIDS-21-30 bin. Then,

- extract the archive ```AIDS_21_30.zip``` from ```Formulations/FORI/AIDS``` folder
- copy into the extracted folder the ```q1_script_execute_mps.sh``` script
- run the script

## With installation
Once our source code has been installed and compiled, move to the ```build``` folder and paste there any of the ```testing_gurobi_*.sh``` script, depending on which bin to replicate.   
Due to the specific Gurobi setting used to answer research question n.2 (cf. our paper for further details), scripts for IMDB, Cora and Pubmed datasets are not provided in this way. You can replicate them by following the Without installation guide.
