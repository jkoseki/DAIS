# DAIS
Dynamical Analysis of Interaction and Structural changes

This is a computational framework for automatically extracting important structural features between two proteins without expertise or experience by simultaneously capturing essential geometric features and their dynamic changes from protein conformational change patterns obtained from molecular dynamics simulations. 

![image](https://user-images.githubusercontent.com/66928602/204697268-0a8c17b4-241c-4833-b211-ee8c30bd223c.png)

## Environment Preparation 
Create a conda environment using conda-TDA.yml. After activating the TDA conda environment, install the bio3d, TDA, readr, data.table, tidyr, stringr, kernlab, tidyverse, dplyr, openxlsx, earth, Rtsne, mclust, gplots, and pheatmap packages in R.

Since the DAIS method can be implemented by executing the R script, the following is not required, but is a procedure to be followed in order to use the uploaded bash script.

1. Create a bin directory directly under the Linux home directory.
2. Save DAIS.bash, AUTO-DAIS.bash, and PDB_for-TDA.bash in the bin directory and grant execution permission (chmod +x).
3. Create a DAIS-Source directory in the bin directory and save the R script in it.

The following assumes the case for parallelized calculations（R-script_parallel computing）.
