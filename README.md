# DAIS
Dynamical Analysis of Interaction and Structural changes

This is a computational framework for automatically extracting important structural features between two proteins without expertise or experience by simultaneously capturing essential geometric features and their dynamic changes from protein conformational change patterns obtained from molecular dynamics simulations. 

Koseki J, et al., Topological data analysis of protein structure and inter/intra-molecular interaction changes attributable to amino acid mutations. Comput Struct Biotechnol J. 2023 May 9;21:2950-2959. doi: 10.1016/j.csbj.2023.05.009.

![image](https://user-images.githubusercontent.com/66928602/204697268-0a8c17b4-241c-4833-b211-ee8c30bd223c.png)

## Environment Preparation 
Create a conda environment using conda-TDA.yml. After activating the TDA conda environment, install the bio3d, TDA, readr, data.table, tidyr, stringr, kernlab, tidyverse, dplyr, openxlsx, earth, Rtsne, mclust, gplots, and pheatmap packages in R.

Since the DAIS method can be implemented by executing the R script, the following is not required, but is a procedure to be followed in order to use the uploaded bash script.

1. Create a bin directory directly under the Linux home directory.
2. Save DAIS.bash, AUTO-DAIS.bash, mdcrd2pdb_for-TDA.bash, and PDB_for-TDA.bash in the bin directory and grant execution permission (chmod +x).
3. Create a DAIS-Source directory in the bin directory and save the R script in it.
4. Set path in the directry, ~/bin.

The following assumes the case for parallelized calculations (R-script_parallel computing).


## How to run DAIS
In order to perform a DAIS calculation, a structure that serves as one Control（WT） and one or more Cases must be prepared in the MD calculation.
If the MD calculation is done using Amber, it is possible to create a PDB for DAIS from prmtop and mdcrd files by using mdcrd2pdb_for-TDA.bash.
On the other hand, if you have used other software to calculate the trajectory, please convert the trajectory into a single PDB that combines the trajectories. For that PDB, you can use PDB_for-TDA.bash to convert it to a PDB for DAIS.

Running DAIS.bash in the execution directory creates all necessary scripts and calculation directories. At that time, the structure in the execution directory will be as follows.

![image](https://github.com/user-attachments/assets/d9b8753f-5f29-49cf-971d-9c1b291d3eeb)

Within each Bash, the number of Cases, the number of PDBs (snapshots), and (only if parallelism is possible) the number of parallelisms must be specified.
- numcase  :  number of Cases 
- numtpdf  :  number of PSBs
- nproc    :  number of parallel CPU cores

When all settings are done, run the bash scripts one by one. Alternatively, you can run AUTO-DAIS.bash to fully automate the DAIS calculation to the end.

## DAIS outputs
Detailed outputs are discharged in various ways, but the easiest way for users to grasp them is to look at the PDB discharged within Pickup-Variation-Data.
Each file name consists of the feature structure site of the Case, the feature structure site of the WT, and the corresponding ID Score value.
Then, by looking at [Case name]-Change-point.csv, you can know the C_alpha of the site where you observed the structural change.
