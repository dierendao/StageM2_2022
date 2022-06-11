# StageM2_2022
These scripts are part of my Master 2 Bioinformatics internship. It took place at the Cancer Research Center of Marseille in the team "epigenetics factors in hematopoeisis".
The subject of my internship is entitled: single cell RNA sequencing (scRNA-seq) analysis of murine promyelocytes transformed by the oncogenic fusion protein PLZF-RARA.

The directories are constituted as follows:
- Create_new-reference:
Contains the scripts that allowed me to create a new reference (mouse genome + trangene (hcG- PLZF-RARA).

- Gene_counting
Contains the scripts for counting the genes after alignment on the new reference

- Integration_scp
 bash script to run the IntegrateCT-AR.R on cluster. Indeed this script takes too much time if it is launched from R studio of the ifb.
 
- Rmd_scr
Contains the set of Rmd scripts that groups my R scripts and they are sorted in the following order:
IntegrateCT-AR.Rmd : All steps from the loading of data from scRNA-Seq, through quality control and normalization, to the integration of untreated (CT) and retinoic acid (RA) treated cells
TransferLabel.Rmd: Label transfer between the data obtained by the host lab (without adding the transgene to the murine genome) and the data I produced. The lab data was used as a reference to predict the clustering of my new data.
AnalysisTransferLabel.Rmd: Analysis and visualization of results such as transgene expression, marker gene identification and statistics associated with differential gene expression ....
UnsupervisedClustering.Rmd: De novo clustering of data and its analysis
PrepareSTREAManalysis.Rmd : Preparation of the data for the study of the cell trajectory by STREAM.
StreamSeurat.Rmd : Concordant analysis of the results of STREAM and Seurat in order to keep the same cluster names and especially to see the distribution of cells on the branches.
AllFunctions : all the functions used in my scripts and that I did not develop.


- Py_scr
contains python file for STREAM code and runing.
