# About this repository
This repository contains R scripts and KNIME workflows associated with the manuscript "XXX".



# KNIME workflows
To run these workflows, the following KNIME extensions need to be installed: Weka Data Mining Integration (3.7), File Handling Nodes, Interactive R Statistics Integration, Virtual Nodes. KNIME uses R for several tasks during the workflow, therefore R also needs to be installed on the same computer, and this installation must include the R packages rJava and Rserve. In addition, before using the workflow for the first time, one needs to open an “R snippet” node in KNIME and execute the lines <i>install.packages(“ggplot2”)</i> and <i>install.packages(“gridExtra”)</i>. To execute code in a single, unconnected R snippet node, highlight the commands and click “Eval Selection”. 

- <b> Knime_RF_p1.1_batch.zip 1.knwf </b> This workflow is designed to be executed locally and uses parallelization to speed up the predictions. All nodes and metanodes are annotated with their purpose. Parameters such as number of decision trees, number of negative training proteins etc. can be modified in the first R snippet node as indicated in the workflow.

# R scripts
- <b> treeClust_OPTICS.R </b> This script determines clusters of tightly co-regulated proteins in ProteomeHD, which were then used as positive training sets for progulonFinder. It first turns ProteomeHD into a distance matrix using the treeClust algorithm, then clusters it using OPTICS.

- <b> Progulon_checkup.R </b> This script takes the progulon calculations coming from the OPTICS clusters and removes results where the model hasn't performed well, thereby reducing the original round of 63 progulons down to 41 final progulons.


# Input files
- <b> Unpolished_progulons.7z </b> Compressed folder containing the individual progulonFinder results (Random Forest scores used stats) for the 63 progulons trained by OPTICS clusters. It's the input for Progulon_checkup.R.



# Input files available elsewhere
- <b> ProteomeHD_v1_1.7z </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios. This file is available in the GitHub repository for the corresponding publication (https://github.com/Rappsilber-Laboratory/ProteomeHD/blob/master/Data/ProteomeHD_v1_1.7z) or as Supplementary Table S1 of the same publication (Kustatscher et al, 2019, <i>The human proteome co-regulation map reveals functional relationships between proteins<i/>)
