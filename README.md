# About this repository
This repository contains R scripts and KNIME workflows associated with the manuscript "XXX".



# KNIME workflows
To run these workflows, the following KNIME extensions need to be installed: Weka Data Mining Integration (3.7), File Handling Nodes, Interactive R Statistics Integration, Virtual Nodes. KNIME uses R for several tasks during the workflow, therefore R also needs to be installed on the same computer, and this installation must include the R packages rJava and Rserve. In addition, before using the workflow for the first time, one needs to open an “R snippet” node in KNIME and execute the lines <i>install.packages(“ggplot2”)</i> and <i>install.packages(“gridExtra”)</i>. To execute code in a single, unconnected R snippet node, highlight the commands and click “Eval Selection”. 

- <b> Knime_RF_p1.1_batch.zip 1.knwf </b> This workflow is designed to be executed locally and uses parallelization to speed up the predictions. All nodes and metanodes are annotated with their purpose. Parameters such as number of decision trees, number of negative training proteins etc. can be modified in the first R snippet node as indicated in the workflow. The input is a table of positive training sets (for example, Clusters_as_training_seeds.csv) and ProteomeHD.

# R scripts
- <b> treeClust_OPTICS.R </b> This script determines clusters of tightly co-regulated proteins in ProteomeHD, which were then used as positive training sets for progulonFinder. It first turns ProteomeHD into a distance matrix using the treeClust algorithm, then clusters it using OPTICS.

- <b> Progulon_checkup.R </b> This script takes the progulon calculations coming from the OPTICS clusters and removes results where the model hasn't performed well, thereby reducing the original round of 63 progulons down to 41 final progulons.

- <b> ProgulonStats.R </b> Calculates how many proteins are assigned to how many progulons (related to Fig. 1f)

- <b> ProgulonOverlap.R </b> Calculates extent of overlap between progulons (related to Fig. 1g)

- <b> ProgulonCor.R </b> Calculates the Spearman correlation between all possible progulon combinations across ProteomeHD. It outputs the mean correlation of all protein pairs made by proteins from two different progulons. Results are written out as ProgulonCor.csv and Fig. 1h.

- <b> ProgulonSNE.R </b> This script creates the scatter, line and t-SNE plots related to Figs. 1b-e. It requires some manual annotation files that are also provided here.


# Input files
- <b> ProteomeHD_v1.7z </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios.

- <b> Clusters_as_training_seeds.csv </b> Table containing the protein IDs for the 63 clusters identified by OPTICS. It is used as input for the offline progulonFinder KNIME workflow (Knime_RF_p1.1_batch.zip 1.knwf), which loops over the columns in this file using each as positive training set.

- <b> Unpolished_progulons.7z </b> Compressed folder containing the individual progulonFinder results (Random Forest scores used stats) for the 63 progulons trained by OPTICS clusters. It's the input for Progulon_checkup.R.

- <b> Progulons.7z </b> Random Forest scores of the 41 progulons in long format (produced by Progulon_checkup.R)

- <b> Progulon_scores.7z </b> Random Forest scores of the 41 progulons in wide format (produced by Progulon_checkup.R)

- <b> Manual_Progulon_Annotation.csv </b> Small table containing the manually annotated names and IDs for the 41 progulons.

- <b> ProgulonCor.csv </b> Small file with the pairwise Spearman correlations between all progulons across ProteomeHD.

- <b> Manual_annotation_PRN02.csv </b> and <b> Manual_annotation_PRN11.csv </b> and <b> Manual_annotation_PRN21.csv </b> These files contain manual annotations of proteins highlighted in Fig. 1c, d, and e, respectively. Required by the ProgulonSNE.R script.

# Input files available elsewhere

