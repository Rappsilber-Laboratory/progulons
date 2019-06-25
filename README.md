# About this repository
This repository contains R scripts and KNIME workflows associated with the manuscript "XXX".



# KNIME workflows
To run these workflows, the following KNIME extensions need to be installed: Weka Data Mining Integration (3.7), File Handling Nodes, Interactive R Statistics Integration, Virtual Nodes. KNIME uses R for several tasks during the workflow, therefore R also needs to be installed on the same computer, and this installation must include the R packages rJava and Rserve. In addition, before using the workflow for the first time, one needs to open an “R snippet” node in KNIME and execute the lines <i>install.packages(“ggplot2”)</i> and <i>install.packages(“gridExtra”)</i>. To execute code in a single, unconnected R snippet node, highlight the commands and click “Eval Selection”. 

- <b> Knime_RF_p1.1_batch.zip 1.knwf </b> This "offline" workflow is designed to be executed locally and uses parallelization to speed up the predictions. All nodes and metanodes are annotated with their purpose. Parameters such as number of decision trees, number of negative training proteins etc. can be modified in the first R snippet node as indicated in the workflow. The input is a table of positive training sets (for example, Clusters_as_training_seeds.csv) and ProteomeHD.

- <b> Knime_RF_online_v4_16A.knwf </b> Essentially the same workflow but designed to be executed online through our web-app. It does not make use of parallelization but includes a part that matches Uniprot IDs to the proteins that are actually present in ProteomeHD.

# R scripts
- <b> treeClust_OPTICS.R </b> This script determines clusters of tightly co-regulated proteins in ProteomeHD, which were then used as positive training sets for progulonFinder. It first turns ProteomeHD into a distance matrix using the treeClust algorithm, then clusters it using OPTICS.

- <b> Progulon_checkup.R </b> This script takes the progulon calculations coming from the OPTICS clusters and removes results where the model hasn't performed well, thereby reducing the original round of 63 progulons down to 41 final progulons.

- <b> ProgulonStats.R </b> Calculates how many proteins are assigned to how many progulons (related to Fig. 1f)

- <b> ProgulonOverlap.R </b> Calculates extent of overlap between progulons (related to Fig. 1g)

- <b> ProgulonCor.R </b> Calculates the Spearman correlation between all possible progulon combinations across ProteomeHD. It outputs the mean correlation of all protein pairs made by proteins from two different progulons. Results are written out as ProgulonCor.csv and Fig. 1h.

- <b> ProgulonSNE.R </b> This script creates the scatter, line and t-SNE plots related to Figs. 1b-e. It requires some manual annotation files that are also provided here.

- <b> ProgulonTopGo.R </b> Uses topGO to calculate GO term enrichment for progulons. Set onotlogy MF, BP, CC and execute sequentially.

- <b> ProgulonTopGo_plotting.R </b> Takes the files created by ProgulonTopGO.R as input and outputs the GO term enrichment plots.

- <b> PRN_genome_loc.R </b> Tests if progulons are enriched for genes from the same chromosome and outputs a volcano plot. Uses genome positions (Human_gene_positions_26Mar18.csv) from ENSEMBL as input.

- <b> PRNs_RNA_prot_Haas.R </b> and <b> PRNs_RNA_prot_Mouse.R </b> and <b> PRNs_RNA_prot_LCL.R </b> These scripts were used to perform all analyses and create all plots related to mRNA and protein comparisons for the breast cancer, mouse and lymphoblastoid cell data, respectively.

- <b> Prn_evol.R </b> Script comparing mouse and human co-regulation within and outside progulons, to assess their potential evolutionary conservation

- <b> Progulon_HL.R </b> This script calculates and compares mRNA and protein half-lifes and turnover kinetics. The mRNA half-lives are from Tani et al (Genome Research, 2012). The script uses their Table S1 as input. The protein half-lives are from McShane et al (Cell, 2016). The script uses their Table S4 as input.

- <b> Replisome.R </b> This script performs the analysis and creates all plots related to the manual prediction of the replisome progulon (the predictions themselves are created by the KNIME-based progulonFinder workflow).

- <b> siRNA_screen_method_comp.R </b> Short script to plot the barcharts found in Supplementary Fig. 8 (comparing SD and RSA scoring methods for the siRNA screen).


# Input files
- <b> ProteomeHD_v1.7z </b> This compressed csv file is ProteomeHD, consisting of 10,323 proteins and 294 SILAC ratios.

- <b> ProteomeHD_v1_1.7z </b> A slightly different version of ProteomeHD used by some scripts. This is median-normalised and the column names have been polished.

- <b> Clusters_as_training_seeds.csv </b> Table containing the protein IDs for the 63 clusters identified by OPTICS. It is used as input for the offline progulonFinder KNIME workflow (Knime_RF_p1.1_batch.zip 1.knwf), which loops over the columns in this file using each as positive training set.

- <b> Unpolished_progulons.7z </b> Compressed folder containing the individual progulonFinder results (Random Forest scores used stats) for the 63 progulons trained by OPTICS clusters. It's the input for Progulon_checkup.R.

- <b> Progulons.7z </b> Random Forest scores of the 41 progulons in long format (produced by Progulon_checkup.R)

- <b> Progulon_scores.7z </b> Random Forest scores of the 41 progulons in wide format (produced by Progulon_checkup.R)

- <b> Manual_Progulon_Annotation.csv </b> Small table containing the manually annotated names and IDs for the 41 progulons.

- <b> ProgulonCor.csv </b> Small file with the pairwise Spearman correlations between all progulons across ProteomeHD.

- <b> Manual_annotation_PRN02.csv </b> and <b> Manual_annotation_PRN11.csv </b> and <b> Manual_annotation_PRN21.csv </b> These files contain manual annotations of proteins highlighted in Fig. 1c, d, and e, respectively. Required by the ProgulonSNE.R script.

- <b> Human_gene_positions_26Mar18.csv </b> ENSEMBL genome annotation used for this manuscript.

- <b> BattleSILAC_PickrellRPKM.7z </b> Compressed csv file containing mRNA and protein abundance changes across lymphoblastoid cell lines. See Kustatscher et al (Mol Syst Biol, 2017) for more details. Original transcriptomics data from Pickrell et al (Nature, 2010), original proteomics data from Battle et al (Science, 2015).

- <b> mouse_SILAC_TPMs_log2_final_min8_features.csv </b> File containing the matched mRNA and protein abundance changes across mouse tissues. See Grabowski et al (MCP, 2018) for more details. Original transcriptomics data from a range of studies, original proteomics data from Geiger et al (MCP, 2013).

- <b> ID_conversion.tab </b> File to convert protein IDs from Lapek et al to ProteomeHD (creating using UniProt's Retrieve tool).

- <b> RefSeq2_to_Uniprot.tab </b> Conversion file used by Progulon_HL>R to match mRNA half-lives from Tani et al to protein half-lives from McShane et al.

- <b> allComplexes.txt </b> File containing CORUM protein complexes. Used by Progulon_HL.R to assess if NEDs are only enriched in progulons because protein complexes are enriched. 

- <b> Rna_Pro_PRN_LCL.csv </b> and <b> Haas_Rna_Pro_PRN.csv </b> and <b> Mouse_Rna_Pro_PRN.csv </b> These three files contain mRNA and protein correlations. They are also available as separate tabs in Supplementary Table 2.

- <b> Replisome_training.csv </b> List of 41 training proteins used for the manual prediction of the replisome progulon (NOTE: This file contains proteinGroups with isoform information. It can be used directly for the offline progulonFinder workflow, but for the online version the isoform annotation needs to be removed).

- <b> QuickGO_DNA_replication_GO0006260.tsv </b> GO annotation for DNA replication obtained from QuickGO.

- <b> Candidates_siRNA_screen.csv </b> IDs of candidates used for the siRNA screen and the corresponding results.

- <b> RF_score_Replisome_all.csv </b> Result of progulonFinder for the manual prediction of the progulon.

- <b> RF_score_Replisome_without_NCC.csv </b> Result of progulonFinder for the manual prediction of the progulon - but without considering the 15 NCC experiments.

- <b> RF_score_Replisome_without15random_1.csv </b> and <b> RF_score_Replisome_without15random_2.csv </b> and <b> RF_score_Replisome_without15random_3.csv </b>  Result of progulonFinder for the manual prediction of the progulon - but without considering the 15 random experiments (three different sets of randomly omitted experiments).

- <b> RF_score_Replisome_NCConly.csv </b>  Result of progulonFinder for the manual prediction of the progulon - considering only the 15 NCC experiments.

- <b> RF_score_Replisome_only15random_1.csv </b> and <b> RF_score_Replisome_only15random_2.csv </b> and <b> RF_score_Replisome_only15random_3.csv </b> Result of progulonFinder for the manual prediction of the progulon - considering only 15 random experiments (three different sets of randomly selected experiments).

