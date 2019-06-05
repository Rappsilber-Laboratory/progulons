# About this repository
This repository contains R scripts and KNIME workflows associated with the manuscript "XXX".



# KNIME workflows
To run these workflows, the following KNIME extensions need to be installed: Weka Data Mining Integration (3.7), File Handling Nodes, Interactive R Statistics Integration, Virtual Nodes. KNIME uses R for several tasks during the workflow, therefore R also needs to be installed on the same computer, and this installation must include the R packages rJava and Rserve. In addition, before using the workflow for the first time, one needs to open an “R snippet” node in KNIME and execute the lines <i>install.packages(“ggplot2”)</i> and <i>install.packages(“gridExtra”)</i>. To execute code in a single, unconnected R snippet node, highlight the commands and click “Eval Selection”.

- <b> Knime_RF_p1.1_batch.zip 1.knwf </b> This workflow is designed to run be executed locally and uses parallelization to speed up the predictions.


