# progulons
R scripts and KNIME workflows associated with our "Progulon" manuscript




To run the workflow, the following KNIME extensions need to be installed: Weka Data Mining Integration (3.7), File Handling Nodes, Interactive R Statistics Integration, Virtual Nodes. KNIME uses R for several tasks during the workflow, therefore R also needs to be installed on the same computer, and this installation must include the R packages rJava and Rserve. Finally, before using the workflow for the first time, the user needs to open  “R snippet” node needs to be opened 


an “R snippet” node. Open the node and type install.packages(“ggplot2”) and install.packages(“gridExtra”). To execute these installation instructions in a loose / unconnected node, just highlight the commands and click “Eval Selection”.
