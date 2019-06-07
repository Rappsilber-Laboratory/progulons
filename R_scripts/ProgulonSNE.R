library(data.table); library(ggplot2); library(viridis); library(treeClust); library(Rtsne); library(GA); library(gridExtra); library(cowplot)


#### Load and prep data ####

# Load the progulon data
prns <- fread("Progulons.csv")

# Function assigning proteins to progulons based on RF score cut-off and feature counts
f_prn_calling <- function(x){                                       # Input is a data.table containing the score and the feature counts
  ifelse( x$Mean_RF_score >= 0.55, "yes" ,                          # Keep all proteins that score >= 0.55
          ifelse( x$Mean_RF_score >= 0.5 & x$Feature_count >= 100 , "yes",  # Keep proteins between 0.5 and 0.55 only if they have enough feature counts
                  "no"))}                                                           # All other proteins are not in the progulon

# Assign proteins to progulons
prns$prot_in_prn <- f_prn_calling(prns)

# Load and prepare ProteomeHD (no normalisation necessary)
ProHD <- read.csv("ProteomeHD_v1.csv", stringsAsFactors=FALSE)
rownames(ProHD) <- ProHD$Majority.protein.IDs       # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]    # Keep only columns with SILAC ratios

# Load the manual annotation data
annot_PRN02 <- fread("Manual_annotation_PRN02.csv")
annot_PRN11 <- fread("Manual_annotation_PRN11.csv")
annot_PRN21 <- fread("Manual_annotation_PRN21.csv")


#### treeClust, tSNE and RF score plots for ATP synthase progulon (PRN21) ####

# Select progulon to be plotted
my_prn <- "PRN21"
my_annot <- annot_PRN21
  
# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(42)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = TRUE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = TRUE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
        annotate(geom = "segment", x=0.50, xend=0.50, y=100, yend=298, size=0.25, linetype="dashed")+
        annotate(geom = "segment", x=0.50, xend=0.55, y=100, yend=100, size=0.25, linetype="dashed")+
        annotate(geom = "segment", x=0.55, xend=0.55, y=2,   yend=100, size=0.25, linetype="dashed")+
        annotate(geom = "segment", x=0.50, xend=0.99, y=298, yend=298, size=0.25, linetype="dashed")+
        annotate(geom = "segment", x=0.55, xend=0.99, y=2,   yend=2,   size=0.25, linetype="dashed")+
        annotate(geom = "segment", x=0.99, xend=0.99, y=2,   yend=298, size=0.25, linetype="dashed")+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta", alpha = 0.8)+
        geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
        scale_x_continuous(limits=c(0,1), expand = c(0,0))+
        scale_y_continuous(limits=c(0,300), expand = c(0,0))+
        xlab("Random forest score")+
        ylab("# available experiments")+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text.x=element_text(size=5), axis.text.y=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)), axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
              axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF
ggsave("PRN21_ATPsyn_RFplot.png", pRF, width=3.7, height=4, units= "cm", dpi = 600)


## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <- c(`ATP synthase` = "#DC3912", `aKG depletion` = "#11069e", `Complex I` = "#FF9900",
             `Complex II` = "#109618", `Complex III` = "#990099", `Complex IV` = "#0099C6",
             `FAO` = "#DD4477", `Uncharacterised` = "#bedf00")

pSNE_21 <- ggplot(SNE, aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_21
ggsave("PRN21_ATPsyn_SNEplot.pdf", pSNE_21, width=8, height=5, units= "cm")


#### treeClust, tSNE and RF score plots for prefoldin progulon (PRN11) ####

# Select progulon to be plotted
my_prn <- "PRN11"
my_annot <- annot_PRN11

# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(42)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = TRUE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = TRUE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                 .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
  annotate(geom = "segment", x=0.50, xend=0.50, y=100, yend=298, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.50, xend=0.55, y=100, yend=100, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.55, xend=0.55, y=2,   yend=100, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.50, xend=0.99, y=298, yend=298, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.55, xend=0.99, y=2,   yend=2,   size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.99, xend=0.99, y=2,   yend=298, size=0.25, linetype="dashed")+
  geom_point(size=0.05, alpha=0.4, colour="grey50")+
  geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta")+
  geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
  #geom_point(size=0.2, data= plot_DT[Used_for_positive_training =="Used"], colour="magenta")+
  scale_x_continuous(limits=c(0,1), expand = c(0,0))+
  scale_y_continuous(limits=c(0,300), expand = c(0,0))+
  xlab("Random Forest score")+
  ylab("Experiments in which protein was quantified")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
        axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)),
        axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF
ggsave("PRN11_Prefoldin_RFplot.pdf", pRF, width=5, height=5, units= "cm")


## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <-  c(`Heat shock proteins` = "#70f0ff",            
              `hnRNPs` = "peru",
              `LSm` = "coral1",
              `Sm` = "lightgoldenrod",
              `Nuclear import/export` = "blue", 
              `Prefoldin` = "#DC3912",                       
              `Prolyl isomerase` = "#f8ff3d",                
              `Proteasome, core` = "#990099",                
              `Proteasome, regulatory` = "#ef6eef",          
              `R2TP/PFDL complex` = "#ef542f",               
              `Ribosome 40S` = "#13a51c",                    
              `Ribosome 60S` = "green",                    
              `Translation factors` = "#00a09b",             
              `TRiC` = "orange",                            
              `Tubulins & tub. chaperones` = "#0087af",      
              `Uncharacterised` = "#bedf00")          


pSNE_11 <- ggplot(SNE[ V1 > -28 ], aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_11
ggsave("PRN11_Prefoldin_SNEplot.pdf", pSNE_11,
       width=10, height=5, units= "cm")
 

  
#### treeClust, tSNE and RF score plots for AP3 progulon (PRN02) ####

# Select progulon to be plotted
my_prn <- "PRN02"
my_annot <- annot_PRN02

# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(41)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = TRUE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = TRUE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                 .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
  annotate(geom = "segment", x=0.50, xend=0.50, y=100, yend=298, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.50, xend=0.55, y=100, yend=100, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.55, xend=0.55, y=2,   yend=100, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.50, xend=0.99, y=298, yend=298, size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.55, xend=0.99, y=2,   yend=2,   size=0.25, linetype="dashed")+
  annotate(geom = "segment", x=0.99, xend=0.99, y=2,   yend=298, size=0.25, linetype="dashed")+
  geom_point(size=0.05, alpha=0.4, colour="grey50")+
  geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta")+
  geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
  #geom_point(size=0.2, data= plot_DT[Used_for_positive_training =="Used"], colour="magenta")+
  scale_x_continuous(limits=c(0,1), expand = c(0,0))+
  scale_y_continuous(limits=c(0,300), expand = c(0,0))+
  xlab("Random Forest score")+
  ylab("Experiments in which protein was quantified")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
        axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)),
        axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF
ggsave("PRN02_AP3_RFplot.pdf", pRF, width=5, height=5, units= "cm")


## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <-  c(`Adaptor protein complex 3` = "#990099",
              `Coatomer` = "lightgoldenrod",
              `COG complex, lobe A` = "blue",
              `COPII coat` = "#DC3912",
              `Exocyst complex` = "#f8ff3d",
              `GIT1/2-PIX complex` = "#990099",
              `Microtubule-based vesicle motility` = "#ef6eef",
              `NRZ vesicle tethering complex` = "green",
              `Other vesicle trafficking` = "#70f0ff",
              `Regulation of actin cytoskeleton` = "orange",
              `Uncharacterised` = "#bedf00")


pSNE_02 <- ggplot(SNE, aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_02
ggsave("PRN02_AP3_SNEplot.pdf", pSNE_02,
       width=10, height=5, units= "cm")


#### Progulon lineplots: seeds + top25 co-regulated proteins in 25 experiments ####

# Select ATP synthase progulon (PRN21) data
PRN21_training_prots <- prns[ Progulon_ID == "PRN21" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
PRN21_progulon_top25 <- prns[ Progulon_ID == "PRN21" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# Select Prefoldin (PRN11) data
PRN11_training_prots <- prns[ Progulon_ID == "PRN11" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
PRN11_progulon_top25 <- prns[ Progulon_ID == "PRN11" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# Select AP3 progulon (PRN02) data
PRN02_training_prots <- prns[ Progulon_ID == "PRN02" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
PRN02_progulon_top25 <- prns[ Progulon_ID == "PRN02" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# These are all the proteins that need to be in the plot
my_prots <- unique( c( PRN21_training_prots, PRN21_progulon_top25, PRN11_training_prots, PRN11_progulon_top25, PRN02_training_prots, PRN02_progulon_top25 ))      

# Fitness function for a genetic algorithm, designed to pick a useful set of experiments to display
# "Useful" is to mean that (a) good coverage and (b) they are experiments were the proteins actually show some (even modest) change
fitness_f <- function(x){ E1 <- ceiling( x[1] )         # Get column index of the currently selected experiments (ratios)
                          E2 <- ceiling( x[2] )         # Need to use ceiling to work with integers
                          E3 <- ceiling( x[3] )
                          E4 <- ceiling( x[4] )
                          E5 <- ceiling( x[5] )
                          E6 <- ceiling( x[6] )
                          E7 <- ceiling( x[7] )
                          E8 <- ceiling( x[8] )
                          E9 <- ceiling( x[9] )
                          E10 <- ceiling( x[10] )
                          E11 <- ceiling( x[11] )  
                          E12 <- ceiling( x[12] )
                          E13 <- ceiling( x[13] )
                          E14 <- ceiling( x[14] )
                          E15 <- ceiling( x[15] )
                          E16 <- ceiling( x[16] )
                          E17 <- ceiling( x[17] )
                          E18 <- ceiling( x[18] )
                          E19 <- ceiling( x[19] )
                          E20 <- ceiling( x[20] )
                          E21 <- ceiling( x[21] )  
                          E22 <- ceiling( x[22] )
                          E23 <- ceiling( x[23] )
                          E24 <- ceiling( x[24] )
                          E25 <- ceiling( x[25] )
                          ratios <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13, E14, E15, E16, E17, E18, E19, E20, E21, E22, E23, E24, E25)
                          
                          # Select the subset of ProteomeHD that covers these proteins and these experiments (ratios)
                          temp_df <- ProHD[ rownames(ProHD) %in% my_prots, ratios]
                          
                          # What's the fraction of missing values across these proteins and experiments?
                          pc_NA <- sum(is.na(temp_df)) / (ncol(temp_df)*nrow(temp_df))
                          
                          # What is the median fold-change in this experiments, and how often does it exceed a certain threshold?
                          PRN21_changes <- apply( temp_df[ rownames(temp_df) %in% PRN21_training_prots,], 2, median, na.rm = TRUE )
                          PRN21_changes <- sum( abs( PRN21_changes ) > 0.5 , na.rm = TRUE)
                          PRN11_changes <- apply( temp_df[ rownames(temp_df) %in% PRN11_training_prots,], 2, median, na.rm = TRUE )
                          PRN11_changes <- sum( abs( PRN11_changes ) > 0.5 , na.rm = TRUE)
                          PRN02_changes <- apply( temp_df[ rownames(temp_df) %in% PRN02_training_prots,], 2, median, na.rm = TRUE )
                          PRN02_changes <- sum( abs( PRN02_changes ) > 0.5 , na.rm = TRUE)
                          
                          # Calculate the output
                          if( length(unique(ratios)) < 25 ){
                            fitness_output <- 0                                                          # If some experiments were picked more than once, the solution is invalid (score zero)
                            } else if( pc_NA > 0.05 ){
                              fitness_output <- 0                                                        # If there are more than 5% missing of values missinge, the display would not be informative (score zero)
                              } else
                                fitness_output <- sum( c(PRN21_changes, PRN11_changes, PRN02_changes))   # In how many experiments is the result above the threshold?
                          
                          # Return the fitness output
                          return(fitness_output) 
                          }                                                                

set.seed(1)
GA <- ga(type="real-valued", fitness = fitness_f, min = rep(1,25), max = rep(294,25))      # Genetic algorithm searching for good combination
my_ratios <- ceiling(GA@solution)          # The 25 ratios identified by GA as good for plotting
my_ratios <- as.integer( my_ratios[1,] )   # Simplified to integer and any duplicate solutions removed     
                          
# Create data frame for plotting
plot_data <- ProHD[ rownames(ProHD) %in% my_prots , my_ratios ]        # The relevant proteins and the relevant ratios
colnames(plot_data) <- paste("exp_", 1:25, sep="")

# Median-normalise
plot_data_exp_medians <- apply(plot_data, 2, median, na.rm=TRUE)       # Median fold-change of these proteins in these experiments
plot_data <- sweep(plot_data, 2, plot_data_exp_medians, FUN="-")       # Set that to zero

# Turn into melted data.table
plot_data$Protein_IDs <- rownames(plot_data)
plot_data <- as.data.table( plot_data )
plot_data <- melt(plot_data, id.vars = "Protein_IDs")

# Make the plot
plot_data[, experiment_number := as.integer( gsub("exp_", "", variable)) ]

pBase <- ggplot(plot_data, aes(x = experiment_number, y = value, group = Protein_IDs))+
         coord_cartesian(ylim = c(-2.5,2.5))+
         xlab("Experiments")+ylab("log2 SILAC ratio")+
         scale_x_continuous( breaks = seq(1,25,1), limits = c(1, 25), expand = c(0,0))+
         scale_y_continuous( breaks = c(-2, 0, 2))+
         theme(legend.position="left", axis.text=element_text(size=5), panel.grid = element_blank(),
                panel.background=element_blank(),
                axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25), axis.ticks.x = element_line(size=0.25), axis.ticks.y = element_line(size=0.25),
                axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

p21 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% PRN21_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% PRN21_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)

p11 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% PRN11_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% PRN11_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)

p02 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% PRN02_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% PRN02_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)

# Save the plots
ggsave("PRN21_ATPsyn_lines.pdf", p21, width=5, height=5, units= "cm")
ggsave("PRN11_prefol_lines.pdf", p11, width=5, height=5, units= "cm")
ggsave("PRN02_AP3_lines.pdf"   , p02, width=5, height=5, units= "cm")


#### Output one combined plot (except RF score) ####
pSNE_21_no_leg <- pSNE_21 + theme( legend.position = "none")
pSNE_11_no_leg <- pSNE_11 + theme( legend.position = "none")
pSNE_02_no_leg <- pSNE_02 + theme( legend.position = "none")

p <- plot_grid(pSNE_21_no_leg, pSNE_11_no_leg, pSNE_02_no_leg, p21, p11, p02,
               ncol = 3, align = "v", rel_heights = c(2.1,1))
p
ggsave("Combined_SNE_and_line.pdf", p, width=18, height=8.3, units= "cm")




