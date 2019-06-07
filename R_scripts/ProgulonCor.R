# Load required libraries
library(data.table); library(ggplot2); library(viridis)

#### Calculate Spearman correlations between all relevant proteins in ProteomeHD ####

# Load and prep ProteomeHD
ProHD <- read.csv("ProteomeHD_v1.csv", stringsAsFactors=FALSE)
rownames(ProHD) <- ProHD$Majority.protein.IDs                   # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]                # Keep only columns with SILAC ratios
feature_count <- apply(ProHD, 1, function(x){ sum(!is.na(x))})  # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 30 ,]                          # Only proteins with >= 30 feature counts will have RF scores

# Transpose and normalise SILAC ratios
SILAC <- t(ProHD)
SILAC_exp_medians <- apply(SILAC, 1, median, na.rm=TRUE)
SILAC <- sweep(SILAC, 1, SILAC_exp_medians, FUN="-")


#### Load progulon data and assign proteins ####

# Load the progulon data
prns <- fread("Progulons.csv")

# Function assigning proteins to progulons based on RF score cut-off and feature counts
f_prn_calling <- function(x){                                       # Input is a data.table containing the score and the feature counts
  ifelse( x$Mean_RF_score >= 0.55, "yes" ,                          # Keep all proteins that score >= 0.55
          ifelse( x$Mean_RF_score >= 0.5 & x$Feature_count >= 100 , "yes",  # Keep proteins between 0.5 and 0.55 only if they have enough feature counts
                  "no"))}                                                           # All other proteins are not in the progulon

# Assign proteins to progulons
prns$prot_in_prn <- f_prn_calling(prns)


#### Calculate Spearman correlations between all progulons (~ 2 hour runtime) ####
all_prn_combinations <- data.table( t( combn( unique(prns$Progulon_ID) , 2 )))

pb <- txtProgressBar(min = 0, max = all_prn_combinations[,.N], style = 3)
for(i in 1:all_prn_combinations[,.N]){ 
  prnA_prots <- prns[ Progulon_ID == all_prn_combinations[i,V1] & prot_in_prn == "yes" , Protein_IDs ]  # Proteins in first prn
  prnB_prots <- prns[ Progulon_ID == all_prn_combinations[i,V2] & prot_in_prn == "yes" , Protein_IDs ]  # Proteins in 2nd prn
  prnA_ratios <- SILAC[, match(prnA_prots, colnames(SILAC))]
  prnB_ratios <- SILAC[, match(prnB_prots, colnames(SILAC))]
  AB_cor <- cor(prnA_ratios, prnB_ratios, use = "pairwise.complete.obs", method = c("spearman"))
  AB_cor <- as.numeric(AB_cor)
  all_prn_combinations[ i ,       N_NA := sum( is.na(AB_cor))        ]
  all_prn_combinations[ i , N_proteins := length(AB_cor)             ]       
  all_prn_combinations[ i ,        RHO := mean(AB_cor, na.rm = TRUE) ]
  all_prn_combinations[ i , N_prnA_prots := length(prnA_prots)       ]
  all_prn_combinations[ i , N_prnB_prots := length(prnB_prots)       ]
  setTxtProgressBar(pb, i)
}

# Control 1: Do we have all data as intended?
all_prn_combinations[ complete.cases(all_prn_combinations) , .N] == all_prn_combinations[,.N]

# Control 2: N_prnA_prots * N_prnB_prots has to equal N_proteins
loop_validation <- all_prn_combinations[, (N_prnA_prots * N_prnB_prots) == N_proteins ]
sum( loop_validation ) == length( loop_validation )

# Control 3: In some cases there were NA in the correlation matrix due to lack of overlapping experiments
# This should be a small fraction of the observations
all_prn_combinations[, missing_values_pc := N_NA / (N_proteins/100) ]
all_prn_combinations[ order(-missing_values_pc) ][ 1:10 ]

#### Write out the results ####
all_prn_combinations <- all_prn_combinations[, .(PRN_A = V1, PRN_B = V2, RHO, N_proteins) ]
fwrite(x = all_prn_combinations, file = "ProgulonCor.csv")


#### Plot a correlation matrix ####

# Load the progulaon annotations
prn_annot <- fread("Manual_Progulon_Annotation.csv")

# Expand data to be able to get a complete matrix, i.e. append duplicates
all_prn_combinations <- rbind(all_prn_combinations[, .(PRN_A, PRN_B, RHO)],
                              all_prn_combinations[, .(PRN_B = PRN_A, PRN_A = PRN_B, RHO)])

# Append functional annotation
all_prn_combinations[, Function_1 := prn_annot[ match( all_prn_combinations[, PRN_A], prn_annot[, Progulon] ), Function ]]
all_prn_combinations[, Function_2 := prn_annot[ match( all_prn_combinations[, PRN_B], prn_annot[, Progulon] ), Function ]]

# Cast into a correlation matrix
cor_mat <- dcast( all_prn_combinations, Function_1 ~ Function_2, value.var = "RHO" )
my_rownames <- cor_mat[, Function_1] 
cor_mat[, Function_1 := NULL ]
cor_mat <- as.data.frame( cor_mat )
rownames(cor_mat) <- my_rownames

# Group progulons by correlation
my_dist <- as.dist( (1-cor_mat)/2 )
my_clust <- hclust(my_dist)
new_prn_order <- rownames(cor_mat)[ my_clust$order ]

# Rearrange progulons (via factor levels) in clustered order
all_prn_combinations[, Function_1 := factor(Function_1, levels = new_prn_order )]
all_prn_combinations[, Function_2 := factor(Function_2, levels = new_prn_order )]

# Create the plot
p <- ggplot(all_prn_combinations, aes(Function_1, Function_2, fill = RHO))+
      geom_tile()+
      scale_fill_viridis()+
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6, angle = 90, hjust=1, vjust = 0.5))

p1 <- p + theme( legend.position = "none")
p2 <- p + theme( legend.position = "right")

# Save the plot
ggsave("ProgulonRho.pdf", p1, width = 10.5, height = 10.5, units = "cm")
ggsave("ProgulonRho_legend.pdf", p2, width = 10.5, height = 10.5, units = "cm")






