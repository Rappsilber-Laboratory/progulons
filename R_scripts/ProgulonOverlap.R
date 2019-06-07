library(ggplot2); library(data.table)

#### Load input data ####

# Load the progulon data
prns <- fread("Progulon_Scores.csv")

# Load the progulaon annotations
prn_annot <- fread("Manual_Progulon_Annotation.csv")


#### Call progulons ####

# Function to "call" progulons by keeping RF score of proteins that are considered part of the progulon
f_prn_calling <- function(x){ ifelse( x >= 0.55, x ,                       # Keep all proteins that score >= 0.55
                              ifelse( x >= 0.5 & Feature_count >= 100 , x, # Keep proteins between 0.5 and 0.55 only if they have enough feature counts
                                      NA))}                                # Set all other proteins to NA

# Apply calling function to progulons
score_columns <- grep("PRN", names(prns), value = TRUE)
Feature_count <- prns[, Feature_count]
prns[, c(score_columns) := lapply(.SD, f_prn_calling), .SDcols = score_columns]


#### Find overlap between progulons ####

# For pairwise comparison, create all progulon combinations
prn_combinations <- as.data.table( t( combn( score_columns, 2 )))

# Calculate number of proteins overlapping between different progulons
overlap <- integer()
for(i in 1:prn_combinations[,.N]){
  prn_pair <- prns[, as.character( prn_combinations[i] ) , with = FALSE ]    # Get the scores for each pair
  overlap[i] <- prn_pair[ complete.cases(prn_pair) , .N ]                    # Count the proteins for which they both have a score (i.e. proteins are in the progulon)
}
prn_combinations[, overlap := overlap ]


#### Get the order in which progulons are plotted in the correlation panel ####

# Load the correlation data
cor_combis <- fread("ProgulonCor.csv")

# Expand data to be able to get a complete matrix, i.e. append duplicates
cor_combis <- rbind(cor_combis[, .(PRN_A, PRN_B, RHO)],
                              cor_combis[, .(PRN_B = PRN_A, PRN_A = PRN_B, RHO)])

# Append functional annotation
cor_combis[, Function_1 := prn_annot[ match( cor_combis[, PRN_A], prn_annot[, Progulon] ), Function ]]
cor_combis[, Function_2 := prn_annot[ match( cor_combis[, PRN_B], prn_annot[, Progulon] ), Function ]]

# Cast into a correlation matrix
cor_mat <- dcast( cor_combis, Function_1 ~ Function_2, value.var = "RHO" )
my_rownames <- cor_mat[, Function_1] 
cor_mat[, Function_1 := NULL ]
cor_mat <- as.data.frame( cor_mat )
rownames(cor_mat) <- my_rownames

# Group progulons by correlation
my_dist <- as.dist( (1-cor_mat)/2 )
my_clust <- hclust(my_dist)
new_prn_order <- rownames(cor_mat)[ my_clust$order ]



#### Plot overlap between progulons ####

# Expand data into a complete "overlap matrix" (minus diagonal)
prn_combinations <- rbind(prn_combinations,
                          prn_combinations[, .(V1 = V2, V2 = V1, overlap) ]) # Append duplicates to complete matrix

# Append functional annotation
prn_combinations[, Function_1 := prn_annot[ match( prn_combinations[, V1], prn_annot[, Progulon] ), Function ]]
prn_combinations[, Function_2 := prn_annot[ match( prn_combinations[, V2], prn_annot[, Progulon] ), Function ]]

# Rearrange progulons (via factor levels) in clustered order
prn_combinations[, Function_1 := factor(Function_1, levels = new_prn_order )]
prn_combinations[, Function_2 := factor(Function_2, levels = new_prn_order )]

# Create the plot
p <- ggplot(prn_combinations, aes(Function_1, Function_2, size = overlap))+
      geom_point()+
      geom_vline(xintercept = seq(0.5,42.5,1), size = 0.25, colour = "grey80")+
      geom_hline(yintercept = seq(0.5,42.5,1), size = 0.25, colour = "grey80")+
      scale_size_continuous(name = "# proteins", limits = c(1,500), range = c(0,4), breaks = c(10,50,100,200))+
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6, angle = 90, hjust=1, vjust = 0.5))

p1 <- p + theme( legend.position = "none")
p2 <- p + theme( legend.position = "right")

# Save the plot
#ggsave("Progulon_overlap.pdf", p1, width = 10.5, height = 10.5, units = "cm")
ggsave("Progulon_overlap_legend.pdf", p2, width = 9.5, height = 7, units = "cm")


















