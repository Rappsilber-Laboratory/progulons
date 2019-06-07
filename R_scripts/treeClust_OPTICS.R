library(treeClust); library(dbscan); library(data.table)

# Define necessary functions
percentNA <- function(x){ sum(is.na(x)) / (nrow(x)*ncol(x)/100) }  # To calculate percentage of missing values
count_features <- function(x){ sum( !is.na(x) ) }                  # To calculate number of SILAC ratios per protein

# Load and prepare ProteomeHD (no normalisation necessary)
ProHD <- read.csv("ProteomeHD_v1.csv", stringsAsFactors=FALSE)
rownames(ProHD) <- ProHD$Majority.protein.IDs       # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]    # Keep only columns with SILAC ratios
percentNA(ProHD)

feature_count <- apply(ProHD, 1, count_features)    # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 95 ,]              # Discard proteins detected in fewer than 95 experiments
percentNA(ProHD)

# Use treeClust to learn a dissimilarity matrix (~ 15 min)
set.seed(42)
tc_dist <- treeClust.dist(ProHD, d.num = 2, verbose = TRUE)

# Perform OPTICS clustering
OPTICS <- optics(tc_dist, eps = 0.47 , minPts = 4)           # Get the clustering order
OPTICS <- extractXi(OPTICS, xi = 0.0031)                     # Call the clusters
DT <- data.table(cl = OPTICS$cluster, ID = rownames(ProHD))  # Assign clusters to protein IDs

# Remove clusters that have less than 4 members
too_small_clusters <- DT[, .N, by = cl][ N < 4, cl]
DT <- DT[ !cl %in% too_small_clusters ]

# Remove proteins that are not assigned to any cluster
DT <- DT[ cl != 0 ]
DT <- DT[ order(cl) ]

# Output some stats
number_of_clusters <- DT[, .N, by = cl ][,.N]
median_cluster_size <- DT[, .N, by = cl ][, median(N)]
mean_cluster_size <- DT[, .N, by = cl ][, mean(N)]
cluster_size_min <- DT[, .N, by = cl ][, min(N)]
cluster_size_max <- DT[, .N, by = cl ][, max(N)]

# Also add protein and gene names
annotation <- read.csv("ProteomeHD_v1.csv", stringsAsFactors=FALSE)
annotation <- annotation[, c("Majority.protein.IDs", "Protein.names", "Gene.names")]
annotation <- as.data.table(annotation)
DT <- merge(annotation, DT, by.x="Majority.protein.IDs", by.y="ID", sort=FALSE)

# Write out clusters as training sets
size_of_largest_cl <- DT[ , .N, by = cl][, max(N)]
training_seeds <- data.table(place_holder_column = 1:size_of_largest_cl)   # Initialise result table
current_cl_names <- DT[, unique(cl) ]                                      # Current numbers have gaps from clusters
                                                                           # that were removed
final_cl_names <- paste("ts", formatC(1:number_of_clusters, width=2, flag="0"), sep="")

for(i in 1:number_of_clusters){
  temp__current_cl_name <- current_cl_names[i]                      # The name/ID of the current cluster
  temp_final_cl_name <- final_cl_names[i]                           # The name/ID the current cluster should have
  IDs <- DT[ cl == temp__current_cl_name , Majority.protein.IDs ]   # Get protein IDs in current cluster
  IDs <- IDs[ 1:size_of_largest_cl ]                                # Fill up with NA to match length of largest cluster
  training_seeds[, c(temp_final_cl_name) := IDs ]                   # Append cluster to result table
}
training_seeds[, place_holder_column := NULL ]  # Remove placeholder as no longer needed
fwrite(training_seeds, "Clusters_as_training_seeds.csv")

# Write out clusters with annotation as supplementary table
DT$Cluster_ID <- final_cl_names[ match( DT$cl, current_cl_names) ]  # Assign final cluster names
DT[, cl := NULL ]
names(DT) <- c("Uniprot_IDs", "Protein_names", "Gene_names", "Cluster_ID")
fwrite(DT, "Annotated_clusters_TableSX.csv")





