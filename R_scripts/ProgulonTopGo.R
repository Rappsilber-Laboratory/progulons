# Read in the necessary libraries
library(plyr); library(data.table)

# Read in ProteomeHD to determine the list of proteins that need to be tested
ProHD <- fread("ProteomeHD_v1.csv")
ProHD_ratios <- ProHD[, .SD, .SDcols = colnames(ProHD) %like% "Ratio"]      # Limit to data columns
feature_count <- apply( ProHD_ratios, 1, function(x){ sum( !is.na(x)) } )   # Calculate number of features per protein
ProHD <- ProHD[ feature_count >= 30 ,]                                      # Keep only proteins that were actually included in the Progulon analysis
ProHD_proteins <- unique( ProHD[, SimpleID_1 := gsub(";.+", "", Majority.protein.IDs)][, gsub("-.+", "", SimpleID_1) ] )  # Unique genes (protein isoforms removed)

# Download the human GO associations file (goa_human.gaf) from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN
# Then read them in and prep as follows:
all_GO <- fread("goa_human.gaf")    # Loading human protein GO associations
all_GO <- all_GO[, .(ID = V2, Qualifier = V4, GO_ID = V5, Aspect = V9)]                # Keep relevant columns only
all_GO <- all_GO[ Qualifier == ""]                                                     # Keep only associations without qualifier
all_GO <- all_GO[ ID %in% ProHD_proteins ]                                             # Keep only proteins that were part of our analysis

# Prepare the data for topGO
GO <- all_GO[ Aspect == "F" ]                                      # Extract the relevant GO aspect
GO <- unique(GO)                                                   # Discard duplicate annotations
geneID2GO <- dlply( GO, "ID", function(x){as.character(x$GO_ID)})  # GO annotation for all proteins in the analysis,
                                                                   # ... for which GO annotation is available
# Load the progulon data
prns <- fread("Progulons.csv")
f_prn_calling <- function(x){                                                 # Input is a data.table containing the score and the feature counts
  ifelse( x$Mean_RF_score >= 0.55, "yes" ,                                    # Keep all proteins that score >= 0.55
          ifelse( x$Mean_RF_score >= 0.5 & x$Feature_count >= 100 , "yes",    # Keep proteins between 0.5 and 0.55 only if they have enough feature counts
                  "no"))}        # Function assigning proteins to progulons based on RF score cut-off and feature counts
prns$prot_in_prn <- f_prn_calling(prns)  # Assign proteins to progulons
prns[, Protein_IDs := gsub(";.+", "", Protein_IDs)][, Protein_IDs := gsub("-.+", "", Protein_IDs)]           # Simplify protein IDs 
prns <- dlply( prns, "Progulon_ID", function(x){  unique( x[ x$prot_in_prn == "yes" , "Protein_IDs" ] ) })   # List of proteins in each progulon

# Turn progulons into allGenes input for topGO
# allGenes_list is a list of progulons, with each progulon being a named factor that can be used for topGO's allGenes parameter.
# The names are all proteins in ProHD for which progulon-membership has been tested AND for which we have GO annotations (the "universe"). 
# The factor levels indicate whether proteins are in a progulon or not
geneNames <- names(geneID2GO)
allGenes_list <- lapply(prns, function(x){  geneList <- factor(as.integer(geneNames %in% x))
                                            names(geneList) <- geneNames
                                            return(geneList) })

# Create pilot topGOdata object (which will be updated with new allGenes in each iteration below)
library(topGO)
GOdata <- new("topGOdata", ontology = "MF", allGenes = allGenes_list[[1]],
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

# Screen every progulon for GO term enrichment
pvalues <- data.frame(GO_term = usedGO(GOdata))
for(i in 1:length(allGenes_list)){
  GOdata <- updateGenes(GOdata, allGenes_list[[i]])                           # Update the pilot GOdata object with the real progulon to be tested
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")     # Calculate enrichment considering graph structure
  current_pvalue <- as.data.frame( score(result))                             # Get the pvalues for enrichment of all GO terms for this gene subset
  colnames(current_pvalue) <- names(allGenes_list)[i]                         # Attach the progulon ID
  pvalues <- merge(pvalues, current_pvalue, by.x="GO_term", by.y="row.names")
  print(i)
}

# Write out result
fwrite(pvalues, "TopGO_result_MF.csv")
















