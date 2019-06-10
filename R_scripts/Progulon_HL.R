library(readxl); library(ggplot2); library(data.table); library(perm); library(gridExtra); library(grid)

## NOTE

# In this script we consider two different sources of mRNA half-lives - Duan et al (Sci Rep 2013) and Tani et al (Genome Research 2012)
# As the plots show, only the data from Tani et al show a nice correlation with the corresponding protein values and Tani et al also claim
# that their novel half-live measurement method is superior to the traditional one. We therefore only considered their data for the 
# final manuscript version

#### Prep Progulon associations ####

# Load the progulon data
prns <- fread("Progulons.csv")

# Function assigning proteins to progulons based on RF score cut-off and feature counts
f_prn_calling <- function(x){                                       # Input is a data.table containing the score and the feature counts
  ifelse( x$Mean_RF_score >= 0.55, "yes" ,                          # Keep all proteins that score >= 0.55
          ifelse( x$Mean_RF_score >= 0.5 & x$Feature_count >= 100 , "yes",  # Keep proteins between 0.5 and 0.55 only if they have enough feature counts
                  "no"))}                                                   # All other proteins are not in the progulon

# Assign proteins to progulons
prns$prot_in_prn <- f_prn_calling(prns)

# Simplify protein IDs
prns[, SimpleID := gsub(";.+", "", Protein_IDs) ][, SimpleID := gsub("-.+", "", SimpleID)]


#### Prep mRNA half-lives from Duan et al (Sci Rep 2013) - LCL cells ####

# Read in mRNA half-lifes (hr) from Table S3 from Duan et al (Sci Rep 2013)
dt <- read_excel(path = "srep01318-s2.xls", sheet = 3, skip = 1)
dt <- as.data.table(dt)
dt <- dt[, names(dt) == "Gene ID" | names(dt) %like% "GM", with=FALSE]   # Get relevant columns

# Write out Refseq mRNA IDs, convert to Uniprot IDs using Retrieve function at www.uniprot.org
write.table(dt$`Gene ID`, "temp_Gene_IDs.csv", row.names = FALSE, col.names = FALSE, sep=",")
converted <- fread("RefSeq_to_Uniprot.tab")   # Read converted results back in 
converted <- converted[, .( UniprotID = Entry, RefSeqID = `yourlist:M2018050948CF0A2DF181CEB7EC2BC48F8F5F3B059B9C8CW`) ]
converted <- converted[ !duplicated(UniprotID) ]  # Remove duplicate assignments 
converted <- converted[ !duplicated(RefSeqID) ]   # Remove duplicate assignments

# Assign Uniprot IDs to half-lives
dt$`Gene ID` <- gsub(" ", "", dt$`Gene ID`)       # Remove a space which is there due to loading from Excel file
dt <- dt[ !duplicated(`Gene ID`) ]                # Remove duplicates
dt$UniprotID <- converted$UniprotID[ match(dt$`Gene ID`, converted$RefSeqID) ]
dt <- dt[ complete.cases(dt) ]                    # Remove genes which haven't been mapped

# Average the half-lifes of replicates
dt[, GM07029 := rowMeans(.SD), .SDcols = c("GM07029-A1", "GM07029-A2", "GM07029-A3") ]
dt[, GM10835 := rowMeans(.SD), .SDcols = c("GM10835A1", "GM10835A2", "GM10835A3") ]
dt[, GM12813 := rowMeans(.SD), .SDcols = c("GM12813-A1", "GM12813-A1dup", "GM12813-A2", "GM12813-A2dup", "GM12813-A3dup") ]
dt[, c("GM07029-A1", "GM07029-A2", "GM07029-A3", "GM10835A1", "GM10835A2", "GM10835A3",
       "GM12813-A1", "GM12813-A1dup", "GM12813-A2", "GM12813-A2dup", "GM12813-A3dup") := NULL ]

# Average the half-lifes across the 7 LCLs
half_life_cols <- grep("GM", names(dt), value = TRUE)
dt <- dt[, .(UniprotID, mean_LCL_mRNA_HL = rowMeans(.SD)), .SDcols = half_life_cols]


#### Prep mRNA half-lives from Tani et al (Genome Research 2012) - HeLa cells ####

# Read in mRNA half-lifes (hr) from Table S1 from Tani et al (Genome Research 2012)
# You can download this file from here: https://genome.cshlp.org/content/suppl/2012/02/14/gr.130559.111.DC1/Tani_Supp_Tables_revised2.xls
dt2 <- read_excel(path = "Tani_Supp_Tables_revised2.xls", sheet = 1, skip = 3)
dt2 <- as.data.table(dt2)
dt2 <- dt2[ !is.na(`t1/2 (h)`) ]           # Keep only mRNAs with a half-life measurement
dt2[, RepName := gsub(",", "", RepName) ]  # Remove commas from mRNA name (this will fuse multi-name entries, but I would discard these anyway)

# Write out Refseq mRNA IDs, convert to Uniprot IDs using Retrieve function at www.uniprot.org
write.table(dt2$RepName, "temp_Gene_IDs2.csv", row.names = FALSE, col.names = FALSE, sep=",")
converted2 <- fread("RefSeq2_to_Uniprot.tab")   # Read converted results back in 
converted2 <- converted2[, .( UniprotID = Entry, RefSeqID = `yourlist:M20180509F725F458AC8690F874DD868E4ED79B88FE4F7BJ`) ]
converted2 <- converted2[ !duplicated(UniprotID) ]  # Remove duplicate assignments 
converted2 <- converted2[ !duplicated(RefSeqID) ]   # Remove duplicate assignments

# Assign Uniprot IDs to half-lives
dt2$UniprotID <- converted2$UniprotID[ match(dt2$RepName, converted2$RefSeqID) ]
dt2 <- dt2[ complete.cases(dt2) ]                  # Remove genes which haven't been mapped

# Rename columns and remove unnessary columns
dt2 <- dt2[, .(UniprotID, HeLa_mRNA_HL = `t1/2 (h)`) ]

# Combine mRNA half-lives into one table
HL <- merge(dt, dt2, all = TRUE)

# Restrict to genes that have been mapped to at least one progulon, in order to avoid any bias
# derived from a potential bias in assigning proteins to progulons
all_prn_prots <- prns[ prot_in_prn == "yes" , unique(SimpleID) ]
HL <- HL[ UniprotID %in% all_prn_prots ]

# Clear workspace
rm( list = ls()[! ls() %in% c("HL", "prns")] )  

# Save the HeLa mRNA HLs for later comparison with RPE1 protein HLs
HeLa_mRNA_HLs <- HL[, .(UniprotID, HeLa_mRNA_HL)]
HeLa_mRNA_HLs <- HeLa_mRNA_HLs[ complete.cases(HeLa_mRNA_HLs) ]


#### Permutation test of mRNA half-lives in LCL cells ####

# Focus on the LCL half-lives
HL_LCL <- HL[ !is.na(mean_LCL_mRNA_HL) , .(UniprotID, mean_LCL_mRNA_HL) ]

# Test if genes within the same progulon have different mRNA half-lives than expected by random chance
HL_LCL_results <- data.table()                                             # Initialise result table
for(i in unique(prns$Progulon_ID)){                                        # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Proteins of the current progulon
     in_prn <- HL_LCL[  UniprotID %in% prn_prots , mean_LCL_mRNA_HL ]      # Their half-lives
 not_in_prn <- HL_LCL[ !UniprotID %in% prn_prots , mean_LCL_mRNA_HL ]      # Half-lives of the non-progulon proteins
temp_pvalue <- permTS(in_prn, not_in_prn, alternative = "two.sided",       # Calculate the pvalue by randomization test
                      method = "exact.mc",                                 # Using 10K Monte-Carlo replications
                      control = permControl(nmc = 10000, setSEED = FALSE,
                                            tsmethod = "abs"))$p.value
     temp_dt <- data.table( Progulon_ID = i,                               # Assemble results into a table
                            mean_mRNA_HL_LCL = mean(in_prn),
                            N_proteins_HL_LCL = length(in_prn),
                            pvalue_mRNA_HL_LCL = temp_pvalue)
     HL_LCL_results <- rbind(HL_LCL_results, temp_dt)                      # Append results per progulons
     print(i)
}
  

#### Permutation test of mRNA half-lives in HeLa cells ####

# Focus on the HeLa half-lives
HL_HeLa <- HL[ !is.na(HeLa_mRNA_HL) , .(UniprotID, HeLa_mRNA_HL) ]

# Test if genes within the same progulon have different mRNA half-lives than expected by random chance
HL_HeLa_results <- data.table()                                            # Initialise result table
for(i in unique(prns$Progulon_ID)){                                        # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Proteins of the current progulon
     in_prn <- HL_HeLa[  UniprotID %in% prn_prots , HeLa_mRNA_HL ]         # Their half-lives
 not_in_prn <- HL_HeLa[ !UniprotID %in% prn_prots , HeLa_mRNA_HL ]         # Half-lives of the non-progulon proteins
temp_pvalue <- permTS(in_prn, not_in_prn, alternative = "two.sided",       # Calculate the pvalue by randomization test
                      method = "exact.mc",                                 # Using 10K Monte-Carlo replications
                      control = permControl(nmc = 10000, setSEED = FALSE,
                                              tsmethod = "abs"))$p.value
  temp_dt <- data.table( Progulon_ID = i,                                  # Assemble results into a table
                         mean_mRNA_HL_HeLa = mean(in_prn),
                         N_proteins_HL_HeLa = length(in_prn),
                         pvalue_mRNA_HL_HeLa = temp_pvalue)
  HL_HeLa_results <- rbind(HL_HeLa_results, temp_dt)                       # Append results per progulons
  print(i)
}


# Merge LCL and HeLa cell mRNA half-lives per progulon
mRNA_HLs <- merge( HL_LCL_results, HL_HeLa_results, by = "Progulon_ID")

# Clear workspace
rm( list = ls()[! ls() %in% c("mRNA_HLs", "prns", "HeLa_mRNA_HLs")] )  


#### Prep protein half-lives from McShane et al (Cell 2016) - RPE1 cells ####

# Read in protein half-lifes (hr) from Table S4 from McShane et al (Cell 2016)
dt <- read_excel(path = "mmc4.xlsx", sheet = 1)
dt <- as.data.table(dt)

# Select and rename relevant columns
HL <- dt[, .(Uniprot_IDs = `Protein IDs (Uniprot)`, 
             protein_HL_RPE_ED = `Half-life (exponential) 1-state-model [h]`,
             protein_delta_score = `Δ-score` ,
             protein_degradation_profile = `Degradation profile`) ]

# Add simplified (non-isoformed protein IDs)
HL[, SimpleID := gsub(";.+", "", Uniprot_IDs) ][, SimpleID := gsub("-.+", "", SimpleID)]

# Remove duplicate genes
HL <- HL[ !duplicated(SimpleID) ]

# Restrict to genes that have been mapped to at least one progulon, in order to avoid any bias
# derived from a potential bias in assigning proteins to progulons
all_prn_prots <- prns[ prot_in_prn == "yes" , unique(SimpleID) ]
HL <- HL[ SimpleID %in% all_prn_prots ]

# Half-lives of > 300 hours are not accurate, so I remove them
HL[ protein_HL_RPE_ED == "> 300",  protein_HL_RPE_ED := NA ]

# Turn HL into numeric vectors
HL[,  protein_HL_RPE_ED := as.numeric(  protein_HL_RPE_ED ) ]

# Clear workspace
rm( list = ls()[! ls() %in% c("HL", "prns", "mRNA_HLs", "HeLa_mRNA_HLs")] )  

# Save protein HLs to compare to mRNA HLs later
RPE1_protein_HLs <- HL[, .(Uniprot_IDs, protein_HL_RPE_ED, protein_degradation_profile)]
RPE1_protein_HLs <- RPE1_protein_HLs[ complete.cases(RPE1_protein_HLs) ]


#### Comparison of mRNA and protein half-lives on a per-gene basis ####

# Simplify protein IDs (remove isoform info)
RPE1_protein_HLs[, SimpleID := gsub("-", "", Uniprot_IDs)]

# Remove any duplicates that may have arisen from that
RPE1_protein_HLs <- RPE1_protein_HLs[ !duplicated(SimpleID) ]

# Merge with mRNA data
RPE1_HeLa_mRNA_prot_HL <- merge(RPE1_protein_HLs, HeLa_mRNA_HLs, by.x = "SimpleID", by.y = "UniprotID")

# Calculate significance of correlation and turn into plot labels
sigRHOest <- RPE1_HeLa_mRNA_prot_HL[, cor.test(protein_HL_RPE_ED, HeLa_mRNA_HL, method = "spearman" )$estimate ] 
sigRHOpva <- RPE1_HeLa_mRNA_prot_HL[, cor.test(protein_HL_RPE_ED, HeLa_mRNA_HL, method = "spearman" )$p.value  ] 
sigRHO <- paste("RHO", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Plot it
pHLcomp <- ggplot(RPE1_HeLa_mRNA_prot_HL, aes(HeLa_mRNA_HL, protein_HL_RPE_ED))+
            geom_point( alpha = 0.5, size = 0.5 )+
            scale_x_continuous( limits = c(0,25), expand = c(0,0), breaks = seq(0,25,5))+
            scale_y_continuous( limits = c(0,300), expand = c(0,0), breaks = seq(0,300,50))+
            xlab("mRNA half-life [h] - HeLa cells")+
            ylab("protein half-life [h] - RPE1 cells")+
            annotate(geom = "text", x = 6, y = 80, size = 2, hjust = 0, label = sigRHO )+      
            theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
                  legend.position = "none", axis.ticks = element_line(size = 0.25))
pHLcomp
ggsave("Indi_gene_HL_mRNA_HeLa_vs_protein.pdf", pHLcomp,
       width = 4.9, height = 4.9, units = "cm")

# Clear workspace
rm( list = ls()[! ls() %in% c("HL", "prns", "mRNA_HLs")] )  

#### Permutation test of protein half-lives in RPE1 cells ####

# Focus on proteins whose (exponentially degraded) half-lives were quantified
HL_ED <- HL[ complete.cases(HL) ]

# Test if genes within the same progulon have different protein half-lives than expected by random chance
HL_RPE_results <- data.table()                                             # Initialise result table
for(i in unique(prns$Progulon_ID)){                                        # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Proteins of the current progulon
     in_prn <- HL_ED[  SimpleID %in% prn_prots , protein_HL_RPE_ED ]       # Their half-lives
 not_in_prn <- HL_ED[ !SimpleID %in% prn_prots , protein_HL_RPE_ED ]       # Half-lives of the non-progulon proteins
temp_pvalue <- permTS(in_prn, not_in_prn, alternative = "two.sided",       # Calculate the pvalue by randomization test
                        method = "exact.mc",                               # Using 10K Monte-Carlo replications
                        control = permControl(nmc = 10000, setSEED = FALSE,
                                              tsmethod = "abs"))$p.value
  temp_dt <- data.table( Progulon_ID = i,                                  # Assemble results into a table
                         mean_protein_HL_RPE = mean(in_prn),
                         N_proteins_HL_RPE = length(in_prn),
                         pvalue_protein_HL_RPE = temp_pvalue)
  HL_RPE_results <- rbind(HL_RPE_results, temp_dt)                        # Append results per progulons
  print(i)
}


#### Annotate protein degradation profiles of progulons #####

NED_results <- data.table()                                                # Initialise result table
for(i in unique(prns$Progulon_ID)){                                        # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Proteins of the current progulon
  med_DS <- HL[ SimpleID %in% prn_prots , median(protein_delta_score) ]    # Quantify extent of deviation from exponential degradation
  N_proteins_degradation_profile <- HL[ SimpleID %in% prn_prots , .N]      # Number of proteins for which degradation profiles were available
  pct_NEDs <- HL[ SimpleID %in% prn_prots ,                                # Percentage of proteins that are non-exponentially degraded
                  sum( protein_degradation_profile == "NED" )/.N*100]  
  
  # Assess statistical signficance of NED enrichment / depletion in progulon by Fisher's Exact test
       N_NED_in_prog <- HL[   SimpleID %in% prn_prots  & protein_degradation_profile == "NED" , .N]
   N_NED_not_in_prog <- HL[ !(SimpleID %in% prn_prots) & protein_degradation_profile == "NED" , .N]
    N_notNED_in_prog <- HL[   SimpleID %in% prn_prots  & protein_degradation_profile != "NED" , .N]
N_notNED_not_in_prog <- HL[ !(SimpleID %in% prn_prots) & protein_degradation_profile != "NED" , .N]
  fisher_dt <- rbind( data.frame( NED = N_NED_in_prog,     not_NED = N_notNED_in_prog ),
                      data.frame( NED = N_NED_not_in_prog, not_NED = N_notNED_not_in_prog ))
  rownames(fisher_dt) <- c("in_prog", "not_in_prog")
  pval_NED_enrichment_or_depletion <- fisher.test(fisher_dt)$p.value
  
  temp_dt <- data.table( Progulon_ID = i,                                  # Assemble results into a table
                         median_protein_delta_score = med_DS,
                         N_proteins_degradation_profile,
                         pct_NEDs,
                         pval_NED_enrichment_or_depletion)
  NED_results <- rbind(NED_results, temp_dt)                               # Append results per progulon
  print(i)
}

# Save the percentage of NEDs in the whole dataset, for the volcano plot p5
pct_NED_overall <- HL[, sum(protein_degradation_profile == "NED")/.N*100]


#### Annotate membership to heteromeric multiprotein complexes ####

# Download and load the CORUM multiprotein complexes (http://mips.helmholtz-muenchen.de/corum/#download)
CORUM <- fread("allComplexes.txt")

# Restrict to human complexes
CORUM <- CORUM[ Organism == "Human" ]

# Get all protein complexes with one subunit per column
CORUM <- CORUM[, tstrsplit(`subunits(UniProt IDs)`, ";") ]

# Remove mono- or homomeric protein complexes
subunit_count <- apply(CORUM, 1, function(x){ length(unique(na.omit(x)))  })
CORUM <- CORUM[ subunit_count > 1 ]

# Get vector of all proteins annotated as complex subunits in CORUM
CORUM <- na.omit( unique( unlist(CORUM) ))

# Remove isoforms and duplicates arising from that
CORUM <- gsub("-", "", CORUM)
CORUM <- unique( CORUM )

# Annotate the percentage of heteromeric protein complex subunits per progulon
CORUM_results <- data.table()                                              # Initialise result table
for(i in unique(prns$Progulon_ID)){                                        # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Proteins of the current progulon
  pct_complex <- sum( prn_prots %in% CORUM ) / length(prn_prots) * 100     # Which percentage of these proteins are subunits of heteromeric protein complexes in CORUM
  CORUM_results <- rbind(CORUM_results,
                         data.table( Progulon_ID = i, pct_complex))        # Append results per progulon
}


#### Merge all data and write out supplementary table ####

# Merge data
DT <- merge( mRNA_HLs, HL_RPE_results )
DT <- merge( DT, NED_results )
DT <- merge( DT, CORUM_results )

# Add manually annotated protein function
annot <- fread("Manual_Progulon_Annotation.csv")
annot <- annot[, .(Progulon_ID = Progulon, Progulon_function = Function)]
DT <- merge(DT, annot)

# Write out supplementary table
fwrite(DT, "Progulon_half_lives.csv")

# Clear workspace
rm( list = ls()[! ls() %in% c("DT", "prns", "pct_NED_overall")] )  


#### Create plot 1: mRNA (HeLa) vs protein progulon stability ####

# Calculate significance of correlation and turn into plot labels
sigRHOest <- DT[, cor.test(mean_mRNA_HL_HeLa, mean_protein_HL_RPE, method = "spearman" )$estimate ] 
sigRHOpva <- DT[, cor.test(mean_mRNA_HL_HeLa, mean_protein_HL_RPE, method = "spearman" )$p.value  ] 
sigRHO <- paste("RHO", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))


p1 <- ggplot(DT, aes(x = mean_mRNA_HL_HeLa, y = mean_protein_HL_RPE))+
      geom_hline(yintercept = median(DT$mean_protein_HL_RPE), size = 0.25, linetype = "dashed", colour = "grey50")+
      geom_vline(xintercept = median(DT$mean_mRNA_HL_HeLa), size = 0.25, linetype = "dashed", colour = "grey50")+
      geom_point( alpha = 0.5, size = 1 )+
      geom_text(data = DT[ Progulon_ID %in% c("PRN12", "PRN04", "PRN07", "PRN30", "PRN33", "PRN16")],
                aes(label = Progulon_function), size = 2, colour = "black", hjust = -0.1)+
      scale_x_continuous( limits = c(5.5,12.5), expand = c(0,0), breaks = seq(0,20,2))+
      scale_y_continuous( limits = c(18,88), expand = c(0,0), breaks = seq(0,100,20))+
      xlab("average mRNA half-life [h]")+
      ylab("average protein half-life [h]")+
      annotate(geom = "text", x = 6, y = 80, size = 2, hjust = 0, label = sigRHO )+      
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
            legend.position = "none", axis.ticks = element_line(size = 0.25))

p1
ggsave("HL_mRNA_HeLa_vs_protein.pdf", p1,
       width = 4.9, height = 4.9, units = "cm")


#### Create plot 2: mRNA (LCL) vs protein progulon stability ####

# Note that Tani et al (HeLa cell data) used a different method to determine mRNA half-lives and claim
# that their method is superior to the standard one used by Duan et al (LCL data). Indeed, I find that 
# there is no correlation between mRNA and protein half-lives when using Duan et al's data, suggesting
# that Tani et al's data are more likely to be correct. Therefore, I only print out plot p2 and p3 as
# reminder of how the data look like when the LCL data are being used (combined plot = p4)

# Calculate significance of correlation and turn into plot labels
sigPCCest <- DT[, cor.test(mean_mRNA_HL_LCL, mean_protein_HL_RPE, method = "pearson" )$estimate ] 
sigPCCpva <- DT[, cor.test(mean_mRNA_HL_LCL, mean_protein_HL_RPE, method = "pearson" )$p.value  ] 
sigPCC <- paste("PCC", round(sigPCCest, 2), ", p value", signif(sigPCCpva, 2))

p2 <- ggplot(DT, aes(x = mean_mRNA_HL_LCL, y = mean_protein_HL_RPE))+
        geom_hline(yintercept = median(DT$mean_protein_HL_RPE), size = 0.25, linetype = "dashed", colour = "grey50")+
        geom_vline(xintercept = median(DT$mean_mRNA_HL_LCL), size = 0.25, linetype = "dashed", colour = "grey50")+
        geom_point( alpha = 0.5, size = 1 )+
        # geom_text(data = DT[ Progulon_ID %in% c("PRN12", "PRN04", "PRN07", "PRN30", "PRN33", "PRN16")],
        #           aes(label = Progulon_function), size = 2, colour = "black", hjust = -0.1)+
        scale_x_continuous( limits = c(4,10), expand = c(0,0), breaks = seq(0,20,2))+
        scale_y_continuous( limits = c(18,88), expand = c(0,0), breaks = seq(0,100,20))+
        xlab("average mRNA half-life [h]")+
        ylab("average protein half-live [h]")+
        annotate(geom = "text", x = 4.5, y = 75, size = 2, hjust = 0, label = sigPCC )+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))


#### Create plot 3: mRNA (LCL) vs mRNA (HeLa) ####

# In addition to LCL mRNA data not correlating with the protein half-lives, there is only a weak correlation
# between half-lives measured by Tani et al and Duan et al. While this could be cell-line specific, it may be
# an effect of the methods as well.

# Calculate significance of correlation and turn into plot labels
sigPCCest <- DT[, cor.test(mean_mRNA_HL_LCL, mean_mRNA_HL_HeLa, method = "pearson" )$estimate ] 
sigPCCpva <- DT[, cor.test(mean_mRNA_HL_LCL, mean_mRNA_HL_HeLa, method = "pearson" )$p.value  ] 
sigPCC <- paste("PCC", round(sigPCCest, 2), ", p value", signif(sigPCCpva, 2))

p3 <- ggplot(DT, aes(x = mean_mRNA_HL_LCL, y = mean_mRNA_HL_HeLa))+
        geom_smooth( method = "lm" )+
        geom_point( alpha = 0.5, size = 1 )+
        xlab("average mRNA half-life in LCLs [h]")+
        ylab("average mRNA half-life in HeLas [h")+
        annotate(geom = "text", x = 4.5, y = 11, size = 2, hjust = 0, label = sigPCC )+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))
      
# Combine p2 and p3 and write out as a common "LCL plot"
p4 <- grid.arrange(p3, p2, nrow = 1)

ggsave("HL_mRNA_LCL_plot.pdf", p4,
       width = 14, height = 7, units = "cm")


#### Create plot 4: progulons enriched / depleted in NEDs ####

# Create a volcano plot
p5 <- ggplot(DT, aes(x = pct_NEDs, y = -log10(pval_NED_enrichment_or_depletion)))+
        geom_vline(xintercept = pct_NED_overall, size = 0.25, linetype = "dashed", colour = "grey50")+
        geom_point(alpha = 0.5, size = 1)+
        geom_point(data = DT[ pval_NED_enrichment_or_depletion < 0.0001 ], size = 1, colour = "red")+
        geom_text( data = DT[ pval_NED_enrichment_or_depletion < 0.0001 ], size = 2, colour = "red",
                   aes( label = Progulon_function ), hjust = 0)+
        scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
        xlab("% NED proteins")+
        ylab("-log10 p-value")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              axis.ticks = element_line(size = 0.25))

p5
ggsave("PRN_NEDs_volcano.pdf", p5,
       width = 4.9, height = 4.9, units = "cm")


#### Create plot 5: NEDs and CORUM multiprotein complexes ####

# Does the enrichment / depletion of NEDs in certain progulons simply reflect the fact that these
# progulons contain more or less subunits of multiprotein complexes?

# Calculate significance of correlation and turn into plot labels
sigPCCest <- DT[, cor.test(pct_NEDs, pct_complex, method = "pearson" )$estimate ] 
sigPCCpva <- DT[, cor.test(pct_NEDs, pct_complex, method = "pearson" )$p.value  ] 
sigRHOest <- DT[, cor.test(pct_NEDs, pct_complex, method = "spearman")$estimate ]
sigRHOpva <- DT[, cor.test(pct_NEDs, pct_complex, method = "spearman")$p.value  ]
sigPCC <- paste("PCC", round(sigPCCest, 2), ", p value", signif(sigPCCpva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Because the correlation is driven mainly by the ribosome and the mitochondrial ribosome,
# which have > 90% complex subunits and > 40% NEDs, I also test it without those two progulons
sigPCCest <- DT[ pct_complex < 90, cor.test(pct_NEDs, pct_complex, method = "pearson" )$estimate ] 
sigPCCpva <- DT[ pct_complex < 90, cor.test(pct_NEDs, pct_complex, method = "pearson" )$p.value  ] 
sigRHOest <- DT[ pct_complex < 90, cor.test(pct_NEDs, pct_complex, method = "spearman")$estimate ]
sigRHOpva <- DT[ pct_complex < 90, cor.test(pct_NEDs, pct_complex, method = "spearman")$p.value  ]
sigPCC2 <- paste("PCC without outliers", round(sigPCCest, 2), ", p value", signif(sigPCCpva, 2))
sigRHO2 <- paste("rho without outliers", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Make the plot
p6 <- ggplot(DT, aes(x = pct_NEDs, y = pct_complex))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_text(data = DT[ pct_complex > 90 ], aes(label = Progulon_function), size = 2, colour = "blue", hjust = 0)+
        scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
        scale_y_continuous( limits = c(0,100), expand = c(0,0), breaks = seq(0,100,20))+
        xlab("% NED proteins")+
        ylab("% subunits of heteromeric protein complexes")+
        annotate(geom = "text", x = 45, y = 50, size = 2, hjust = 0, label = sigPCC )+      
        annotate(geom = "text", x = 45, y = 40, size = 2, hjust = 0, label = sigRHO )+    
        annotate(geom = "text", x = 45, y = 30, size = 2, hjust = 0, label = sigPCC2)+      
        annotate(geom = "text", x = 45, y = 20, size = 2, hjust = 0, label = sigRHO2)+   
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p6 
ggsave("NEDs_vs_CORUM.pdf", p6,
       width = 4.9, height = 4.9, units = "cm")


#### NEDs vs (co)expression ####

# To compare the percentage of NEDs per progulon with (co)expression behaviour,  load the three different coexpression datasets 
# we have available for that. NOTE: These are three separate tables here, but they correspond to the three different tabs in
# the manuscript's Supplementary Table 2
  coexpr_dt_LCL <- fread("Rna_Pro_PRN_LCL.csv")
 coexpr_dt_Haas <- fread("Haas_Rna_Pro_PRN.csv")
coexpr_dt_Mouse <- fread("Mouse_Rna_Pro_PRN.csv")

# Select revelant columns only
  coexpr_dt_LCL <- coexpr_dt_LCL[,   .(Progulon_ID, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC)]
 coexpr_dt_Haas <- coexpr_dt_Haas[,  .(Progulon_ID, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC)]
coexpr_dt_Mouse <- coexpr_dt_Mouse[, .(Progulon_ID, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC)]

# Merge with the progulon HL data
  coexpr_dt_LCL <- merge(DT, coexpr_dt_LCL)
 coexpr_dt_Haas <- merge(DT, coexpr_dt_Haas)
coexpr_dt_Mouse <- merge(DT, coexpr_dt_Mouse)

# Create a function to calculate correlations and their significance, with and without the two ribosome PRN outliers,
# both PCC and RHO
my_cor_function <- function(x){ 
  PCC_RNARNA <- x[, cor.test(pct_NEDs, med_RNA_RNA_PCC, method = "pearson") ]
  PCC_propro <- x[, cor.test(pct_NEDs, med_pro_pro_PCC, method = "pearson") ]
  PCC_RNApro <- x[, cor.test(pct_NEDs, med_RNA_pro_PCC, method = "pearson") ] 
  RHO_RNARNA <- x[, cor.test(pct_NEDs, med_RNA_RNA_PCC, method = "spearman") ]
  RHO_propro <- x[, cor.test(pct_NEDs, med_pro_pro_PCC, method = "spearman") ]
  RHO_RNApro <- x[, cor.test(pct_NEDs, med_RNA_pro_PCC, method = "spearman") ]
  PCC_RNARNA_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_RNA_RNA_PCC, method = "pearson") ]
  PCC_propro_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_pro_pro_PCC, method = "pearson") ]
  PCC_RNApro_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_RNA_pro_PCC, method = "pearson") ] 
  RHO_RNARNA_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_RNA_RNA_PCC, method = "spearman") ]
  RHO_propro_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_pro_pro_PCC, method = "spearman") ]
  RHO_RNApro_noOut <- x[ pct_NEDs < 40 , cor.test(pct_NEDs, med_RNA_pro_PCC, method = "spearman") ]
                 
  cor_dt <- rbind(               
  data.table( Type = "PCC", Format = "RNARNA", Outlier = "With", Value = PCC_RNARNA$estimate, p_value = PCC_RNARNA$p.value),
  data.table( Type = "PCC", Format = "propro", Outlier = "With", Value = PCC_propro$estimate, p_value = PCC_propro$p.value),
  data.table( Type = "PCC", Format = "RNApro", Outlier = "With", Value = PCC_RNApro$estimate, p_value = PCC_RNApro$p.value),
  data.table( Type = "RHO", Format = "RNARNA", Outlier = "With", Value = RHO_RNARNA$estimate, p_value = RHO_RNARNA$p.value),
  data.table( Type = "RHO", Format = "propro", Outlier = "With", Value = RHO_propro$estimate, p_value = RHO_propro$p.value),
  data.table( Type = "RHO", Format = "RNApro", Outlier = "With", Value = RHO_RNApro$estimate, p_value = RHO_RNApro$p.value),
  data.table( Type = "PCC", Format = "RNARNA", Outlier = "Without", Value = PCC_RNARNA_noOut$estimate, p_value = PCC_RNARNA_noOut$p.value),
  data.table( Type = "PCC", Format = "propro", Outlier = "Without", Value = PCC_propro_noOut$estimate, p_value = PCC_propro_noOut$p.value),
  data.table( Type = "PCC", Format = "RNApro", Outlier = "Without", Value = PCC_RNApro_noOut$estimate, p_value = PCC_RNApro_noOut$p.value),
  data.table( Type = "RHO", Format = "RNARNA", Outlier = "Without", Value = RHO_RNARNA_noOut$estimate, p_value = RHO_RNARNA_noOut$p.value),
  data.table( Type = "RHO", Format = "propro", Outlier = "Without", Value = RHO_propro_noOut$estimate, p_value = RHO_propro_noOut$p.value),
  data.table( Type = "RHO", Format = "RNApro", Outlier = "Without", Value = RHO_RNApro_noOut$estimate, p_value = RHO_RNApro_noOut$p.value))
  
  return(cor_dt)
  
  }

# Apply function to calculate correlations and p-values for each of the three datasets
  cor_coexpr_dt_LCL <- my_cor_function( coexpr_dt_LCL )
 cor_coexpr_dt_Haas <- my_cor_function( coexpr_dt_Haas )
cor_coexpr_dt_Mouse <- my_cor_function( coexpr_dt_Mouse )

# Append dataset name
cor_coexpr_dt_LCL[, Dataset := "LCL" ]
cor_coexpr_dt_Haas[, Dataset := "Haas" ]
cor_coexpr_dt_Mouse[, Dataset := "Mouse" ]

# Combine into one table
my_cors <- rbind(cor_coexpr_dt_LCL, cor_coexpr_dt_Haas, cor_coexpr_dt_Mouse)


#### Create plot 6: NEDs vs cor_coexpression scatterplots #### 

# Create LCL plotting dataset
LCL_plot_dt <- melt( coexpr_dt_LCL[, .(Progulon_ID, Progulon_function, pct_NEDs, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC) ],
                     measure.vars = c("med_RNA_RNA_PCC", "med_pro_pro_PCC", "med_RNA_pro_PCC"))

# Append LCL correlations and p-values
LCL_cors <- my_cors[ Dataset == "LCL" ]
LCL_cors[ Format == "RNARNA" , Format := "med_RNA_RNA_PCC" ]
LCL_cors[ Format == "propro" , Format := "med_pro_pro_PCC" ]
LCL_cors[ Format == "RNApro" , Format := "med_RNA_pro_PCC" ]
LCL_cors[, Outlier := ifelse( Outlier == "With", "", "(w/o outliers)") ]
LCL_cors[, Label := paste(Type, round(Value, 2), ", p", signif(p_value, 2), Outlier)]
LCL_cors <- LCL_cors[, .(Label = paste(Label, collapse = "\n")), Format]
LCL_plot_dt <- merge( LCL_plot_dt , LCL_cors, by.x = "variable", by.y = "Format")
LCL_plot_dt[ variable == "med_RNA_RNA_PCC" , variable := "mRNA - mRNA coexpression"]
LCL_plot_dt[ variable == "med_pro_pro_PCC" , variable := "protein - protein coexpression"]
LCL_plot_dt[ variable == "med_RNA_pro_PCC" , variable := "mRNA - protein correlation"]
LCL_plot_dt[, variable := factor(variable, levels = c("mRNA - mRNA coexpression", "protein - protein coexpression", "mRNA - protein correlation"))]

p7a <- ggplot( LCL_plot_dt, aes( x = pct_NEDs, y = value ))+
        facet_wrap(~variable)+
        geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
        geom_smooth( data = LCL_plot_dt[ pct_NEDs < 40 ], method = "lm", se = FALSE, size = 0.25, 
                     colour = "grey50", fullrange = TRUE)+
        geom_point( alpha = 0.5, size = 1 )+
        geom_point( data = LCL_plot_dt[ pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
        geom_text(data = LCL_plot_dt[ pct_NEDs > 40 ], aes(label = Progulon_function), size = 2, colour = "dodgerblue2", hjust = -0.1)+
        scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
        scale_y_continuous( limits = c(-0.2,1), expand = c(0,0), breaks = seq(-0.2,1,0.2))+
        xlab("% NED proteins")+
        ylab("median correlation [PCC]")+
        geom_text( data = unique( LCL_plot_dt[, .(variable, Label)]), aes( label = Label), x = 10, y = 0.8, size = 2, hjust = 0)+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              axis.ticks = element_line(size = 0.25), strip.text = element_text(size = 7),
              strip.background = element_rect( fill = "paleturquoise", colour = "black", size = 0.25))


# Create Haas plotting dataset
Haas_plot_dt <- melt( coexpr_dt_Haas[, .(Progulon_ID, Progulon_function, pct_NEDs, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC) ],
                      measure.vars = c("med_RNA_RNA_PCC", "med_pro_pro_PCC", "med_RNA_pro_PCC"))

# Append Haas correlations and p-values
Haas_cors <- my_cors[ Dataset == "Haas" ]
Haas_cors[ Format == "RNARNA" , Format := "med_RNA_RNA_PCC" ]
Haas_cors[ Format == "propro" , Format := "med_pro_pro_PCC" ]
Haas_cors[ Format == "RNApro" , Format := "med_RNA_pro_PCC" ]
Haas_cors[, Outlier := ifelse( Outlier == "With", "", "(w/o outliers)") ]
Haas_cors[, Label := paste(Type, round(Value, 2), ", p", signif(p_value, 2), Outlier)]
Haas_cors <- Haas_cors[, .(Label = paste(Label, collapse = "\n")), Format]
Haas_plot_dt <- merge( Haas_plot_dt , Haas_cors, by.x = "variable", by.y = "Format")
Haas_plot_dt[ variable == "med_RNA_RNA_PCC" , variable := "mRNA - mRNA coexpression"]
Haas_plot_dt[ variable == "med_pro_pro_PCC" , variable := "protein - protein coexpression"]
Haas_plot_dt[ variable == "med_RNA_pro_PCC" , variable := "mRNA - protein correlation"]
Haas_plot_dt[, variable := factor(variable, levels = c("mRNA - mRNA coexpression", "protein - protein coexpression", "mRNA - protein correlation"))]

p7b <- ggplot( Haas_plot_dt, aes( x = pct_NEDs, y = value ))+
        facet_wrap(~variable)+
        geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
        geom_smooth( data = Haas_plot_dt[ pct_NEDs < 40 ], method = "lm", se = FALSE, size = 0.25, 
                     colour = "grey50", fullrange = TRUE)+
        geom_point( alpha = 0.5, size = 1 )+
        geom_point( data = Haas_plot_dt[ pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
        geom_text(data = Haas_plot_dt[ pct_NEDs > 40 ], aes(label = Progulon_function), size = 2, colour = "dodgerblue2", hjust = -0.1)+
        scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
        scale_y_continuous( limits = c(-0.2,1), expand = c(0,0), breaks = seq(-0.2,1,0.2))+
        xlab("% NED proteins")+
        ylab("median correlation [PCC]")+
        geom_text( data = unique( Haas_plot_dt[, .(variable, Label)]), aes( label = Label), x = 10, y = 0.8, size = 2, hjust = 0)+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              axis.ticks = element_line(size = 0.25), strip.text = element_text(size = 7),
              strip.background = element_rect( fill = "paleturquoise", colour = "black", size = 0.25))


# Create Mouse plotting dataset
Mouse_plot_dt <- melt( coexpr_dt_Mouse[, .(Progulon_ID, Progulon_function, pct_NEDs, med_RNA_RNA_PCC, med_pro_pro_PCC, med_RNA_pro_PCC) ],
                       measure.vars = c("med_RNA_RNA_PCC", "med_pro_pro_PCC", "med_RNA_pro_PCC"))
 
# Append Mouse correlations and p-values
Mouse_cors <- my_cors[ Dataset == "Mouse" ]
Mouse_cors[ Format == "RNARNA" , Format := "med_RNA_RNA_PCC" ]
Mouse_cors[ Format == "propro" , Format := "med_pro_pro_PCC" ]
Mouse_cors[ Format == "RNApro" , Format := "med_RNA_pro_PCC" ]
Mouse_cors[, Outlier := ifelse( Outlier == "With", "", "(w/o outliers)") ]
Mouse_cors[, Label := paste(Type, round(Value, 2), ", p", signif(p_value, 2), Outlier)]
Mouse_cors <- Mouse_cors[, .(Label = paste(Label, collapse = "\n")), Format]
Mouse_plot_dt <- merge( Mouse_plot_dt , Mouse_cors, by.x = "variable", by.y = "Format")
Mouse_plot_dt[ variable == "med_RNA_RNA_PCC" , variable := "mRNA - mRNA coexpression"]
Mouse_plot_dt[ variable == "med_pro_pro_PCC" , variable := "protein - protein coexpression"]
Mouse_plot_dt[ variable == "med_RNA_pro_PCC" , variable := "mRNA - protein correlation"]
Mouse_plot_dt[, variable := factor(variable, levels = c("mRNA - mRNA coexpression", "protein - protein coexpression", "mRNA - protein correlation"))]

p7c <- ggplot( Mouse_plot_dt, aes( x = pct_NEDs, y = value ))+
        facet_wrap(~variable)+
        geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
        geom_smooth( data = Mouse_plot_dt[ pct_NEDs < 40 ], method = "lm", se = FALSE, size = 0.25, 
                     colour = "grey50", fullrange = TRUE)+
        geom_point( alpha = 0.5, size = 1 )+
        geom_point( data = Mouse_plot_dt[ pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
        geom_text(data = Mouse_plot_dt[ pct_NEDs > 40 ], aes(label = Progulon_function), size = 2, colour = "dodgerblue2", hjust = -0.1)+
        scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
        scale_y_continuous( limits = c(-0.2,1), expand = c(0,0), breaks = seq(-0.2,1,0.2))+
        xlab("% NED proteins")+
        ylab("median correlation [PCC]")+
        geom_text( data = unique( Mouse_plot_dt[, .(variable, Label)]), aes( label = Label), x = 10, y = 0.8, size = 2, hjust = 0)+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              axis.ticks = element_line(size = 0.25), strip.text = element_text(size = 7),
              strip.background = element_rect( fill = "paleturquoise", colour = "black", size = 0.25))


# Combine the plots and print as one
p7 <- grid.arrange(p7a, p7b, p7c, nrow = 3)
ggsave("NEDs_vs_coexpression_scatter.pdf", p7,
       width = 15, height = 17, units = "cm")


#### Create plot 7: NEDs vs cor_coexpression correlation tile plot #### 

# Create ordered factors to impose plotting order
my_cors[, Format := factor(Format, levels = c("RNARNA", "propro", "RNApro"))]
my_cors[, Dataset := factor(Dataset, levels = c("Mouse", "Haas", "LCL"))]


p8 <- ggplot(my_cors, aes(x = Format, y = Dataset, fill = Value))+
       facet_grid(Outlier ~ Type)+
       geom_tile()+
       geom_text(data = my_cors[ p_value < 0.05 & p_value > 0.01 ], label = "*", size = 2)+
       geom_text(data = my_cors[ p_value < 0.01 & p_value > 0.001 ], label = "**", size = 2)+
       geom_text(data = my_cors[ p_value < 0.001 ], label = "***", size = 2)+
       scale_fill_gradient2( limits = c(-1, 1), name = "correlation")+
       theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.title.y = element_text(size = 6), axis.title.x = element_blank(),
             axis.text = element_text(size = 5, colour = "black"), axis.ticks = element_blank(),
             strip.background = element_rect(colour = "black", size = 0.25),
             strip.text = element_text(size = 7), legend.title = element_text(size = 7),
             legend.text = element_text(size = 6), legend.key.width = unit(0.3, "cm"))

p8
ggsave("NEDs_vs_coexpression_corr.pdf", p8,
       width = 7, height = 4.5, units = "cm")


#### Create plot 8: The main figure composite plot ####

# The p5 volcano plot is ready as is
# But the HeLa scatterplots need to be broken down into three separate plots first

p7b1 <- ggplot( Haas_plot_dt[ variable == "mRNA - mRNA coexpression"], aes( x = pct_NEDs, y = value ))+
          geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
          geom_point( alpha = 0.5, size = 1 )+
          geom_point( data = Haas_plot_dt[ variable == "mRNA - mRNA coexpression" & pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
          scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
          scale_y_continuous( limits = c(0,0.601), expand = c(0,0), breaks = seq(0,1,0.2))+
          xlab("% NED proteins")+
          ylab("median mRNA - mRNA coexpression [PCC]")+
          geom_text( data = unique( Haas_plot_dt[ variable == "mRNA - mRNA coexpression", .(variable, Label)]),
                     aes( label = Label), x = 10, y = 0.5, size = 2, hjust = 0)+      
          theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
                axis.ticks = element_line(size = 0.25))

p7b2 <- ggplot( Haas_plot_dt[ variable == "protein - protein coexpression"], aes( x = pct_NEDs, y = value ))+
          geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
          geom_point( alpha = 0.5, size = 1 )+
          geom_point( data = Haas_plot_dt[ variable == "protein - protein coexpression" & pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
          scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
          scale_y_continuous( limits = c(0,0.8), expand = c(0,0), breaks = seq(-0.2,1,0.2))+
          xlab("% NED proteins")+
          ylab("median protein - protein coexpression [PCC]")+
          geom_text( data = unique( Haas_plot_dt[ variable == "protein - protein coexpression", .(variable, Label)]),
                     aes( label = Label), x = 10, y = 0.6, size = 2, hjust = 0)+      
          theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
                axis.ticks = element_line(size = 0.25))

p7b3 <- ggplot( Haas_plot_dt[ variable == "mRNA - protein correlation"], aes( x = pct_NEDs, y = value ))+
          geom_smooth(method = "lm", se = FALSE, size = 0.25, colour = "dodgerblue2", fullrange = TRUE)+
          geom_point( alpha = 0.5, size = 1 )+
          geom_point( data = Haas_plot_dt[ variable == "mRNA - protein correlation" & pct_NEDs > 40 ], size = 1 , colour = "dodgerblue2")+
          scale_x_continuous( limits = c(0,75), expand = c(0,0), breaks = seq(0,100,25))+
          scale_y_continuous( limits = c(0,0.8), expand = c(0,0), breaks = seq(-0.2,1,0.2))+
          xlab("% NED proteins")+
          ylab("median mRNA - protein correlation [PCC]")+
          geom_text( data = unique( Haas_plot_dt[ variable == "mRNA - protein correlation", .(variable, Label)]),
                     aes( label = Label), x = 50, y = 0.6, size = 2, hjust = 0)+      
          theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
                axis.ticks = element_line(size = 0.25))

# To output a properly aligned plot grid, combine the gtables first
g5 <- ggplotGrob(p5)
g7b1 <- ggplotGrob(p7b1)
g7b2 <- ggplotGrob(p7b2)
g7b3 <- ggplotGrob(p7b3)
g <- cbind( rbind(g5, g7b1, size = "first"), rbind(g7b2, g7b3, size = "first"), size = "first")
grid.newpage()
grid.draw(g)

ggsave("Haas_NED_composite.pdf", g,
       width = 6.7, height = 6.5, units = "cm")


# In addition, output a plot where the y-axis of the volcano plot is limited to a -log10 pvalue of 10
# I will split and combine these plots in Inkscape, so I can create a discontinuous y-axis and avoid 
# lots of white space in the figure
p5b <- p5 + scale_y_continuous( limits = c(0,13), breaks = seq(0,15,5))
g5b <- ggplotGrob(p5b)
gb <- cbind( rbind(g5b, g7b1, size = "first"), rbind(g7b2, g7b3, size = "first"), size = "first")
grid.newpage()
grid.draw(gb)

ggsave("Haas_NED_composite_B.pdf", gb,
       width = 6.7, height = 6.5, units = "cm")



