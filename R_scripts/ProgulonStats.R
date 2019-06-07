library(ggplot2); library(data.table)

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


#### Get some stats ####

# How many proteins are, on mean / median, in each progulon?
proteins_in_progulon <- prns[ prot_in_prn == "yes", .N, by = Progulon_ID]
mean(   proteins_in_progulon[,N] )
median( proteins_in_progulon[,N] )
min(    proteins_in_progulon[,N] )
max(    proteins_in_progulon[,N] )

# How many training proteins have been used on mean/median?
training_prots_per_prog <- prns[ Used_for_positive_training == "Used", .N, by = Progulon_ID]
mean(   training_prots_per_prog[,N] )
median( training_prots_per_prog[,N] )
min(    training_prots_per_prog[,N] )
max(    training_prots_per_prog[,N] )

# How many proteins have been assigned to how many progulons?
progs_per_protein <- prns[ prot_in_prn == "yes", .N, by = Protein_IDs]
progs_per_protein[,.N]                                      # Proteins found in any progulon
progs_per_protein[,.N] / prns[,length(unique(Protein_IDs))] # Fraction of all analysed proteins
mean(   progs_per_protein[,N] )
median( progs_per_protein[,N] )
min(    progs_per_protein[,N] )
max(    progs_per_protein[,N] )


#### Plot progulon membership ####
progs_per_protein[, category := ifelse(N <= 3, N, "> 3")]
progs_per_protein[, category := factor( category, levels = c("> 3", "3", "2", "1"))]

p <- ggplot(progs_per_protein, aes(x = 1, fill = category))+
      geom_bar()+
      ylab("# proteins")+
      scale_fill_manual(values = c("#7FDBFF", "#39CCCC", "#0074D9", "#001f3f"))+
      theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=7), 
            axis.text.y = element_text(size=6), axis.text.x = element_blank(), axis.line.y = element_line(size=0.25),
            axis.ticks.x = element_blank())

ggsave("Progs_per_protein.pdf", p, width = 5, height = 5, units = "cm")
