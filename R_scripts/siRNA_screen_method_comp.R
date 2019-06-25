## This script creates the plots comparing siRNA screen scoring methods

# Load the necessary libraries
library(data.table); library(ggplot2)

# Note: Set working directory to the script file location 

# Load list of proteins that we followed up
siRNA_results <- fread("Candidates_siRNA_screen.csv")


## SD-based plot

# Get column order for SD-based scoring
my_levels_SD <- siRNA_results[, sum(SD_score, na.rm = TRUE), by = Candidate ][ order(-V1) , Candidate ]

# Set column order
siRNA_results[, Candidate := factor( Candidate, levels = my_levels_SD )]

pScreenSD <- ggplot( siRNA_results, aes(x = Candidate, y = SD_score, fill = Process ))+
              geom_bar( stat = "identity" )+
              geom_hline( yintercept = 13/2, linetype = "dashed", size = 0.25, colour = "grey50")+
              geom_hline( yintercept = 13/3, linetype = "dashed", size = 0.25, colour = "grey50")+
              scale_y_continuous( limits = c(0,13), breaks = seq(0,13,1))+
              ylab("combined siRNA screen score (SD-based)")+
              scale_fill_manual( values = c("chartreuse3", "yellow2", "mediumslateblue"))+
              theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                    axis.text.y=element_text(size=5), axis.text.x=element_text(size=5, angle = 90, hjust = 1), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
                    plot.background = element_blank(), axis.title.y=element_text(size=6), axis.title.x = element_blank(),
                    legend.title = element_blank(), legend.text = element_text(size=6), legend.position = "top",
                    panel.grid.major.y = element_line(colour = "grey50", size = 0.25, linetype = "dotted"))

# Plot result
ggsave("Screen_SD_bar.pdf", pScreenSD, width=15, height=7, units= "cm")


## RSA-based plot

# Get column order for SD-based scoring
my_levels_RSA <- siRNA_results[, sum(RSA_score, na.rm = TRUE), by = Candidate ][ order(-V1) , Candidate ]

# Set column order
siRNA_results[, Candidate := factor( Candidate, levels = my_levels_RSA )]

pScreenRSA <- ggplot( siRNA_results, aes(x = Candidate, y = RSA_score, fill = Process ))+
              geom_bar( stat = "identity" )+
              geom_hline( yintercept = 13/2, linetype = "dashed", size = 0.25, colour = "grey50")+
              geom_hline( yintercept = 13/3, linetype = "dashed", size = 0.25, colour = "grey50")+
              scale_y_continuous( limits = c(0,13), breaks = seq(0,13,1))+
              ylab("combined siRNA screen score (RSA-based)")+
              scale_fill_manual( values = c("chartreuse3", "yellow2", "mediumslateblue"))+
              theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                    axis.text.y=element_text(size=5), axis.text.x=element_text(size=5, angle = 90, hjust = 1), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
                    plot.background = element_blank(), axis.title.y=element_text(size=6), axis.title.x = element_blank(),
                    legend.title = element_blank(), legend.text = element_text(size=6), legend.position = "top",
                    panel.grid.major.y = element_line(colour = "grey50", size = 0.25, linetype = "dotted"))

# Plot result
ggsave("Screen_RSA_bar.pdf", pScreenRSA, width=15, height=7, units= "cm")


## Output high-confidence validation overlap
siRNA_results[, lapply(.SD, sum, na.rm = T), .SDcols = c("SD_score", "RSA_score"), by = Candidate][
  SD_score > (13/2) & RSA_score > (13/2) ]



