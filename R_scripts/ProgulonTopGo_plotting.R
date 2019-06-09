library(data.table); library(ggplot2); library(viridis)

# Load the data
BP <- fread("TopGO_result_BP.csv")
MF <- fread("TopGO_result_MF.csv")
CC <- fread("TopGO_result_CC.csv")

p_value_cutoff <- 1e-06

# Number of GO terms that are significantly enriched, per progulon
BP <- melt(BP, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
BP <- BP[, .(BP_terms = sum(pvalue < p_value_cutoff)), by = prn ]
MF <- melt(MF, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
MF <- MF[, .(MF_terms = sum(pvalue < p_value_cutoff)), by = prn ]
CC <- melt(CC, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
CC <- CC[, .(CC_terms = sum(pvalue < p_value_cutoff)), by = prn ]

GO <- merge(BP, MF)
GO <- merge(GO, CC)

#### Re-order progulons as in the correlation figure panel ####

# Load the correlation data
all_prn_combinations <- fread("ProgulonCor.csv")

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

# Add annotated names and reorder progulons 
GO <- merge(GO, prn_annot[, .(prn = Progulon, Function)])
GO <- GO[, .(Function, BP_terms, MF_terms, CC_terms)]
GO <- melt(GO, id.vars = "Function")
GO[, Function := factor(Function, levels = new_prn_order)]
GO[, variable := factor(variable, levels = c("MF_terms", "CC_terms", "BP_terms"))]

#### Create the GO plots ####
p <- ggplot(GO, aes(x = Function, y = value, fill = variable))+
      geom_bar(stat="identity")+
      coord_flip()+
      scale_fill_viridis(discrete = TRUE, option = "D", direction = -1)+
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size=0.25), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6))

p1 <- p + theme( legend.position = "none")
p2 <- p + theme( legend.position = "top")

# Save the plot
ggsave("GO_terms.pdf", p1, width = 6, height = 10.5, units = "cm")
ggsave("GO_terms_legend.pdf", p2, width = 10.5, height = 10.5, units = "cm")



