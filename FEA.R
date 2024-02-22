###Functional enrichment analysis
#1. Get inputs
#2. check organism (default is human)
#3.convert gene names from input2 to species specific orthologs (in many to one, use first)
#4. make gene_list
#5 run gseGO
#6. overlap input2 orthologs and pull common ones from input1
#7. run sgeGo on these


library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(biomaRt)
library(enrichplot)
library(ggupset)

#-------------------------------------------------------------------------------#
##################################### INPUTS ####################################
#-------------------------------------------------------------------------------#

#take user args for degenes and geneset files
# Get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Check if at least 3 arguments are provided
if(length(args) < 3) {
  cat("Error: At least 3 arguments are required.\n")
  quit(status = 1)
}

# Assign arguments to variables, with default values for optional arguments
input1 <- args[1]

#next iteration*
#input2 <- args[2]

#Human and mouse for now
species <- args[2] 

outdir <- args[3]

#next iteration*
#currently produce both.
#filter <- ifelse(length(args) >= 5, args[5], "Default selected. Inputs wont be filtered.")

#List of supported species
species_list <- c("mouse", "human")

###############################################################################
#############Functions below.##################################################
###############################################################################

#add handling of other species when > mouse/human.
is_third_arg_in_list <- function(args, my_list) {
  if (length(args) < 3) {
    stop("Not enough arguments provided")
  }
  
  third_arg <- args[2]
  return(third_arg %in% my_list)
}
###############################################################################
#Function to run gene enrichment analysis. 
#Inputs:
#1. genelist
#Output: output of gseGO
###############################################################################

gseGO_analysis <-function(gene_list, species) {
  if(species == 'mouse'){
    db=org.Mm.eg.db
  } else {
    db=org.Hs.eg.db
  }
  #gseout=paste0(outdir,"/", "gseGO_output.txt")
  gse_out <- gseGO(geneList=gene_list ,
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 1, 
                   verbose = TRUE, 
                   OrgDb = db, 
                   eps = 0,
                   pAdjustMethod = "none")
  return(gse_out)
}

gseGO_analysis_pvalue <-function(gene_list, species) {
  if(species == 'mouse'){
    db=org.Mm.eg.db
  } else {
    db=org.Hs.eg.db
  }
  #gseout=paste0(outdir,"/", "gseGO_output.txt")
  gse_out <- gseGO(geneList=gene_list ,
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = db,
                   eps = 0,
                   pAdjustMethod = "none")
  return(gse_out)
}

###############################################################################
#Function to run functional enrichment analysis. 
#Inputs:
#1. genelist
#Output: output of enrichGO
###############################################################################

enrichGO_analysis <-function(gene_list,universe,species) {
  ###Output1
  #gseout=paste0(outdir,"/", "gseGO_output.txt")
  if(species == 'mouse'){
    db='org.Mm.eg.db'
  } else {
    db='org.Hs.eg.db'
  }
  #Question: Should I be giving a list of significant genes as first input or all?
  go_enrich <- enrichGO(gene = names(gene_list),
                        universe = names(universe),
                        OrgDb = db, 
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = "ALL",
                        pAdjustMethod = "none",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
  
  return(go_enrich)
}


###############################################################################
#Function to plot functional enrichment analysis. 
#Inputs:
#1. enrichGO output
#2 Output directory
#Output: output of enrichGO
###############################################################################

plots_FEA <-function(go, sname, outdir) {
  #Upset plot
  upset<- upsetplot(go)
  outn<-paste0(outdir,"/",sname,"_upset_plot.pdf")
  pdf(file=outn) # or other device
  print(upset)
  dev.off()
  #Bar plot
  bar<-barplot(go, 
               drop = TRUE, 
               showCategory = 10, 
               title = "GO Biological Pathways",
               font.size = 8)
  outb<-paste0(outdir,"/",sname,"_bar_plot.pdf")
  pdf(file=outb)
  print(bar)
  dev.off()
  #Dot plot
  #outd<-paste0(outdir,"/",sname,"_dot_plot.pdf")
  #dot<-dotplot(go, showCategory=20, split=".sign",font.size = 8) + facet_wrap(.~.sign)
  #pdf(file=outd)
  #print(dot)
  #dev.off()
  
}

###############################################################################
#Function to plot functional enrichment analysis. 
#Inputs:
#1. enrichGO output
#2 Output directory
#Output: output of enrichGO
###############################################################################

plots_GSEA <-function(go, sname, outdir) {
  #Dot plot
  outd<-paste0(outdir,"/",sname,"_dot_plot.pdf")
  dot<-dotplot(go, showCategory=20, split=".sign",font.size = 8) + facet_wrap(.~.sign)
  pdf(file=outd)
  print(dot)
  dev.off()
  #ridgeplot(gse) + labs(x = "enrichment distribution")
  #https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
}
###############################################################################
###############################################################################
###############################################################################

#check if passed species is in list supoorted. Useful for later when more species.
speciesL <- is_third_arg_in_list(args, species_list)
if(speciesL == TRUE){
  species = args[2]
}

print(species)
#temp--> only for local testing
#input1 = "/Users/apaala.chatterjee/ptran/PTran/functional_enrichment_input/DE_unfilter_PGRT/RPGT_vs_Control.de_genes.txt"
#input2 = "/Users/apaala.chatterjee/ptran/PTran/functional_enrichment_input/SCLC_SUBTYPE_MARKER.txt"
#outdir = "/Users/apaala.chatterjee/ptran/testout"
#species = "mouse"

#-------------------------------------------------------------------------------#
################################## Start Processing #############################
#-------------------------------------------------------------------------------#

#Read DESeq2 results (unfiltered)
input1f = read.table(input1, sep = '\t', header = TRUE)

# Define Universe--------------------------------------------------

# we want the log2 fold change 
original_gene_list <- unlist(input1f[5])

# name the vector
names(original_gene_list) <- unlist(input1f[1])

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
universe_gene_list = sort(gene_list, decreasing = TRUE)
################################################################
# From significant results, we want to filter on log2fold change- AMOL Q
#
# Extract significant results (padj < 0.05) *later?

# filter on min log2fold change (log2FoldChange <-1 / >1)
sig_genes_df = subset(input1f, input1f[5] > 1 | input1f[5]< -1)

#Set vector values
sig_genes <- unlist(sig_genes_df[5])

# Name the vector with feature IDs
names(sig_genes) <- unlist(sig_genes_df[1])

# omit NA values
sig_genes <- na.omit(sig_genes)

# I need to sort it
sig_genes_list = sort(sig_genes, decreasing = TRUE)

#-------------------------------------------------------------------------------#
################################## Start FEA ####################################
#-------------------------------------------------------------------------------#

#For universe
universe_enrich <- enrichGO_analysis(universe_gene_list, universe_gene_list, species)
print(head(universe_enrich))
#write to file
out1=paste0(outdir,'/FEA_universe_results.txt')
write.table(universe_enrich, file = out1, sep = "\t")

#Make plots
plots_FEA(universe_enrich,"universe", outdir)

 #For significant ones
sig_enrich <- enrichGO_analysis(sig_genes_list, universe_gene_list, species)

#write to file
out1=paste0(outdir,'/FEA_significant_results.txt')
write.table(sig_enrich, file = out1, sep = "\t")

#Make plots
plots_FEA(sig_enrich,"significant", outdir)

#-------------------------------------------------------------------------------#
################################## Start GSE ####################################
#-------------------------------------------------------------------------------#

#Run gseGO universe
go_enrich <- gseGO_analysis(universe_gene_list,species)

#write to file
out1=paste0(outdir,'/GSEA_universe_results.txt')
write.table(go_enrich, file = out1, sep = "\t")

#Make plots
plots_GSEA(go_enrich,"universe_gsea", outdir)

#Run gseGO universe
go_enrich_sig <- gseGO_analysis(sig_genes_list,species)

#write to file
out1=paste0(outdir,'/GSE_significant_results.txt')
write.table(go_enrich_sig, file = out1, sep = "\t")

#Make plots
plots_GSEA(go_enrich_sig,"significant_gsea", outdir)


#Run gseGO pvalue
go_enrich_sig_pval <- gseGO_analysis_pvalue(sig_genes_list,species)

#write to file
out1=paste0(outdir,'/GSE_significant_results_pvalue.txt')
write.table(go_enrich_sig_pval, file = out1, sep = "\t")

#Make plots
plots_GSEA(go_enrich_sig_pval,"significant_gsea_pvalue", outdir)


