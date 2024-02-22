

###Functional enrichment analysis
#1. Get inputs
#2. check organism (default is human)
#3.convert gene names from input2 to species specific orthologs (in many to one, use first)
#4. make gene_list
#5 run gseGO
#6. overlap input2 orthologs and pull common ones from input1
#7. run sgeGo on these


Usage
Rscript FEA.R path/to/deseq/out species path/to/outdir
