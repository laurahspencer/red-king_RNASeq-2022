
input="WGCNA-module_lightgreen.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotation="WGCNA-genes_for-GOMWU.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotation, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=.99,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=3,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.15, #, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE
           #Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# ----------- Plotting results

#quartz()
results=gomwuPlot(input,goAnnotation,goDivision,
                  absValue=0.001,  #-log(0.05,10) # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
print(results[[1]])


# ------- extracting representative GOs
# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
print(bestGOs)

(read_delim("WGCNA-genes_for-GOMWU.tab", delim = "\t") %>%  filter(grepl(go.interest, gene_ontology_i_ds))) %>% View()
# Looking for my GO term and associated genes of interest in main data file 

# Double check that my GO term (and associated genes) is in the background gene list and input data 
(go.interest <- "GO:0006313") # GO term of interest, transposition, DNA-mediated, which in the Biological Process category, http://www.informatics.jax.org/vocab/gene_ontology/GO:0006313
(go.genes <- (read_delim("WGCNA-genes_for-GOMWU.tab", delim = "\t") %>%  filter(grepl(go.interest, gene_ontology_i_ds)))$id) #these genes in background gene list are associated with GO term of interest 
((go.genes.sign <- read_delim("WGCNA-module_lightgreen.csv", delim = ",") %>% filter(sign!=0) %>% filter(id %in% go.genes))) #these genes in input file are associated with GO of interest

read_delim("WGCNA-module_lightgreen.csv", delim = ",") %>% filter(sign!=0) %>% View()

#read_delim("WGCNA-genes_for-GOMWU.tab", delim = "\t") %>% select(gene_ontology_i_ds) %>% unlist()

(read_delim("WGCNA-genes_for-GOMWU.tab", delim = "\t") %>%  filter(grepl(go.interest, gene_ontology_i_ds))) %>% nrow()
(read_delim("WGCNA-genes_for-GOMWU.tab", delim = "\t")) %>% nrow()

nrow(read_delim("WGCNA-module_lightgreen.csv", delim = ",") %>% filter(sign!=0) %>% filter(id %in% go.genes))
nrow(read_delim("WGCNA-module_lightgreen.csv", delim = ",") %>% filter(sign!=0))


# Are all genes associated with GO term of interest in both my inputs, background and gene set? 
all(go.genes.sign %in% go.genes)

# Read in the "main data table" that contains analyzed GO terms and filter for GO ID and genes of interest
read_delim("BP_WGCNA-module_lightgreen.csv", delim = "\t") %>% filter(grepl(go.interest, term)) # GO term of interest is not in the "main data file" - why not? 
read_delim("BP_WGCNA-module_lightgreen.csv", delim = "\t") %>% filter(seq %in% go.genes.sign$id) # genes associated with that go term are also not in the "main data file" - why are they gone? 

# Are there any GO terms that contain the word "transposition" in the output? All offspring of my GO term of interest contain "transposition"
(read_delim("BP_WGCNA-module_lightgreen.csv", delim = "\t") %>% filter(grepl("transpos", name))) #yes, four genes associated 
(read_delim("BP_WGCNA-module_lightgreen.csv", delim = "\t") %>% filter(grepl("transpos", name)))$seq %in% go.genes.sign$id # but none of the genes assigned to my GO term of interest are those associated with the offspring!

# so weird!