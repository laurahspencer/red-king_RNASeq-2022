#===========================================================================#
# Script created by Janna Willoughby
# Script modified by Avril Harder
# Script further modified by Laura Spencer 
#
# This script generates metacoder trees for a given list of GO terms;
# See code_notes.pdf for additional info. on functionality
#
# Update March 17, 2021: script was previously only run in R v3 and some 
# required packages may not be compatible with R v4.

# Laura's updates:
# Run on RStudio Version 3.6.3 using RStudio Cloud
#===========================================================================#

list.of.packages <- c("tidyveres", "stringr", "ggplot2", "igraph", "scales", "taxize", "seqinr", "reshape2", "zoo", "traits", "RColorBrewer", "RCurl", "ape", "reshape", "lazyeval", "dplyr", "magrittr", "readr", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
# Install all new packages
if(length(new.packages)) install.packages(new.packages)
# Load all required libraries
all.packages <- c(list.of.packages, bioconductor.packages)
lapply(all.packages, FUN = function(X) {
  do.call("require", list(X))
})

phangorn <- "https://cran.r-project.org/src/contrib/Archive/phangorn/phangorn_2.5.5.tar.gz"
install.packages(phangorn, repos=NULL, type="source") ## version 1.3 is the only version compatible with the code below
install.packages("taxize")

mc <- "https://cran.r-project.org/src/contrib/Archive/metacoder/metacoder_0.1.3.tar.gz"
install.packages(mc, repos=NULL, type="source") ## version 1.3 is the only version compatible with the code below

# Install BiocManager if needed
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioconductor.packages <- c("GO.db", "org.Hs.eg.db", "airway", "AnnotationDbi")
new.bioc.packages <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[, "Package"])]
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

library(GO.db)
library(org.Hs.eg.db)
library(airway)
library(AnnotationDbi)
library(metacoder)

## Custom function to create transparent colors
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
## END



# set module of interest
module <- "orange" # lightsteelblue1  orange  paleturquoise

setwd("/Users/Avril/Desktop/")

##### READ IN ALL MODULE GENE:GO INFORMATION #####
## get list of all GO IDs and corresponding terms
terms <- Term(GOTERM)
terms <- as.data.frame(terms)
terms <- cbind(rownames(terms), terms)
rownames(terms) <- c()
colnames(terms) <- c("go.id","term")

# ## get list of gene:GO term relationships
# ##
# ## gene.name       go.id
# ##
# rep.mods <- read.csv("/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/wgcna/go_for_janna/rep4_module_gene_go_info.csv")
# ## keep info for module of interest
# rep.mods <- rep.mods[which(rep.mods$module==module),]
# rep.mods$go.id <- as.character(rep.mods$go.id)
# rep.mods <- rep.mods[,c(2,4)]


##### GET GO TERMS UNIQUE TO MODULE OF INTEREST #####
##
## ID       term       count
##
uni.terms <- read.csv(paste0(module,"_unique_GO_terms.txt"))
load(file = "BP.metacoder") 
load(file = "MF.metacoder")
load(file = "CC.metacoder")



# Ambient vs. Moderate
uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "amb-mod") %>% filter(delta.rank>0) %>% arrange(pval) #UPREGULATED
#uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "amb-mod") %>% filter(delta.rank<0) %>% arrange(pval) #DOWNREGULATED

#Ambient vs. Low
uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "amb-low") %>% filter(delta.rank>0) %>% arrange(pval) #UPREGULATED IN OA
#uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "amb-low") %>% filter(delta.rank<0) %>% arrange(pval) #DOWNREGULATED IN OA

# Don't specify upreg/downreg here, do so in the network 
uni.terms <- BP.metacoder %>% filter(p.adj<0.01) %>% filter(gene_set == "amb-mod") %>% arrange(pval) #UPREGULATED IN OA
uni.terms <- BP.metacoder %>% filter(p.adj<0.01) %>% filter(gene_set == "amb-low") %>% arrange(pval) #UPREGULATED IN OA

# Don't specify upreg/downreg here, do so in the network 
uni.terms <- MF.metacoder %>% filter(p.adj<0.01) %>% filter(gene_set == "amb-mod") %>% arrange(pval) #UPREGULATED IN OA
uni.terms <- MF.metacoder %>% filter(p.adj<0.01) %>% filter(gene_set == "amb-low") %>% arrange(pval) #UPREGULATED IN OA

# #Moderate vs. Low
# uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "mod-low") %>% filter(delta.rank>0) %>% arrange(pval) #UPREGULATED IN OA
# uni.terms <- BP.metacoder %>% filter(p.adj<0.05) %>% filter(gene_set == "mod-low") %>% filter(delta.rank<0) %>% arrange(pval) #DOWNREGULATED IN OA



network.DEGs <- function(contrast, GO.results, GO.type, GO.parents, alpha) {
  
  # for testing purposes
  modules = c("amb-mod")
  GO.results = BP.metacoder
  GO.type = GOBPPARENTS
  alpha = 0.05
  
  # Define input dataframe from the desired GO type (BP, MF, or CC)
  uni.terms <- GO.results %>% filter(p.adj<alpha) %>% filter(gene_set %in% contrast) %>% arrange(pval)
  
  ##### GET LIST OF INTERESTING GO TERMS (filt.go) THAT ARE ALSO UNIQUE TO MODULE #####
  # n <- 30 ## can limit the # of GO terms you want to include, just to see how it runs
  n <- nrow(uni.terms)
  final.go <- head(uni.terms, n=n)
  # final.go <- uni.terms[which(uni.terms$ID %in% filt.go),]
  # rm(list=ls()[! ls() %in% c("final.go","module")])
  
  ##### PASS FINAL GO LIST TO METACODER TO PLOT TREES #####
  #setwd("~/Desktop/")
  
  bpdata <- as.data.frame(final.go$ID)
  colnames(bpdata) = c("GO.ID")
  bpdata$GO.ID <- as.character(bpdata$GO.ID)
  
  #function needed to parse GO IDs
  term_class = function(x, current=x, all_paths = FALSE, type = GO.type, verbose = TRUE, valid_relationships = c("is_a")) {
    # Get immediate children of current taxon  #type = GOMFPARENTS, 
    parents = tryCatch({
      possible_parents <- as.list(type[x[1]])[[1]] #this line doesn't function?
      if (! is.null(valid_relationships)) {
        possible_parents <- possible_parents[names(possible_parents) %in% valid_relationships]
      }
      names(AnnotationDbi::Term(possible_parents))
    }, error = function(e) {
      c()
    })
    
    # only go down one path if desired
    if (! all_paths) {
      parents <- parents[1]
    }
    parents <- parents[parents != "all"]
    
    if (is.null(parents)) {
      return(c())
    } else if (length(parents) == 0) {
      cat(length(x))
      return(paste0(collapse = "|", AnnotationDbi::Term(x)))
    } else {
      next_x <- lapply(parents, function(y) c(y, x))
      
      # Run this function on them to get their output
      child_output <- lapply(next_x, term_class, all_paths = all_paths, type = type)
      output <- unlist(child_output, recursive = FALSE)
      
      return(output)
    }
  }
  
  #modify, pull go term relationships
  bpterms = lapply(bpdata$GO.ID, term_class, all_paths = FALSE, type = GO.type) #this line can take forever
  bpres   = data.frame(class=unlist(bpterms))
  
  #write/write data (annoyingly the only way I can get it to work, but at least it works)
  write.table(bpres, "temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
  bpres=read.table("temp.csv", header=TRUE, sep=",")
  
  ## bpdata contains the GO IDs for the terminal nodes
  ## bpres and bpdata contain the same information, different formats = path information for each terminal node (GO Terms)
  
  #parse GO data
  data = parse_taxonomy_table("temp.csv",
                              taxon_col = c("class" = -1),
                              other_col_type = "obs_info",
                              sep=",",
                              class_sep = "\\|")
  
  # parse GO data (parse_taxonomy_table replaced by taxa::parse_tax_data)
  # data <- parse_tax_data(bpres,
  #                        class_sep = "\\|")
  
  
  #create figure
  tempdata = filter_taxa(data, n_supertaxa <= 500) #filters terms that are WAAAAY out from the middle
  ## n_supertaxa sets # of nodes that can be in a single path from center to terminal nodes (cuts from terminal end, not internal nodes)
  
  # Laura's modifications - add some data
  tempdata$taxon_data <- tempdata$taxon_data %>% mutate(name=str_remove(name, '"')) %>%
    left_join(uni.terms %>% select(term, nseqs, pval, p.adj, delta.rank), by=c("name"="term")) %>%
    mutate(nseqs = replace_na(nseqs, 1)) %>% #for GO terms that were filled in by metacoder, replace nseqs with 1 for size mapping
    
    # Create a new column indicating upregulation=1 or downregulation=0 for coloring nodes
    mutate(change=case_when(
      delta.rank > 0 ~ 1,
      delta.rank < 0 ~ 0)) %>% 
    
    # Create a new column that retains only GO names for enriched GO terms AND the central parent nodes
    mutate(name2=p.adj) %>% 
    mutate(name2 = ifelse(supertaxon_ids == 1, name,
                          ifelse(is.na(p.adj), NA, name))) 
  
  # remotes::install_github("mjdufort/miscHelpers")
  # require(miscHelpers)
  
  # Create a vector of colors designating upregulated (green) or downregulated (red) GO categories 
  node.colors <- tempdata$taxon_data$change %>% 
    values2colors(col.start = "#CC3333", col.end = "#009933", breaks = 2) %>%
    adjustcolor(alpha.f = 0.65) #set transparency
  
  # pdf(paste("output/",module,"meta_labels.pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=FALSE)
  #grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
  
  par(bg=NA)
  set.seed(100) ## seeds: orange = 35       lightsteelblue1 = 75       paleturquoise = 100
  heat_tree(tempdata, node_label = tempdata$taxon_data$name2,
            node_size = tempdata$taxon_data$nseqs,
            # node_size_trans = "log10",
            node_size_range = c(0.01, 0.04),
            # node_label_size_trans = "log10",
            node_label_size_range = c(0.009, 0.009),
            # edge_size_trans = "log10",
            edge_size_range = c(0.004, 0.004),
            #node_color = tempdata$taxon_data$p.adj %>% values2colors(col.start = "red", col.end = "blue", breaks = 15),
            
            # This node color indicates upregulation vs. downregulation 
            node_color = node.colors,
            # node_color_trans = "linear",
            #node_color_range = quantative_palette(),
            #node_color_interval = c(-4, 4),
            # edge_color_trans = "linear",
            # edge_color_range = diverging_palette(),
            # edge_color_interval =  c(-4, 4),
            node_label_max = 500,
            # node_color_axis_label = "Factor change",
            node_size_axis_label = "Number of genes",
            layout = "da", initial_layout = "re",
            overlap_avoidance=2, make_legend=FALSE
  )
  # dev.off()
}

# OPTIONS: 
# contrast: "amb-mod", "amb-low", "mod-low"
# GO.results: BP.metacoder, MF.metacoder, CC.metacoder
# GO.type: "GOBPPARENTS", "GOMFPARENTS", "GOCCPARENTS"
# alpha= any value 0-1 designating level of significance 

network.DEGs(contrast = "amb-mod", GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.05)
network.DEGs(contrast = "amb-low", GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.05)
network.DEGs(contrast = "mod-low", GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.05)

network.DEGs(contrast = "amb-mod", GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.05)
network.DEGs(contrast = "amb-low", GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.05)
network.DEGs(contrast = "mod-low", GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.05)

network.DEGs(contrast = "amb-mod", GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.05)
network.DEGs(contrast = "amb-low", GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.05)
network.DEGs(contrast = "mod-low", GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.05)


### =================== USING WGCNA MODULE RESULS =================== ### 

network.modules <- function(modules, GO.results, GO.type, GO.parents, alpha) {
  
  # # for testing purposes
  # modules = c("pink", "royalblue", "darkviolet", "firebrick4", "lightcyan", "magenta")
  # GO.results = BP.metacoder
  # GO.type = GOBPPARENTS
  # alpha = 0.05
  
  # Define input dataframe from the desired GO type (BP, MF, or CC)
  uni.terms <- GO.results %>% filter(p.adj<alpha) %>% filter(gene_set %in% modules) %>% arrange(pval)
  
  ##### GET LIST OF INTERESTING GO TERMS (filt.go) THAT ARE ALSO UNIQUE TO MODULE #####
  # n <- 30 ## can limit the # of GO terms you want to include, just to see how it runs
  n <- nrow(uni.terms)
  final.go <- head(uni.terms, n=n)
  # final.go <- uni.terms[which(uni.terms$ID %in% filt.go),]
  # rm(list=ls()[! ls() %in% c("final.go","module")])
  
  ##### PASS FINAL GO LIST TO METACODER TO PLOT TREES #####
  #setwd("~/Desktop/")
  
  bpdata <- as.data.frame(final.go$ID)
  colnames(bpdata) = c("GO.ID")
  bpdata$GO.ID <- as.character(bpdata$GO.ID)
  
  #function needed to parse GO IDs
  term_class = function(x, current=x, all_paths = FALSE, type = GO.type, verbose = TRUE, valid_relationships = c("is_a")) {
    # Get immediate children of current taxon  #type = GOMFPARENTS, 
    parents = tryCatch({
      possible_parents <- as.list(type[x[1]])[[1]] #this line doesn't function?
      if (! is.null(valid_relationships)) {
        possible_parents <- possible_parents[names(possible_parents) %in% valid_relationships]
      }
      names(AnnotationDbi::Term(possible_parents))
    }, error = function(e) {
      c()
    })
    
    # only go down one path if desired
    if (! all_paths) {
      parents <- parents[1]
    }
    parents <- parents[parents != "all"]
    
    if (is.null(parents)) {
      return(c())
    } else if (length(parents) == 0) {
      cat(length(x))
      return(paste0(collapse = "|", AnnotationDbi::Term(x)))
    } else {
      next_x <- lapply(parents, function(y) c(y, x))
      
      # Run this function on them to get their output
      child_output <- lapply(next_x, term_class, all_paths = all_paths, type = type)
      output <- unlist(child_output, recursive = FALSE)
      
      return(output)
    }
  }
  
  #modify, pull go term relationships
  bpterms = lapply(bpdata$GO.ID, term_class, all_paths = FALSE, type = GO.type) #this line can take forever
  bpres   = data.frame(class=unlist(bpterms))
  
  #write/write data (annoyingly the only way I can get it to work, but at least it works)
  write.table(bpres, "temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
  bpres=read.table("temp.csv", header=TRUE, sep=",")
  
  ## bpdata contains the GO IDs for the terminal nodes
  ## bpres and bpdata contain the same information, different formats = path information for each terminal node (GO Terms)
  
  #parse GO data
  data = parse_taxonomy_table("temp.csv",
                              taxon_col = c("class" = -1),
                              other_col_type = "obs_info",
                              sep=",",
                              class_sep = "\\|")
  
  # parse GO data (parse_taxonomy_table replaced by taxa::parse_tax_data)
  # data <- parse_tax_data(bpres,
  #                        class_sep = "\\|")
  
  
  #create figure
  tempdata = filter_taxa(data, n_supertaxa <= 500) #filters terms that are WAAAAY out from the middle
  ## n_supertaxa sets # of nodes that can be in a single path from center to terminal nodes (cuts from terminal end, not internal nodes)
  
  # Laura's modifications - add some data
  tempdata$taxon_data <- tempdata$taxon_data %>% mutate(name=str_remove(name, '"')) %>%
    left_join(uni.terms %>% select(gene_set_up, term, nseqs, pval, p.adj), by=c("name"="term")) %>%
    
    # Multiple modules have the same enriched GO terms which can't be mapped, so here I break those ties using teh lowest p-adj
    group_by(taxon_ids, supertaxon_ids, name) %>%    
    summarise(p.adj=min(p.adj)) %>% distinct() %>%  #for rows with duplicate grouping variabes, select one with lowest p-value 
    left_join(uni.terms %>% select(gene_set_up, term, nseqs, p.adj), by=c("name"="term", "p.adj")) %>% #re-add data 
    ungroup() %>% group_by(across(c(-gene_set_up))) %>% sample_n(1) %>%  # Still some duplicate terms (since there were dup. p-values), select one at random
    
    mutate(nseqs = replace_na(nseqs, 1)) %>% #for GO terms that were filled in by metacoder, replace nseqs with 1 for size mapping
    mutate(gene_set_up  = replace_na(gene_set_up , "gray")) %>%
    mutate(gene_set_up = str_replace(gene_set_up, "ivory", "wheat")) %>% #replace colors that are too light to see
    
    # Create a new column that retains only GO names for enriched GO terms AND the central parent nodes
    mutate(name2=p.adj) %>% 
    mutate(name2 = ifelse(supertaxon_ids == 1, name,
                          ifelse(is.na(p.adj), NA, name))) 
  
  # Create a vector of colors designating upregulated (green) or downregulated (red) GO categories 
  node.colors <- tempdata$taxon_data$gene_set_up %>% 
    adjustcolor(alpha.f = 0.65) #set transparency
  
  
  
  # pdf(paste("output/",module,"meta_labels.pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=FALSE)
  #grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
  
  par(bg=NA)
  set.seed(100) ## seeds: orange = 35       lightsteelblue1 = 75       paleturquoise = 100
  heat_tree(tempdata, node_label = tempdata$taxon_data$name2,
            node_size = tempdata$taxon_data$nseqs, # Use number of genes to set node size
            #node_size_trans = "log10",
            node_size_range = c(0.01, 0.04),
            # node_label_size_trans = "log10",
            node_label_size_range = c(0.009, 0.009),
            # edge_size_trans = "log10",
            edge_size_range = c(0.004, 0.004),
            #node_color = tempdata$taxon_data$p.adj %>% values2colors(col.start = "red", col.end = "blue", breaks = 15),
            
            # This node color indicates module membership
            node_color = node.colors,
            # node_color_trans = "linear",
            #node_color_range = quantative_palette(),
            #node_color_interval = c(-4, 4),
            # edge_color_trans = "linear",
            # edge_color_range = diverging_palette(),
            # edge_color_interval =  c(-4, 4),
            node_label_max = 500,
            # node_color_axis_label = "Factor change",
            node_size_axis_label = "Number of genes",
            layout = "da", initial_layout = "re",
            overlap_avoidance=2, make_legend=FALSE
  )
  # dev.off()
}

### ===== CONSTRUCT NETWORKS OF ENRICHED GO TERMS FROM PCO2-ASSOCIATED MODULES === #### 

## BIOLOGICAL PROCESSES
# Genes that decrease with pCO2 - DOWNREGULATED 
network.modules(modules = c("pink", "royalblue", "darkviolet", "firebrick4", "lightcyan", "magenta"), 
                GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.05)

# Genes that decrease with pCO2 - UPREGULATED 
network.modules(modules = c("blue2", "lightsteelblue1", "ivory", "lightgreen", "purple", "green"), 
                GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.1)

# Genes that go up in moderate OA then down in severe OA 
network.modules(modules = c("coral1", "plum1"), 
                GO.results = BP.metacoder, GO.type = GOBPPARENTS, alpha = 0.1)



## MOLECULAR FUNCTIONS
# Genes that decrease with pCO2 - DOWNREGULATED 
network.modules(modules = c("pink", "royalblue", "darkviolet", "firebrick4", "lightcyan", "magenta"), 
                GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.01)

# Genes that decrease with pCO2 - UPREGULATED 
network.modules(modules = c("blue2", "lightsteelblue1", "ivory", "lightgreen", "purple", "green"), 
                GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.1)

# Genes that go up in moderate OA then down in severe OA 
network.modules(modules = c("coral1", "plum1"), 
                GO.results = MF.metacoder, GO.type = GOMFPARENTS, alpha = 0.1)

## CELULAR COMPONENTS
# Genes that decrease with pCO2 - DOWNREGULATED 
network.modules(modules = c("pink", "royalblue", "darkviolet", "firebrick4", "lightcyan", "magenta"), 
                GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.05)

# Genes that decrease with pCO2 - UPREGULATED 
network.modules(modules = c("blue2", "lightsteelblue1", "ivory", "lightgreen", "purple", "green"), 
                GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.1)

# Genes that go up in moderate OA then down in severe OA 
network.modules(modules = c("coral1", "plum1"), 
                GO.results = CC.metacoder, GO.type = GOCCPARENTS, alpha = 0.1)


## run this version to find best seed #s for each module's tree - have to highlight and run manually :/
# par(bg=NA)
# sede <- sample(1:100,1)
# pdf(paste0("/Users/Avril/Desktop/",module,"/",sede,".pdf"), width=5, height=5, useDingbats=FALSE, onefile=F)
# set.seed(sede) #grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
# heat_tree(tempdata, #node_label = tempdata$taxon_data$name,
#           # node_size = colsandsize$ngenes,
#           # node_size_trans = "log10",
#           node_size_range = c(0.01, 0.01),
#           # node_label_size_trans = "log10",
#           node_label_size_range = c(0.01, 0.01),
#           # edge_size_trans = "log10",
#           edge_size_range = c(0.004, 0.004),
#           node_color = module,
#           # node_color_trans = "linear",
#           # node_color_range = diverging_palette(),
#           # node_color_interval = c(-4, 4),
#           # edge_color_trans = "linear",
#           # edge_color_range = diverging_palette(),
#           # edge_color_interval =  c(-4, 4),
#           node_label_max = 500,
#           # node_color_axis_label = "Factor change",
#           # node_size_axis_label = "Number of genes",
#           layout = "da", initial_layout = "re"
# )
# dev.off()

sessionInfo()
