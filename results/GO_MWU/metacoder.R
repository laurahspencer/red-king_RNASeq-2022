#===========================================================================#
# Script created by Janna Willoughby
# Script modified by Avril Harder
#
# This script generates metacoder trees for a given list of GO terms;
# See code_notes.pdf for additional info. on functionality
#
# Update March 17, 2021: script was previously only run in R v3 and some 
# required packages may not be compatible with R v4.
#===========================================================================#

mc <- "https://cran.r-project.org/src/contrib/Archive/metacoder/metacoder_0.1.3.tar.gz"
install.packages(mc, repos=NULL, type="source") ## version 1.3 is the only version compatible with the code below

library(GO.db)
library(metacoder)
library(org.Hs.eg.db)
library(airway)
library(AnnotationDbi)
library(metacoder)

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
uni.terms <- read.csv(paste0(module,"_unique_GO_terms.csv"))


##### GET LIST OF INTERESTING GO TERMS (filt.go) THAT ARE ALSO UNIQUE TO MODULE #####
n <- 30 ## can limit the # of GO terms you want to include, just to see how it runs
final.go <- head(uni.terms, n=n)
# final.go <- uni.terms[which(uni.terms$ID %in% filt.go),]
# rm(list=ls()[! ls() %in% c("final.go","module")])

##### PASS FINAL GO LIST TO METACODER TO PLOT TREES #####
setwd("~/Desktop/")

bpdata <- as.data.frame(final.go$ID)
colnames(bpdata) = c("GO.ID")
bpdata$GO.ID <- as.character(bpdata$GO.ID)

#function needed to parse GO IDs
term_class = function(x, current=x, all_paths = FALSE, type = GOBPPARENTS, verbose = TRUE, valid_relationships = c("is_a")) {
  # Get immediate children of current taxon
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
bpterms = lapply(bpdata$GO.ID, term_class, all_paths = FALSE, type = GOBPPARENTS) #this line can take forever
bpres   = data.frame(class=unlist(bpterms))

#write/write data (annoyingly the only way I can get it to work, but at least it works)
write.table(bpres, "data/temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
bpres=read.table("data/temp.csv", header=TRUE, sep=",")

## bpdata contains the GO IDs for the terminal nodes
## bpres and bpdata contain the same information, different formats = path information for each terminal node (GO Terms)

#parse GO data
data = parse_taxonomy_table("data/temp.csv",
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

# pdf(paste("output/",module,"meta_labels.pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=FALSE)
#grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
par(bg=NA)
set.seed(100) ## seeds: orange = 35       lightsteelblue1 = 75       paleturquoise = 100
heat_tree(tempdata, node_label = tempdata$taxon_data$name,
          # node_size = colsandsize$ngenes,
          # node_size_trans = "log10",
          node_size_range = c(0.01, 0.01),
          # node_label_size_trans = "log10",
          node_label_size_range = c(0.01, 0.01),
          # edge_size_trans = "log10",
          edge_size_range = c(0.004, 0.004),
          node_color = module,
          # node_color_trans = "linear",
          # node_color_range = diverging_palette(),
          # node_color_interval = c(-4, 4),
          # edge_color_trans = "linear",
          # edge_color_range = diverging_palette(),
          # edge_color_interval =  c(-4, 4),
          node_label_max = 500,
          # node_color_axis_label = "Factor change",
          # node_size_axis_label = "Number of genes",
          layout = "da", initial_layout = "re"
)
# dev.off()



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

