#FILE :    generate_network.R
#AUTHOR :  Karsten Suhre
#DATE :   10 January 2021
#         12 Jan 2021 - update IDs to be more informative
#         14 Jan 2021 - add STAT phenotype assocs
#         15 Jan 2021 - add eQTL meQTL, eQTMs
#         21 Jan 2021 - add GWAS catalogue data
#         29 July 2022 - comment out writing of nodes and edges, as this conflicts with writing rights on the server
#                        adapt for running in a shiny server:
#                        https://www.r-bloggers.com/2021/06/running-shiny-server-in-docker/
#                        add passing a query string to the server
#         1 August 2022 - add HTML code to the about page 
#                       - add an examples page
#         4 August 2022 - fixed a weird bug when clicking again on a selected table entry
#                       - added filter option to tables
#                       - fixed a problem that prevented the display of GWAS traits that contain single quotes, like Crohn's
#
#
#PURPOSE: reformat info from Supplementary Tables to generate a network
#MODIF:
#
# BUGS:    need to update the trait annotations
#
# https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html
# # minimal example
# nodes <- data.frame(id = 1:3)
# edges <- data.frame(from = c(1,2), to = c(1,3))
# visNetwork(nodes, edges, width = "100%")

# books
# https://www.jessesadler.com/post/network-analysis-with-r/

# https://www.rdocumentation.org/packages/tidygraph/versions/1.2.0

# https://datastorm-open.github.io/visNetwork/
# https://datastorm-open.github.io/visNetwork/image_icon.html

# https://rstudio.github.io/shinydashboard/structure.html

# https://www.urlencoder.io/ to encode special character in the use cases

rm(list=ls())
library("readxl")
library("tidyverse")
library("visNetwork")
library("shiny")
library("shinydashboard")
library("DT")

infile = "DNA_RNA_METH_supptables_v5.xlsx"
cat("reading data from", infile, "\n")
GWAS = read_excel(infile, sheet = "GWAS")
RWAS = read_excel(infile, sheet = "RWAS")
EWAS = read_excel(infile, sheet = "EWAS")
MBH = read_excel(infile, sheet = "MBH")
GGM = read_excel(infile, sheet = "GGM")
STAT = read_excel(infile, sheet = "STAT")
GENO = read_excel(infile, sheet = "GENO")
CATA = read_excel(infile, sheet = "CATA")
ANNO = read_excel(infile, sheet = "ANNO")

# fix a problem with quotes in trait names (this blocks using them in URL queries)
cat("removing ' from CATA$TRAIT\n")
CATA$TRAIT = gsub("'", " ", CATA$TRAIT)
ANNO$TRAIT = gsub("'", " ", ANNO$TRAIT)
ANNO$TRAITID = gsub("'", " ", ANNO$TRAITID)
ANNO$SHORTNAME = gsub("'", " ", ANNO$SHORTNAME)

GWAS_edges = data.frame(
  id = GWAS$ID,
  from = GWAS$LOCUS,
  to = GWAS$TRAITID,
  type = rep("GWAS", dim(GWAS)[1]),
  weight = abs(GWAS$BETA),
  sign = sign(GWAS$BETA),
  pvalue = GWAS$PVALUE
)

RWAS_edges = data.frame(
  id = RWAS$ID,
  from = RWAS$`RNA TRANSCRIPT`,
  to = RWAS$TRAITID,
  type = rep("RWAS", dim(RWAS)[1]),
  weight = abs(RWAS$BETA),
  sign = sign(RWAS$BETA),
  pvalue = RWAS$PVALUE
)

EWAS_edges = data.frame(
  id = EWAS$ID,
  from = EWAS$CPGID,
  to = EWAS$TRAITID,
  type = rep("EWAS", dim(EWAS)[1]),
  weight = abs(EWAS$BETA),
  sign = sign(EWAS$BETA),
  pvalue = EWAS$PVALUE
)

MBH_edges = data.frame(
  id = MBH$ID,
  from = MBH$TRAITID1,
  to = MBH$TRAITID2,
  type = rep("MBH", dim(MBH)[1]),
  weight = abs(MBH$`Spearman RHO`),
  sign = sign(MBH$`Spearman RHO`),
  pvalue = MBH$PVALUE
)

GGM_edges = data.frame(
  id = GGM$ID,
  from = GGM$TRAITID1,
  to = GGM$TRAITID2,
  type = rep("GGM", dim(GGM)[1]),
  weight = abs(GGM$PR),
  sign = sign(GGM$PR),
  pvalue = GGM$PVALUE
)

GENO_edges = data.frame(
  id = GENO$ID,
  from = GENO$TRAITID1,
  to = GENO$TRAITID2,
  type = GENO$TYPE,
  weight = abs(GENO$BETA),
  sign = sign(GENO$BETA),
  pvalue = GENO$PVALUE
)

CATA_edges = data.frame(
  id = CATA$ID,
  from = CATA$LOCUS,
  to = CATA$TRAIT,
  type = rep("CATA", dim(CATA)[1]),
  weight = rep(1, dim(CATA)[1]),
  sign = rep(1, dim(CATA)[1]),
  pvalue = 10^(-CATA$LOGPVALUE)
)

all_edges = rbind(GWAS_edges, RWAS_edges, EWAS_edges, MBH_edges, GGM_edges, GENO_edges, CATA_edges)

# get all nodes
all_nodes = data.frame(TRAITID = unique(sort(c(all_edges$from, all_edges$to))))

# get node annotations
all_nodes = left_join(all_nodes, ANNO)

# add phenotype associations

# keep only traits that are already in the model
cat("processing STAT, keeping only traits that are in the network already\n")
TRAITID_ignore = setdiff(STAT$TRAITID, all_nodes$TRAITID)
TRAITID_keep = intersect(all_nodes$TRAITID, STAT$TRAITID)
cat("keeping ", length(TRAITID_keep), " traits, dropping ", length(TRAITID_ignore), "\n")
STATX = left_join(data.frame(TRAITID = TRAITID_keep), STAT)

phenolist = names(STATX) %>% grep("PVALUE_", .) %>% names(STATX)[.] %>% sub("PVALUE_","", .)

# remove the somaPCs
phenolist = phenolist[-grep("SOMA", phenolist)]

cat("these phenotypes will be added to the nodes list:\n")
print(phenolist)
all_nodes %>% tail()

for (pp in phenolist) {
  df.PHENO = data.frame(
    TRAITID = paste("STAT:", pp),
    SORT = 0,
    TRAIT = pp,
    PLAT = "STAT",
    SHORTNAME = paste("STAT:", pp)
  )
  all_nodes = rbind(all_nodes, df.PHENO)
}

pbonf = 0.05/ length(STAT$TRAITID)
cat("adding edges for STAT, using Bonferroni cutoff P=", pbonf, "\n")

for (pp in phenolist) {
  paste("PVALUE_", pp, sep="")
  ix = which(STATX[[paste("PVALUE_", pp, sep="")]] < pbonf)
  cat("Adding", length(ix), "associations with p <", pbonf, "for", pp, "\n")
  STAT_edges = data.frame(
    id = paste("oASSO_", seq(length(ix)), sep=""),
    from = paste("STAT:", pp),
    to = STATX$TRAITID[ix],
    type = rep("STAT", length(ix)),
    weight = abs(STATX[[paste("BETA_", pp, sep="")]][ix]),
    sign = sign(STATX[[paste("BETA_", pp, sep="")]][ix]),
    pvalue = STATX[[paste("PVALUE_", pp, sep="")]][ix]
  )
  all_edges = rbind(all_edges, STAT_edges)
}

# replace all_edges$to and $from ids with all_nodes$SHORTNAME
TRAITID2SHORTNAME = all_nodes$SHORTNAME
names(TRAITID2SHORTNAME) = all_nodes$TRAITID

SHORTNAME2TRAITID = all_nodes$TRAITID
names(SHORTNAME2TRAITID) = all_nodes$SHORTNAME

all_edges$from = TRAITID2SHORTNAME[all_edges$from]
all_edges$to = TRAITID2SHORTNAME[all_edges$to]
all_nodes$id = all_nodes$SHORTNAME
names(all_nodes)[which(names(all_nodes) == "PLAT")] = "plat"

# round the p-values and weights
all_edges$pvalue = signif(all_edges$pvalue, digits = 2)
all_edges$weight = signif(all_edges$weight, digits = 2)

# write nodes and edges to file
cat("there are", dim(all_edges)[1], "edges\n")
# write.table(file = "MultiomicsEdges.tsv", all_edges, col.names = T, row.names = F, sep = "\t")
cat("there are", dim(all_nodes)[1], "nodes\n")
# write.table(file = "MultiomicsNodes.tsv", all_nodes, col.names = T, row.names = F, sep = "\t")

#########################################################
# a data structure for a network
#########################################################

fullnet = list(
  edges = all_edges,
  nodes = all_nodes
)

##################################################################
# function: neighbors - extract all neighbors of a given node list
# input: a node list
# output: a node list
##################################################################
neighbors = function(nodes, network) {
  # test input:
  # network = fullnet  
  # nodes = c("cg19693031", "3485-28_2")
  ix = lapply(nodes, function(x){union(which(network$edges$from == x), which(network$edges$to == x))}) %>% unlist() %>% unique()
  d = union(network$edges$from[ix], network$edges$to[ix])
  d
}

##################################################################
# function: maxneighbors - extract all nodes connected to a given node list
# input: a node list
# output: a node list
# limit: stop if more than limit nodes were found
##################################################################
maxneighbors = function(nodes, network, limit = 0) {
  # test input:
  # network = fullnet  
  # nodes = c("cg19693031", "3485-28_2")
  if (limit == 0) {limit = 1E99}
  nnodes = length(nodes)
  nnodeslast = 0
  nodeslast = nodes
  while ((nnodes > nnodeslast) & (nnodes<=limit)) {
    nodeslast = nodes
    nodes = neighbors(nodes, network)
    nnodeslast = nnodes
    nnodes = length(nodes)
  }
  nodeslast
}

##################################################################
# function: maxneighbors_noSTAT - extract all nodes connected to a given node list
#           but stop growing STAT: nodes
# input: a node list
# output: a node list
# limit: stop if more than limit nodes were found
##################################################################
maxneighbors_noSTAT = function(nodes, network, limit = 0) {
  # test input:
  # network = fullnet  
  # nodes = c("cg19693031", "3485-28_2")
  if (limit == 0) {limit = 1E99}
  nodesin = nodes
  nnodes = length(nodes)
  nnodeslast = 0
  nodeslast = nodes
  while ((nnodes > nnodeslast) & (nnodes<=limit)) {
    nodeslast = nodes
    nnodeslast = nnodes
    nogrow_nodes = nodes[grep("STAT: |GWAS: ", nodes)]
    grow_nodes   = nodes[grep("STAT: |GWAS: ", nodes, invert = TRUE)]
    nodes = c(nogrow_nodes, neighbors(grow_nodes, network))
    nnodes = length(nodes)
  }
  unique(c(nodesin, nodeslast))
}

##################################################################
# function: nodes2network - connect a node list (extract all nodes between them)
# input: a node list
# output: a network
##################################################################
nodes2network = function(nodes, network) {
  
  # test input:
  # network = fullnet  
  # nodes = maxneighbors("cg19693031", fullnet, limit = 100)
  
  ixfrom = lapply(nodes, function(x){which(network$edges$from == x)}) %>% unlist() %>% unique()
  ixto = lapply(nodes, function(x){which(network$edges$to == x)}) %>% unlist() %>% unique()
  ix = intersect(ixfrom, ixto)
  
  iy = lapply(nodes, function(x){which(network$nodes$id == x)}) %>% unlist() %>% unique()
  list ( edges = network$edges[ix,], nodes = network$nodes[iy,])
}

##################################################################
# function: network2node - extract all nodes from a network
# input: a network
# output: a node list
##################################################################
network2nodes = function(network) {
  
  # test input:
  # network = maxneighbors(c("1","2"), fullnet) %>%  nodes2network(fullnet)
  
  nodes = union(network$edges$from, network$edges$to) %>% unlist() %>% unique()
  
  iy = lapply(nodes, function(x){which(network$nodes$id == x)}) %>% unlist() %>% unique()
  
  network$nodes$id[iy]
}

##################################################################
# function: summarize_network - print info about a network
# input: a network
##################################################################
summary_network = function(network) {
  
  cat("The following variables are defined for the edges:\n")
  names(network$edges) %>% print()
  
  cat("The following variables are defined for the nodes:\n")
  names(network$nodes) %>% print()
  
  for (i in names(network$edges)) {
    cat(">>> network$edges$",i, sep="", "\n")
    head(network$edges[[i]]) %>% print()
    print("...")
    tail(network$edges[[i]]) %>% print()
  }
  
  for (i in names(network$nodes)) {
    cat(">>> netnwork$nodes$",i, sep="", "\n")
    head(network$nodes[[i]]) %>% print()
    print("...")
    tail(network$nodes[[i]]) %>% print()
  }
  
  # check whether the node list and the nodes attached top the edges map
  nodes = network2nodes(network)
  diff1 = setdiff(nodes, network$nodes$id)
  if (length(diff1) > 0) {
    cat("WARNING: ", length(diff1), " nodes are in the node list, but not in the edge list\n")
    diff1 %>% head(5) %>% print()
  }
  diff2 = setdiff(network$nodes$id, nodes)
  if (length(diff2) > 0) {
    cat("WARNING: ", length(diff2), " nodes are in the edge list, but not in the node list\n")
    diff2 %>% head(5) %>% print()
  }
  
  cat("The network has", length(network$edges$id), "edges and", length(network$nodes$id), "nodes\n")
  
}

##########################################################################
##########################################################################
# 
# # reduce the network to exclude STAT
# traitnetnodes = fullnet %>% network2nodes() 
# traitnetnodes = traitnetnodes[grep("STAT: ", traitnetnodes, invert = TRUE)]
# traitnet = nodes2network(traitnetnodes, fullnet) 
# 
# # get the STATnodes
# STATnetnodes = fullnet %>% network2nodes() 
# STATnetnodes = STATnetnodes[grep("STAT: ", STATnetnodes, invert = FALSE)]

##########################################################################
##########################################################################
# "cg19693031:chr1:144152909 TXNIP" is TXNIP
# 3485-28_2 is B2M
##########################################################################
if (FALSE) { # test code

  summary_network(fullnet)

  TXNIP = fullnet$nodes$id[grep("cg19693031", fullnet$nodes$id)]
  neighbors(TXNIP, fullnet)

  LEPR = fullnet$nodes$id[grep("Leptin receptor", fullnet$nodes$id)]
  neighbors(LEPR, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = 20)
  neighbors(LEPR, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = 200)
  
  neighbors("STAT: AGE", fullnet)
  neighbors("SOMA: C9 : Complement component C9", fullnet)
  neighbors("some error", fullnet)
  neighbors(c(), fullnet)
  
  maxneighbors(c("1","2"), fullnet)
  maxneighbors("cg19693031", fullnet, limit = 100)
  
  neighbors(c(TXNIP), fullnet) 
  neighbors(c(TXNIP), fullnet) %>% neighbors(fullnet) 
  neighbors(c(TXNIP), fullnet) %>% neighbors(fullnet) %>% neighbors(fullnet) 
  
  subnet = maxneighbors(c(TXNIP), fullnet, limit = 100) %>%  nodes2network(fullnet)
  summary_network(subnet)

  subnet = maxneighbors_noSTAT(c(TXNIP), fullnet, limit = 100) %>%  nodes2network(fullnet)
  summary_network(subnet)
  
  neighbors(TXNIP, subnet)
  neighbors(TXNIP, fullnet)
  
}
##########################################################################
##########################################################################

# add colors for platforms to network
# to be used as: PLATcols[fullnet$nodes$plat]
PLATlist = unique(fullnet$nodes$plat)

PLATcols = rep("#999999",length(PLATlist))
names(PLATcols) = PLATlist

PLATshapes = rep("star",length(PLATlist))
# PLATshapes = rep("icon",length(PLATlist))
names(PLATshapes) = PLATlist

# define PLAT colors manually
PLATcols["DNA"] = "#23bbee"
PLATcols["SOMA"] = "#a62281"
PLATcols["BRAIN"] = "#f2921f"
PLATcols["BM"] = "#ffc815"
PLATcols["LD"] = "#ffc815"
PLATcols["CPG"] = "#145da9"
PLATcols["CLIN"] = "#a0b6a8"
PLATcols["CM"] = "#57ba47"
PLATcols["IgA"] = "#e41d30"
PLATcols["RNA"] = "#5c2d83"
PLATcols["HD4"] = "#57ba47"
PLATcols["IgG"] = "#e41d30"
PLATcols["miRNA"] = "#5c2d83"
PLATcols["OLINK"] = "#a62281"
PLATcols["PGP"] = "#e41d30"
PLATcols["PM"] = "#57ba47"
PLATcols["SM"] = "#57ba47"
PLATcols["UM"] = "#57ba47"

PLATcols["STAT"] = "#EEEEEE"
PLATcols["GWAS"] = "#EEEEEE"

# define PLAT colors manually
PLATshapes["UM"] = "square"
PLATshapes["CM"] = "square"
PLATshapes["SM"] = "triangle"
PLATshapes["DNA"] = "diamond"
PLATshapes["RNA"] = "diamond"
PLATshapes["CPG"] = "diamond"

PLATshapes["STAT"] = "circle"
PLATshapes["GWAS"] = "dot"

fullnet$nodes$color = PLATcols[fullnet$nodes$plat]

fullnet$nodes$group = fullnet$nodes$plat

fullnet$nodes$shape = PLATshapes[fullnet$nodes$plat]

fullnet$edges$title = paste(fullnet$edges$type, ": ", fullnet$edges$id, 
                            ", p=", fullnet$edges$pvalue,
                            ", beta=", 
                            ifelse(fullnet$edges$sign >0, "", "-"), fullnet$edges$weight,
                            sep = "")
fullnet$edges$color = ifelse(fullnet$edges$sign > 0, "blue", "red")

# list of ids to choose from
plat_list = c("ALL", fullnet$nodes$plat %>% unique %>% sort())

# limit the number of nodes here
max_nodes_list = c(1,20,40,60,80,100,150,200)

# define a shiny server
server <- function(input, output, session) {

  # storage for reactiveValues
  storage <- reactiveValues(
    focus = "CLIN: HbA1c (%)",
    plat = "CLIN"
  )

  
  # read the URL parameter from session$clientData$url_search
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['focus']])) {
      storage$focus = query[['focus']]
    }
    if (!is.null(query[['maxnodes']])) {
      updateSelectInput(session, "maxnodes", selected = query[['maxnodes']])    
    }
  })
  
  # select trait pair from association tables  
  observeEvent(input$GWAStable_rows_selected,{storage$asso_selected = GWAS$ID[as.numeric(input$GWAStable_rows_selected)]})
  observeEvent(input$EWAStable_rows_selected,{storage$asso_selected = EWAS$ID[as.numeric(input$EWAStable_rows_selected)]})
  observeEvent(input$RWAStable_rows_selected,{storage$asso_selected = RWAS$ID[as.numeric(input$RWAStable_rows_selected)]})
  observeEvent(input$MBHtable_rows_selected,{storage$asso_selected = MBH$ID[as.numeric(input$MBHtable_rows_selected)]})
  observeEvent(input$GGMtable_rows_selected,{storage$asso_selected = GGM$ID[as.numeric(input$GGMtable_rows_selected)]})
  observeEvent(input$GENOtable_rows_selected,{storage$asso_selected = GENO$ID[as.numeric(input$GENOtable_rows_selected)]})
  observeEvent(input$STATtable_rows_selected,{storage$asso_selected = STAT$ID[as.numeric(input$STATtable_rows_selected)]})
  observeEvent(input$CATAtable_rows_selected,{storage$asso_selected = CATA$ID[as.numeric(input$CATAtable_rows_selected)]})
  
  # set focus if pull-down menu changes
  observeEvent(storage$asso_selected,
               {
                 message("storage$asso_selected:", storage$asso_selected)

                 if (substr(storage$asso_selected,1,4) == "STAT") {
                   # special case for STAT
                   ix = which(STAT$ID == storage$asso_selected)   
                   storage$selected_TRAIT1 = TRAITID2SHORTNAME[STAT$TRAITID[ix]]
                   storage$selected_TRAIT2 = TRAITID2SHORTNAME[STAT$TRAITID[ix]]
                   if (is.na(storage$selected_TRAIT1)){
                     storage$selected_TRAIT1 = ""
                     storage$selected_TRAIT2 = ""
                   }
                 } else {
                   ix = which(all_edges$id == storage$asso_selected)
                   storage$selected_TRAIT1 = all_edges$from[ix]
                   storage$selected_TRAIT2 = all_edges$to[ix]
                 }
                 message("storage$selected_TRAIT1:", storage$selected_TRAIT1)
                 message("storage$selected_TRAIT2:", storage$selected_TRAIT2)
               }
  )             
               
  
  # set focus if pull-down menu changes
  observeEvent(input$trait,
               {
                 if (input$trait != "") {
                   storage$focus <- input$trait
                   # message("focus set by dropdown to ", storage$focus)
                 }
               })
  
  # set focus if network is clicked
  observeEvent(input$network_selected,
               {
                 if (input$network_selected != "") {
                   storage$focus <- input$network_selected 
                   message("focus set by network click to ", storage$focus)
                 }
               })

  # observe click onto SubmitPair
  observeEvent(input$submitInfo, {
    updateTabItems(session = session, inputId = "tabs", selected = "Tab1")
    storage$focus = "special:pair"
    message("focus set by table click to ", storage$focus)
  })
  
  # text output of focus
  output$focus <- renderText({
    paste("Click here to focus on", storage$focus)
  })

  # text output of asso_selected
  output$asso_selected <- renderText({
    if (length(storage$selected_TRAIT1) == 0) {
      paste("Select an association, then click here to view the network context")
    } else if (storage$selected_TRAIT1 != "") {
      paste("Click here to focus network on association ", storage$asso_selected)
    } else {
      paste("Sorry - no data for association ", storage$asso_selected)
    }   
  })
  
  output$SelectTrait = renderUI({
    
    if (input$plat == "ALL") {
      traitlist = c(fullnet$nodes$id %>% sort())
    } else {
      ix = which(fullnet$nodes$plat == input$plat)
      traitlist = c(fullnet$nodes$id[ix] %>% sort())
    }
    
    tagList(
      selectInput(inputId = "trait", label = "select trait", 
                  choices = traitlist, selected = "")
    )
  })
  
  output$network <- renderVisNetwork({
    
    maxnodes = input$maxnodes
    if (length(maxnodes) == 0) {maxnodes = max_nodes_list[2]}
    maxnodes = as.numeric(maxnodes)
    message("limiting to maxnodes = ", maxnodes)
    message("GoButton pressed: ", input$GoButton)
    message("submitInfo pressed: ", input$submitInfo)
    message("storage$focus:", isolate(storage$focus))

    # to fix some weird bug
    zwi = isolate(storage$focus)
    zwi = zwi[1]
    
    if (length(zwi) == 0) {
      message("storage$focus empty, returning")
      act_trait = ""
      return()
    } else if (zwi == "") {
      message("storage$focus empty, returning")
      act_trait = ""
      return()
    } else if (zwi == "special:pair") {
      message("selected an edge")
      act_trait = c(isolate(storage$selected_TRAIT1), isolate(storage$selected_TRAIT2))
      # storage$focus = storage$selected_TRAIT1 # set focus to the first trait of the edge
      storage$focus = c(isolate(storage$selected_TRAIT1), isolate(storage$selected_TRAIT2))
    } else {
      act_trait = isolate(storage$focus)
    }

    message("getting subnet for act_trait: '", act_trait, "'")
    
    # get a network that is tractable
    # subnet = neighbors(act_trait, fullnet) %>% maxneighbors(fullnet, limit = maxnodes) %>% nodes2network(fullnet)
    subnet = neighbors(act_trait, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = maxnodes) %>% nodes2network(fullnet)
    # summary_network(subnet)
    
    # store in a reactive variable
    storage$subnet = subnet
    
    edges = data.frame(from = subnet$edges$from, 
                       to = subnet$edges$to,
                       title = subnet$edges$title,
                       color = subnet$edges$color)
    nodes = data.frame(id = subnet$nodes$id, 
                       label = subnet$nodes$id,
                       title = subnet$nodes$id,
                       color = subnet$nodes$color,
                       shape = subnet$nodes$shape)
    #group = subnet$nodes$group
    
    message("rendering subnet with ", length(nodes$id), " nodes and ", length(edges$from), " edges for act_trait: ", act_trait)
    
    visNetwork(nodes, edges, height = "500px", width = "100%") %>% 
      visNodes(shadow = list(enabled = TRUE, size = 10), 
               scaling = list(min=10,max=30)) %>%
      visLayout(randomSeed = 4711) %>% 
      visOptions(nodesIdSelection = list(enabled = TRUE, style = 'width: 1px; height: 1px;')) %>% 
      visPhysics(stabilization = FALSE) %>% 
      visEdges(smooth = TRUE)
    
  })
 
  # export the current edges as a table
  output$EXPORTtable = DT::renderDataTable({
    
    message("UPDATING EXPORTtable")
    zwi = storage$subnet$edges

    out = data.frame(ID = zwi$id,
                     TYPE = zwi$type,
                     TRAITID1 = SHORTNAME2TRAITID[zwi$from],
                     TRAIT1 = zwi$from,
                     BETA = zwi$sign*zwi$weight,
                     PVALUE = zwi$pvalue,
                     TRAITID2 = SHORTNAME2TRAITID[zwi$to],
                     TRAIT2 = zwi$to
                     )
    out
    }, filter = 'top', rownames = FALSE, selection = 'single', server = FALSE, extensions = 'Buttons', options = list(pageLength = 10, scrollX = TRUE, dom = 'Blfrtip', 
                                                                                                      buttons = 
                                                                                                        list("copy", list(
                                                                                                          extend = "collection"
                                                                                                          , buttons = c("csv", "excel", "pdf")
                                                                                                          , text = "Download"
                                                                                                        ) ) # end of buttons customization
                                                                                                      ))
  
   
  ############# mytable
  output$GWAStable = DT::renderDataTable({GWAS}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$EWAStable = DT::renderDataTable({EWAS}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$RWAStable = DT::renderDataTable({RWAS}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$MBHtable = DT::renderDataTable({MBH}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$GGMtable = DT::renderDataTable({GGM}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$GENOtable = DT::renderDataTable({GENO}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$CATAtable = DT::renderDataTable({CATA}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  output$STATtable = DT::renderDataTable({STAT}, filter = 'top', rownames = FALSE, selection = 'single', options = list(pageLength = 10, scrollX = TRUE))
  
}

#####################################################################
# UI dashboard
#####################################################################
ui <- dashboardPage(skin="red",
  
  ##################### TITLE ####################
  dashboardHeader(title = "QMDiab comics"),
  
  ##################### SIDEBAR ##################
  dashboardSidebar( width = 150,
    sidebarMenu( id = "tabs",
      menuItem("Network", tabName = "Tab1", icon = icon("dashboard")),
      menuItem("Tables", tabName = "Tab2", icon = icon("th")),
      menuItem("Export", tabName = "Tab3", icon = icon("th")),
      menuItem("HowTo", tabName = "Tab4", icon = icon("th")),
      menuItem("Use-cases", tabName = "Tab5", icon = icon("th")),
      menuItem("About", tabName = "Tab6", icon = icon("th"))
    )
  ),
  
  ##################### BODY #####################
  dashboardBody(
    tabItems(
      
      ##################### TAB1 #####################
      tabItem( tabName = "Tab1",
        fluidPage(
          # Boxes need to be put in a row (or column)
          fluidRow(
            column( width = 2,
              # controls above the network
              selectInput(inputId = "plat", label = "select platform", choices = plat_list, selected = plat_list[1])
            ),
            column( width = 6,
              uiOutput("SelectTrait")
            ),
            column( width = 2,
                    selectInput(inputId = "maxnodes", label = "max nodes", 
                                choices = max_nodes_list, selected = max_nodes_list[2])
            )
          ),
          fluidRow(
            column( width = 9,
              # the network itself
              actionButton("GoButton", textOutput("focus")),
              visNetworkOutput("network")
            ),
            column(width = 3,
                   img(src="legend_grey.jpg", width = 150)
            )
          )
        )
      ),
      ##################### TAB2 #####################
      tabItem( tabName = "Tab2",
        fluidPage(
          fluidRow(
            ###################################
            
            h3("Data used in the network"),
            # jump to network with selection
            actionButton(inputId = "submitInfo", label = textOutput("asso_selected")),
            p("."),
            tabBox( width = 32,
              # title = "Data used in the network",
              tabPanel("GWAS", DT::dataTableOutput("GWAStable")),
              tabPanel("EWAS", DT::dataTableOutput("EWAStable")),  
              tabPanel("RWAS", DT::dataTableOutput("RWAStable")),  
              tabPanel("MBH", DT::dataTableOutput("MBHtable")),  
              tabPanel("GGM", DT::dataTableOutput("GGMtable")),  
              tabPanel("GENO", DT::dataTableOutput("GENOtable")),  
              tabPanel("CATA", DT::dataTableOutput("CATAtable")),
              tabPanel("STAT", DT::dataTableOutput("STATtable"))
            )      
            ###################################
          )
        )
      ),
    
      ##################### TAB3 #####################
      tabItem( tabName = "Tab3",
        fluidPage(
          fluidRow(
            column( width = 12,
              h3("Export"),
              p("Download or view the network data of the active view"),
              DT::dataTableOutput("EXPORTtable")
            )
            ###################################
          )
        )
      ),
      ##################### TAB4 #####################
      tabItem( tabName = "Tab4",
        fluidPage(
          fluidRow(
            column( width = 12,
                h3("HowTo"),
                includeHTML("howto.html")
            )
            ###################################
          )
        )
      ),
      ##################### TAB5 #####################
      tabItem( tabName = "Tab5",
               fluidPage(
                 fluidRow(
                   column( width = 12,
                           h3("Use-cases"),
                           includeHTML("usecases.html")
                   )
                   ###################################
                 )
               )
      ),
      ##################### TAB6 #####################
      tabItem( tabName = "Tab6",
               fluidPage(
                 fluidRow(
                   column( width = 12,
                           h3("About"),
                           img(src="about.jpg", width = "100%")
                   )
                   ###################################
                 )
               )
      )
      ####################################################
    )
  )
)


#####################################################################
#####################################################################


#####################################################################
#####################################################################

# run the shiny app
shinyApp(ui = ui, server = server)