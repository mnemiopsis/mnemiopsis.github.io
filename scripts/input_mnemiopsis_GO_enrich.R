#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
    message('usage: path/to/input_smed_GO_enrich.R input1 input2 input3
    Runs GO term enrichment for a list of genes.
    Output:  A table of enrichment, either "all_go_enrichment.xls", or "go_slim_enrichment.xls"
    Args:
    input1            File with genes to be enriched, one gene per line.
    input2            Either "all_go" or "go_slim".
    input3            Optional output prefix.
    -h, --help        to print help messages')
    q('no')
}


library(topGO)
library(reshape)
library(plyr)
library(GO.db)
library(dplyr)


run.get.go.terms <- function(index, list.terms){
    # Runs get.*.GO but checks that it is runable. If the get.*.GO function is not runable, and error message is returned.
    geneList             <- factor(as.integer(GOterms$Gene %in% list.terms[[index]]))
    names(geneList)      <- GOterms$Gene
    lev <- levels(geneList)
    if(length(lev) < 2){
        message(paste0("Pvalue too stringent for ", dx[index, "contrasts"], ". Choose a higher pvalue to yield results for this contrast."))
    } else {
        terms <- get.GOterms(geneList[[index]], 10, c("MF", "CC", "BP"))
        invisible(terms)
    }
}

get.bh <- function(go.id, df){
    # Returns the benjamini-hochberg correction so it can later be appended to the final GOterms.
    #
    # Args:
    #   go.id: a go.id to find the benjamini-hochberg correction for.
    #   df: the data.frame holding the bh correction values.
    bh     <- grep(go.id, rownames(df))
    bh.val <- df[bh, "BH.correction"]
    invisible(bh.val)
}


create.file <- function(enriched.terms, name){
    # Creats a file holding the casted table.
    #
    # Args:
    #   prefix: a string of what should be the file identifier.
    #   casted.data: a resultant table from the function "cast.table".
    file.create(name)
    zz   <- file(name, "w")
    write.table(enriched.terms, file=zz, sep = "\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    close(zz)
}


Get.go.terms <- function(GOterms, list.terms, ontol, xx){
    # Runs the base line GO term analysis as provided through topGO.
    ## contrast             <- dx[index,"contrasts"]
    #terms                <- list.terms[[index]]
    #genes                <- table[terms, "smed_id"]
    geneList             <- factor(as.integer(GOterms$id %in% list.terms))
    names(geneList)      <- GOterms$id
    #contrast <- dx[index,"contrasts"]
    #geneList <- geneList[[index]]
    GO                   <- new("topGOdata", ontology=ontol, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = readmap)
    resultFisher         <- runTest(GO, algorithm = "classic", statistic = "fisher")
    pVal                 <- data.frame(pval=signif(score(resultFisher), 6),BH.correction=seq(1, length(score(resultFisher))))
    pVal                 <- pVal[order(pVal$pval),]
    pVal$BH.correction   <- signif(p.adjust(pVal$pval, method="BH"), 6)
    pVal.sub             <- pVal[pVal$BH.correction<0.01,]
    pVal.sub$go_id       <- rownames(pVal.sub)
    pVal.try <- cbind(pVal.sub, Term=sapply(pVal.sub$go_id, FUN=function(n){Term(xx[[n]])}), Ontol=sapply(pVal.sub$go_id, FUN=function(n){Ontology(xx[[n]])}))
    if (nrow(pVal.try) >0){
        df2                  <- pVal.try[,c("go_id", "Term", "Ontol", "pval", "BH.correction")]
        out <- list(df2, GO, resultFisher)
        invisible(out)
} else {
    message("No enrichment")
}
}

retrieve.my.genes <- function(final, readmap.inverse, terms){
    ## retrieve the genes for the GO terms that are in the readmapping file
    ## final: the final df to add the genes to.
    ## readmap.inverse: the readmapping to find the genes in
    ## terms: the list of GO terms
    ids     <- unlist(lapply(final$go_id, toString))
    a.small <- readmap.inverse[ids]
    my.terms <- lapply(a.small, function(x){
        new.terms <- na.omit(match(x, terms))
        b         <- terms[new.terms]
        invisible(b)}
                       )
    invisible(my.terms)
}

lookup.smed.zeros <- function(off.gos, all.readmap){
    ## lookup SMED genes present in the analysis for the offspring terms
    ##      off.gos: GO terms with zero SMED genes associated
    ##      all.readmap: the map of transcripts to GO terms for the analysis
    
    resul <- na.omit(lapply(off.gos, function(x){all.readmap[[x]]}))
    names(resul) <- names(off.gos)
    invisible(resul)
}


get.parent.length <- function(zeros.dx, Ontol){
    ## return the SMED genes associated with a parent GO term
    ##     zeros.dx: dataframe of the parent GO terms with zero values
    ##     Ontol: One of either "CC", "MF", "BP"
    require(GO.db)
    require(dplyr)
    if (nrow(zeros.dx) < 1){
        message(paste0(Ontol, " has no zero gene length terms for this analysis"))
    } else {
        if(Ontol == "CC"){## initialize ontology hierarchies
            ontol.off <- as.list(GOCCOFFSPRING)
        }
        if (Ontol == "BP"){
            ontol.off <- as.list(GOBPOFFSPRING)
        }
        if(Ontol == "MF"){
            ontol.off <- as.list(GOMFOFFSPRING)
        }
        ## lookup the offspring GO terms
        ontol.zero.off   <- unlist2(lapply(zeros.dx$go_id, function(x){ontol.off[x]}))
        ## lookup SMED genes for the offspring terms
        ontol.smed       <- lookup.smed.zeros(ontol.zero.off, readmap.inverse)
        ## pull out unique names of the GO terms
        ontol.names.uniq <- unique(names(ontol.smed))
        ## pull out the number of SMED genes associated with the parent GO terms
        ontol.ids.genes        <- lapply(ontol.names.uniq, function(x){unique(unlist(ontol.smed[grep(x, names(ontol.smed))]))})
        ontol.ids.len          <- unlist2(lapply(ontol.ids.genes, length))
        Genes                  <- unlist2(lapply(ontol.ids.genes, paste, collapse=", "))
        dff                    <- data.frame(GO=ontol.names.uniq, Length=ontol.ids.len, Genes=Genes, stringsAsFactors=FALSE)
        zeros.dx               <- zeros.dx[order(zeros.dx$go_id), ]
        dff                    <- dff[order(dff$GO), ]
        zeros.dx$Branch_Num_Genes     <- dff$Length
        zeros.dx$Branch_Genes         <- dff$Genes
        invisible(zeros.dx)
    }
}

create.hierarchy <- function(go.results, go.df, ontol, go.terms, out){
    if (nrow(go.df) == 0){
        message(paste0('Hierarchy for ', ontol, " cannot be assigned."))
    } else {
        ## creates ontology hierarchy 
        pvals <- go.df$BH.correction
        names(pvals) <- rownames(go.df)
        if (length(pvals) > 50){
            numNodes <- 50
        } else {
            numNodes <- length(pvals)
        }
        pdf(file=paste0(out, "/", go.terms, "_", ontol, "_hierarchy.pdf"), width=10, height=10)
        showSigOfNodes(go.results, pvals, firstSigNodes = numNodes, useInfo = "all")
        dev.off()
    }
}




##args <- c("genes.txt", "all_go")
if(length(args) == 3){
    outfile <- paste0(args[3], "_", args[2], "_enrichment.xls")
} else {
    outfile <- paste0(args[2], "_enrichment.xls")
}


if(args[2] == "all_go"){
    GOterms <- read.table("/home/klg/projects/Mnemiopsis/website/output/gene2go.txt", quote="", stringsAsFactors=FALSE, sep="\t", header=TRUE)
    colnames(GOterms) <- c("id", "term")
    readmap           <- readMappings("/home/klg/projects/Mnemiopsis/website/output/gene2go.txt", sep = "\t", IDsep = ";")
} else {
    GOterms <- read.table("/n/projects/klg/RNAseq_pipeline/smed/gene2GOslim.txt", quote="", stringsAsFactors=FALSE, sep="\t", header=TRUE)
    #colnames(GOterms) <- c("id", "term")
    readmap           <- readMappings("/n/projects/klg/RNAseq_pipeline/smed/gene2GOslim.txt", sep = "\t", IDsep = ",")
}


## edit this to put in custom read mappings
## GOterms <- read.table(args[2], quote="", stringsAsFactors=FALSE, sep="\t", header=TRUE)
## colnames(GOterms) <- c("id", "GO_terms")
## readmap           <- readMappings(args[2], sep = "\t", IDsep = ",")


terms <- scan(args[1], "character")



## read in the entire gene2allGo annotation and only keep the terms that are present in the chosen significant genes.
all.readmap <- readMappings("/home/klg/projects/Mnemiopsis/website/output/gene2go.txt", sep = "\t", IDsep = ";")
all.readmap <- all.readmap[terms]
readmap.inverse    <- inverseList(all.readmap)


## access the GO ontology hierarchy
xx <- as.list(GOTERM)


## run the enrichment analysis
bp <- Get.go.terms(GOterms, terms, "BP", xx)
cc <- Get.go.terms(GOterms, terms, "CC", xx)
mf <- Get.go.terms(GOterms, terms, "MF", xx)

bp.df <- bp[[1]]
mf.df <- mf[[1]]
cc.df <- cc[[1]]

bp.go <- bp[[2]]
mf.go <- mf[[2]]
cc.go <- cc[[2]]

bp.fish <- bp[[3]]
mf.fish <- mf[[3]]
cc.fish <- cc[[3]]


## bring to together all of the enrichments
final           <- rbind(mf.df, cc.df, bp.df)
rownames(final) <- NULL


## get the SMED genes associated with the enrichment and attach them to the enrichment data.frame
my.genes        <- retrieve.my.genes(final, readmap.inverse, terms)
my.genes        <- lapply(my.genes, unique)
len.my.genes    <- lapply(my.genes, length)
my.genes        <- lapply(my.genes, toString)
final$Num_Genes <- unlist(len.my.genes)
final$Genes     <- unlist(my.genes)



all.my.genes        <- readmap.inverse[final$go_id]
all.my.genes        <- lapply(all.my.genes, unique)
len.all.my.genes    <- lapply(all.my.genes, length)
all.my.genes        <- lapply(all.my.genes, toString)
final$Num_Leaf_Genes <- unlist(len.all.my.genes)
final$Leaf_Genes     <- unlist(all.my.genes)




## first filter the final df of the zero values for parent genes to be added back later

final.t                        <- filter(final, Num_Genes != 0)
final.t$Branch_Num_Genes     <- NA
final.t$Branch_Genes         <- NA


## next filter out any parent terms that have zero SMED Genes currently annotated
zero <- filter(final, Num_Genes == 0)




## initialize ontologies for parent GO terms that are not annotated to SMED genes
zero.cc <- filter(zero, Ontol == "CC")
zero.bp <- filter(zero, Ontol == "BP")
zero.mf <- filter(zero, Ontol == "MF")




## find the number of SMED genes associated with the parent ID
cc.final <- get.parent.length(zero.cc, "CC")
bp.final <- get.parent.length(zero.bp, "BP")
mf.final <- get.parent.length(zero.mf, "MF")




## Bind together everything and write out the final table
final.done <- rbind(final.t, cc.final, bp.final, mf.final)


if (nrow(final.done)>0){
    create.file(final.done, outfile)
} else {
    message("No enrichment")
}


## out.dir <- dirname(args[1])

## create.hierarchy(cc.go, cc.df, "CC", args[2], out.dir)
## create.hierarchy(mf.go, mf.df, "MF", args[2], out.dir)
## create.hierarchy(bp.go, bp.df, "BP", args[2], out.dir)

save.image(file = paste0(args[2], ".Rdata"))
