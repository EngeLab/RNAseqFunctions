
#' jet.colors
#'
#' @name jet.colors
#' @rdname jet.colors
#' @author Martin Enge
#' @examples
#'
#' jet.colors()
#'
#'
NULL
#' @export
#' @importFrom grDevices colorRampPalette

jet.colors <- function() {
    colorRampPalette(
    c(
    "#00007F", "blue", "#007FFF",
    "cyan", "#7FFF7F", "yellow",
    "#FF7F00", "red", "#7F0000"
    )
    )
}


#' norm.log.counts
#'
#' @name norm.log.counts
#' @rdname norm.log.counts
#' @param counts Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

norm.log.counts <- function(counts) {
    norm.fact <- colSums(counts)
    counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
    counts.log <- log2(counts.norm)
    counts.log
}

#' load.count.data
#'
#' @name load.count.data
#' @rdname load.count.data
#' @param fname Enter description.
#' @param mincount Enter description.
#' @param omit.bad.genes Enter description.
#' @param omit.bad.cells Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

load.count.data <- function(fname, mincount = 4e5, omit.bad.genes=TRUE, omit.bad.cells=TRUE) {
    counts <- read.table(fname, header=T, sep="\t")
    
    gene.names <- counts[[1]]
    counts <- as.matrix(counts[,-1])
    rownames(counts) <- gene.names
    
    ercc.counts <- counts[grepl("^ERCC\\-[0-9]*$", gene.names),]
    counts <- counts[!grepl("^ERCC\\-[0-9]*$", gene.names),]
    counts[is.na(counts)] <- 0
    last3.counts <- counts[(dim(counts)[1]-5):dim(counts)[1],]
    counts <- counts[1:(dim(counts)[1]-5),]
    
    
    # omit genes with no counts
    if(omit.bad.genes) {
        counts <- counts[rowSums(counts)>0,];
        
    }
    # omit cells with very poor coverage
    hist(colSums(counts), breaks=100)
    abline(v=mincount, col="red")
    ercc.counts <- ercc.counts[,colSums(counts)>mincount];
    last3.counts <- last3.counts[,colSums(counts)>mincount];
    counts <- counts[,colSums(counts)>mincount];
    
    list(counts=counts, ercc.counts=ercc.counts, last3.counts=last3.counts)
}

my.plot.callback <- function(x) {
    plot(x)
}

#' run.tsne
#'
#' @name run.tsne
#' @rdname run.tsne
#' @param my.dist Enter description.
#' @param plot.callback Enter description.
#' @param k Enter description.
#' @param max_iter Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

run.tsne <- function(my.dist, plot.callback=my.plot.callback, k=3, max_iter=5000, ...) {
    require(tsne)
    my.tsne <- tsne(my.dist, k=k, epoch_callback=plot.callback, initial_dims=50, max_iter=max_iter, ...)
    rownames(my.tsne) <- rownames(my.dist)
    my.tsne
}

#' geneset.colors
#'
#' @name geneset.colors
#' @rdname geneset.colors
#' @param gene Enter description.
#' @param counts.log Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

geneset.colors <- function(gene, counts.log) {
    if(!is.na(match(gene[1], rownames(counts.log)))) {
        #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
        #        panc.genes.max <- panc.genes.count
        panc.g.col1 <- counts.log[gene[1],]
        panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
        panc.g.col2 <- rep(0, dim(counts.log)[2])
        panc.g.col3 <- rep(0, dim(counts.log)[2])
        if(length(gene) > 1) {
            cat("two colors!")
            panc.g.col2 <- counts.log[gene[2],]
            panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
        }
        if(length(gene) > 2) {
            cat("three colors!")
            panc.g.col3 <- counts.log[gene[3],]
            panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
        }
        #        sum(panc.genes.max <= 1)
        #        panc.genes.max[panc.genes.max <= 1] <- 1
        #        panc.genes.max[panc.genes.max == 0] <- 0.0001
        #        panc.g.col <- as.numeric(log2(panc.genes.max))
        col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
        #        cat(col)
        return(col)
    }
    return(NA)
}

#' vecs2rgb
#'
#' @name vecs2rgb
#' @rdname vecs2rgb
#' @param r Enter description.
#' @param g Enter description.
#' @param b Enter description.
#' @param NA.ret Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

vecs2rgb <- function(r, g=NULL, b=NULL, NA.ret = NA) {
    r[is.na(r)] <- NA.ret
    r <- (r-min(r))/(max(r)-min(r))
    G <- rep(0, length(r))
    B <- rep(0, length(r))
    if(!is.null(g)) {
        G[is.na(G)] <- NA.ret
        G <- (g-min(g))/(max(g)-min(g))
    }
    if(!is.null(b)) {
        B[is.na(B)] <- NA.ret
        B <- (b-min(b))/(max(b)-min(b))
    }
    rgb(r, G, B)
}

#' plot.coreg
#'
#' @name plot.coreg
#' @rdname plot.coreg
#' @param geneset Enter description.
#' @param counts.log Enter description.
#' @param my.tsne Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

plot.coreg <- function(geneset, counts.log, my.tsne) {
    require(rgl)
    for(gene in geneset) {
        #    invisible(readline(prompt="Press [enter] to continue"))
        if(!is.na(match(gene[1], rownames(counts.log))) & !is.na(match(gene[2], rownames(counts.log)))) {
            #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
            #        panc.genes.max <- panc.genes.count
            panc.g.col1 <- counts.log[gene[1],]
            panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
            panc.g.col2 <- rep(0, dim(counts.log)[1])
            panc.g.col3 <- rep(0, dim(counts.log)[1])
            if(length(gene) > 1) {
                cat("two colors!")
                panc.g.col2 <- counts.log[gene[2],]
                panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
            }
            if(length(gene) > 2) {
                cat("three colors!")
                panc.g.col3 <- counts.log[gene[3],]
                panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
            }
            #        sum(panc.genes.max <= 1)
            #        panc.genes.max[panc.genes.max <= 1] <- 1
            #        panc.genes.max[panc.genes.max == 0] <- 0.0001
            #        panc.g.col <- as.numeric(log2(panc.genes.max))
            col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
            #        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
            plot3d(my.tsne, col=col, size=10, main=gene)
        }
        else {
            cat(gene, " not found\n");
        }
    }
}

#' get.cutoff.lognorm
#'
#' Select a cutoff for bad cells based on actin expression.
#'
#' @name get.cutoff.lognorm
#' @rdname get.cutoff.lognorm
#' @param my.counts.log Enter description.
#' @param quantile.cut Enter description.
#' @param gene.name Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

get.cutoff.lognorm <- function(my.counts.log, quantile.cut=0.001, gene.name='ACTB') {
    cl.act <- my.counts.log[gene.name,]
    cl.act.m <- median(cl.act)
    cl.act.sd <- sqrt(sum((cl.act[cl.act > cl.act.m] - cl.act.m)^2)/(sum(cl.act  > cl.act.m)-1))
    my.cut <- qnorm(p=quantile.cut, mean=cl.act.m, sd=cl.act.sd)
    my.cut
}

#' pick.cluster
#'
#' @name pick.cluster
#' @rdname pick.cluster
#' @param layout Enter description.
#' @param classification Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

pick.cluster <- function(layout, classification) {
    plot(layout, col=classification, pch=19, cex=0.8)
    #    invisible(readline(prompt="Press [enter] to pick clusters"))
    selected.clusts <- unique(classification[identify(layout, labels=classification)])
    #    invisible(readline(prompt="Press [enter] to show selected cells"))
    selected.cells <- !is.na(match(classification, selected.clusts))
    plot(layout, col=selected.cells+1, pch=19, cex=0.8)
    selected.cells
}

#' pick.clusters
#'
#' @name pick.clusters
#' @rdname pick.clusters
#' @param layout Enter description.
#' @param classification Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

pick.clusters <- function(layout, classification) {
    plot(layout, col=classification, pch=19, cex=0.8)
    #    invisible(readline(prompt="Press [enter] to pick clusters"))
    selected.clusts <- unique(classification[identify(layout, labels=classification)])
    selected.cells <- rep("undef", length=length(classification))
    a <- invisible(readline(prompt="Write name of cluster, [Enter] to stop"))
    while(a != "") {
        selected.cells[!is.na(match(classification, selected.clusts))] <- a
        plot(layout, col=as.integer(as.factor(selected.cells))+1, pch=19, cex=0.8)
        invisible(readline(prompt="Press [enter] to pick clusters"))
        plot(layout, col=classification, pch=19, cex=0.8)
        selected.clusts <- unique(classification[identify(layout, labels=classification)])
        a <- invisible(readline(prompt="Write name of cluster, [Enter] to stop"))
    }
    plot(layout, col=as.integer(as.factor(selected.cells))+1, pch=19, cex=0.8)
    selected.cells
}

#' sel.by.feature
#'
#' @name sel.by.feature
#' @rdname sel.by.feature
#' @param counts Enter description.
#' @param annot Enter description.
#' @param feature Enter description.
#' @param featureVals Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

sel.by.feature <- function(counts, annot, feature, featureVals) {
    o <- !is.na(match(annot[,feature], featureVals))
    o.n <- annot$rmangled.name[o]
    o <- !is.na(match(colnames(counts), o.n))
}

#' revcomp
#'
#' @name revcomp
#' @rdname revcomp
#' @param dna Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

revcomp <- function(dna) {
    sapply(strsplit(chartr("ATGC", "TACG", dna), NULL), function(x) {paste(rev(x), collapse='')})
}

#' col.from.targets
#'
#' @name col.from.targets
#' @rdname col.from.targets
#' @param targets vector, ngenes long
#' @param values matrix, ngenes x ncells long
#' @return char vector: ncells long color list
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

col.from.targets <- function(targets, values) {
    if(is.matrix(values)) {
        targets <- targets[1:(dim(values)[1])]
        v <- t(apply(values, 1, function(x) {(x-min(x))/(max(x)-min(x))}))
        fractions <- apply(v, 2, function(x) {x/sum(x)})
        fractions[is.na(fractions)] <- 1.0/(dim(fractions)[1])
        targets.rgb <- col2rgb(targets)
        res <- vector("character", length=length(targets))
        for(i in 1:(dim(values)[2])) {
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(x) {(255-x) * v[,i]}))
            mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(x) {x * fractions[,i]})))
            #        cat(fractions)
            #        cat(i)
            res[i] <- rgb(red=mytarget.rgb['red']/256, green=mytarget.rgb['green']/256, blue=mytarget.rgb['blue']/256)
        }
        return(res)
    } else {
        v <- (values-min(values))/(max(values)-min(values))
        #        fractions <- apply(v, 2, function(values) {values/sum(values)})
        #        fractions[is.nan(fractions)] <- 1.0/(dim(fractions)[1])
        targets.rgb <- col2rgb(targets)
        res <- vector("character", length=length(targets))
        for(i in 1:length(values)) {
            mytarget.rgb <- 255-t(apply(targets.rgb, 1, function(values) {(255-values) * v[i]}))
            #            mytarget.rgb <- rowSums(t(apply(mytarget.rgb, 1, function(values) {values * fractions[,i]})))
            res[i] <- rgb(red=mytarget.rgb[1]/256, green=mytarget.rgb[2]/256, blue=mytarget.rgb[3]/256)
        }
        return(res)
    }
}

#' plot.nice
#'
#' @name plot.nice
#' @rdname plot.nice
#' @param layout Enter description.
#' @param gene.vals Enter description.
#' @param pal Enter description.
#' @param maxVal Enter description.
#' @param cex Enter description.
#' @param rim.modifier Enter description.
#' @param legend.pos Enter description.
#' @param rim.col Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export

plot.nice <- function(layout, gene.vals, genes, pal=NULL, maxVal=NULL, cex=1, rim.modifier=0.85, legend.pos="bottomright", rim.col="black") {
    if(is.null(pal)) {
        require(RColorBrewer)
        pal <- brewer.pal(9, "Set1")
    }
    targets <- pal[1:length(genes)]
    values <- gene.vals[genes,]
    if(!is.null(maxVal)) {
        if(maxVal < 1) {
            values <- t(apply(values, 1, function(x) {
                x[x > quantile(x[x != min(x)], maxVal)] <- quantile(x[x != min(x)],maxVal)
                x
            }))
        } else {
            values <- apply(values, 2, function(x) {x[x > maxVal] <- maxVal})
        }
        #        values[values > max] <- max
    }
    cols <- col.from.targets(targets, values)
    plot(layout, col=rim.col, pch=19, cex=cex)
    points(layout, col=cols, pch=19, cex=cex*rim.modifier)
    legend(legend.pos, legend=genes, fill=pal)
}

#' plot.nice.whist
#'
#' @name plot.nice.whist
#' @rdname plot.nice.whist
#' @param coords Enter description.
#' @param gene.vals Enter description.
#' @param genes Enter description.
#' @param pal Enter description.
#' @param cex Enter description.
#' @author Martin Enge
#' @examples
#'
#'
#'
#'
NULL
#' @export
#' @importFrom RColorBrewer brewer.pal

plot.nice.whist <- function(coords, gene.vals, genes, pal=NULL, cex=1) {
    if(is.null(pal)) {
        #require(RColorBrewer)
        pal <- brewer.pal(9, "Set1")
    }
    ngenes <- length(genes)
    ll <- rep(1, length.out=ngenes*2)
    ll[seq(1:ngenes)*2] <- (1:ngenes)+1
    zones=matrix(ll, ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5))
    targets <- pal[1:length(genes) ]
    values <- gene.vals[genes,]
    cols <- col.from.targets(targets, values)
    plot(coords, col="black", pch=19, cex=cex)
    points(coords, col=cols, pch=19, cex=cex*0.85)
    legend("bottomright", legend=genes, fill=pal)
    sapply(genes, function(gene) {hist(gene.vals[gene,], main=gene)})
}

