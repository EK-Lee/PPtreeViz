#' Tree with LDA.
#' 
#' Find tree structure using linear discriminant(LD) in each split.
#' @usage LDA.Tree(origclass, origdata, weight = TRUE, ...) 
#' @param origclass class information
#' @param origdata data  without class information
#' @param weight weight flag in LDA index
#' @param ... arguments to be passed to methods
#' @return Tree.Struct Tree structure of PPtree result
#' @return projbest.node 1-dim optimal projections of each split node
#' @return splitCutoff.node cutoff values of each split node 
#' @return origclass original class 
#' @return origdata original data
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- LDA.Tree(iris[,5],iris[,1:4])
#' Tree.result
LDA.Tree<-function(origclass, origdata, weight = TRUE, ...) 
{
    origdata <- as.matrix(origdata)
    Find.proj <- function(origclass, origdata, weight, ...) {
        n <- nrow(origdata)
        p <- ncol(origdata)
        g <- table(origclass)
        g.name <- as.numeric(factor(names(g)))
        G <- length(g)
        a.proj.best <- LDAopt(as.numeric(as.factor(origclass)),origdata,weight,q=1)$projbest
        proj.data <- as.matrix(origdata) %*% a.proj.best
        sign <- sign(a.proj.best[abs(a.proj.best) == max(abs(a.proj.best))])   
        index <- (1:p) * (abs(a.proj.best) == max(abs(a.proj.best)))
        index <- index[index > 0]
        if (G == 2) {
            class <- origclass
        }
        else {
            m <- tapply(c(proj.data), origclass, mean)
            sd <- tapply(c(proj.data), origclass, sd)
            sd.sort <- sort.list(sd)
            m.list <- sort.list(m)
            m.sort <- sort(m)
            m.name <- as.numeric(factor(names(m.sort)))
            G <- length(m)
            dist <- 0
            split <- 0
            for (i in 1:(G - 1)) {
                if (m[m.list[i + 1]] - m[m.list[i]] > dist) {
                  split <- i
                  dist <- m[m.list[i + 1]] - m[m.list[i]]
                }
            }
            class <- rep(0, n)
            for (i in 1:split) class <- class + 
                               (as.numeric(as.factor(origclass)) == m.name[i])
            class <- 2 - class
            g <- table(class)
            g.name <- as.numeric(names(g))
            G <- length(g)
            n <- nrow(origdata)
            a.proj.best <- LDAopt(as.numeric(factor(class)),origdata,weight,q=1)$projbest
            if (sign != sign(a.proj.best[index])) 
                a.proj.best <- -a.proj.best
            proj.data <- as.matrix(origdata) %*% a.proj.best
        }
        m.LR <- tapply(proj.data, class, mean)
        temp.list<-sort.list(m.LR)
        m.LR<-m.LR[temp.list]
        sd.LR <- tapply(proj.data, class, sd)[temp.list]
        IQR.LR <- tapply(proj.data, class, IQR)[temp.list]
        median.LR <- tapply(proj.data, class, median)[temp.list]
        n.LR <- table(class)[temp.list]
        
        c1 <- (m.LR[1] + m.LR[2])/2
        c2 <- (m.LR[1] * n.LR[2] + m.LR[2] * n.LR[1])/sum(n.LR)
        c3 <- (m.LR[1] * sd.LR[2] + m.LR[2] * sd.LR[1])/sum(sd.LR)
        c4 <- (m.LR[1] * sd.LR[2]/sqrt(n.LR[2]) + 
                 m.LR[2] * sd.LR[1]/sqrt(n.LR[1]))/
                 (sd.LR[1]/sqrt(n.LR[1])+sd.LR[2]/sqrt(n.LR[2]))
        c5 <- (m.LR[1] * IQR.LR[2] + m.LR[2] * IQR.LR[1])/sum(IQR.LR)      
        c6 <- (m.LR[1] * (IQR.LR[2]/sqrt(n.LR[2])) + 
                 m.LR[2] * (IQR.LR[1]/sqrt(n.LR[1])))/
               ((IQR.LR[1]/sqrt(n.LR[1]))+(IQR.LR[2]/sqrt(n.LR[2])))
        C <- c(c1, c2, c3, c4,c5,c6)
        
        Index <-LDAindex(as.numeric(as.factor(class)),proj.data,weight)
        Alpha <- t(a.proj.best)
        IOindexR <- NULL
        IOindexL <- NULL
        sort.LR <- as.numeric(names(sort(m.LR)))
        IOindexL <- class == sort.LR[1]
        IOindexR <- class == sort.LR[2]
        list(Index=Index,Alpha = Alpha, C = C, IOindexL = IOindexL, 
            IOindexR = IOindexR)
    }
    Tree.construct <- function(origclass, origdata, Tree.Struct,  id, 
                               rep, rep1, rep2, projbest.node, splitCutoff.node, ...) {
        origclass <- as.integer(origclass)
        n <- nrow(origdata)
        g <- table(origclass)
        G <- length(g)
        if (length(Tree.Struct) == 0) {
            Tree.Struct <- matrix(1:(2 * G - 1), ncol = 1)
            Tree.Struct <- cbind(Tree.Struct, 0, 0, 0, 0)
        }
        if (G == 1) {
            Tree.Struct[id, 3] <- as.numeric(names(g))
            list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
                splitCutoff.node = splitCutoff.node, rep = rep, rep1 = rep1, rep2 = rep2)
        } else {
            Tree.Struct[id, 2] <- rep1
            rep1 <- rep1 + 1
            Tree.Struct[id, 3] <- rep1
            rep1 <- rep1 + 1
            Tree.Struct[id, 4] <- rep2
            rep2 <- rep2 + 1
            a <- Find.proj(origclass,origdata,weight)
            splitCutoff.node <- rbind(splitCutoff.node, a$C)
            Tree.Struct[id, 5] <- a$Index
            projbest.node <- rbind(projbest.node, a$Alpha)
            t.class <- origclass
            t.data <- origdata
            t.class <- t.class * a$IOindexL
            t.n <- length(t.class[t.class == 0])
            t.index <- sort.list(t.class)
            t.index <- sort(t.index[-(1:t.n)])
            t.class <- t.class[t.index]
            t.data <- origdata[t.index, ]
            b <- Tree.construct(t.class, t.data, Tree.Struct, Tree.Struct[id, 2], 
                                rep, rep1, rep2, projbest.node, splitCutoff.node)
            Tree.Struct <- b$Tree.Struct
            projbest.node <- b$projbest.node
            splitCutoff.node <- b$splitCutoff.node
            rep <- b$rep
            rep1 <- b$rep1
            rep2 <- b$rep2
            t.class <- origclass
            t.data <- origdata
            t.class <- (t.class * a$IOindexR)
            t.n <- length(t.class[t.class == 0])
            t.index <- sort.list(t.class)
            t.index <- sort(t.index[-(1:t.n)])
            t.class <- t.class[t.index]
            t.data <- origdata[t.index, ]
            n <- nrow(t.data)
            G <- length(table(t.class))
            b <- Tree.construct(t.class, t.data, Tree.Struct, 
                           Tree.Struct[id, 3], rep, rep1, rep2, projbest.node,splitCutoff.node)
            Tree.Struct <- b$Tree.Struct
            projbest.node <- b$projbest.node
            splitCutoff.node <- b$splitCutoff.node
            rep <- b$rep
            rep1 <- b$rep1
            rep2 <- b$rep2
        }
        list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
            splitCutoff.node = splitCutoff.node, rep = rep, rep1 = rep1, rep2 = rep2)
    }
    splitCutoff.node <- NULL
    projbest.node <- NULL
    Tree.Struct <- NULL
    id <- 1
    rep1 <- 2
    rep2 <- 1
    rep <- 1

    Tree.final <- Tree.construct(origclass, origdata, Tree.Struct, 
        id, rep, rep1, rep2, projbest.node, splitCutoff.node)
    Tree.Struct <- Tree.final$Tree.Struct
    colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID","Coef.ID","Index")    
    projbest.node <- Tree.final$projbest.node
    splitCutoff.node <- Tree.final$splitCutoff.node
    treeobj<-list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
        splitCutoff.node = splitCutoff.node,origclass = origclass,origdata= origdata)
    class(treeobj)<-append(class(treeobj),"PPtree")
    return(treeobj)
}