#testvoronoi.R
source("voronoi.R")

test.vni <- function() {

    # generate 10 points in the 100 x 100 square   
    points <- vni.genPoints(10,100,100); points

    # generate bounding points, edges, and cells
    bounding <- vni.genBounding(1,1); bounding

    # test edge intersection detection
    o <- vni.makePoint(0, 0); o
    p <- vni.makePoint(0, 1); p
    q <- vni.makePoint(2, 1); q
    r <- vni.makePoint(2, 2); r
    s <- vni.makePoint(0, 2); s
    t <- vni.makePoint(3, 3); t

    E_OR <- vni.makeEdge(o$x, o$y, r$x, r$y); E_OR
    E_PQ <- vni.makeEdge(p$x, p$y, q$x, q$y); E_PQ
    
    vni.ixEdge(E_OR, E_PQ)
    
    E_SR <- vni.makeEdge(s$x, s$y, r$x, r$y); E_SR
    E_OQ <- vni.makeEdge(o$x, o$y, q$x, q$y); E_OQ
    
    vni.ixEdge(E_SR, E_OQ)
    vni.ixEdge(E_OQ, E_OR)
    vni.ixEdge(E_SR, E_OR)
    
    # check multiple edges
    vni.ixEdges(E_PQ, list(E_SR, E_OR)) 
    
    # midpoints
    vni.midpointL2(o,p) 
    vni.midpointL2(o,q) 
    vni.midpointL2(o,r) 
    
    # perpendicular bisector
    pb_OR <- vni.getPBL2(E_OR,3,3); pb_OR
    pb_ix_OR <- vni.ixEdge(E_OR,pb_OR); pb_ix_OR
    pb_OR <- vni.getPB(E_OR,"L2",3,3); pb_OR
    pb_ix_OR <- vni.ixEdges(E_OR,pb_OR); pb_ix_OR
    
    # proximal detection
    E_RT <- vni.makeEdge(r$x, r$y, t$x, t$y); E_RT
    E_PS <- vni.makeEdge(p$x, p$y, s$x, s$y); E_PS
    pb_OQ <- vni.getPB(E_OQ,geometry = "L2",100,100); pb_OQ
    vni.sameSide(o, E_RT, pb_OQ)
    vni.sameSide(o, E_PS, pb_OQ)
    vni.sameSide(o, E_SR, pb_OQ)
    
    # edge equality
    vni.sameEdge(E_OR,E_RT)
    vni.sameEdge(E_PS,E_SR)
    E_RO <- vni.makeEdge(r$x, r$y, o$x, o$y)
    vni.sameEdge(E_OR, E_RO)
    
    # find indices of equivalent edges
    Some_Edges <- list(E_OR,E_PQ,E_SR,E_RO,E_OQ,E_RT) # 1 and 4 equivalent
    vni.findEdges(E_OR,Some_Edges) # should be c(1, 4)
    
    # delete indices from a list
    length(Some_Edges)
    Fewer_Edges <- vni.removeEdgeMatches(E_OR, Some_Edges); length(Fewer_Edges)
    
    # testing out plot
    extractLabelToList <- function(l,label='x') {
        X <- list()
        counter <- 0
        for (p in l) {
            counter <- counter + 1
            if (length(p) > 0) {
                X[[length(X)+1]] <- p[[label]]
            }
        }
        return(X)
    }
    extractLabel <- function(l,label='x') {
        X <- c()
        for (p in l) {
            if (length(p) > 0) {
                X <- c(X,p[[label]])
            }
        }
        return(X)
    }
    
    points <- list(o,p,q,r,s,t); points
    X <-extractLabel(points,"x"); X
    Y <-extractLabel(points,"y"); Y
    plot(X,Y)
    AllEdges <- list(E_OQ,E_OR,E_PQ,E_PS,E_RT,E_SR)
    edge_to_debug <- NULL
    plotEdges <- function(E, color="black") {
        if (length(E) > 0) {
            edge_to_debug <<- E
            startPoints <- extractLabelToList(E,1)
            endPoints <- extractLabelToList(E,2)
            startX <- extractLabel(startPoints,"x")
            startY <- extractLabel(startPoints,"y")
            endX <- extractLabel(endPoints,"x")
            endY <- extractLabel(endPoints,"y")
            segments(startX,startY,endX,endY, col=color)
        } else {
            print("Empty edge list")
        }
    }
    plotEdges(AllEdges,"pink")
    
    width <- 500
    height <- 500
    # points <- list(q,r,t)
    n <- 50
    points <- vni.genPoints(n,0.99*width,0.99*height) 
    geometry <- "L2"
   
    # testing extract label 
    X <-extractLabel(points,"x"); X
    Y <-extractLabel(points,"y"); Y
    plot(X,Y)
    
    voronoi_cells_l2 <- vni.genCells(points,geometry,width,height)
    cells <- voronoi_cells_l2$cells
    extractLabelToList(cells,"point")
    status <- "" 
    plotCells <- function(cells, cex = 0.5, pch = 16) {
        
        extractPointsXY <- function(cells) {
            points <- extractLabelToList(cells, "point")
            X <- extractLabel(points,"x")
            Y <- extractLabel(points,"y")
            return(list(X=X,Y=Y))
        }
        extractEdges <- function(cells) {
            edgesLists <- extractLabelToList(cells, "edges")
            allEdges <- list()
            for (l in edgesLists) {
                for (e in l) {
                    allEdges[[length(allEdges)+1]] <- e
                }
            }
            print(list(edgecount=length(allEdges)))
            return(allEdges)
        }
        
        status <<- "extracting points"
        pointsXY <- extractPointsXY(cells)
        status <<- "plotting points"
        
        plotTitle <- paste("Voronoi Diagram,",n,"random points")
        
        par(cex=cex,pch=pch)
        plot(pointsXY$X,pointsXY$Y, main=plotTitle,xlab="x",ylab="y")
        status <<- "extracting edges and plotting edges"
        allEdges <- extractEdges(cells)
        plotEdges(allEdges)
        
        print(list(edges=length(allEdges)))
        return(allEdges)
    } 
    #cells <- Cells
    allEdges <- plotCells(cells)
    # plotEdges(allEdges)
    
    # plotEdges(voronoi_cells_l2$edges)
    # ed <- voronoi_cells_l2$edges[[1]]        
    # segments(c(ed[[1]]$x), c(ed[[1]]$y), c(ed[[2]]$x), c(ed[[2]]$y), col="pink")
}

test.vni.plot <- function() {
    plotSize <- 600 
    defaultCellsCount <- 160
    defaultGeometry <- "L2"
    
    points <- NULL
    cells <- NULL
    cellsCount <- defaultCellsCount
    geometry <- defaultGeometry
    height <- 0.8 * plotSize
    width <- 0.8 * plotSize
    
    regeneratePoints <- function() {
      points <<- vni.genPoints(cellsCount, 0.99*width, 0.99*height)
        cells <<- vni.genCells(points, geometry, width, height)          
    }
    
    redrawPlot <- function(){
        regeneratePoints()
        if (length(cells) > 0) {
          print(paste0("Cells OK!"))
        }
        if (length(cells) > 0) {
            print(paste0("Cells still OK!",length(cells)))
        }
        allEdges <- vni.plotCells(cells, width, height)
    }
   
    redrawPlot() 
}
