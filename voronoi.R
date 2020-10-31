#voronoi.R

# generate a set of random points
library(pracma)
library(rlist)

vni.makePoint <- function(x,y) {
# x: numeric
# y: numeric
    return(list(x=x,y=y))
}

vni.validatePoint <- function(p) {
    x <- p$x
    y <- p$y
    if (is.null(x) || is.null(y)) {
        print("Invalid point:")
        print(c(x,y))
        return(FALSE)
    }
    return(TRUE)
}

vni.genPoints <- function(n=1, xmax, ymax) {
# n: numeric (an integer)
# x: numeric
# y: numeric 
# @return A list of random points
    xvals <- runif(n)*xmax
    yvals <- runif(n)*ymax
    P <- list()
    for (i in 1:length((xvals))) {
        q <- vni.makePoint(round(xvals[i],0), round(yvals[i],0))
        P[[length(P) + 1]] <- q
    }
    return(P)
}

vni.makeEdge <- function(x1, y1, x2, y2) {
    p1 <- vni.makePoint(x1,y1)
    p2 <- vni.makePoint(x2,y2)
    result <- list(p1,p2)
    if (vni.validateEdge(result)) {}
    return(result)
}

vni.validateEdge <- function(e) {
    p1 <- e[[1]]
    p2 <- e[[2]]
    if (!vni.validatePoint(p1) || !vni.validatePoint(p2)) {
        stop(list(msg="invalid edge!",edge=e))
    }
    return(TRUE)
}

# define basic operations on the 'cell' data type
vni.makeCell <- function(point, edges) {
# A cell is a point and a collection (list) of edges
    result <- list(point=point,edges=edges)
    return(result)
}

vni.validateCell <- function(c) {
    p <- c$point
    validPoint <- vni.validatePoint(p)
    el <- c$edges
    validEdges <- function(el) {
        for (e in el) {
            if (!vni.validateEdge(e)) {
                return(FALSE)
            } 
        }
        return(TRUE)
    }
    if (!validPoint || !validEdges(el)) {
        print("Invalid cell!")
        return(FALSE)
    }
    return(TRUE)
}


vni.makeCells <- function(P,E) {
    result <- list()
    if (length(P) == 0 || length(P) != length(E)) {
        return(NULL)
    }
    for (i in 1:length(P)) {
        cell <- vni.makeCell(P[[i]],E[[i]])
        validCell <- vni.validateCell(cell) 
        
        if (is.null(cell) || !validCell) {
            print("problem cell:")
            print(cell)
            stop("encountered an error making a batch of cells") 
        }
        result[[length(result)+1]] <- cell
    }
    return(result)
}





vni.ixEdge <- function(e, e_0) {
# algorithm adapted from this page:
# http://geomalgorithms.com/a05-_intersect-1.html
# returns point of intersection if edge e intersects e0
    e_p <- c(e[[1]]$x, e[[1]]$y)
    e_q <- c(e[[2]]$x, e[[2]]$y)
    e_v <- e_q - e_p # get vector by fixing one endpoint
    
    e_0_p <- c(e_0[[1]]$x, e_0[[1]]$y)
    e_0_q <- c(e_0[[2]]$x, e_0[[2]]$y)
    e_0_v <- e_0_q - e_0_p # get vector by fixing one endpoint
    if (length(e_v) != 2 || length(e_0_v) != 2) {
        print("e_p")
        print(e_p)
        print("e_q")
        print(e_q) 
        print("e_v")
        print(e_v) 
        print("e_0_p")
        print(e_0_p)
        print("e_0_q")
        print(e_0_q) 
        print("e_0_v")
        print(e_0_v) 
        stop("Problem with finding intersection")
    }
    
    
    
    w <- e_p - e_0_p
    
    perp <- function(u,v) {
        return(det(cbind(u,v)))
    }
    
    D <- perp(e_v, e_0_v)
    
    # test if parallel
    if (abs(round(D,3)) == 0) { # parallel
        # check for collinearity   
        if (perp(e_v,w) != 0 || perp(e_0_v, w) != 0) {
            # not collinear
            return(NULL) 
        }
        # get the squared lengths
        dev <- dot(e_v,e_v)
        de0v <- dot(e_0_v, e_0_v)
        if (dev == 0 && de0v == 0) {
            # both are points
            if (e_p == e_0_p) {
                # the same point 
                return(vni.makePoint(e_p[1], e_p[2]))
            } 
            # else different points
            return(NULL)
        }
        if (dev == 0) {
            # e is a point 
            # should check if it lies on e0
            # todo
            return(NULL)
        }
        if (de0v == 0) {
            # e0 is a point 
            # should check if it lies on e0
            # todo
        }
        # else collinear segments
        # todo
    }
    # skew and may intersect
    # get slope products
    s <- perp(e_0_v,w) / D 
    if (s < 0 || s > 1) {
        # no intersection
        return(NULL)
    }
    t <- perp(e_v, w) / D
    if (t < 0 || t > 1) {
        # no intersection
        return(NULL)
    }
    ix <- e_p + s*e_v
    return(vni.makePoint(ix[1], ix[2]))
}

vni.ixEdges <- function(e0, edges) {
    # returns the first intersection of edge e0 and the list of edges
    # NULL if none found
    for (i in 1:length(edges)) {
        ix <- vni.ixEdge(edges[[i]], e0)
        if (is.list(ix)) {
            return(ix) 
        }
    }
    return(NULL)
}


vni.getPBL1 <- function(e) {
 
}

vni.rotateR2 <- function(v, a=pi) {
    # rotates a two dimensional vector clockwise by angle a (in radians)
    R <- rbind(c(cos(a), -sin(a)),
               c(sin(a), cos(a)))
    result <- round(R%*%v,6)
    return(result[,1])
}

vni.midpointL2 <- function(p, q) {
    rx <- (p$x + q$x)/2
    ry <- (p$y + q$y)/2 
    return(c(rx,ry))
}

vni.getPBL2 <- function(e, width=400, height=400) {
    # we need to know the width and height so we can scale to infinity
    # get one vector corresponding to the edge
    v <- c(e[[2]]$x-e[[1]]$x,e[[2]]$y-e[[1]]$y)
    # get perpendicular vector
    w <- vni.rotateR2(v, pi/2)
    # get midpoint of the segment
    m <- vni.midpointL2(e[[1]],e[[2]])
    # perpendicular bisector is vertical
    if (w[1] == 0) {
        return(vni.makeEdge(m[1], height*2, m[1], -height*2))
    }
    # perpendicular bisector is horizontal
    if (w[2] == 0) {
        return(vni.makeEdge(width*2, m[2], -width*2,m[2]))
    }
    # default: pointing to the right
    x0 <- -3*width
    x1 <- 3*width
    y <- function(x,m,w) {
        return(w[2]/w[1]*(x-m[1])+m[2])
    }
    return(vni.makeEdge(x0, y(x0,m,w), x1, y(x1,m,w)))
}

vni.getPB <- function(e, geometry="L2", width=400, height=400) {
    # given an edge, returns a bisector as a LIST of edges
    if (geometry == "L2") {
        result <- list(vni.getPBL2(e,width, height))
        return(result)
    }
}

vni.sameSide <- function(c, e, bound) {
# c: a list with components x and y representing a point
# e: an edge made of two points like c
# bound: a list of edges forming a boundary
# this sort of inspired by the ray casting algorithm
#
# looking at rays from point c to each of the endpoints of e
    e_1 <- e[[1]]
    e_2 <- e[[2]]
    v1 <- vni.makeEdge(c$x, c$y, e_1$x, e_1$y)
    v2 <- vni.makeEdge(c$x, c$y, e_2$x, e_2$y)
    ix1 <- vni.ixEdges(v1, bound)
    ix2 <- vni.ixEdges(v2, bound)
    if (is.null(ix1) && is.null(ix2)) {
    # if there are no intersections, the edge is on the same side as the point
        return(TRUE)
    }
    # if both have an intersection with the boundary, we are on the other side
    # if only one has an intersection, the edge crosses the boundary
    # in both cases, not on the same side
    return(FALSE)
}

vni.sameEdge <- function(e1, e2) {
    # uses the triangle inequality to see if this is the same edge:
    # 
    # fix the first endpoint of e1
    # obtain the vectors pointing to the other endpoints of e2
    # the sum of their lengths must be exactly equal to the length of e1
    # 
    # e1: the first edge
    e1_a <- c(e1[[1]]$x, e1[[1]]$y) # first endpoint of e1
    e1_b <- c(e1[[2]]$x, e1[[2]]$y) # second endpoint of e1
    v_e1 <- e1_a - e1_b
    len2_e1 <- dot(v_e1,v_e1)
    #
    # e2: the second edge
    # first endpoint of e2
    e2_a <- c(e2[[1]]$x, e2[[1]]$y)
    # second endpoint of e2
    e2_b <- c(e2[[2]]$x, e2[[2]]$y)
    #
    # vectors from the first endpoint of e1 pointing to the endpoints of e2
    aa <- e2_a - e1_a # from e1_a to e2_a
    ab <- e2_b - e1_a # from e1_a to e2_b
    # why can't we just use vector addition? e1 and e2 might be parallel
    # each other and have the same length!
    # compare sum of squared lengths
    if (dot(aa,aa) + dot(ab,ab) == len2_e1) {
        return(TRUE)
    }
    return(FALSE)
}

vni.findEdges <- function(e0, edges) {
# find all matches of e0 in edges and return their indices
    matches = c()
    for (e in 1:length(edges)) {
        if (vni.sameEdge(e0,edges[[e]])) {
            matches[length(matches)+1] <- e        
        }
    }
    return(matches)
}


vni.removeEdgeMatches <- function(e0, edges) {
   matches <- vni.findEdges(e0, edges) 
   if (length(matches) > 0) {
       return(list.remove(edges,range=matches))
   }
   return(edges) 
}

# bounding points at "infinity"
vni.genBounding <- function(width, height) {
    
    P1 <- list(x=-width, y=-height)
    P2 <- list(x=2*width, y=-height)
    P3 <- list(x=2*width, y=2*height)
    P4 <- list(x=-width,y=2*height)
    
    # edges for cells of bounding points   
    E1 <- list(
        vni.makeEdge(0.5*width,0.5*height,0.5*width,-10*height),
        vni.makeEdge(0.5*width,-10*height, -10*width,0.5*height),
        vni.makeEdge(-10*width,0.5*height, 0.5*width,0.5*height)
    )
    
    E2 <- list(
        vni.makeEdge(0.5*width,0.5*height,0.5*width,-10*height),
        vni.makeEdge(0.5*width,-10*height,10*width,0.5*height),
        vni.makeEdge(10*width,0.5*height,0.5*width,0.5*height)
    )
    
    E3 <- list(
        vni.makeEdge(0.5*width,0.5*height,0.5*width,10*height),
        vni.makeEdge(0.5*width,10*height,10*width,0.5*height),
        vni.makeEdge(10*width,0.5*height,0.5*width,0.5*height)
    )
    
    E4 <- list(
        vni.makeEdge(0.5*width,0.5*height, 0.5*width,10*height),
        vni.makeEdge(0.5*width,10*height, -10*width,0.5*height),
        vni.makeEdge(-10*width,0.5*height, 0.5*width,0.5*height)
    )
    
    
    bound_points <- list(P1,P2,P3,P4)
    bound_edges <- list(E1,E2,E3,E4)
    
    for (el in bound_edges) {
        for (e in el) {
            if (!vni.validateEdge(e)) {
                stop("for edge: ",e," in list ", "el") 
            }
        }
    }
    
    bound_cells <- vni.makeCells(bound_points, bound_edges)
    return(list(points=bound_points, edges=bound_edges, cells=bound_cells))
}

vni.genCells <- function(points, geometry="L2", width=400, height=400) {
    # this procedure adapted from
    # https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html 
    Sites <- points # generator point for the Voronoi cell
    # initialize cells of bounding points "at infinity" (outside the finite space)
    bounding <- vni.genBounding(width = width, height = height) 
    Cells <- bounding$cells  # treat as immutable
    siteIndex <- 0
    status <- ""
    for (s in Sites) {
        # === debug
        # status <- "start of sites"
        # siteIndex <- siteIndex + 1
        # print(list(site=siteIndex)) # debug
        # s <- Sites[[siteIndex]] # for debug
        # create a new cell for site
        
        cell_point <- s
        cell_edges <- list() # to add edges
        CellsUpdate <- list()
        # cellIndex <- 0 # debug
        for (c in Cells) {
            # === debug ===
            # status <- "start of cells"
            # cellIndex <- cellIndex + 1
            # print(list(cell=cellIndex)) # debug
            # c <- Cells[[cellIndex]] # for debug
            
            s2 <- c$point
            # get edge connecting cell sites
            seg <- vni.makeEdge(cell_point$x, cell_point$y, s2$x, s2$y)
            # find pb, the perpendicular bisector of the points
            pb <- vni.getPB(seg, geometry, width, height)
            # create a place to store intersections with pb (critical points)
            ix <- list() # intersections
            c_edges <- list() # edges to keep in the cell
            # edgeCounter <- 0
            for (e in c$edges) {
                # debug ==
                # status <- "start of edges"
                # edgeCounter <- edgeCounter + 1
                # print(list(edge=edgeCounter)) 
                # e <- c$edges[[edgeCounter]] # debug
                
                # validate edge
                if (!vni.validateEdge(e)) {
                    print("invalid edge debug")
                    print("[e]")
                    print(e)
                    print("in cell with site:")
                    print(s2)
                    print("cell index")
                    print(length(CellsUpdate)+1)
                    stop("encountered invalid edge") 
                }
                # test spatial relationship between e and pb 
                #
                # if e intersects pb (we will get an intersection point,
                # which will be a list with $x $y)
                # status <<- "detecting intersections"
                ix_pb <- vni.ixEdges(e, pb)
                if (is.list(ix_pb) && vni.validatePoint(ix_pb)) {
                    # clip e to keep the part on the far side of pb
                    # - we can tell which endpoint is on the other side
                    #   by making a "degenerate" edge with the same endpoint
                    e1 <- e[[1]]
                    e2 <- e[[2]]
                    makeDegenEdge <- function(p) {
                        return(vni.makeEdge(p$x,p$y,p$x,p$y))    
                    }
                    e1e <- makeDegenEdge(e[[1]])
                    clippedEdge <- NULL
                    # status <- "edge intersection"
                    if (vni.sameSide(cell_point,e1e,pb)) {
                        # if e1 is on the same side as s (the site of the cell)
                        # then use e2 as the endpoint for the clipped edge
                        clippedEdge <- vni.makeEdge(e2$x,e2$y,ix_pb$x,ix_pb$y) 
                    } else {
                        # otherwise  e1 is on the other side, so use that
                        clippedEdge <- vni.makeEdge(e1$x,e1$y,ix_pb$x,ix_pb$y) 
                    }
                    # keep this edge
                    c_edges[[length(c_edges)+1]] <- clippedEdge
                    # point of intersection with pb
                    ix[[length(ix)+1]] <- ix_pb
                } else if (vni.sameSide(cell_point, e, pb)) {
                    # if e is entirely on the proximal side of pb relative to s
                    # don't save this edge in c
                    # find this edge in Edges and delete it 
                    # status <- "edge: entirely same side, need to remove"
                    # Edges <- vni.removeEdgeMatches(e, Edges)
                } else {
                    # do nothing to the edge if it is distal to s and pb
                    # keep this edge 
                    # validate edge
                    c_edges[[length(c_edges)+1]] <- e
                }
            }
            # status <- "looking at critical points"
            # should be 0 or 2 critical points
            if (length(ix) == 2) {
                # if 2, create a new edge
                e_prime <- vni.makeEdge(ix[[1]]$x,ix[[1]]$y,ix[[2]]$x,ix[[2]]$y)
                # print(list(msg="New edge for cell ", cell=cellIndex, edge_num=(length(cell_edges)+1)))
                # and add to:
                # - cell edges
                cell_edges[[length(cell_edges)+1]] <- e_prime
                # - c
                c_edges[[length(c_edges)+1]] <- e_prime 
            }
            # update c before we move to the next one
            # status <- "updating cell"
            c_update <- vni.makeCell(c$point, c_edges)
            if (!vni.validateCell(c_update)) {
                print("invalid cell after update")
                print("old cell")
                print(c)
                print("new cell")
                print(c_update)
                stop("problem with updated cell")
            }
            # add the updated cell back into the overall Cells collection
            CellsUpdate[[length(CellsUpdate)+1]] <- c_update
            # status <<- "end of cells loop"
        } # for c in Cells
        # add new cell to Cells
        CellsUpdate[[length(CellsUpdate)+1]] <- vni.makeCell(cell_point, cell_edges)
        # update Cells
        Cells <- CellsUpdate
    }    
    # now to clip all edges to the bounding rectangle

    return(Cells[5:length(Cells)])
}

# plotting functions
# ------------------

### data helpers

vni.extractLabelToList <- function(l,label='x') {
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
    
vni.extractLabel <- function(l,label='x') {
    X <- c()
    for (p in l) {
        if (length(p) > 0) {
            X <- c(X,p[[label]])
        }
    }
    return(X)
}
    
vni.plotEdges <- function(E, color="black") {
    if (length(E) > 0) {
        edge_to_debug <<- E
        startPoints <- vni.extractLabelToList(E,1)
        endPoints <- vni.extractLabelToList(E,2)
        startX <- vni.extractLabel(startPoints,"x")
        startY <- vni.extractLabel(startPoints,"y")
        endX <- vni.extractLabel(endPoints,"x")
        endY <- vni.extractLabel(endPoints,"y")
        segments(startX,startY,endX,endY, col=color)
    } else {
        print("warning - vni.plotEdges: got empty edge list")
    }
}
    
vni.plotCells <- function(cells, width, height, cex = 0.5, pch = 16) {
    extractPointsXY <- function(cells) {
        points <- vni.extractLabelToList(cells, "point")
        X <- vni.extractLabel(points,"x")
        Y <- vni.extractLabel(points,"y")
        return(list(X=X,Y=Y))
    }
    extractEdges <- function(cells) {
        edgesLists <- vni.extractLabelToList(cells, "edges")
        allEdges <- list()
        for (l in edgesLists) {
            for (e in l) {
                allEdges[[length(allEdges)+1]] <- e
            }
        }
        return(allEdges)
    }
    
    print(cells)
    
    pointsXY <- extractPointsXY(cells)
    n <- length(pointsXY$X)
    plotTitle <- paste("Voronoi Diagram,",n,"random points")
    print(paste("Plotting",plotTitle))
    
    par(cex=cex,pch=pch)
    plot(pointsXY$X,pointsXY$Y,
         main=plotTitle,
         col="red",
         xlab="x",
         ylab="y",
         xlim=c(0,width),
         ylim=c(0,height)
    )
    allEdges <- extractEdges(cells)
    vni.plotEdges(allEdges)
    return(allEdges)
} 
