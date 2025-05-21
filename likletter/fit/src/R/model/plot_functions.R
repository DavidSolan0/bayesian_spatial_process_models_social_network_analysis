adjacency.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
        JJ <- dim(mat)[1]
        colorscale <- c("white", rev(heat.colors(100)))
        if(missing(labs))     labs <- 1:JJ
        if(missing(col.axis)) col.axis <- rep("black", JJ)
        if(missing(cex.axis)) cex.axis <- 0.5
        if(missing(tick))     tick <- TRUE
        ## adjacency matrix
        image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
              col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
        for(j in 1:JJ){
                axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
                axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
        }
        box()
        if(show.grid) grid(nx=JJ, ny=JJ)
}

heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis,...)
{ 
        JJ <- dim(mat)[1]
        colorscale <- c("white", rev(heat.colors(100)))
        if(missing(labs))     labs <- 1:JJ
        if(missing(col.axis)) col.axis <- rep("black", JJ)
        if(missing(cex.axis)) cex.axis <- 0.5
        if(missing(tick))     tick <- TRUE
        ## adjacency matrix
        image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
              col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
        for(j in 1:JJ){
                axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
                axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
        }
        box()
        if(show.grid) grid(nx = JJ, ny = JJ)
}

chain_plot <- function(x, main, truth) 
{
        windows(width = 7, height = 3.5)
        par(mfrow = c(1,1), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
        plot(x = 1:length(x), y = c(x), type = "l", cex.axis = 0.6, col = "gray", xlab = "", ylab = "", main = main)
        abline(h = mean(x), lty = 2, col = 1)
        abline(h = quantile(x, c(.025,.975)), lty = 2, col = 4)
        if (!missing(truth)) abline(h = truth, lty = 1, col = 2, lwd = 2)
}

hist_plot <- function(x, main, truth) 
{ 
        windows(width = 3.5, height = 3.5)
        hist(x, freq = F, nclass = 20, col = "gray95", cex.axis = 0.6, border = "gray", xlab = "", ylab = "", main = main)
        # lines(density(x, adjust = 1.5))
        abline(v = mean(x), lty = 2, col = 1)
        abline(v = quantile(x, c(.025,.975)), lty = 2, col = 4)
        if (!missing(truth)) abline(v = truth, lty = 1, col = 2, lwd = 2)
}

data_gen2 <- function(n, mu, X, g, plotit = T,seed=NULL,
                      font_size=40,main='Graph',layout='layout_with_kk',
                      gravitation=150) 
{
        ### simulates a data set following a LSSP model as in Linkletter 2007
        ### p = 2 covariates
        # LSSP scores
        z <- NULL
        for (i in 1:n)
                z[i] <- g(X[i,1], X[i,2])
        # network simulation
        Y <- Theta <- matrix(0, n , n) 
        for (i in 1:(n-1)) {
                for (j in (i+1):n) {
                        Theta[i,j] <- Theta[j,i] <- exp(pnorm(mu - abs(z[i] - z[j]), log.p = T))
                        Y[i,j] <- Y[j,i] <- rbinom(1, 1, Theta[i,j])
                }
        }
        if (plotit) {
                # graph & adjacency
                suppressMessages(suppressWarnings(library(igraph)))
                suppressMessages(suppressWarnings(library(visNetwork)))
                suppressMessages(suppressWarnings(library(plotly)))
                windows(width = 10.5, height = 3.5)
                par(mfrow = c(1,3), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
                heat.plot0(Theta, main = "True probabilities")
                adjacency.plot(Y , main = "Adjacency matrix")
                plot(x = graph.adjacency(adjmatrix = Y, mode = "undirected", diag = F), 
                     vertex.color = "black", edge.color = "firebrick1", vertex.label = NA, 
                     vertex.size = 4.5, layout = layout_with_kk, main = "Graph")
                
                net.data <- graph_from_adjacency_matrix(Y,diag = F,mode='undirected')
                graph.net <- toVisNetworkData(net.data)
                nodes = graph.net$nodes; edges =  graph.net$edges
                graph = visNetwork(nodes, edges, main = main) %>%
                        visIgraphLayout(layout = layout) %>%
                        visNodes(font = list(size = font_size)) %>%
                        visIgraphLayout(randomSeed = seed) %>%
                        visPhysics(solver = "forceAtlas2Based", 
                                   forceAtlas2Based = list(gravitationalConstant = -gravitation))
                
                
                # surface
                x <- seq(from = 0, to = 1, length = 20)
                y <- seq(from = 0, to = 1, length = 20)
                f <- outer(x, y, g)
                nrz <- nrow(f)
                ncz <- ncol(f)
                nbcol <- 100
                jet.colors <- colorRampPalette(c("blue", "green"))
                color <- jet.colors(nbcol)
                zfacet <- (f[-1, -1] + f[-1, -ncz] + f[-nrz, -1] + f[-nrz, -ncz])/4
                windows(width = 7.5, height = 5)
                par(mfrow = c(2,3), mar = c(3,3,2,1), mgp = c(1.6,0.6,0), oma = c(0,0,0,0))
                for (theta in c(30,60,120,150,210,240))
                        persp(x = x, y = y, z = f, theta = theta, phi = 30, col = color[cut(zfacet, nbcol)], 
                              axes = T, ticktype = "simple", cex.axis = 0.5, zlab = "z")
                title(main = "True Surface", outer = T, line = -1)
                
                fig <- plot_ly(x =~ x, y =~ y, z =~ f) %>% add_surface(
                        contours = list(
                                z = list(
                                        show=TRUE,
                                        usecolormap=TRUE,
                                        highlightcolor="#ff0000",
                                        project=list(z=TRUE))))
                
                fig <- fig %>% layout(
                        scene = list(
                                camera=list(
                                        eye = list(x=1.87, y=0.88, z=-0.64))))
                
                
        }
        # return
        list(z = z, Y = Y, Theta = Theta,graph=graph,surface = fig)
}
