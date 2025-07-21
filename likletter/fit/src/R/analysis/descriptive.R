rm(list = ls())

loc <- "likletter/"
model <- "LSSP"
dataset <- "synthetic"
prefix <- paste0(model, "_", dataset)
path_outs <- "results/likletter/"
seed <- 1

base::source(paste0(loc, "utils.R"))
Rcpp::sourceCpp(paste0(loc, "fit/src/cpp/cfunctions.cpp"))

### Data generation
### p: no. covariates
### n: no. actors
### m: no. knots (locations)
### X: covariates matrix (n x p)
### g: surface

set.seed(seed)
p <- 2
n <- 40
m <- 10 * p
mu <- -.5
X <- matrix(runif(n * p), n, p)
g <- function(x, y) 1.5 * exp(x^2) * x^2

dat <- data_gen2(n, mu, X, g)
z <- dat$z # scores
Theta <- dat$Theta # probabilities
Yadj <- dat$Y # adjacency matrix

#Y <- as.matrix(Yadj[lower.tri(Yadj)])
n_sams <- 10000

net.data <- graph_from_adjacency_matrix(Yadj, diag = FALSE, mode = 'undirected')
graph.net <- toVisNetworkData(net.data)
visNetwork(
  nodes = graph.net$nodes,
  edges = graph.net$edges,
  main = 'Simulated Network'
) %>%
  visIgraphLayout(layout = "layout_with_kk") %>%
  visNodes(font = list(size = 50)) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50)
  )

# DESCRIPTIVE ANALYSIS
Adj.obj = Yadj
ecount(net.data)
vcount(net.data) # Number of edges and vertices

#* Degree analysis
d.net = degree(net.data)
fig <- plot_ly(
  x = d.net,
  type = "histogram",
  marker = list(
    color = 'rgb(158,202,225)',
    line = list(color = 'rgb(8,48,107)', width = 1.5)
  )
) %>%
  layout(
    title = "",
    xaxis = list(title = " Grado del vértice"),
    yaxis = list(title = "Frecuencia"),
    font = list(size = 14)
  )
fig

a.nn.deg <- knn(net.data, V(net.data))$knn
data = data.frame(x = d.net, y = a.nn.deg) %>% na.omit()
fig <- plot_ly(
  data = data,
  x = ~x,
  y = ~y,
  marker = list(
    size = 10,
    color = 'rgba(255, 182, 193, .9)',
    line = list(color = 'rgba(152, 0, 0, .8)', width = 2)
  )
) %>%
  layout(
    xaxis = list(title = " Grado del vértice"),
    yaxis = list(title = "Grado Promedio de los Vecinos"),
    font = list(size = 14)
  )
fig

# Degree
net.data. <- net.data
colores <- rep('#97C2FC', vcount(net.data.))
colores[top(degree(net.data.))] = c('tomato', 'orange', 'yellow')
V(net.data.)$color = colores

graph.net <- toVisNetworkData(net.data.)
nodes = graph.net$nodes
edges = graph.net$edges
nodes$label[-which(nodes$label %in% c('USA', 'IRQ', 'JOR'))] = NA

visNetwork(nodes, edges, main = 'Grado') %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visNodes(font = list(size = 80)) %>%
  visEdges(color = list(opacity = 0.4))

# Closeness
net.data. <- decompose(net.data)[[1]]
colores <- rep('#97C2FC', vcount(net.data.))
colores[top(closeness(net.data.))] = c('tomato', 'orange', 'yellow')
V(net.data.)$color = colores

graph.net <- toVisNetworkData(net.data.)
nodes = graph.net$nodes
edges = graph.net$edges
nodes$label[-which(nodes$label %in% c('USA', 'IRQ', 'JOR'))] = NA

visNetwork(nodes, edges, main = 'Closeness') %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visNodes(font = list(size = 80)) %>%
  visEdges(color = list(opacity = 0.4))

# Eigenvalue
net.data. <- net.data
colores <- rep('#97C2FC', vcount(net.data.))
colores[top(eigen_centrality(net.data.)$vector)] = c(
  'tomato',
  'orange',
  'yellow'
)
V(net.data.)$color = colores

graph.net <- toVisNetworkData(net.data.)
nodes = graph.net$nodes
edges = graph.net$edges
nodes$label[-which(nodes$label %in% c('USA', 'IRQ', 'JOR'))] = NA

visNetwork(nodes, edges, main = 'Eigen Value') %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visNodes(font = list(size = 80)) %>%
  visEdges(color = list(opacity = 0.4))

# Betweenness
net.data. <- net.data
colores <- rep('#97C2FC', vcount(net.data.))
colores[top(betweenness(net.data.))] = c('tomato', 'orange', 'yellow')
V(net.data.)$color = colores

graph.net <- toVisNetworkData(net.data.)
nodes = graph.net$nodes
edges = graph.net$edges
nodes$label[-which(nodes$label %in% c('USA', 'IRQ', 'JOR'))] = NA

visNetwork(nodes, edges, main = 'Betweenness') %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visNodes(font = list(size = 80)) %>%
  visEdges(color = list(opacity = 0.4))

#* LARGEST COMPONENT
net.data.gc <- decompose(net.data)[[1]]

graph.net.gc <- toVisNetworkData(net.data.gc)
visNetwork(
  nodes = graph.net.gc$nodes,
  edges = graph.net.gc$edges,
  main = 'Largest component of simulated data'
) %>%
  visIgraphLayout(layout = "layout_with_kk") %>%
  visNodes(font = list(size = 50)) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -100)
  )

#* Some measures
# Total graph
mean_distance(net.data)
diameter(net.data)
transitivity(net.data)
vertex_connectivity(net.data)
edge_density(net.data)

# Largest component
mean_distance(net.data.gc)
diameter(net.data.gc)
transitivity(net.data.gc)
vertex_connectivity(net.data.gc)
edge_density(net.data.gc)

(net.data.cut.vertices <- articulation_points(net.data.gc))
length(net.data.cut.vertices)

#* Largest component partition
kc <- cluster_fast_greedy(net.data.gc)
length(kc)
sizes(kc)
group = membership(kc)
cluster = data.frame(id = as.numeric(V(net.data.gc)), group = as.numeric(group))
plot(kc, net.data.gc)

graph.net.gc <- toVisNetworkData(net.data.gc)
nodes = graph.net.gc$nodes
edges = graph.net.gc$edges
nodes = merge(nodes, cluster)
visNetwork(nodes, edges) %>%
  visNodes(font = list(size = 80)) %>%
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -100)
  )

#* Assortativity analysis
assortativity_degree(net.data)
assortativity_degree(net.data.gc)
