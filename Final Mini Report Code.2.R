'5. Application to Bacilus subtilis gene expression data'



##########################################################################


'Load and manipulate data'


set.seed(1)


load('Documents/Year 4/High Dimensional Statistics Practicals/riboflavin.data.rda')
# Load our data


riboflavin <- data
# Store our data


dim(riboflavin)
# Consists of 71 observations and 40088 genes


riboflavin <- riboflavin[,1:50]
# We transform our data to only consist of the first 50 genes


flavin <- scale(riboflavin)
# Standardise the data


S <- cov(flavin)
# Store the empirical covariance matrix



##########################################################################


'Perform Graphical Lasso'


install.packages('glasso')
library(glasso)


rho <- c(0.01,0.1,0.2)
# Multiple values of rho stored


flavin_glasso_1 <- glasso::glasso(S,rho[1])
# Glasso with rho = 0.01


flavin_glasso_2 <- glasso::glasso(S,rho[2])
# Glasso with rho = 0.1


flavin_glasso_3 <- glasso::glasso(S,rho[3])
# Glasso with rho = 0.2


flavin_inv_1 <- flavin_glasso_1$wi
# Incovariance matrix for rho = 0.01


flavin_inv_2 <- flavin_glasso_2$wi
# Incovariance matrix for rho = 0.1


flavin_inv_3 <- flavin_glasso_3$wi
# Incovariance matrix for rho = 0.2


colnames(flavin_inv_1) <- rownames(flavin_inv_1) <- colnames(flavin)
colnames(flavin_inv_2) <- rownames(flavin_inv_2) <- colnames(flavin)
colnames(flavin_inv_3) <- rownames(flavin_inv_3) <- colnames(flavin)
# Inserting column and row names 



##########################################################################


'Heat Maps'


colour_key <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
# Define colour key


# For rho = 0.1


flavin_vis_2 <- flavin_inv_2- diag(diag(flavin_inv_2))


flavin_vis_2 <- abs(flavin_vis_2) / max(abs(flavin_vis_2)) # mapping into [0,1]
diag(flavin_vis_2) <- -1 # to make the diagonal stand out visually


flavin_vis_2_long <- melt(flavin_vis_2)


ggplot(flavin_vis_2_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colour_key,
    limits = c(-1, 1),
    name = " Colour Key"
  )  + labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )


# Appendix #


# For rho = 0.01


flavin_vis_1 <- flavin_inv_1- diag(diag(flavin_inv_1))


flavin_vis_1 <- abs(flavin_vis_1) / max(abs(flavin_vis_1)) # mapping into [0,1]
diag(flavin_vis_1) <- -1 # to make the diagonal stand out visually


flavin_vis_1_long <- melt(flavin_vis_1)


ggplot(flavin_vis_1_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colour_key,
    limits = c(-1, 1),
    name = " Colour Key"
  )  + labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

# For rho = 0.2


flavin_vis_3 <- flavin_inv_3- diag(diag(flavin_inv_3))


flavin_vis_3 <- abs(flavin_vis_3) / max(abs(flavin_vis_3)) # mapping into [0,1]
diag(flavin_vis_3) <- -1 # to make the diagonal stand out visually


flavin_vis_3_long <- melt(flavin_vis_3)


ggplot(flavin_vis_3_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colour_key,
    limits = c(-1, 1),
    name = " Colour Key"
  )  + labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )



##########################################################################


'Undirected Graphs'


install.packages('igraph')
library(igraph)
# Load igraph


flavin_graph_1 <- flavin_inv_1[5:15,5:15]
# Only considering 10 genes for rho = 0.01



flavin_graph_2 <- flavin_inv_2[5:15,5:15]
# Only considering 10 genes for rho = 0.1


flavin_graph_3 <- flavin_inv_3[5:15,5:15]
# Only considering 10 genes for rho = 0.2


igraph_options(vertex.size = 25,
               vertex.color = "orange",
               vertex.frame.color = 'black',
               vertex.label.cex = 0.9,
               vertex.label.color = "black",
               edge.width = 1,
               edge.color = "darkgrey")
# Setting igraph options


dev.off()
par(mfrow = c(1,3))
# Adjust format 


# For rho = 0.01


network_initial_1 <- graph.adjacency(abs(flavin_graph_1), weighted=TRUE,
                                     mode="undirected", diag=FALSE)
# Undirected graph


network_layout_1 <- layout_in_circle(network_initial_1)
# Specifies where the nodes go


plot.igraph(network_initial_1, layout=network_layout_1)
# Plots final graph


# For rho = 0.1


network_initial_2 <- graph.adjacency(abs(flavin_graph_2), weighted=TRUE,
                                     mode="undirected", diag=FALSE)
# Undirected graph


network_layout_2 <- layout_in_circle(network_initial_2)
# Specifies where the nodes go


plot.igraph(network_initial_2, layout=network_layout_2)
# Plots final graph


# For rho = 0.2


network_initial_3 <- graph.adjacency(abs(flavin_graph_3), weighted=TRUE,
                                     mode="undirected", diag=FALSE)
# Undirected graph


network_layout_3 <- layout_in_circle(network_initial_3)
# Specifies where the nodes go


plot.igraph(network_initial_3, layout=network_layout_3)
# Plots final graph


##########################################################################


'CVglasso'


install.packages('CVglasso')
library(CVglasso)
# Load packages


lam_seq <- seq(0.01, 0.5, length.out = 100)
# Specify lamda sequence we want it to analyse

CV_glasso <- CVglasso::CVglasso(X = flavin, S = S,
                                K = 10, crit.cv = 'loglik',
                                trace = 'none',
                                lam.min.ratio = 1e-5,
                                lam = lam_seq)
# Perform 10 fold cross validation


CV_glasso$Tuning[2]
# Optimal Lambda we found is 0.1436364


plot(CV_glasso, xlab = 'log10(rho)')
# Plot of Cross-Validation Errors


glasso_CV <- glasso::glasso(S, rho = CV_glasso$Tuning[2])
# Perform glasso


glasso_inv_CV <- glasso_CV$wi


###
'Heat map'


colour_key <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
# Define colour key


colnames(glasso_inv_CV) <- rownames(glasso_inv_CV) <- colnames(flavin)

glasso_vis_CV <- glasso_inv_CV- diag(diag(glasso_inv_CV))


glasso_vis_CV <- abs(glasso_vis_CV) / max(abs(glasso_vis_CV)) # mapping into [0,1]
diag(glasso_vis_CV) <- -1 # to make the diagonal stand out visually


glasso_vis_CV_long <- melt(glasso_vis_CV)


ggplot(glasso_vis_CV_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colour_key,
    limits = c(-1, 1),
    name = " Colour Key"
  )  + labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )
# Heatmap for optimal Rho


###
'Undirected graph'


igraph_options(vertex.size = 25,
               vertex.color = "orange",
               vertex.frame.color = 'black',
               vertex.label.cex = 0.9,
               vertex.label.color = "black",
               edge.width = 1,
               edge.color = "darkgrey")
# Setting igraph options


glasso_CV_graph <- glasso_inv_CV[5:15,5:15]


network_initial <- graph.adjacency(abs(glasso_CV_graph), weighted=TRUE,
                                     mode="undirected", diag=FALSE)
# Undirected graph


network_layout <- layout_in_circle(network_initial)
# Specifies where the nodes go


plot.igraph(network_initial, layout=network_layout)
# Plots final graph

