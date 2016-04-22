###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
###########################################
load("data/functions.RData")

# real data example 2 : facebook
########### Facebook with Phi ############
circle <- read.csv("data/0.circles", sep ="\t", header=FALSE)
edges <- read.csv("data/0.edges", sep =" ", header=FALSE) 
egofeat <- read.csv("data/0.egofeat", sep = " ", header=FALSE)
feat <- read.csv("data/0.feat", sep = " ", header=FALSE)
feat <- feat[,-1]

#77 gender;anonymized feature 77
#78 gender;anonymized feature 78

gender1 <- feat[,78]
gender2 <- feat[,79]


## Generate Network(this network includes six nodes)
n <- 347
A <- matrix(0, nrow = n, ncol = n)

for(i in 1:nrow(edges)){
  e1 <- edges[i,1]
  e2 <- edges[i,2]
  A[e1, e2] <- 1
  A[e2, e1] <- 1
}

fbnetwork <- graph.adjacency(A, mode = c("undirected"),
                             weighted = NULL)

fbnetwork1<-induced_subgraph(fbnetwork
                             , V(fbnetwork)[clusters(fbnetwork)$membership == 1])


### outcome matrix ### (does not contain six nodes)
V(fbnetwork1)$gen <- gender1[clusters(fbnetwork)$membership == 1]
V(fbnetwork1)$gen2 <- gender2[clusters(fbnetwork)$membership == 1]
missing <- ifelse(V(fbnetwork1)$gen + V(fbnetwork1)$gen2 == 1, 0, 1)

outcome <- V(fbnetwork1)$gen[-which(missing == 1)]

net <- fbnetwork1

average.path.length(net)

np <- 500 ; ndelta = 5

overall_prob <- rep(NA, 2)

for(i in 1:length(overall_prob)){
  overall_prob[i] <- sum(outcome == i) / length(outcome) 
}

A_hat <- get.adjacency(net)
A_hat <- as.matrix(A_hat)
which(missing == 1)
A_hat <- A_hat[-which(missing==1) , -which(missing==1)] # 318 * 318

np <- 500
fb_newPhi <- rep(NA, np)
fb_moran <- rep(NA, np)
outcome <- outcome + 1

fb_moran <- MoranI(A_hat, outcome)[4]

s0 <- sum(A_hat + t(A_hat)) / 2
s1 <- sum((A_hat + t(A_hat))^2) / 2
s2 <- sum((apply(A_hat,1,sum) + apply(A_hat, 2, sum))^2)

fb_matrix <- matrix(NA, ncol = 3, nrow = 500)
fb_matrix[1,] <- weight.matrix(A_hat)


fb_newPhi <- newphi(A_hat, outcome)
fb_newPhi <- fb_newPhi - mean.newphi(s0, outcome)
fb_newPhi <- fb_newPhi / sqrt(m2.newphi(s0, s1, s2, outcome) - (mean.newphi(s0,outcome))^2 )


B <- get.adjacency(fbnetwork1)


##################### different distance matrix ###########
## adjacency matrix by distance(pairwise)
distance.matrix <- shortest.paths(fbnetwork1, v = V(fbnetwork1), to= V(fbnetwork1))
distance.matrix <- distance.matrix[-which(missing == 1), -which(missing == 1)]

## histogram of distance matrix
pdf("figures/fbdist.pdf")
par(mfrow = c(1,1))
hist(distance.matrix, breaks = rep(0:11, each = 2) + c(-0.4, 0.4), 
     main = "Distribution of distance", 
     xlab = "shortest length between two users in Facebook",
     ylab = "frequency", col = "dodgerblue",
     xaxt = "n")
axis(1, at = 0:10, labels = c(0:10))
dev.off()

weight.matrix2 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
weight.matrix3 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
weight.matrix4 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
weight.matrix5 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
weight.matrix6 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
weight.matrix7 <- matrix(0, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
fb_reci1 <- rep(0, 4)
fb_reci2 <- rep(0, 4)
fb_reci3 <- rep(0, 4)
fb_reci4 <- rep(0, 4)
fb_reci5 <- rep(0, 4)
fb_reci6 <- rep(0, 4)
fb_reci7 <- rep(0, 4)


for(k in 1:2){
  weight.matrix2 <- ifelse(distance.matrix == k, 1/k, weight.matrix2)
}

for(k in 1:3){
  weight.matrix3 <- ifelse(distance.matrix == k, 1/k, weight.matrix3)
}

for(k in 1:4){
  weight.matrix4 <- ifelse(distance.matrix == k, 1/k, weight.matrix4)
}

for(k in 1:5){
  weight.matrix5 <- ifelse(distance.matrix == k, 1/k, weight.matrix5)
}

for(k in 1:6){
  weight.matrix6 <- ifelse(distance.matrix == k, 1/k, weight.matrix6)
}

for(k in 1:7){
  weight.matrix7 <- ifelse(distance.matrix == k, 1/k, weight.matrix7)
}

fb_reci1 <- MoranI(A_hat, outcome)
fb_reci2 <- MoranI(weight.matrix2, outcome)
fb_reci3 <- MoranI(weight.matrix3, outcome)
fb_reci4 <- MoranI(weight.matrix4, outcome)
fb_reci5 <- MoranI(weight.matrix5, outcome)
fb_reci6 <- MoranI(weight.matrix6, outcome)
fb_reci7 <- MoranI(weight.matrix7, outcome)

distance.matrix12 <- ifelse(distance.matrix <= 2 & distance.matrix >0,1,0)
distance.matrix13 <- ifelse(distance.matrix <= 3 & distance.matrix >0,1,0)
distance.matrix14 <- ifelse(distance.matrix <= 4 & distance.matrix >0,1,0)
distance.matrix15 <- ifelse(distance.matrix <= 5 & distance.matrix >0,1,0)
distance.matrix16 <- ifelse(distance.matrix <= 6 & distance.matrix >0,1,0)
distance.matrix17 <- ifelse(distance.matrix <= 7 & distance.matrix >0,1,0)

fb_moran1 <- MoranI(A_hat, outcome)
fb_moran12 <- MoranI(distance.matrix12, outcome)
fb_moran13 <- MoranI(distance.matrix13, outcome)
fb_moran14 <- MoranI(distance.matrix14, outcome)
fb_moran15 <- MoranI(distance.matrix15, outcome)
fb_moran16 <- MoranI(distance.matrix16, outcome)
fb_moran17 <- MoranI(distance.matrix17, outcome)



########### Tables ###############
fb_summary <- matrix(NA, ncol = 4, nrow = 7)
colnames(fb_summary) <- c("Distance", "number of neighbors", "Moran's I", "p-value")
fb_summary[,1] <- c("<=1", "<=2", "<=3",
                    "<=4", "<=5", "<=6", "<=7")
fb_summary[,2] <- c(sum(distance.matrix), sum(distance.matrix12),
                    sum(distance.matrix13), sum(distance.matrix14),
                    sum(distance.matrix15), sum(distance.matrix16),
                    sum(distance.matrix17))
fb_summary[,2] <- as.integer(fb_summary[,2])
fb_summary[,3] <- c(fb_moran1[4], fb_moran12[4], fb_moran13[4],
                    fb_moran14[4], fb_moran15[4], fb_moran16[4],
                    fb_moran17[4])
fb_summary[,4] <- c(1-pnorm(fb_moran1[4]), 1-pnorm(fb_moran12[4]),
                    1-pnorm(fb_moran13[4]), 1-pnorm(fb_moran14[4]), 
                    1-pnorm(fb_moran15[4]), 1-pnorm(fb_moran16[4]),
                    1-pnorm(fb_moran17[4]))
fb_summary <- as.data.frame(fb_summary)
# export a table
print(xtable(fb_summary, digits = 2, row.names = FALSE))

#########################
fb_weight <- matrix(NA, ncol = 4, nrow = 7)
colnames(fb_weight) <- c("Distance", "number of neighbors", "Moran's I", "p-value")
fb_weight[,1] <- c("<=1", "<=2", "<=3",
                   "<=4", "<=5", "<=6", "<=7")
fb_weight[,2] <- c(sum(distance.matrix1), sum(distance.matrix12),
                   sum(distance.matrix13), sum(distance.matrix14),
                   sum(distance.matrix15), sum(distance.matrix16),
                   sum(distance.matrix17))
fb_weight[,3] <- c(fb_moran1[4], fb_reci2[4], fb_reci3[4],
                   fb_reci4[4], fb_reci5[4], fb_reci6[4],
                   fb_reci7[4])
fb_weight[,4] <- c(1-pnorm(fb_moran1[4]), 1-pnorm(fb_reci2[4]),
                   1-pnorm(fb_reci3[4]), 1-pnorm(fb_reci4[4]), 
                   1-pnorm(fb_reci5[4]), 1-pnorm(fb_reci6[4]),
                   1-pnorm(fb_reci7[4]))
fb_weight <- as.data.frame(fb_weight)
print(xtable(fb_weight, digits = 2, row.names = FALSE))


## join-count analysis only for delta = 1
install.packages("spdep")
library(spdep)

AA <- mat2listw(A_hat)
AA <- AA$neighbours
fbjoin <- joincount.multi(as.factor(outcome), nb2listw(AA, style = "B"))
# just mention it
fbjoin <- as.data.frame(fbjoin)
print(xtable(fbjoin, digits = 4, row.names = FALSE))



############ picture of graph ###################

# igraph with 318 nodes
fbnetwork2 <- graph.adjacency(A_hat, mode = "undirected")
V(fbnetwork2)$gen <- outcome

# visualizing network data

pdf("../WriteUp/figures/fbgraph.pdf")
par(mfrow = c(1,1))
igraph.options(vertex.size = 3, vertex.label = NA, 
               edge.arrow.size = 0.5)
V(fbnetwork1)$color <- "grey"
V(fbnetwork1)[V(fbnetwork1)$gen == 1]$color <- "dodgerblue"
V(fbnetwork1)[V(fbnetwork1)$gen2 == 1]$color <- "hotpink"
plot(fbnetwork1, layout = layout.fruchterman.reingold)
legend("topleft", c("gender 1", "gender 2", "unknown"),
       col = c("dodgerblue", "hotpink", "grey"), pch = 19 )

dev.off()





