###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
library(foreign)
###########################################
load("data/functions.RData")

household <- read.dta("data/indian/Data/Demographics_and_Outcomes/household_characteristics.dta")
for(i in 1:12){
  assign(paste0("participation", i), household[household$village == i, ]$hhSurveyed)
}
for(i in 14:21){
  assign(paste0("participation", i), household[household$village == i, ]$hhSurveyed)
}
for(i in 23:77){
  assign(paste0("participation", i), household[household$village == i, ]$hhSurveyed)
}

a<-c()
for(i in 1:12){
  a[i] <- sum(get(paste0("participation",i)) == 1)/ length(get(paste0("participation",i)))
}

for(i in 14:21){
  a[i] <- sum(get(paste0("participation",i)) == 1)/ length(get(paste0("participation",i)))
}

for(i in 23:77){
  a[i] <- sum(get(paste0("participation",i)) == 1)/ length(get(paste0("participation",i)))
}


## participation rates
pdf("figures/partirate.pdf")
hist(a[!is.na(a)], main = "participation rate for 75 villages",
     breaks = seq(0.3, 0.6, 0.01), col = "orange", 
     xlab = "participation rates")
dev.off()

## only include households involved in microfinance study
dd<- c(1,2,3,4,6,9,12,15,19,20,21,
       23,24,25, 29, 31,32,33,36,39,
       42,43,45,46,47,48,50,51,52,55,57,
       59,60,62,64,65,67,68,70,71,72,73,75)



### adjacency matrix ###
for(i in 1:12){
  dta <- read.csv(paste0("data/indian/Data/Network_Data/Adjacency_Matrices/adj_allVillageRelationships_HH_vilno_", i, ".csv" ), sep = ",",
                  header = FALSE)
  assign(paste0("adj", i), dta)
}

for(i in 14:21){
  dta <- read.csv(paste0("data/indian/Data/Network_Data/Adjacency_Matrices/adj_allVillageRelationships_HH_vilno_", i, ".csv" ), sep = ",",
                  header = FALSE)
  assign(paste0("adj", i), dta)
}

for(i in 23:77){
  dta <- read.csv(paste0("data/indian/Data/Network_Data/Adjacency_Matrices/adj_allVillageRelationships_HH_vilno_", i, ".csv" ), sep = ",",
                  header = FALSE)
  assign(paste0("adj", i), dta)
}

## Moran's I statistics 
parti.moran <- rep(NA, 77)

for(i in 1:12){
  parti.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("participation",i)) )[4]
}

for(i in 14:21){
  parti.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("participation",i)) )[4]
}

for(i in 23:77){
  parti.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("participation",i)) )[4]
}

pdf("figures/include.pdf")
par(mfrow= c(1,1))
hist(parti.moran[dd], main = "Microfinance network of 43 included villages",
     xlab = "Moran's I", col = "orange",
     breaks = seq(-6, 15, 0.5))
dev.off()


net50 <- graph.adjacency(as.matrix(adj50))
is.connected(net50) # it is not connected
V(net50)$participation <- participation50
V(net50)$leader <- household[household$village == 50, ]$leader

igraph.options(vertex.size = 5, vertex.label = NA, 
               edge.arrow.size = 0.1)
V(net50)[leader == 1]$shape <- "square"
V(net50)[leader == 0]$shape <- "circle"
V(net50)[participation == 1]$color <- "firebrick1"
V(net50)[participation == 0]$color <- "lightgreen"

subadj50 <- adj50[clusters(net50)$membership == 1, clusters(net50)$membership == 1]
subnet50 <- graph.adjacency(as.matrix(subadj50), mode = "undirected")
V(subnet50)$participation <- participation50[clusters(net50)$membership == 1] 


### snowball sampling 
snow <- c()
starter <- sample(V(subnet50)$name, 1)
current <- c()
current[1] <- starter
count <- 1
snow[1] <- current[1]
samn <- 100
set.seed(123)
while(count < samn){
  nnode <- length(current)
  for(i in 1:nnode){
    
    ngh <- neighbors(subnet50, current[i])
    snow <- c(snow, V(subnet50)$name[ngh])
    snow <- unique(snow)
    
  }
  tmp_sample <- snow[(count+1):length(snow)]
  if ( samn < length(snow)){
    need <- samn - count
    tmp_sample <- sample(tmp_sample, need)
    snow[(count+1):samn] <- tmp_sample
    snow <- snow[-c((samn + 1): length(snow))] 
  }
  current <- tmp_sample
  count <- length(snow)
}

# exclude 'snow' nodes from a network
snowsubnet50 <- subnet50 - vertices(V(subnet50)$name %in% snow) 
V(snowsubnet50)$participation <- V(subnet50)$participation[-which(V(subnet50)$name %in% snow)] 

V(snowsubnet50)[participation == 1]$color <- "firebrick1"
V(snowsubnet50)[participation == 0]$color <- "lightgreen"



####### do the same for random network !!
net25 <- graph.adjacency(as.matrix(adj25))
is.connected(net25) # it is not connected
V(net25)$participation <- participation25

subadj25 <- adj25[clusters(net25)$membership == 1, clusters(net25)$membership == 1]
subnet25 <- graph.adjacency(as.matrix(subadj25), mode = "undirected")
V(subnet25)$participation <- participation25[clusters(net25)$membership == 1] 


### snowball sampling 
snow <- c()
starter <- sample(V(subnet25)$name, 1)
current <- c()
current[1] <- starter
count <- 1
snow[1] <- current[1]
samn <- 100
set.seed(123)
while(count < samn){
  nnode <- length(current)
  for(i in 1:nnode){
    
    ngh <- neighbors(subnet25, current[i])
    snow <- c(snow, V(subnet25)$name[ngh])
    snow <- unique(snow)
    
  }
  tmp_sample <- snow[(count+1):length(snow)]
  if ( samn < length(snow)){
    need <- samn - count
    tmp_sample <- sample(tmp_sample, need)
    snow[(count+1):samn] <- tmp_sample
    snow <- snow[-c((samn + 1): length(snow))] 
  }
  current <- tmp_sample
  count <- length(snow)
}

# exclude 'snow' nodes from a network
snowsubnet25 <- subnet25 - vertices(V(subnet25)$name %in% snow) 
V(snowsubnet25)$participation <- V(subnet25)$participation[-which(V(subnet25)$name %in% snow)] 

V(snowsubnet25)[participation == 1]$color <- "firebrick1"
V(snowsubnet25)[participation == 0]$color <- "lightgreen"


pdf("figures/indiannetwork_max.pdf")
igraph.options(vertex.size = 7, edge.arrow.size = 0.1,
               vertex.label = NULL)
par(mfrow = c(1,1))
plot(decompose(snowsubnet50)[[1]],
     layout = layout.fruchterman.reingold,
     vertex.label = "",
     main = "Social Network of Village 50")
dev.off()

pdf("figures/indiannetwork_random.pdf")
igraph.options(vertex.size = 7, edge.arrow.size = 0.1,
               vertex.label = NULL)
par(mfrow = c(1,1))
plot(decompose(snowsubnet25)[[2]],
     layout = layout.fruchterman.reingold,
     vertex.label = "",
     main = "Social Network of Village 25")
dev.off()


#################################################################
############ Room number ################
############# (household level) room number and networks ######
for(i in 1:12){
  assign(paste0("room", i), household[household$village == i, ]$room_no)
}
for(i in 14:21){
  assign(paste0("room", i), household[household$village == i, ]$room_no)
}
for(i in 23:77){
  assign(paste0("room", i), household[household$village == i, ]$room_no)
}

room.moran <- rep(NA, 77)

for(i in 1:12){
  room.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("room",i)) )[4]
}

for(i in 14:21){
  room.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("room",i)) )[4]
}

for(i in 23:77){
  room.moran[i] <- MoranI(as.matrix(get(paste0("adj", i))), get(paste0("room",i)) )[4]
}


pdf("figures/roommoran.pdf")
par(mfrow = c(1,1))
hist(room.moran[!is.na(room.moran)], 
     col = "skyblue", 
     breaks =seq(-6, 15, 0.5), main = "Bed number in 75 villages",
     xlab = "Standardized Moran's I")
dev.off()



########## frequency table room 25
################### Network Graph #####################
net25 <- graph.adjacency(as.matrix(adj25), mode = "undirected")
subadj25 <- adj25[clusters(net25)$membership == 1, clusters(net25)$membership == 1 ]
subnet25 <- graph.adjacency(as.matrix(subadj25), mode = "undirected")
V(subnet25)$room <- room25[clusters(net25)$membership == 1]

V(subnet25)$color[V(subnet25)$room == 1]  <- "lightblue1"
V(subnet25)$color[V(subnet25)$room == 2]  <- "skyblue2"
V(subnet25)$color[V(subnet25)$room == 3]  <- "skyblue3"
V(subnet25)$color[V(subnet25)$room == 4]  <- "skyblue4"
V(subnet25)$color[V(subnet25)$room == 5]  <- "steelblue"
V(subnet25)$color[V(subnet25)$room == 6]  <- "steelblue3"
V(subnet25)$color[V(subnet25)$room == 7]  <- "steelblue4"
V(subnet25)$color[V(subnet25)$room == 8]  <- "royalblue3"
V(subnet25)$color[V(subnet25)$room == 9]  <- "royalblue4"
V(subnet25)$color[V(subnet25)$room == 10]  <- "mediumblue"
V(subnet25)$color[V(subnet25)$room == 19]  <- "navyblue"


igraph.options(vertex.size = 5, edge.arrow.size = 0.1,
               vertex.label = NULL)
pdf("figures/room25.pdf")
par(mfrow = c(1,1))
plot(subnet25, layout = layout.fruchterman.reingold,
     vertex.label = "", main = "Social Network in Village 25")
dev.off()



############ subnet25
snow <- c()
starter <- sample(V(subnet25)$name, 1)
current <- c()
current[1] <- starter
count <- 1
snow[1] <- current[1]
samn <- 130
while(count < samn){
  nnode <- length(current)
  for(i in 1:nnode){
    
    ngh <- neighbors(subnet25, current[i])
    snow <- c(snow, V(subnet25)$name[ngh])
    snow <- unique(snow)
    
  }
  tmp_sample <- snow[(count+1):length(snow)]
  if ( samn < length(snow)){
    need <- samn - count
    tmp_sample <- sample(tmp_sample, need)
    snow[(count+1):samn] <- tmp_sample
    snow <- snow[-c((samn + 1): length(snow))] 
  }
  current <- tmp_sample
  count <- length(snow)
}

# exclude 'snow' nodes from a network
snowsubnet25 <- subnet25 - vertices(V(subnet25)$name %in% snow) 
V(snowsubnet25)$room <- V(subnet25)$room[-which(V(subnet25)$name %in% snow)] 

V(snowsubnet25)$color[V(snowsubnet25)$room == 1]  <- "lightblue1"
V(snowsubnet25)$color[V(snowsubnet25)$room == 2]  <- "skyblue2"
V(snowsubnet25)$color[V(snowsubnet25)$room == 3]  <- "skyblue3"
V(snowsubnet25)$color[V(snowsubnet25)$room == 4]  <- "skyblue4"
V(snowsubnet25)$color[V(snowsubnet25)$room == 5]  <- "steelblue"
V(snowsubnet25)$color[V(snowsubnet25)$room == 6]  <- "steelblue3"
V(snowsubnet25)$color[V(snowsubnet25)$room == 7]  <- "steelblue4"
V(snowsubnet25)$color[V(snowsubnet25)$room == 8]  <- "royalblue3"
V(snowsubnet25)$color[V(snowsubnet25)$room == 9]  <- "royalblue4"
V(snowsubnet25)$color[V(snowsubnet25)$room == 10]  <- "mediumblue"
V(snowsubnet25)$color[V(snowsubnet25)$room == 19]  <- "navyblue"


pdf("figures/room25_simple.pdf")
par(mfrow = c(1,1))
plot(decompose(snowsubnet25)[[1]], # 94 nodes
     layout = layout.fruchterman.reingold,
     main = "Subgraph for Village 25",
     vertex.label = "")
dev.off()



###### frequency table of village 25 ###########
freq_table <- matrix(NA, ncol = 10, nrow = 2)
freq_table[1,] <- as.integer(names(table(V(subnet25)$room)))
freq_table[2,] <- table(V(subnet25)$room)
rownames(freq_table) <- c("number of rooms", "frequency")
freq_table <- as.data.frame(freq_table)
print(xtable(freq_table , digits = 4, row.names = FALSE))



