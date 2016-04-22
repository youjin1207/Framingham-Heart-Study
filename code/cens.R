###########################################
library(igraph)
library(MASS)
library(mvrtn)
library(MultinomialCI)
library(xtable)
library(XML)
###########################################
load("data/functions.RData")
###########################################
data <- read.graph(file = "data/cens_networks_anonymized.xml",
                   format = "graphml")


## biggest component in networks
subdata <- decompose(data)[[1]]

table(V(subdata)$pos) # 9 categories
table(V(subdata)$aff) # too diverse
table(V(subdata)$dep) # 23 categories


# edgea attributes : coauthorship / communication edges / acquanintanceship
subdataadj <- get.adjacency(subdata)
subdataadj <- as.matrix(subdataadj)
subdataadj1 <- ifelse(subdataadj >= 1, 1, 0)
subdataadj2 <- ifelse(subdataadj >= 2, 1, 0)
subdataadj3 <- ifelse(subdataadj >= 3, 1, 0)


position <- as.integer(as.factor(V(subdata)$pos))
department <- as.integer(as.factor(V(subdata)$dep))

### standardize phi 
weight.matrix1 <- weight.matrix(subdataadj1)
weight.matrix2 <- weight.matrix(subdataadj2)
weight.matrix3 <- weight.matrix(subdataadj3)

a <- m2.newphi(weight.matrix1[1], weight.matrix1[2],
               weight.matrix1[3], position) 
b <- mean.newphi(weight.matrix1[1], position)
stdphi <- newphi(subdataadj1, position) - b
stdphi <- stdphi / sqrt(a - b^2)



### standardize phi 
weight.matrix1 <- weight.matrix(subdataadj1)
weight.matrix2 <- weight.matrix(subdataadj2)
weight.matrix3 <- weight.matrix(subdataadj3)

a <- m2.newphi(weight.matrix1[1], weight.matrix1[2],
               weight.matrix1[3], department) 
b <- mean.newphi(weight.matrix1[1], department)
stdphi <- newphi(subdataadj1, department) - b
stdphi <- stdphi / sqrt(a - b^2)

############### visualize ##############
############ position ########
V(subdata)$color[V(subdata)$pos == "AssistProf"]  <- "orange3"
V(subdata)$color[V(subdata)$pos == "AssocProf"] <- "orangered"
V(subdata)$color[V(subdata)$pos == "Lecturer"] <- "orchid"
V(subdata)$color[V(subdata)$pos == "Masters"] <- "palegreen"
V(subdata)$color[V(subdata)$pos == "Phd"] <- "dodgerblue"
V(subdata)$color[V(subdata)$pos == "Postdoc"] <- "cyan"
V(subdata)$color[V(subdata)$pos == "Professor"] <- "lightgreen"
V(subdata)$color[V(subdata)$pos == "Researcher"] <- "yellow"
V(subdata)$color[V(subdata)$pos == "Undergrad"] <- "honeydew"


pdf("figures/position_full.pdf")
par(mfrow = c(1,1))
igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
plot(subdata, layout = layout.fruchterman.reingold,
     main = "Collaboration Network and Position", vertex.label="")
dev.off()


######### snowball sampling 
V(subdata)$name <- 1:385
snow <- c()
starter <- sample(V(subdata)$name, 1)
current <- c()
current[1] <- starter
count <- 1
snow[1] <- current[1]
samn <- 250
while(count < samn){
  nnode <- length(current)
  for(i in 1:nnode){
    
    ngh <- neighbors(subdata, current[i])
    snow <- c(snow, V(subdata)$name[ngh])
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
snowsubdata <- subdata - vertices(V(subdata)$name %in% snow) 
V(snowsubdata)$pos <- V(subdata)$pos[-which(V(subdata)$name %in% snow)] 

V(snowsubdata)$color[V(snowsubdata)$pos == "AssistProf"]  <- "orange3"
V(snowsubdata)$color[V(snowsubdata)$pos == "AssocProf"] <- "orangered"
V(snowsubdata)$color[V(snowsubdata)$pos == "Lecturer"] <- "orchid"
V(snowsubdata)$color[V(snowsubdata)$pos == "Masters"] <- "palegreen"
V(snowsubdata)$color[V(snowsubdata)$pos == "Phd"] <- "dodgerblue"
V(snowsubdata)$color[V(snowsubdata)$pos == "Postdoc"] <- "cyan"
V(snowsubdata)$color[V(snowsubdata)$pos == "Professor"] <- "lightgreen"
V(snowsubdata)$color[V(snowsubdata)$pos == "Researcher"] <- "yellow"
V(snowsubdata)$color[V(snowsubdata)$pos == "Undergrad"] <- "honeydew"



pdf("figures/position_simple.pdf")
par(mfrow = c(1,1))
igraph.options(vertex.size = 5, edge.arrow.size = 0.1,
               vertex.label = NULL)
plot(decompose(snowsubdata)[[1]], # 100 nodes
     layout = layout.fruchterman.reingold,
     vertex.label = "",
     main = "Collaboration Network and Position")
dev.off()


####################### Department ########################
V(subdata)$color[V(subdata)$dep == "Biology"]  <- "orange3"
V(subdata)$color[V(subdata)$dep == "Botany"] <- "orangered"
V(subdata)$color[V(subdata)$dep == "ChemicalEng"] <- "orchid"
V(subdata)$color[V(subdata)$dep == "Civil&EnvEng"] <- "palegreen"
V(subdata)$color[V(subdata)$dep == "CivilEng"] <- "palevioletred"
V(subdata)$color[V(subdata)$dep == "CS"] <- "yellow"
V(subdata)$color[V(subdata)$dep == "Earth&Space"] <- "royalblue"
V(subdata)$color[V(subdata)$dep == "Ecology"] <- "skyblue"
V(subdata)$color[V(subdata)$dep == "Education"] <- "purple"
V(subdata)$color[V(subdata)$dep == "EE"]  <- "aquamarine4"
V(subdata)$color[V(subdata)$dep == "Engineering"] <- "beige"
V(subdata)$color[V(subdata)$dep == "Environment"] <- "blue"
V(subdata)$color[V(subdata)$dep == "Film"] <- "brown"
V(subdata)$color[V(subdata)$dep == "Geology"] <- "cyan"
V(subdata)$color[V(subdata)$dep == "IS"] <- "firebrick"
V(subdata)$color[V(subdata)$dep == "Law"] <- "gold"
V(subdata)$color[V(subdata)$dep == "Linguistics"] <- "deeppink3"
V(subdata)$color[V(subdata)$dep == "Marine"] <- "khaki"
V(subdata)$color[V(subdata)$dep == "Mathematics"] <- "ivory"
V(subdata)$color[V(subdata)$dep == "MechEng"] <- "lightsalmon"
V(subdata)$color[V(subdata)$dep == "Meteorology"] <- "lightcyan"
V(subdata)$color[V(subdata)$dep == "Statistics"] <- "dodgerblue"
V(subdata)$color[V(subdata)$dep == "UrbanPlanning"] <- "coral"
E(subdata)$color <- "azure3"

pdf("figures/dep_full.pdf")
par(mfrow = c(1,1))
igraph.options(vertex.size = 4, edge.arrow.size = 0.1,
               vertex.label = NULL)
plot(subdata, layout = layout.fruchterman.reingold,
     main = "Collaboration Network and Department", vertex.label="")
dev.off()



V(snowsubdata)$dep <- V(subdata)$dep [-which(V(subdata)$name %in% snow)] 

V(snowsubdata)$color[V(snowsubdata)$dep == "Biology"]  <- "orange3"
V(snowsubdata)$color[V(snowsubdata)$dep == "Botany"] <- "orangered"
V(snowsubdata)$color[V(snowsubdata)$dep == "ChemicalEng"] <- "orchid"
V(snowsubdata)$color[V(snowsubdata)$dep == "Civil&EnvEng"] <- "palegreen"
V(snowsubdata)$color[V(snowsubdata)$dep == "CivilEng"] <- "palevioletred"
V(snowsubdata)$color[V(snowsubdata)$dep == "CS"] <- "yellow"
V(snowsubdata)$color[V(snowsubdata)$dep == "Earth&Space"] <- "royalblue"
V(snowsubdata)$color[V(snowsubdata)$dep == "Ecology"] <- "skyblue"
V(snowsubdata)$color[V(snowsubdata)$dep == "Education"] <- "purple"
V(snowsubdata)$color[V(snowsubdata)$dep == "EE"]  <- "aquamarine4"
V(snowsubdata)$color[V(snowsubdata)$dep == "Engineering"] <- "beige"
V(snowsubdata)$color[V(snowsubdata)$dep == "Environment"] <- "blue"
V(snowsubdata)$color[V(snowsubdata)$dep == "Film"] <- "brown"
V(snowsubdata)$color[V(snowsubdata)$dep == "Geology"] <- "cyan"
V(snowsubdata)$color[V(snowsubdata)$dep == "IS"] <- "firebrick"
V(snowsubdata)$color[V(snowsubdata)$dep == "Law"] <- "gold"
V(snowsubdata)$color[V(snowsubdata)$dep == "Linguistics"] <- "deeppink3"
V(snowsubdata)$color[V(snowsubdata)$dep == "Marine"] <- "khaki"
V(snowsubdata)$color[V(snowsubdata)$dep == "Mathematics"] <- "ivory"
V(snowsubdata)$color[V(snowsubdata)$dep == "MechEng"] <- "lightsalmon"
V(snowsubdata)$color[V(snowsubdata)$dep == "Meteorology"] <- "lightcyan"
V(snowsubdata)$color[V(snowsubdata)$dep == "Statistics"] <- "dodgerblue"
V(snowsubdata)$color[V(snowsubdata)$dep == "UrbanPlanning"] <- "coral"
E(snowsubdata)$color <- "azure3"

pdf("figures/dep_simple.pdf")
par(mfrow = c(1,1))
igraph.options(vertex.size = 5, edge.arrow.size = 0.1,
               vertex.label = NULL)
plot(decompose(snowsubdata)[[1]], # 100 nodes
     layout = layout.fruchterman.reingold,
     vertex.label = "",
     main = "Collaboration Network and Department")
dev.off()



######### frequency table ############
### position
freq_pos <- matrix(NA, ncol = 9, nrow = 1)
colnames(freq_pos) <-  names(table(V(subdata)$pos))
freq_pos[1,] <- as.integer(table(V(subdata)$pos))
rownames(freq_pos) <-  "frequency"
freq_pos <- as.data.frame(freq_pos)
print(xtable(freq_pos , digits = 4, row.names = FALSE))



### departments
freq_dep <- matrix(NA, ncol = 8, nrow = 6)
freq_dep[1,] <- names(table(V(subdata)$dep))[1:8]
freq_dep[2,] <-  as.integer(table(V(subdata)$dep))[1:8]
freq_dep[3,] <- names(table(V(subdata)$dep))[9:16]
freq_dep[4,] <-  as.integer(table(V(subdata)$dep))[9:16]
freq_dep[5,1:7] <- names(table(V(subdata)$dep))[17:23]
freq_dep[6,1:7] <-  as.integer(table(V(subdata)$dep))[17:23]


freq_dep <- as.data.frame(freq_dep)
print(xtable(freq_dep , digits = 4, row.names = FALSE))

