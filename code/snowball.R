library(igraph)

## snowball sampling here is in the contexts of sampling depending on latent variable (X)
## so that the observations in the sample are likely to dependent each other via X.

snowball_sampling = function(G, samn){
  ### input  
  # G : population igraph (having N subjects)
  # samn : a targeted sample size (n)
  
  if (vcount(G) < samn){
    # exit if the population number is less than the sample size
    return("Population size is not enough for snowball sampling")
  }
  
  snow = c()
  V(G)$name = c(1:length(V(G)))
  # starter : intiator
  starter = sample(1:length(V(G)),1)
  current = c()
  current[1] = V(G)$name[starter]
  count = 1
  
  snow[1] = current[1] 
  
  while (count < samn ){

    nnode = length(current) # the number of subjects in the current stage
      for (i in 1:nnode){
        ngh = neighbors(G, current[i]) # vertex index
        snow = c(snow, V(G)$name[ngh])
        snow = unique(snow) 
    }  
    
    tmp_sample = snow[(count+1):length(snow)]
    

    if (samn < length(snow)){ # if we reach more than the targeted sample size
      need = samn - count # number of subjects needed
      tmp_sample = sample(tmp_sample ,need) 
      snow[(count+1):samn] = tmp_sample
      snow = snow[-c((samn+1):length(snow))]  
    }
    
    current = tmp_sample  
    count = length(snow)

  
  }
  
  if (count == samn){
    return(snow)
  } else {
    return("somthing goes wrong.")
  }
  
  
}