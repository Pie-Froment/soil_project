
# additionnal analyses for the presentation.
# source("analysis/cleaning.data.R")
# source("analysis/statistical.analysis.R")

# peace to get Paracoudata:
library(EcoFoG)
paracou = Paracou2df(
  WHERE = NULL,
  UID = NULL,
  PWD = NULL,
  Driver = "SQL Server Native Client 10.0"
)

paracou$ssp = paste(paracou$Genus,paracou$Species, sep = "_")
length(unique(paracou$ssp))
length(unique(paracou$Genus))


# Calculate bate diversity between forest and agroforestry plot:----
#create the metacommunity object: 
tsbf_la %>% 
  group_by(site)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop")%>%
  select(-site)%>%t()%>% data.frame() -> tsbf_meta


tsbf_meta = MetaCommunity(tsbf_meta)
BetaDiversity(tsbf_meta, q = 4, Correction = "None")

# new function to get the gamma vector from the two new samples community:
betadiv = function(MC, q = 1, Correction = "None"){
  return(BetaDiversity(MetaCommunity(MC), q = q,
                       Correction = Correction)$Total)
}

# Function to compute betadiversity for new sampled communities:
newcom2 = function(community1,community2, q.seq){
  sampleco1 = table(sample(x = names(community1),
                                    prob = as.ProbaVector(community1),
                                    size =  sum(community1),
                                    replace = T))
  sampleco2 = table(sample(x = names(community2),
                                     prob = as.ProbaVector(community2),
                                     size =  sum(community2),
                                     replace = T))
  gammaco <- data.frame(matrix(ncol = 29,nrow = 2,data = 0))
  names(gammaco) = taxon_focus
  
  gammaco[1,intersect(names(gammaco), names(sampleco1))] <-
    gammaco[1,intersect(names(gammaco), names(sampleco1))]+
    sampleco1[intersect(names(gammaco), names(sampleco1))]
  
  gammaco[2,intersect(names(gammaco), names(sampleco2))] <-
    gammaco[2,intersect(names(gammaco), names(sampleco2))]+
    sampleco2[intersect(names(gammaco), names(sampleco2))]
  gammaco <- data.frame(t(gammaco))
  return(unlist(lapply(q.seq, function(q.seq) betadiv(MC = gammaco, q = q.seq, 
                                                      Correction = "None"))))
}


# Function to create a confidence interval for the betadiversity:
conf.int2 = function(community1,community2,repp = 10,problow = 0.025,
                    probhigh = 0.975,
                    q.seq = seq(0,2,0.1)){
  vectors = lapply(1:repp,function(x)newcom2(community1,community2, q.seq))# draw repp commmunities
  matv = t(do.call(cbind,vectors))
  newvector = split(matv, col(matv))
  
  
  betaobs <- data.frame(matrix(ncol = 29,nrow = 2,data = 0))
  names(betaobs) = taxon_focus
  
  betaobs[1,intersect(names(betaobs), names(community1))] <-
    betaobs[1,intersect(names(betaobs), names(community1))]+
    community1[intersect(names(betaobs), names(community1))]
  betaobs[2,intersect(names(betaobs), names(community2))] <-
    betaobs[2,intersect(names(betaobs), names(community2))]+
    community2[intersect(names(betaobs), names(community2))]
  betaobs = data.frame(t(betaobs))

  return(list(lowconf = sapply(newvector, quantile,probs = problow),
              highconf= sapply(newvector, quantile,probs = probhigh),
              obs = unlist(lapply(q.seq, 
                                  function(q.seq) betadiv(MC = betaobs, 
                                                          q = q.seq,
                                                          Correction = "None"))) ,
              x = q.seq))
}


bootlist2 = conf.int2(community1 = datalist$cocoa, community2 = datalist$forest,
                      repp = 100, q.seq = seq(0,4,0.2))


ggplot()+geom_line(mapping = aes(x = bootlist2$x, y = bootlist2$obs),
                   colour = "black", linewidth = 0.9)+
  theme_bw()+ labs(x = "Order of Diversity",y = "Beta Diversity")+
  geom_ribbon(mapping = aes(x = bootlist2$x, ymin = bootlist2$lowconf,
                            ymax= bootlist2$highconf), 
              alpha = 0.3,colour ="lightgray")


# plot(CommunityProfile(Diversity, colSums(tsbf_la[,taxon_focus]),
#                  Correction = "None"))
# ccc = CommunityProfile(Diversity, colSums(tsbf_la[,taxon_focus]),
#                  Correction = "None")
# ccc1 = CommunityProfile(Diversity, datalist$cocoa,
#                  Correction = "None")
# ccc2 = CommunityProfile(Diversity, datalist$forest,
#                  Correction = "None")
# plot(x = ccc1$x, y = ccc$y/((ccc1$y+ccc2$y)/2) )
# 
# CommunityProfile(BetaDiversity, tsbf_meta,q.seq = seq(0,4,0.2), 
#                  Correction = "None")

