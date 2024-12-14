# Soil project
# Statistical analyses



# Libraries: ----

library(lme4) # for mixed models (in case)
library(lmtest) # to check linear models assumptions
library(car)# for type II Anovas

library(FactoMineR) # To carry out multivariate analyses (PCA, CA)
library(factoextra)
library(vegan) # for PCoA using Bray-Curtis distances
library(ade4) # for multivariate analyses (PcoA)
library(ape) # for multivariate analyses (PcoA)

library(entropart) # for diversity analyses
library(gtools) # to try the order() function that gives wrong order.
library(ggpubr)

# setting the seed:
set.seed(76)

# Importing data: ----

# this code will be embedded within the main Rmarkdown. Run the line below if 
# you are running it independently:
# source("analysis/cleaning.data.R")


# for the analyses, we get rid of the yellow and green eggs taxon (probably
# fertilisers) and of unknonwn taxons + taxa with 0 observations.
# We create a taxon_focus taxa list which we will use for our analyses to focus
# only on taxa of interests:
taxon_focus = lexicon$code
#  taxa with 0 observations:
selec= names(colSums(tsbf[,taxon_focus])[colSums(tsbf[,taxon_focus])==0])
taxon_focus = taxon_focus[!(taxon_focus%in%c("yl.egg","gr.egg","unk", selec))]
rm(selec)


#1: TSBF diversity analysis: ----


## A: Correspondence analysis: ----

CA_tsbf = CA(tsbfPA[,taxon_focus], graph = F, ncp = 6)

fviz_eig(CA_tsbf)
fviz_ca_biplot(CA_tsbf, axes = c(1,2),
        col.row=tsbfPA$site)
fviz_ca_biplot(CA_tsbf, axes = c(1,2),
               col.row=tsbfPA$layer)
fviz_ca_row(CA_tsbf, col.row = tsbfPA$layer,
            axes = c(1,2))

fviz_ca_row(CA_tsbf, col.row = tsbfPA$site,
            axes = c(1,2))


## B: PCoA/Bray Curtis distances: ----
# We use bray curtis distances as it integrates both the composition and the 
# relative abundance of each species in its distance.

# Computing the distance matrix: 
bcdist = vegdist(tsbf[,taxon_focus], "bray")

# Bray-Curtis distances are non-metric: can convert them to a metric by using 
# their Square root: 
bcdist_sqrt = sqrt(bcdist)
# Checking the distortion of the square root transformation compared to initial
# Bray-curtis distances:
plot(bcdist~bcdist_sqrt, xlim = c(0,1), ylim = c(0,1))
cor(bcdist,bcdist_sqrt)
is.euclid(bcdist_sqrt) # it is euclidean, so we can use a PcoA to represent the
# distances.

# Computing the PcoA:
pcoa_tsbf = pcoa(bcdist_sqrt)

# Checking the eigen values: 
ggplot()+ geom_bar(mapping = aes(x = order(pcoa_tsbf$values$Relative_eig, 
                                           decreasing = T),
                                 y = pcoa_tsbf$values$Relative_eig), 
                   stat = "identity")+
  xlab("Axes")+ylab("Proportion of variance")+
  theme_minimal()
# Here, we see that only 25% of the information is retained is the first 2 axes.

# Plotting the first two axes (colour = Area):
ggplot()+ geom_point(data= pcoa_tsbf$vectors, 
                     mapping = aes(x = Axis.1, y = Axis.2, 
                                   colour = tsbf$site),
                     shape = 16)+
  theme_bw()
# on axes 3 and 4:
# ggplot()+ geom_point(data= pcoa_tsbf$vectors, 
#                      mapping = aes(x = Axis.3, y = Axis.4, 
#                                    colour = tsbf$site),
#                      shape = 16)+
#   theme_bw()

# again but we colour by layers:
ggplot()+ geom_point(data= pcoa_tsbf$vectors, 
                     mapping = aes(x = Axis.1, y = Axis.2, 
                                   colour = tsbf$layer),
                     shape = 16)+
  theme_bw()
# Here, we see that the points clearly seperate based on the layer, but not so
# much based on the site. 

# Computing a Permanova on Bray-Curtis distances:

# check whether the dispersion in bray-Curtis distance (sqrt) varies between 
# the groups (site, layer):
anova(betadisper(bcdist_sqrt,tsbf$site)) 
anova(betadisper(bcdist_sqrt,tsbf$layer))
# results: homogeneous dispersion for both sites and layers in bray-curtis 
# distances. 

adonis2(bcdist_sqrt~layer, data = tsbf, permutations = 999)
adonis2(bcdist_sqrt~site, data = tsbf, permutations = 999)
adonis2(bcdist_sqrt~site*layer, data = tsbf, permutations = 999)


## UPDATE 22/11/2024: ----
# After discussing with Irene, choice to summarise data by
# sampling point (only 10 points remaining, 5 in each site) to get rid of
# the dispersion caused by the different layer (we know that the layers have
# different composition).

# Here we create a new dataset where we summed all taxa abundance in the 
# different layers in each tsbf sites (for the accumulation curve)
tsbf %>% select(c("site","replicat",all_of(taxon_focus)))%>%
  group_by(site,replicat)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") ->
  tsbf_co 

# we repeat the steps in B to do a PCOA based on bray-curtis measures:
bcdist2 = vegdist(tsbf_co[,taxon_focus],"bray")
bcdist_sqrt2 = sqrt(bcdist2)
plot(bcdist2~bcdist_sqrt2, xlim = c(0,1), ylim = c(0,1))
cor(bcdist2,bcdist_sqrt2)
is.euclid(bcdist_sqrt2) 

# We test whether these beta diversity measures differentiate forest from cocoa
# sites using a Permanova:
# testing for overdispersion first:
anova(betadisper(bcdist_sqrt2, tsbf_co$site)) # no overdispersion.
finalTSBFbetatest = adonis2(bcdist_sqrt2~site, data = tsbf_co, permutations = 999) 
# no differences in centroids.

# graphical representation using a PcoA:
# Computing the PcoA:
pcoa_tsbf2 = pcoa(bcdist_sqrt2)

# Checking the eigen values: 
eigplot = ggplot()+ geom_bar(mapping = aes(x = order(pcoa_tsbf2$values$Relative_eig, 
                                           decreasing = T),
                                 y = pcoa_tsbf2$values$Relative_eig), 
                   stat = "identity")+
  labs(x = "Axes", y= "Variance")+
  theme_minimal()+theme(plot.background = element_rect(fill = "white"),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank())
# Again quite bad, only 45% of information retained in the first two axes

# Plotting the first two axes (colour = Area):
axeplot = ggplot()+ geom_point(data= pcoa_tsbf2$vectors, 
                     mapping = aes(x = Axis.1, y = Axis.2, 
                                   colour = tsbf_co$site),
                     shape = 16)+
  scale_color_manual(values = c("#e31a1c","darkgreen"))+
  theme_bw()+labs(color = "Sites")+
  stat_ellipse(data= pcoa_tsbf2$vectors, 
               mapping = aes(x = Axis.1, y = Axis.2, 
                             colour = tsbf_co$site),level= 0.95,
               alpha = 0.1, geom = "polygon")

axeplot2 = ggplot()+ geom_point(data= pcoa_tsbf2$vectors, 
                     mapping = aes(x = Axis.3, y = Axis.4, 
                                   colour = tsbf_co$site),
                     shape = 16)+
  scale_color_manual(values = c("#e31a1c","darkgreen"))+
  theme_bw()+labs(color = "Sites")+
  stat_ellipse(data= pcoa_tsbf2$vectors, 
               mapping = aes(x = Axis.3, y = Axis.4, 
                             colour = tsbf_co$site),level= 0.95,
               alpha = 0.1, geom = "polygon")+
  theme(legend.position = "none")

eigplot = ggplotGrob(eigplot)
axeplot = axeplot+annotation_custom(grob = eigplot,
                          xmin = -1, xmax = -0.25, ymin = -0.97, ymax = -0.22)
finalTSBFbetaplot = ggarrange(axeplot, axeplot2, ncol = 2, common.legend = T,
                           legend = "right")

# repeating this steps within the different layers (we have already tested for
# distortion and euclideanity):
# Here no differences between sites, we put the loop in comment

# layers = c("L","10","20")
# 
# for (i in layers){
#   distmat = sqrt(vegdist(tsbf[tsbf$layer == i,taxon_focus]))# dist matrix
#   anova(betadisper(distmat, tsbf[tsbf$layer == i,"site"])) # no overdispersion.
#   adonis2(bcdist_sqrt2~site, data = tsbf[tsbf$layer == i,], permutations = 999)
# 
#   pcoa_tsbf2 = pcoa(distmat) # pcoa
#   #eigen values: 
#   ggplot()+ geom_bar(mapping = aes(x = order(pcoa_tsbf2$values$Relative_eig, 
#                                              decreasing = T),
#                                    y = pcoa_tsbf2$values$Relative_eig), 
#                      stat = "identity")+
#     xlab("Axes")+ylab("Proportion of variance")+
#     theme_minimal()
#     # Plotting the first two axes (colour = Area):
#     ggplot()+ geom_point(data= pcoa_tsbf2$vectors, 
#                          mapping = aes(x = Axis.1, y = Axis.2, 
#                                        colour = tsbf_co$site),
#                          shape = 16)+
#     scale_color_manual(values = c("#e31a1c","darkgreen"))+
#     theme_bw()+labs(color = "Sites")
#   
# }

## C: Species abundances analyses: ----

### Abundances between sites ----

# Tests to compare for median abundances of 
# the different orders between the forest and the cocoa plots: 

# We will us the tsbf_co matrix for our analyse, which summed abundances of all
# orders at the different layers (litter, 0-10cm and 0-20cm) for each TSBF 
# digging site.

# We will restrain our analyses for taxa counted at least 10 times overall:

taxa_test = 
  names(colSums(tsbf_co[,taxon_focus])[colSums(tsbf_co[,taxon_focus])>10])

testres = data.frame(taxa = taxon_focus, w = NA, pw = NA,
                     Fstat = NA, pF = NA, df = NA, shaptest = NA,
                     bptest = NA, dwtest = NA)
for(i in taxa_test){
  # non-parametric wilcox.test
  testwil = wilcox.test(unlist(tsbf_co[tsbf_co$site == "Forest",i]),
                        unlist(tsbf_co[tsbf_co$site == "Cocoa",i]))
  testres$w[testres$taxa == i] = testwil$statistic
  testres$pw[testres$taxa == i] = testwil$p.value
  
  # We also look at analyses of variances: 
  testan = lm(get(i)~site, data = tsbf_co)
  testres$Fstat[testres$taxa == i] = Anova(testan)$`F value`[1]
  testres$df[testres$taxa == i] = paste(Anova(testan)$Df[1],",",
                                        Anova(testan)$Df[2])
  testres$pF[testres$taxa == i] = Anova(testan)$`Pr(>F)`[1]
  testres$shaptest[testres$taxa == i] = shapiro.test(resid(testan))$p.value
  testres$bptest[testres$taxa == i] = bptest(testan)$p.value
  testres$dwtest[testres$taxa == i] = dwtest(testan)$p.value
  
}
rm(testan, i , testwil)
# Looking at testres, we see that there is no big differences of significances
# between wilcox rank test and anovas. Only differences is that the normality is
# breached for some taxa in the anova. Therefore, I choose the wilcox test to
# compare abundances of each taxon between the two sites.

### Pyramid graph ----
# Creation of a nice graph showing the pyramid of age abundances. Show also the 
# significance obtained from the wilcox test inside "testres".
# data.frame:
tsbf %>% select(c("site",all_of(taxon_focus)))%>%
  group_by(site)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") %>%
  pivot_longer(cols = all_of(taxon_focus),
               names_to = "taxa",
               values_to = "abundance") -> tsbfplot

# We add a column with the taxon order:
col = colSums(tsbf[,taxon_focus])
col = sort(col)
newvec = 1:29
# very dirty corde but it works (did not manage to use order())
tsbfplot$rank=newvec[match(tsbfplot$taxa,names(col))]
rm(newvec,col)

# change abundances for Cocoa to negative abundance for the plotting:
tsbfplot$abundance = ifelse(tsbfplot$site == "Cocoa",
                            -tsbfplot$abundance,
                            tsbfplot$abundance)
tsbfplot %>% 
  left_join(testres, by = "taxa") %>%
  select(site,taxa,abundance,rank,pw) -> tsbfplot

# Computing labels position to show significances of wilcox tests:
labelss = unlist(tsbfplot[tsbfplot$site == "Cocoa", "abundance"])

sign = c(unlist(ifelse(tsbfplot[tsbfplot$site=="Cocoa","pw"] > 0 &
                       tsbfplot[tsbfplot$site=="Cocoa","pw"]<0.05 & 
                       !is.na(tsbfplot[tsbfplot$site=="Cocoa","pw"]),
                     "*","")))

#Plotting abundances:
finalTSBFabundanceplot = ggplot()+
  geom_bar(data= tsbfplot, mapping = aes(x = as.numeric(rank),
                                         y = abundance, fill = site),
           stat = "identity", position = "identity")+
  scale_fill_manual(values = c("#e31a1c","darkgreen"))+
  coord_flip()+
  theme_minimal()+
  scale_y_continuous(breaks = c(seq(-300,300,100)),
                     labels = c(seq(300,0,-100),seq(100,300,100)),
                     limits = c(-300,300))+
  scale_x_continuous(breaks = c(1:29),
                     labels = unique(tsbfplot$taxa[order(tsbfplot$rank)]))+
  ylab("Abundance")+xlab("Taxa (order)")+
  geom_text(mapping = aes(x = tsbfplot$rank[tsbfplot$site=="Cocoa"],
                          y = labelss-10, 
                          label = sign), fontface = "bold", size = 6)+
  theme(panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(angle = 0))



### trophic preferences: ----

# creating a vector for indexing columns of interest: 
trophic_focus = c("phytophagous","saprophagous","predator",
                  "rhizophagous","xylophagous")

# creating an empty table to fill the values:
tsbf %>% select(c("ID","site","plot","replicat","layer")) %>%
  mutate(phytophagous = 0, saprophagous = 0,predator = 0,
         rhizophagous = 0, xylophagous = 0) ->tsbftrophic

# filling the values (to be vectorised): 
# when an order has more than one trophic preference, we assume that each 
# individual only has one (trophic preferences are mutually exclusive). 
# We assume that each trophic preferences are as common in this order.
# I equally split the order abundances equally in the order trophic preferences
for (i in seq_along(tsbf$replicat)){
  vec = tsbf[i,taxon_focus]
  for(p in seq_along(vec)){
    index = which(lexicon$code ==names(vec)[p])
    for(j in trophic_focus){
      if(lexicon[index,names(lexicon)==j]==1){
        tsbftrophic[tsbftrophic$ID==tsbf$ID[i],j]=
          tsbftrophic[tsbftrophic$ID==tsbf$ID[i],j]+
          tsbf[i,names(vec)[p]]/sum(lexicon[index,trophic_focus])
      }
    }
  }
}
# checking that the total number of indiv is equal between tsbftrophic and
# tsbf:
sum(colSums(tsbftrophic[,trophic_focus])) == sum(colSums(tsbf[,taxon_focus]))

# This loop is considering that each individual in orders with multiple trophic 
# preferences has every trophic preferences
# for (i in seq_along(tsbf$replicat)){
#   vec = tsbf[i,taxon_focus]
#   for(p in seq_along(vec)){
#     index = which(lexicon$code ==names(vec)[p])
#     for(j in trophic_focus){
#       if(lexicon[index,names(lexicon)==j]==1){
#         tsbftrophic[tsbftrophic$ID==tsbf$ID[i],j]=
#           tsbftrophic[tsbftrophic$ID==tsbf$ID[i],j]+
#           tsbf[i,names(vec)[p]]
#       }
#     }
#   }
# }

# Investigating the associations between trophic groups and the two habitats
# by doing a correspondence analysis:
pca_trophic = CA(tsbftrophic[,trophic_focus],
                  graph = F)
eigplot = fviz_eig(pca_trophic, ncp = 10, # it worked well
         linecolor = NA, title = "")

pca_trophicp = data.frame(pca_trophic$row$coord)

caplot = fviz_ca_col(pca_trophic, axes =c(1,2),
             col.col = "black", repel = F, alpha = "cos2",
            title = "", labelsize = 5, labelfont = "bold")+
  geom_point(data = pca_trophicp, 
             mapping = aes(x = Dim.1,y = Dim.2,
                           color = tsbftrophic$site,
                           shape =  tsbftrophic$layer), size = 2.5,
             stroke = 1.4, alpha = 0.7) +
  scale_shape_manual(values = c(1,2,3))+
  scale_color_manual(values = c("#e31a1c","darkgreen"))+
  labs(color = "Sites", shape = "Layers", 
       alpha = "cos2 (columns)")

finalTSBFtrpohicplot = ggarrange(caplot, 
          eigplot,
          widths = c(16,7),
          ncol = 2, labels = c("A","B"))

# to check rhizophagous projection: 
# fviz_ca_biplot(pca_trophic, axes = c(1,4))

# mm = cbind(pca_trophic$row$coord, tsbftrophic$site,
#       tsbftrophic$layer)
pca_trophic$col$cos2


## D: Species-accumulation curve ----

# Create the accumulation curve using the specaccum function (method = random):
# here we use the tsbf_co data, where abundances in the different layers have
# have been summed.
Cocoa = specaccum(tsbf_co[tsbf_co$site=="Cocoa",taxon_focus],
          method = "random", permutations = 100,
          ci = 0.95)
Forest = specaccum(tsbf_co[tsbf_co$site=="Forest",taxon_focus],
          method = "random", permutations = 100,
          ci = 0.95)

# We plot these curves in a very tedious way:

finalTSBFaccumplot= ggplot()+ # Cocoa line:
  geom_line(mapping = aes(x = unlist(Cocoa$sites), y = unlist(Cocoa$richness)),
            linewidth = 1.2, col = "#e31a1c", show.legend = T)+
  geom_ribbon(mapping = aes(ymin = unlist(Cocoa$richness) - unlist(Cocoa$sd), 
                            ymax = unlist(Cocoa$richness) + unlist(Cocoa$sd),
                            x = unlist(Cocoa$sites)), colour ="lightgray",
              alpha = 0.3, show.legend = T)+
  theme_bw()+# forest line: 
  geom_line(mapping = aes(x = unlist(Forest$sites), y = unlist(Forest$richness)),
            linewidth = 1.2, col = "darkgreen", show.legend = T)+
  geom_ribbon(mapping = aes(ymin = unlist(Forest$richness) - unlist(Forest$sd), 
                            ymax = unlist(Forest$richness) + unlist(Forest$sd),
                            x = unlist(Forest$sites)), colour ="lightgray",
              alpha = 0.3, show.legend = T)+
  xlab("Sites")+ylab("Richness (order)")+
  geom_line(mapping = aes(x = c(1,1), y = c(1,2), 
                          # geom_line just to make a legend appear
                           colour = c("#e31a1c", "darkgreen")))+
  scale_color_manual(values = c("#e31a1c", "darkgreen"),# modifying the legend
                     labels = c("Cocoa","Forest"))+
  ylim(c(9,26))+ labs(color = "Sites")



## E: Diversity profiles: ----

# For diversity profiles, we want to differentiate between the different
# layers (could be interesting), as well as between the sites:

# Creating a new dataset summing data by layers and site:

tsbf %>% select(c("site","layer",all_of(taxon_focus)))%>%
  group_by(site,layer)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") ->
  tsbf_la 


# Plotting the diversity profiles using the whittaker plots:

datalist = list(cocoa = colSums(tsbf_la[tsbf_la$site == "Cocoa",taxon_focus]),
                forest = colSums(tsbf_la[tsbf_la$site == "Forest",taxon_focus]),
                Lco =colSums(tsbf_la[tsbf_la$layer == "L"&
                                       tsbf_la$site=="Cocoa",taxon_focus]),
                Lfo = colSums(tsbf_la[tsbf_la$layer == "L"&
                                        tsbf_la$site=="Forest",taxon_focus]),
                tenco = colSums(tsbf_la[tsbf_la$layer == "10"&
                                        tsbf_la$site=="Cocoa",taxon_focus]),
                tenfo = colSums(tsbf_la[tsbf_la$layer == "10"&
                                         tsbf_la$site=="Forest",taxon_focus]),
                twentyco = colSums(tsbf_la[tsbf_la$layer == "20"&
                                             tsbf_la$site=="Cocoa",taxon_focus]),
                twentyfo = colSums(tsbf_la[tsbf_la$layer == "20"&
                                             tsbf_la$site=="Forest",taxon_focus]))

#Whittaker plots in the different sub-grouping (exploratory):
  # All:
  plot(as.AbdVector(colSums(tsbf_la[,taxon_focus])))
# Cocoa:
as.AbdVector(datalist$cocoa)%>%
  plot()
# Forest:
as.AbdVector(datalist$forest)%>%
  plot()
# maybe do a ggplot combining the two whitakker's plots together. 

# 
# divcocoa = CommunityProfile(Diversity,datalist$cocoa,
#                  Correction = "None")
# divforest=CommunityProfile(Diversity,datalist$forest,
#                  Correction = "None")


# Finding a way to compute confidence interval by bootstrapping:

# function to create a new taxa distribution based on the frequency of the ones
# we obtained using the TSBF method:
newcom = function(community, q.seq){
  sampleco = as.vector(table(sample(x = names(community),
                          prob = as.ProbaVector(community),
                          size =  sum(community),
                          replace = T)))
  return(CommunityProfile(Diversity,sampleco,
                          Correction = "None", q.seq = q.seq)$y)
}
# to check if it works: newcom(datalist$cocoa)

# new functions to repeat community draws repp times (100 times, 200 etc.)
conf.int = function(community,repp = 10,problow = 0.025,probhigh = 0.975,
                    q.seq = seq(0,2,0.1)){
  vectors = lapply(1:repp,function(x)newcom(community, q.seq))# draw repp commmunities
  matv = t(do.call(cbind,vectors))
  newvector = split(matv, col(matv))
  obscom= CommunityProfile(Diversity,community, Correction = "None", q.seq=q.seq)
  return(list(lowconf = sapply(newvector, quantile,probs = problow),
              highconf= sapply(newvector, quantile,probs = probhigh),
              obs = obscom$y,
              x = q.seq))
}

# we compute confidence interval for all the layers+overall sites in datalist: 
bootlist = lapply(datalist,conf.int,repp = 500, q.seq = seq(0,4,0.2))

# bootlist = lapply(datalist,conf.int,repp = 500, q.seq = seq(0,16,0.8))

# function to do the graph quickly:
plotfun = function(bootcocoa, bootforest, title){
  plotdiv= ggplot()+geom_line(mapping = aes(x = bootcocoa$x, y = bootcocoa$obs),
                              colour = "#e31a1c", linewidth = 0.9)+
    theme_bw()+ labs(x = "Order of Diversity",y = "Diversity")+
    geom_ribbon(mapping = aes(x = bootcocoa$x, ymin = bootcocoa$lowconf,
                              ymax= bootcocoa$highconf), 
                alpha = 0.3,colour ="lightgray")+
    geom_line(mapping = aes(x = bootforest$x, y = bootforest$obs), 
              colour = "darkgreen", linewidth = 0.9)+
    geom_ribbon(mapping = aes(x = bootforest$x, ymin = bootforest$lowconf,
                              ymax= bootforest$highconf), 
                alpha = 0.3,colour ="lightgray")+
    # here just to make a legend appear: 
    geom_line(mapping = aes(x = c(1000,1000), y = c(1000,2000), 
                            # geom_line just to make a legend appear
                            colour = c("#e31a1c", "darkgreen")))+
    scale_color_manual(values = c("#e31a1c", "darkgreen"),# modifying the legend
                       labels = c("Cocoa","Forest"))+
    ylim(c(min(c(bootforest$lowconf,bootcocoa$lowconf))-1,
           c(max(c(bootforest$highconf,bootcocoa$highconf))+1)))+ 
    xlim(c(0,max(bootcocoa$x)))+labs(color = "Sites", title = title)
  return(plotdiv)
}

plot1 = plotfun(bootlist$cocoa,bootlist$forest, title = "All layers")
plot2 = plotfun(bootlist$Lco,bootlist$Lfo, title = "Leaf litter")
plot3 = plotfun(bootlist$tenco,bootlist$tenfo, title = "0-10 cm")
plot4 = plotfun(bootlist$twentyco,bootlist$twentyfo, title = "10-20 cm")

finalTSBFdivplot =ggarrange(plot1, plot2, plot3, plot4, 
          ncol = 2, nrow = 2,                
          labels = c("A", "B", "C", "D"), common.legend = T, legend = "right")

#2: Wood decomposition rate analysis: ----

# here, we will model the wood mass loss (initial mass-final mass) as a function
# of mesh size (effect of mesofauna vs microfauna), site (forest effect 
# vs Cocoa effect) and their interaction (differential effect of mesofauna in the
# different sites)

# Since we have two sticks by bags and these sticks are complete pseudoreplicates 
# we merge the two sticks measurements into one average weight loss per bag.

wood %>% 
  pivot_wider(names_from =mesh, values_from = proploss,
              values_fill = list(proploss = 0))%>%
  group_by(site,plot, replicat)%>%
  summarise(tight.mesh = max(tight.mesh),
         wide.mesh = max(wide.mesh))%>%
  mutate(faunaproploss = wide.mesh-tight.mesh) -> woodfauna

# Plotting the data: 
# Plotting the percentage of weight loss versus site (mesh = colour)
ggplot(data = wood)+ 
  geom_point(mapping = aes(x = jitter(as.numeric(as.factor(site)), 0.1), 
                           y = (startwht-endwht)/startwht, colour = mesh))+
  theme_bw()+
  xlim(0.5,2.5)

# same plot but a box plot:
ggplot(data = wood)+ 
  geom_boxplot(mapping = aes(x = site, 
                           y = (startwht-endwht)/startwht, colour = mesh))+
  theme_bw()



## A: first model ----
# (here, a simple linear model with an interaction):
# Notes: since the response variable is a percentage (bounded between 0 and 100),
# maybe we have to adopt a generalised linear model (with a beta distribution).
# UPDATE: we abandoned this model to focus on modelling the loss caused by fauna
# (losses of wide mesh-losses of tight mesh)for each locations:

M_wood = lm((startwht-endwht)/startwht~site*mesh, data= wood)

# Checking model assumptions: 
par(mfrow = c(2,2))
plot(M_wood)
par(mfrow = c(1,1))

hist(resid(M_wood))# distribution of residuals

shapiro.test(resid(M_wood)) # normality
bptest(M_wood) # Homogeneity
dwtest(M_wood) # independence

Anova(M_wood) # type 2 anova
summary(M_wood) # summary

## B: second model ----
# Following meeting with Irene, (22/11/2024):
# we focus on predicting the mass loss proportion caused by fauna (proportion of
# mass lost in wide-meshed bags- proportion of mass lost in tight-mesh bags):

# quick plot: 
ggplot(data = woodfauna)+
  geom_boxplot(mapping = aes(x = site, y = faunaproploss))

M_wood = lm(faunaproploss~site, data = woodfauna)

hist(resid(M_wood))# distribution of residuals: not normal

# Checking model assumptions: 
par(mfrow = c(2,2))
plot(M_wood)
par(mfrow = c(1,1))

shapiro.test(resid(M_wood)) # normality: there is non-normality of the residuals
bptest(M_wood) # Homogeneity
dwtest(M_wood) # independence

Anova(M_wood) # type 2 anova
summary(M_wood) # summary

# We move to a comparison of median using a Wilcoxon text (non-parametric):
cacaomed = woodfauna$faunaproploss[woodfauna$site=="Cacao"]
woodmed = woodfauna$faunaproploss[woodfauna$site=="Forest"]

finalwoodtest = wilcox.test(cacaomed,woodmed, conf.int = T)
# tendency for median to be different, but they are not (over 0.05).



finalwoodplot =ggplot()+
  geom_point(data = woodfauna,
             mapping = aes(x = jitter(as.numeric(as.factor(site)), 0.2), 
                           y = faunaproploss))+
  scale_x_continuous(breaks = c(1,2),labels = c("Cocoa","Forest"),
                     limits = c(0.5,2.5))+
  xlab("Site")+ylab("Differences in proportion of mass lost")+
  geom_segment(mapping = aes(x = c(rep(0.75,3),rep(1.75,3)), 
                          y = c(median(cacaomed),
                                quantile(cacaomed,probs = c(0.25,0.75)),
                                median(woodmed),
                                quantile(woodmed,probs = c(0.25,0.75))),
                          xend = c(rep(1.25,3),rep(2.25,3)),
                          yend = c(median(cacaomed),
                                quantile(cacaomed,probs = c(0.25,0.75)),
                                median(woodmed),
                                quantile(woodmed,probs = c(0.25,0.75))),
                          linetype = factor(c(1,3,3,1,3,3))), col = "red")+
  theme_bw()+theme(legend.position = "none")


