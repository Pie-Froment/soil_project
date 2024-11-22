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
library(gtools) # because of this *** order() function that gives wrong order.

# Importing data: ----

# this code will be embedded within the main Rmarkdown. Run the line below if 
# you are running it independently:
# source("analysis/cleaning.data.R")


#1: TSBF diversity analysis: ----

# list of trophic levels we have so far:
unique(unlist(str_split(lexicon$trophic, "-")))

# for the analyses, we get rid of the yellow and green eggs taxon (probably
# fertilisers) and of unknonwn taxons + taxa with 0 observations:
taxon_focus = lexicon$code
#  taxa with 0 observations:
selec= names(colSums(tsbf[,taxon_focus])[colSums(tsbf[,taxon_focus])==0])
taxon_focus = taxon_focus[!(taxon_focus%in%c("yl.egg","gr.egg","unk", selec))]
rm(selec)

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

## C: Species abundances plot: ----
# Creation of a nice graph showing the pyramid of age abundances instead:
# data.frame:
tsbf %>% select(c("site",all_of(taxon_focus)))%>%
  group_by(site)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") %>%
  pivot_longer(cols = all_of(taxon_focus),
               names_to = "taxa",
               values_to = "abundance") -> tsbfplot

# We add a column with the taxon order:
col = colSums(tsbf[,taxon_focus])
col = sort(col, decreasing = T)
newvec = 1:29
# very dirty corde but it works (did not manage to use order())
tsbfplot$rank=newvec[match(tsbfplot$taxa,names(col))]
rm(newvec,col)

#Plotting abundances:
ggplot(data= tsbfplot, mapping = aes(x = rev(as.numeric(rank)),
                                             y = abundance, fill = site))+
  geom_bar(stat = "identity", position = "identity")+
  scale_fill_manual(values = c("#e31a1c","darkgreen"))+
  coord_flip()+
  theme_minimal()+
  theme(panel.grid.major.y= element_blank(),
        panel.grid.minor.y= element_blank(),
        axis.text.y=element_blank())

## D: Species-accumulation curve ----

# Here we create a new dataset where we summed all taxa abundance in the 
# different layers in each tsbf sites (for the accumulation curve)
tsbf %>% select(c("site","replicat",all_of(taxon_focus)))%>%
  group_by(site,replicat)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") ->
  tsbf_co 


# Create the accumulation curve using the specaccum function (method = random):
Cocoa = specaccum(tsbf_co[tsbf_co$site=="Cocoa",taxon_focus],
          method = "random", permutations = 100,
          ci = 0.95)
Forest = specaccum(tsbf_co[tsbf_co$site=="Forest",taxon_focus],
          method = "random", permutations = 100,
          ci = 0.95)

# We plot these curves in a very tedious way:

ggplot()+ # Cocoa line:
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

# deleting variables we dont need anymore: 
rm(Cocoa,Forest)

## E: Diversity profiles: ----

# For diversity profiles, we want to differentiate between the different
# layers (could be interesting), as well as between the sites:

# Creating a new dataset summing data by layers and site:

tsbf %>% select(c("site","layer",all_of(taxon_focus)))%>%
  group_by(site,layer)%>%
  summarise(across(all_of(taxon_focus), sum, na.rm = T), .groups = "drop") ->
  tsbf_la 

# Whittaker plots in the different sub-grouping (exploratory):
# All:
plot(as.AbdVector(colSums(tsbf_la[,taxon_focus])))
# Cocoa:
as.AbdVector(colSums(tsbf_la[tsbf_la$site == "Cocoa",taxon_focus]))%>%
  plot()
# Forest:
as.AbdVector(colSums(tsbf_la[tsbf_la$site == "Forest",taxon_focus]))%>%
  plot()

# maybe do a ggplot combining the two whitakker's plots together. 




#2: Wood decomposition rate analysis: ----

# here, we will model the wood mass loss (initial mass-final mass) as a function
# of mesh size (effect of mesofauna vs microfauna), site (forest effect 
# vs Cocoa effect) and their interaction (differential effect of mesofauna in the
# different sites)

# Since we have two sticks by bags and these sticks are complete pseudoreplicates 
# we merge the two sticks measurements into one average weight loss per bag.

# merging the two sticks measurements per bag into a single one:
wood %>% group_by(ID)%>%
  reframe(plot = unique(plot), replicat = unique(replicat), 
          mesh = unique(mesh),site = unique(site),
          startwht= mean(startwht),
            endwht = mean(endwht)) -> wood

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

## a: first model ----
# (here, a simple linear model with an interaction):
# Notes: since the response variable is a percentage (bounded between 0 and 100),
# maybe we have to adopt a generalised linear model (with a beta distribution). 
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

