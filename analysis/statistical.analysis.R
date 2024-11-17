# Soil project
# Statistical analyses



# Libraries: ----

library(lme4) # for mixed models (in case)
library(lmtest) # to check linear models assumptions.
library(car)# for type II Anovas.

library(FactoMineR) # To carry out multivariate analyses (PCA, CA, PCoA).
library(factoextra)
library(vegan)# for PCoA using Bray-Curtis distances.

# Importing data: ----

# this code will be embedded within the main Rmarkdown. Run the line below if 
# you are running it independently:
# source("analysis/cleaning.data.R")


#1: TSBF diversity analysis: ----

# list of trophic levels we have so far:
unique(unlist(str_split(lexicon$trophic, "-")))

# for the analyses, we get rid of the yellow and green eggs taxon (probably
# fertilisers) and of unknonwn taxons:
taxon_focus = lexicon$code
taxon_focus = taxon_focus[!(taxon_focus%in%c("yl.egg","gr.egg","unk"))]

CA_tsbf = CA(tsbfPA[,taxon_focus], graph = F, ncp = 6)

fviz_eig(CA_tsbf)
fviz_ca_biplot(CA_tsbf, axes = c(1,2),
        col.row=tsbfPA$site)
fviz_ca_biplot(CA_tsbf, axes = c(1,2),
               col.row=tsbfPA$layer)



#2: Wood decomposition rate analysis: ----

# here, we will model the wood mass loss (initial mass-final mass) as a function
# of mesh size (effect of mesofauna vs microfauna), site (forest effect 
# vs cacao effect) and their interaction (differential effect of mesofauna in the
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
# maybe we have to adopt a generalised linear model. 
M_wood = lm((startwht-endwht)/startwht~site*mesh, data= wood)

# Checking model assumptions: 
par(mfrow = c(2,2))
plot(M_wood)
par(mfrow = c(1,1))

hist(resid(M_wood))

shapiro.test(resid(M_wood)) # normality
bptest(M_wood) # Homogeneity
dwtest(M_wood) # independence

Anova(M_wood) # type 2 anova
summary(M_wood) # summary

