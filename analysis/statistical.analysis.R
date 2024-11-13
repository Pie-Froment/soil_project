# Soil project
# Statistical analyses



# Libraries: ----

library(lme4)
library(lmtest)
library(car)

# Importing data: ----

source("analysis/cleaning.data.R")


#1: TSBF diversity analysis: ----




#2: Wood decomposition rate analysis: ----

# here, we will model the wood mass loss (initial mass-final mass) as a function
# of mesh size (effect of mesofauna vs microfauna), site (forest effect 
# vs cacao effect) and their interaction (differential effect of mesofauna in the
# different sites)

# Since we have two sticks by bags and these sticks are complete pseudoreplicates 
# we merge the two sticks measurements into one average weight loss per bag.

wood %>% group_by(ID)%>%
  reframe(plot = unique(plot), replicat = unique(replicat), 
          mesh = unique(mesh),site = unique(site),
          startwht= mean(startwht),
            endwht = mean(endwht)) -> wood

# Plotting the data: 

ggplot(data = wood)+ 
  geom_point(mapping = aes(x = jitter(as.numeric(as.factor(site)), 0.1), 
                           y = startwht-endwht, colour = mesh))+
  theme_bw()+
  xlim(0.5,2.5)
 
ggplot(data = wood)+ 
  geom_boxplot(mapping = aes(x = site, 
                           y = startwht-endwht, colour = mesh))+
  theme_bw()

M_wood = lm(startwht-endwht~site*mesh, data= wood)
par(mfrow = c(2,2))
plot(M_wood)
par(mfrow = c(1,1))

shapiro.test(resid(M_wood))
bptest(M_wood)
dwtest(M_wood)

summary(M_wood)
Anova(M_wood)
