# Making the method figure
# Pierrot Froment
# I have used this code to create part of the map included in the report.
# I do not include this code in the Rmarkdown workflow as I could not create
# the map entirely on R (too complicated). I therefore saved the output as an
# image and will use that image (on gimp)+ a map from QGIS for the output!

# 16/11/2024

# Library: ----
library(terra)

library(sf)
library(raster)
library(tidyverse)

# Importing files: ----

guiana = st_read("images/method_images/geoBoundaries-GUF-ADM3_simplified.shp")

guiana2 = st_as_sf(data.frame("French Guiana",st_union(guiana$geometry)))
names(guiana2)=c("shapeName","geometry")

guiana3 = st_read("images/method_images/guf_admbnda_adm2_ign.shp")

brazil = st_read("images/method_images/geoBoundaries-BRA-ADM0_simplified.shp")
suriname = st_read("images/method_images/geoBoundaries-SUR-ADM0_simplified.shp")

# Merging brazil, suriname and guiana shapefiles together:
# note: the function below returns an error. To correct if we have time.
southamerica = rbind(brazil[,c("shapeName","geometry")],
                     suriname[,c("shapeName","geometry")],
                     guiana2[,c("shapeName","geometry")])



# cropping this shapefiles to only have the area around french guiana:
southamerica = st_crop(southamerica,xmin = -55.5, 
                       ymax = 6.276353876558921, ymin = 1.6,
                       xmax = -50.9)

# Working out with the data: ----

# First map of French Guiana, brazil and Suriname:
map1 = ggplot()+geom_sf(data = southamerica, mapping = aes(fill = shapeName))+
  scale_fill_manual(values=c("lemonchiffon2","palegreen4","lemonchiffon"))+
  xlab("Longitude")+ylab("Latitude")+
  geom_text(mapping=aes(x = -51.9, y = 2.4, label = "Brazil", 
                        size= 1.5, fontface = "bold", angle = 45),
            col = "grey40")+
  geom_text(mapping=aes(x = -53.2, y = 3.6, label = "French Guiana", 
                        size= 1.5, fontface = "bold", angle = 45))+
  geom_text(mapping=aes(x = -55, y = 4.3, label = "Suriname", 
                        size= 1.5, fontface = "bold", angle = 45),
            col = "grey40")+
  geom_rect(mapping = aes(xmin = -53.0, xmax =-52.88,ymin = 5.233,
                          ymax = 5.41), color = "red", linewidth = 0.8)+
  theme_minimal()+
  theme(legend.position = "none",
                       panel.grid.major = element_blank(),
                       axis.ticks = element_blank())

# +geom_point(mapping = aes(x = -52.954156568810006, y = 5.3342159621687815,
#                            col= "red"))

# colours = colorRampPalette(c("lightgreen","darkgreen"))(22)

# Map of the study sites: near Sinamary etc...
# ggplot()+
#   geom_sf(data = guiana3, mapping = aes(fill = ADM2_REF))+
#   theme_minimal()+
#   scale_fill_manual(values=colours)+
#   theme(legend.position = "none",
#         panel.grid.major = element_blank(),
#         axis.text = element_blank())+
#   geom_point(mapping = aes(x = 283445, y = 589951,
#                            col= "red"))

# the four cardinal points I need in R for my map:
# 5.4046501159429585, -53.02823097360837
# 5.236247913203531, -52.8847024261405
# sinamary: 5.374140303492019, -52.95569086988812

# saving the map:

# ggsave("images/french_guiana.png", plot = map1, dpi = 300, width = 560,
#        height = 380, units = "px")




