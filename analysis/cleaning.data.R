# Pierrot Froment
# Importing and cleaning out data
# small code chunk to run each time to clean the 3 initial datasets


# Libraries: ----

library(tidyverse)

# Importing datasets: ----

tsbf = read.csv("data/tsbf.csv", header = T)
wood = read.csv("data/wood.csv", header = T)
lamina = read.csv("data/lamina.csv", header = T)

# Cleaning the data: ----

## TSBF ----

# removing the points in the column names
names(tsbf)=gsub("[.]", "",names(tsbf))

# creating a lexicon and change order names to code names:

lexicon = names(tsbf)[which(names(tsbf)=="snails"):which(names(tsbf)=="unknown")]

code = c("sna","slu","plan","lum.ol","tub.ol","hir","di.la","cur.la","eru.la",
         "mel.la","car.la","cuc.la","orth","blat","col","derm","het","collem",
         "thys","dipr","prot","hym","isopt","ara","sco","ps.sco","opi","aca",
         "isopd","sym","chil","dipd","neu","lep.la","amb","yl.egg","gr.egg",
         "unk")
# are each code unique: 
# length(unique(code)) == length(code)

lexicon = data.frame(cbind(lexicon,code))
names(lexicon) = c("order","code")
rm(code)

# replacing order name by their code:
names(tsbf)=ifelse(is.na(lexicon$code[match(names(tsbf),lexicon$order)]),
       names(tsbf),
       lexicon$code[match(names(tsbf),lexicon$order)])


# Uniformisation of factors levels across the different data.frames. Here,
# I choose site (2 levels, Forest-cacao),plot (the 2 plots in the cacao and
# forest site), replicat (the five points in each plots). I also correct for 
# dates' formats.

tsbf %>% rename(plot = site,
                site = parcellePCPF, layer = coucheL_10_20)%>%
  mutate(site = recode(site, "SAF"="Cacao","PF"="Forest"),
         date = recode(date, "17/09/2024"=17092024, "18/09/2024"=18092024,
                       "17092024" = 17092024),
         layer = case_when(layer == "Oct-20"~"20",
                           layer == "0-10"~"10",
                           .default = layer),
         ID = paste(site,plot,replicat,layer, sep = "-"))->tsbf



## wood: ----
wood %>% rename_with(~c("ID","site","plot","replicat","mesh","stickreplicat",
                  "startwht","endwht"))%>%
  mutate(ID = paste(site,plot,replicat,mesh,sep = "-"))-> wood


## lamina: ----
lamina %>% rename_with(~c("ID","site","plot","replicat","laminalab","hole"))%>%
  mutate(site = case_when(site == "PC"~"Cacao",
                          site == "PF"~"Forest"),
         ID= paste(site,plot, replicat,laminalab,hole, sep = "-")) -> lamina



