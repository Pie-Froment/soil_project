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

# replacing Oct-20 by 10-20 in the couche column.
names(tsbf)[names(tsbf)=="coucheL_10_20"]="couche"
tsbf$couche[tsbf$couche == "Oct-20"]="10-20"


## wood: ----
names(wood)=c("code","site","plot","replicat","mesh","surreplicat","startwht","endwht")

## lamina: ----




