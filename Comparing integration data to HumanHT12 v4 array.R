
# GOAL: find out which genes from the integrations studies were NOT
#on the Human HT12 v4 array used in  the tenofovir study



library(dplyr)
library(HGNChelper)
#read in array data, file sent from LMF
arrayData<-read.table("exvivo_TNF_HVE_allProbes.txt", sep="\t",
                      col.names=c("TargetID","ProbeID","ENTREZ_GENE","DEFINITION"))


#remove first line which has the column names repeated
arrayData<-arrayData[2:47324,]

#the column with the symbols is called TargetID. There are repeats b/c there are
#multiple probes per gene

arrayTargetIDs<-unique(as.character(arrayData$TargetID))
#HGNC checking
arrayHgncCheck<-checkGeneSymbols(arrayTargetIDs,
                            unmapped.as.na = FALSE)
arrayTargetIDs<-arrayHgncCheck$Suggested.Symbol

##################### MALDERELLI COMPARISON ##############################
#read in the Malderelli data
load("MalderelliData.formatted.likeSCRIData.Rda")

#using HGNChelper to check symbols and provide corrections
hgncCheck<-checkGeneSymbols(MalderelliData.formatted.likeSCRIData$Gene,
                            unmapped.as.na = FALSE)


#keep the correct names

#new gene list for comparison has all correct symbols
MalderelliIDs<-unique(hgncCheck$Suggested.Symbol)

#Which symbols are in Mald but not the array
MalderelliNotArray<-MalderelliIDs[!(MalderelliIDs %in% arrayTargetIDs)]

length(MalderelliNotArray)
#[1] 13


save(MalderelliNotArray,file = "MalderelliNotOnArray.Rda")

#what if I don't hgnc check them??

notCheckedArray<-unique(as.character(arrayData$TargetID))

notCheckedMald<-unique(MalderelliData.formatted.likeSCRIData$Gene)

notCheckedOverlap<-notCheckedMald[!(notCheckedMald %in% notCheckedArray)]
length(notCheckedOverlap)
#[1] 93
#Check causes MORE overlap so there are some genes that are the same from
#each raw data set, but we can't tell that until the names are standardized 
#by HGNC


############################ WANG COMPARISON ##############################

load("found_genes_wang.Rda")

WangIDs<-unique(found_genes_wang)

WangHgncCheck<-checkGeneSymbols(WangIDs,unmapped.as.na = FALSE)

WangIDs<-WangHgncCheck$Suggested.Symbol

WangNotArray<-WangIDs[!(WangIDs %in% arrayTargetIDs)]

save(WangNotArray,file="WangNotOnArray.Rda")


#################### COHN COMPARISON ##################################
CohnData<-read.csv("Cohn et al Integration list.csv")

#select just the symbol.isoform column
symbols.isoforms<-dplyr::select(CohnData, Symbol.Isoform)

#symbol and isoform are separated by |
#need to escape the |
symbols.isoforms<-as.data.frame(str_split_fixed(symbols.isoforms$Symbol.Isoform,"\\|",2))


CohnHgncCheck<-checkGeneSymbols(symbols.isoforms$V1,
                            unmapped.as.na = FALSE)

CohnIDs<-unique(CohnHgncCheck$Suggested.Symbol)

CohnNotArray<-CohnIDs[!(CohnIDs %in% arrayTargetIDs)]

save(CohnNotArray, file="CohnNotOnArray.Rda")
