
# GOAL: find out which genes from the integrations studies were NOT
#on the in the lists of probes(rectal and ex vivo) that were considered
#detectable on the array used in  the tenofovir study.



library(dplyr)
library(HGNChelper)
################################# EX VIVO #################################

#read in data, files sent from LMF
detectableExVivo<-read.table("exvivo_TNF_HVEexc;udedNoExpression.txt", sep="\t")

#remove first line which has the column names repeated
detectableExVivo<-detectableExVivo[2:3188,]
colnames(detectableExVivo)<-c("TargetID","ProbeID")
#the column with the symbols is called TargetID. There are repeats b/c there are
#multiple probes per gene

detectableExVivoTargetIDs<-unique(as.character(detectableExVivo$TargetID))
#HGNC checking
detectableExVivoHgncCheck<-checkGeneSymbols(detectableExVivoTargetIDs,
                            unmapped.as.na = FALSE)

detectableExVivoTargetIDs<-(detectableExVivoHgncCheck$Suggested.Symbol)


################################ RECTAL ####################################

#read in rectal data

detectableRectal<-read.table("MTN007_rectal excludedNoExpression.txt",sep="\t")

detectableRectal<-detectableRectal[2:13567,]

names(detectableRectal)[1]<-"TargetID"

detectableRectalTargetIDs<-unique(as.character(detectableRectal$TargetID))
#HGNC checking
detectableRectalHgncCheck<-checkGeneSymbols(detectableRectalTargetIDs,
                                            unmapped.as.na = FALSE)

detectableRectalTargetIDs<-(detectableRectalHgncCheck$Suggested.Symbol)

##################### MALDERELLI COMPARISON  ##############################
#read in the Malderelli data
load("MalderelliData.formatted.likeSCRIData.Rda")

#using HGNChelper to check symbols and provide corrections
hgncCheck<-checkGeneSymbols(MalderelliData.formatted.likeSCRIData$Gene,
                            unmapped.as.na = FALSE)


#keep the correct names

#new gene list for comparison has all correct symbols
MalderelliIDs<-unique(hgncCheck$Suggested.Symbol)


# EX VIVO
#Which symbols are in Mald but not the ex  vivo
MalderelliNotExVivo<-MalderelliIDs[!(MalderelliIDs %in% detectableExVivoTargetIDs)]

save(MalderelliNotExVivo,file = "MalderelliNotExVivo.Rda")

#RECTAL

MalderelliNotRectal<-MalderelliIDs[!(MalderelliIDs %in% detectableRectalTargetIDs)]

save(MalderelliNotRectal,file="MalderelliNotRectal.Rda")

############################ WANG COMPARISON ##############################
#EX VIVO
load("found_genes_wang.Rda")

WangIDs<-unique(found_genes_wang)

WangHgncCheck<-checkGeneSymbols(WangIDs,unmapped.as.na = FALSE)

WangIDs<-WangHgncCheck$Suggested.Symbol

WangNotExVivo<-WangIDs[!(WangIDs %in% detectableExVivoTargetIDs)]

save(WangNotExVivo,file="WangNotExVivo.Rda")

#RECTAL

WangNotRectal<-WangIDs[!(WangIDs %in% detectableRectalTargetIDs)]

save(WangNotRectal,file="WangNotRectal.Rda")


#################### COHN  COMPARISON ##################################
library(stringr)
CohnData<-read.csv("Cohn et al Integration list.csv")

#select just the symbol.isoform column
symbols.isoforms<-dplyr::select(CohnData, Symbol.Isoform)

#symbol and isoform are separated by |
#need to escape the |
symbols.isoforms<-as.data.frame(str_split_fixed(symbols.isoforms$Symbol.Isoform,"\\|",2))


CohnHgncCheck<-checkGeneSymbols(symbols.isoforms$V1,
                            unmapped.as.na = FALSE)

CohnIDs<-unique(CohnHgncCheck$Suggested.Symbol)
# EX VIVO
CohnNotExVivo<-CohnIDs[!(CohnIDs %in% detectableExVivoTargetIDs)]

save(CohnNotExVivo, file="CohnNotExVivo.Rda")

# RECTAL

CohnNotRectal<-CohnIDs[!(CohnIDs %in% detectableRectalTargetIDs)]

save(CohnNotRectal, file="CohnNotRectal.Rda")
##################### SCRI COMPARISON ######################################

SCRIgenes<-read.table("SCRIGenesInGO.tsv", sep="\t")


#need to split up the go id and the gene symbol
SCRIgenes$V1<-as.character(SCRIgenes$V1)
SCRIgenes<-as.data.frame(str_split_fixed(SCRIgenes$V1," ",2))

SCRIsymbols<-as.data.frame(SCRIgenes$V2)

names(SCRIsymbols)[1]<-"TargetID"

SCRIhgncCheck<-checkGeneSymbols(SCRIsymbols$TargetID,
                                unmapped.as.na = FALSE)
SCRIIDs<-unique(SCRIhgncCheck$Suggested.Symbol)


# EX VIVO
SCRINotExVivo<-SCRIIDs[!(SCRIIDs %in% detectableExVivoTargetIDs)]

save(SCRINotExVivo,file="SCRINotExVivo.Rda")

# RECTAL

SCRINotRectal<-SCRIIDs[!(SCRIIDs %in% detectableRectalTargetIDs)]

save(SCRINotRectal,file="SCRINotRectal.Rda")

