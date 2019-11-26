#I report this script for your ease, if you want to reproduce the analyses on the vaccination cohorts, this is the example for the influenza vaccination cohort (Herati, R. S. et al. Successive annual influenza vaccination induces a recurrent oligoclonotypic memory response in circulating T follicular helper cells. Sci Immunol 2, doi:10.1126/sciimmunol.aag2152 (2017)).


rm(list=ls())

library(rms)
library(dplyr)
library(data.table)
library(base)


setwd("/...Flu_jab/working_files")
#merge the sorted subsets for Day 0 and Day 7

files <- list.files(pattern = "*.tsv", recursive= FALSE)
#files.bas<-list.files(pattern = "*Day0.tsv", recursive= FALSE)
#files.post<-list.files(pattern = "*Day7.tsv", recursive= FALSE)

FS1415106_0<-list.files(pattern = glob2rx("FS1415106_ICOS*Day0.tsv"), recursive=FALSE)
t_FS1415106_0<-lapply(FS1415106_0, read.table, header=TRUE, sep="\t")


combined<-do.call(rbind, t_FS1415106_0)
#FS1415106_BL<-rbindlist(FS1415106_0, use.names=TRUE, fill=FALSE, )

#check the number of rows
describe(combined[,1])
describe(t_FS1415106_0[[1]][,1])
describe(t_FS1415106_0[[2]][,1])

#try to build a loop to group automatically
samples<-c("FS1314999","FS1415100","FS1415106")
for (i in samples) {
	filename<-paste(i, "_ICOS*Day0.tsv", sep="")
	x_0<-list.files(pattern = glob2rx(filename), recursive=FALSE)
	print(x_0)
	t_0<-lapply(x_0, read.table, header=TRUE, sep="\t")
	combined_x<-do.call(rbind, t_0)
	form<-sprintf("%s_BL.tsv", i)
	write.table(combined_x, file=form, quote=FALSE, sep="\t", col.names=NA)
}




setwd("/...Flu_jab/raw_data")
samples<-c("FS1314999","FS1415100","FS1415106","FS1415108","FS1415117","FS1415999","FS1516100","FS1516101","FS1516102","FS1516106","FS1516108","FS1516110","FS1516111","FS1516112","FS1516113","FS1516114","FS1516117","FS1516999")
for (i in samples) {
	setwd("/...Flu_jab/raw_data")
	filename<-paste(i, "_ICOS*Day0.tsv", sep="")
	x_0<-list.files(pattern = glob2rx(filename), recursive=FALSE)
	print(x_0)
	t_0<-lapply(x_0, read.table, header=TRUE, sep="\t")
	combined_x<-do.call(rbind, t_0)
	form<-sprintf("%s_BL_TCRB.tsv", i)
	setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/TCR/TCR_validation/TCR_viral/Flu_jab/working_files")
	write.table(combined_x, file=form, quote=FALSE, sep="\t", col.names=NA)
}

for (i in samples) {
	setwd("/...Flu_jab/raw_data")
	filename<-paste(i, "_ICOS*Day7.tsv", sep="")
	x_0<-list.files(pattern = glob2rx(filename), recursive=FALSE)
	print(x_0)
	t_0<-lapply(x_0, read.table, header=TRUE, sep="\t")
	combined_x<-do.call(rbind, t_0)
	form<-sprintf("%s_D7_TCRB.tsv", i)
	setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/TCR/TCR_validation/TCR_viral/Flu_jab/working_files")
	write.table(combined_x, file=form, quote=FALSE, sep="\t", col.names=NA)
}

#manually delete empty files and not coupled files


#calculate Renyi
#I used the script from Spreafico, R. et al. A circulating reservoir of pathogenic-like CD4+ T cells shares a genetic and phenotypic signature with the inflamed synovial micro-environment. Ann Rheum Dis 75, 459-465, doi:10.1136/annrheumdis-2014-206226 (2016). Please cite that paper, they have the credit for writing the script!

 
 
rm(list=ls())
#rename column name with altered characters 
setwd("/.../Flu_jab/working_files")
#Gini coefficient calculation failed, column 4 and 5 names were messed up
#
files <- list.files(pattern = "*.tsv", recursive= FALSE)
mock<-read.table("FS1415106_BL_TCRB.tsv", header=TRUE, sep="\t")
head(mock)
head(mock[,4])
names(mock)[4]<-"count (templates/reads)"

filez<-lapply(files, read.table, header=TRUE, sep="\t")
names(filez)[4]<-"count (templates/reads)"
names(filez)[4]
head(filez[[1]])

lapply(filez,head())

for (i in files) {
	file<-read.table(i, header=TRUE, sep="\t")
	names(file)[4]<-"count (templates/reads)"
	names(file)[5]<-"frequencyCount (%)"
	write.table(file, file=i, quote=FALSE, sep="\t", col.names=NA)
		}






 #calculate Gini
setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/TCR/TCR_validation/TCR_viral/Flu_jab/working_files")
library(LymphoSeq)
file.list<-readImmunoSeq(path=".")
warnings()
#clonality(file.list=file.list)
x<-clonality(file.list=file.list)
write.csv(x, file = "validation_FluJab_PBMC.csv")


 
 
#put the indexes together
rm(list=ls())

setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/TCR/TCR_validation/TCR_viral/Flu_jab/working_files/results")

file_list<-list.files(path=".", pattern=".txt")

dataset<-data.frame()
for (file in file_list){
		temp_dataset<-read.delim(file, header=TRUE, sep="")
		dataset<-rbind(dataset, temp_dataset)
}

write.csv(file="Renyi_draft.csv", dataset)
y<-read.table(file= "Renyi_draft.csv", header=TRUE, sep=",")
describe(y)
j<-data.frame(matrix(nrow=32))
j$cycle<-y$Renyi
j$Renyi_index<-y$Index
head(j)
j[,1]<-j$code

names<-as.data.frame(file_list)
names$file_list<-gsub(".txt","", file_list)
write.csv(file="Renyi_names.csv", file_list)

j$names[j$cycle=="BL"]<-names$file_list
j$names[j$cycle=="D7"]<-names$file_list
head(j)
write.csv(file="Renyi_indexes.csv", j)

j$names<-paste(j$names,j$cycle, sep="_")
j$names<-paste(j$names,"TCRB", sep="_")

setwd("/.../TCR_viral/Flu_jab/working_files")


k<-read.table(file= "validation_FluJab_PBMC.csv", header=TRUE, sep=",")
x<-merge(j, k, by.x="names", by.y="samples")
write.csv(file="FluJab_indexes_complete_PBMC.csv", x)


#######


