rm(list=ls())

library(rms)
library(mudata)
library(LymphoSeq)

setwd("/...")#set the working directory where you have saved your files

#calculate Gini coefficient using LymphoSeq package.

file.list<-readImmunoSeq(path=".")
xx<-clonality(file.list=file.list)
write.csv(xx, file = "TCR_clonality.csv")


#copy cfDNA and tot PBMC files into TCR4 subfolders


#this is to identify the shared TCR sequences between the peripheral pools and the biopsy for patient #12; the CDR3 sequences have been identified by ImRep using #12 biopsy RNA-Seq. The identity of the peripheral sample ImmunoSeq files is in the Tcell_paper_samples_summary.xls file.

aa_t<-c("CASSPQGNLYYEQYF","CSAHGLASPYEQYF","CSARGPEETQYF","CASGSGTVYEQFF", "CSAPPGNTGELFF","CSARDWLAGDTGELFF","CSLPAENIQYF", "CASSVGRSSYNEQFF", "CASSRASQPYEQYF","CASSLKITGELFF", "CSYNEQFF","CSVAKGTGGSEQYF","CASSLQGANYEQYF","CASREGTASTDTQYF","CATSDYVGPDTQYF","CSAREGSKNIQYF")

rm(list=ls())

file.list <-readImmunoSeq(path='.', recursive=FALSE)
names(file.list)
bio_circ<-searchSeq(list= file.list, sequence=aa_t, type="aminoAcid", match="partial", editDistance=0)
write.csv(bio_circ, "bio_circ.csv")

productive.aa<-productiveSeq(file.list= file.list, aggregate="aminoAcid")
top.freq<-topFreq(productive.aa=productive.aa, percent=0.001)
sequence.matrix<-seqMatrix(productive.aa=productive.aa, sequences=top.freq$aminoAcid)

x.limits<-c("3220_BL_CF", "3220_BL_PBMC", "C003220_T0_D_TCRB","C003220_T0_H_TCRB","C003220_T0_L_TCRB", "3220_C2_CF","3220_C2_PBMC", "C003220_W3_D_TCRB","C003220_W3_H_TCRB","C003220_W3_L_TCRB","3220_C4_PBMC", "3220_C4_CF")
sequence.matrix<-sequence.matrix[,c("aminoAcid", x.limits)]

cloneTrack<-cloneTrack(sequence.matrix=sequence.matrix, productive.aa=productive.aa, track=aa_t, unassigned=FALSE)
cloneTrack

cloneTrack+ggplot2::scale_y_log10()
quartz.save("cloneTrack_C003220_bio_and_subsets.pdf", type="pdf")
dev.off()

#example of how to draw the Bhattacharyya matrices used for Fig 3b. It may be useful to create a subfolder for each patient, in this example the subfolder C003220 contained the ImmunoSeq tsv files for CD8+ sorted subsets of patient C003220.

setwd(".../subsets/C003220")

file.list_3220 <-readImmunoSeq(path='.', recursive=FALSE)
names(file.list_3220)
productive.aa<-productiveSeq(file.list= file.list_3220, aggregate="aminoAcid")
A<-similarityMatrix(productive.seqs= productive.aa)
write.csv(A, "similarity_matrix_C003220.csv")
B<-bhattacharyyaMatrix(productive.seqs= productive.aa)
write.csv(B, "bhattacharyya_matrix_C003220.csv")
P<-pairwisePlot(matrix=B) 
P + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5), axis.text.y = ggplot2::element_text(size = 5))
quartz.save("bhattacharyya_matrix_C003220.pdf", type="PDF")



