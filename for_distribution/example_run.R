#CYCLeR is a tool uses BAM file as an inpit and produces GTF, FASTA and/or flat file with circRNA annotation  
#needed packages
library("tidyverse")
library("SGSeq")
library(GenomicAlignments)
library(Rsamtools)
library("stringr")
library(igraph)
library(DEXSeq)
library(polyester)
#those packages need to be selected based on the model organim
library("BSgenome.Dmelanogaster.UCSC.dm6")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
#the script requires: samples table(accouding to user manual), BSJ files, BAM file paths
###############################################################
#first we need a workign set of BSJs to start the work
#if you have a set of BSJ that you are confident in, you can skip this step 
###############################################################
#load the BSJ files
BSJ_files_prefix<-"/home/sstefan/data/test_sim/CIRI/ciri_"
bam_file_prefix<-"/home/sstefan/data/test_sim/"
filenames<-c("sample1_75","sample2_75","sample3_75","sample4_75")
BSJ_files_ciri<-paste0(BSJ_files_prefix,filenames)
bam_files<-paste0(bam_file_prefix,filenames,"_sorted.bam")
#mark the samples control and enriched or bare the consequences  
sample_table<-data.frame(filenames,c("control","control","enriched","enriched"),bam_files,stringsAsFactors = F)
colnames(sample_table)<-c("sample_name","treatment","file_bam")
si<- DataFrame(sample_table[,c("sample_name","file_bam")])
si$file_bam <-BamFileList(si$file_bam, asMates = T)
#this holds all the needed info of the bam files for downstream processing 
sc <- getBamInfo(si)
sample_table$lib_size<-sc@listData$lib_size
sample_table$read_len<-sc@listData$read_length
# the pipeline needs a set of BSJ; anyone can provide a set of BSJ; I can suggest selecting the set of BSJ this way
#we need to process the raw BSJ files 
#the processing is just a quick check of enrichment
ciri_table<-parse.files(BSJ_files_ciri,file_path="","CIRI")
colnames(ciri_table)<-c("circ_id",sample_table$sample_name)
ciri_bsjs<-process.BSJs(ciri_table,sample_table)
# i would suggest combine the output of pipelines using different mapping tools
#CE_files<-gsub("ciri_","",BSJ_files_ciri)
BSJ_files_prefix_CE<-"/home/sstefan/data/test_sim/CE/CE_"
ce_table<-parse.files(sample_table$sample_name,BSJ_files_prefix_CE,"CE")
colnames(ce_table)<-c("circ_id",sample_table$sample_name)
ce_bsjs<-process.BSJs(ce_table,sample_table)
#we need to unify the results from the BSJ identification and counting 
table_circ<-combine.two.BSJ.tables(ce_bsjs,ciri_bsjs)
#further downstream we need just the mean values for enriched samples  
table_circ<-table_circ[,c("chr","start","end","meanRR")]
colnames(table_circ)<-c("chr","start","end","count")
#combine
BSJ_set<-union(ciri_bsjs$circ_id,ce_bsjs$circ_id)
BSJ_set<-BSJ_set[!grepl("caffold",BSJ_set)]#just in case
BSJ_set<-BSJ_set[!grepl("mitochondrion",BSJ_set)]
#converting the BSJ set into a GRanges object
BSJ_gr<-make.BSJ.gr(BSJ_set)




###################################################
#With this we have finished the preprocessing of the input 
#now we need to process the BAM files
####################################################
#get the gene/transcript info 
#restoreSeqlevels(txdb)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
restoreSeqlevels(txdb)
txdb <- keepSeqlevels(txdb, c("chr2L","chr2R","chr3R","chr3L","chr4","chrX","chrY"))
seqlevelsStyle(txdb) <- "Ensembl"
gene_ranges <- genes(txdb)
txf <- convertToTxFeatures(txdb)
#asnnotation as sg-object
sgf <- convertToSGFeatures(txf)
###################################################
# we extract the regions involving circles ; either whole genes for annotated regions or the small region around the unannotated gene circles 
samtools_prefix<-"/home/sstefan/software/samtools-1.10/bin/"
trimmed_bams<-filter.bam(BSJ_gr,sample_table,samtools_prefix)
##################################################
#time to predict the features from the bam files 
#we need to the change the sample info object to accoun for the new files, but we need to keep the rest of the information e.g. library size
sc@listData[["file_bam"]]<-trimmed_bams
#sc@listData[["sample_name"]]<-sample_table_ce$sample_name
sgfc_pred <- analyzeFeatures(sc, min_junction_count=4, beta =0.1 , min_n_sample=2,cores=1,verbose=F)
sgfc_pred <- SGSeq::annotate(sgfc_pred, txf)
#sgfc_pred_an<- analyzeFeatures(sc, features=txf ,min_junction_count=4, beta =0.2 , min_n_sample=2,cores=1,verbose=F)
#you can plot the features 
df <- plotFeatures(sgfc_pred,  geneID = "131",assay = "counts",color_novel = "red", include = "both",tx_view=F,Rowv=NA, square=T)





########################################################################
#CircRNAs splice graph creation
########################################################################
#extract BSJ-corrected splice graphs  (sg)
full_sg<-overlap.SG.BSJ(sgfc_pred,BSJ_gr) #includes linear and circular features 
#get the correct genome for sequence info 
bs_genome=Dmelanogaster
# we have made new feature set so we need to recount 
full_fc<-recount.features(full_sg,sample_table)#fc==feature counts
#removing super low coverage features
full_sg<-full_sg[rowSums(as.data.frame(full_fc[,sample_table$treatment=="enriched"]))>15]
full_fc<-full_fc[rowSums(as.data.frame(full_fc[,sample_table$treatment=="enriched"]))>15,]
circ_sg<-full_sg[full_sg%over%BSJ_gr] #includes features within BSJ enclosed region 
lin_sg<-full_sg[full_sg%outside%BSJ_gr] #includes features outside of BSJ enclosed region
#annotate the BSJ with the corresponding geneIDs
BSJ_sg<-make.BSJ.sg(circ_sg,BSJ_gr)
#RPKM calculation for exons
full_rpkm<-RPKM.calc(full_fc, full_sg, BSJ_gr, bs_genome=bs_genome , sample_table=sample_table, feature_type = "e", gc_correction = F)
lin_rpkm<-full_rpkm[full_sg%outside%BSJ_gr,]
#extracting circ specific counts 
circ_fc_adj<-full_rpkm[full_sg%over%BSJ_gr,]
depleted_exons<-find.depleted.features(circ_fc_adj,sample_table,circ_sg)
#making sure that the circ edge exons remian in the mix; they coudl be depleted in case of very low levels of the circle
edge_features<-union(full_sg@featureID[start(full_sg)%in%start(BSJ_gr)],full_sg@featureID[end(full_sg)%in%end(BSJ_gr)])
depleted_exons<-setdiff(depleted_exons,edge_features)
circ_exons<-circ_sg[!circ_sg@featureID%in%depleted_exons]
circ_exons_counts<-circ_fc_adj[!circ_sg@featureID%in%depleted_exons,]
#########################################################################
#now for junctions 
#we need to normalize the junction read counts to the exon counts 
count_matrix<-as.data.frame(counts(sgfc_pred))
count_matrix <- apply (count_matrix, c (1, 2), function (x) {(as.integer(x))})
############################################
sg_gr<-rowRanges(sgfc_pred)
sg_gr_j<-sg_gr[sg_gr@type=="J"]
circ_sg_j<-sg_gr_j[unique(queryHits(findOverlaps(sg_gr_j,BSJ_gr,type = "within")))]
count_matrix_j<-count_matrix[circ_sg_j@featureID,]
#get the relative sequences of around a junction 
seqs_j<-paste0(seqs[match(start(circ_sg_j),end(full_sg))],seqs[match(end(circ_sg_j),start(full_sg))])
#calculate the scaled read count for junction 
junc_rpkm<-RPKM.calc(count_matrix=count_matrix_j, sg=circ_sg_j, bsj_granges =  BSJ_gr, sample_table =  sample_table, feature_type = "j", gc_correction = F)
#note that 
deplted_j<-find.depleted.features(junc_rpkm,sample_table,circ_sg_j)
circ_junc<-circ_sg_j[!circ_sg_j@featureID%in%deplted_j]
circ_junc_counts<-junc_rpkm[!circ_sg_j@featureID%in%deplted_j,]
circ_junc_counts[circ_junc_counts==0]<-1
colnames(circ_junc_counts)<-sample_table$sample_name




########################################################################
#transcript assembly and output
########################################################################
qics_out_sim1_75<-transcripts.per.sample("sample3_75")
qics_out_sim2_75<-transcripts.per.sample("sample4_75")

qics_out_sim_75<-merge.qics(qics_out_sim1_75,qics_out_sim2_75)

qics_out_sim_75$chr<-sgfc_pred@rowRanges@seqnames@values[qics_out_sim_75$chr]
qics_out_sim_75$circID<-paste(1:length(qics_out_sim_75$circID),qics_out_sim_75$chr,qics_out_sim_75$start,qics_out_sim_75$end, sep = "_")
qics_out_sim75_fa<-DNAStringSet(qics_out_sim_75$seq)
names(qics_out_sim75_fa)<-qics_out_sim_75$circID

#if you have a known set of circRNA fasta the CYCLeR output can be combined with it
fasta_circ<-readDNAStringSet("...")
final_ref_fa<-merge.fasta(qics_out_sim75_fa,fasta_circ)
writeXStringSet(qics_out_sim75_fa,'...')



