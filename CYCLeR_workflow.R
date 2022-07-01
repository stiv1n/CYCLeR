library(CYCLeR)
#load BSJ files
bam_file_prefix<-system.file("extdata", package = "CYCLeR")
filenames<-c("sample1_75","sample2_75","sample3_75","sample4_75") 
BSJ_files_ciri<-paste0(bam_file_prefix,"/",filenames) 
bam_files<-paste0(bam_file_prefix,"/",filenames,".bam") 
#mark the samples control and enriched or bare the consequences   
sample_table<-data.frame(filenames,c("control","control","enriched","enriched"),bam_files,stringsAsFactors = F) 
colnames(sample_table)<-c("sample_name","treatment","file_bam")
si<- DataFrame(sample_table[,c("sample_name","file_bam")])
si$file_bam <-BamFileList(si$file_bam, asMates = F) 
#this holds all the needed info of the bam files for downstream processing
sc <- getBamInfo(si) 
sample_table$lib_size<-sc@listData$lib_size 
sample_table$read_len<-sc@listData$read_length


BSJ_files_prefix<-paste0(system.file("extdata", package = "CYCLeR"),"/ciri_")
ciri_table<-parse_files(sample_table$sample_name,BSJ_files_prefix,"CIRI")
colnames(ciri_table)<-c("circ_id", "sample1_75","sample2_75","sample3_75","sample4_75")
ciri_bsjs<-process_BSJs(ciri_table,sample_table)
# i would suggest combine the output of pipelines using different mapping tools 
BSJ_files_prefix_CE<-paste0(system.file("extdata", package = "CYCLeR"),"/CE_") 
ce_table<-parse_files(sample_table$sample_name,BSJ_files_prefix_CE,"CE") 
colnames(ce_table)<-c("circ_id", "sample1_75","sample2_75","sample3_75","sample4_75")
ce_bsjs<-process_BSJs(ce_table,sample_table) 
#we need to unify the results from the BSJ identification and counting 
table_circ<-combine_two_BSJ_tables(ce_bsjs,ciri_bsjs) 
#further downstream we need just the mean values for enriched samples  
table_circ<-table_circ[,c("chr","start","end","meanRR")]
colnames(table_circ)<-c("chr","start","end","count") 
#combine 
BSJ_set<-union(ciri_bsjs$circ_id,ce_bsjs$circ_id)
BSJ_set<-BSJ_set[!grepl("caffold",BSJ_set)]
#just in case
BSJ_set<-BSJ_set[!grepl("mitochondrion",BSJ_set)] 
############################################################### 
#converting the BSJ set into a GRanges object 
BSJ_gr<-make_BSJ_gr(BSJ_set)
###################################################
samtools_prefix<-""
trimmed_bams<-filter_bam(BSJ_gr,sample_table,samtools_prefix)
sc@listData[["file_bam"]]<-trimmed_bams
####################################################
#get the gene/transcript info; we heavily suggest users to familiarize themselves with the TxDb packages
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#restoreSeqlevels(txdb)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
txdb <- keepSeqlevels(txdb, c("chr2L","chr2R","chr3R","chr3L","chr4","chrX","chrY"))
seqlevelsStyle(txdb) <- "Ensembl"
gene_ranges <- genes(txdb)
txf <- convertToTxFeatures(txdb)
#asnnotation as sg-object
sgf <- convertToSGFeatures(txf)
###################################################
sgfc_pred <- analyzeFeatures(sc, min_junction_count=2, beta =0.1 , min_n_sample=1,cores=1,verbose=F)
sgfc_pred <- SGSeq::annotate(sgfc_pred, txf)
#extract BSJ-corrected splice graphs (sg)
full_sg<-overlap_SG_BSJ(sgfc_pred,BSJ_gr,sgf) #includes linear and circular features
# we have made new feature set so we need to recount the exons
full_fc<-recount_features(full_sg,sample_table)#fc==feature counts
# time to prepare the circular splice graph
#get the correct genome for sequence info
#requires the appropriate BSgenome library 
library(BSgenome.Dmelanogaster.UCSC.dm6)
bs_genome=Dmelanogaster
circ_sgfc<-prep_circular_sg(full_sg, full_fc,sgfc_pred, bs_genome, BSJ_gr, th=15)
qics_out1<-transcripts_per_sample(sgfc=circ_sgfc,BSJ_gr = BSJ_gr,"sample3_75")
qics_out2<-transcripts_per_sample(sgfc=circ_sgfc,BSJ_gr = BSJ_gr,"sample4_75")
qics_out_final<-merge_qics(qics_out1,qics_out2,sgfc_pred)

gtf.table<-prep_output_gtf(qics_out_final,circ_sgfc)
write.table(qics_out_final[,-9],file = "dm_circles.txt", sep = "\t",row.names = F, col.names = T,quote=F)
qics_out_fa<-DNAStringSet(qics_out_final$seq)
names(qics_out_fa)<-qics_out_final$circID

#prepping the circRNA sequences for quantification  
extended_seq<-paste0(qics_out_final$seq,substr(qics_out_final$seq,1,30),strrep("N",mean(sc@listData$frag_length[sample_table$treatment=="enriched"])))
qics_out_fa_extended<-DNAStringSet(extended_seq)
names(qics_out_fa_extended)<-qics_out_final$circID
writeXStringSet(qics_out_fa_extended,'circles_seq_extended_padded.fa')
#if you have a known set of circRNA in FASTA format the CYCLeR output can be combined with it
fasta_circ<-readDNAStringSet("...")
final_ref_fa<-merge_fasta(qics_out_fa,fasta_circ)
writeXStringSet(final_ref_fa,'...')
#the same function can be used for merging with known linear annotation for the quantification step
fasta_lin<-readDNAStringSet("...")
final_ref_fa<-merge_fasta(qics_out_fa_extended,fasta_lin)
writeXStringSet(final_ref_fa,'...')
