# This is the CYCLeR software, comprising the files CYCLeR_functions.R, CYCLeR_manual.html, example_run.R and license_information.txt, is copyrighted 2020 by Stefan R. Stefanov and
# Irmtraud M.Meyer (irmtraud.meyer@cantab.net) and distributed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License, see
# 
# https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
# 
# For a short summary of what this license entrails, see
# 
# https://creativecommons.org/licenses/by-nc-sa/4.0/
# 
# In a nutshell:
# 
# You are free to:
# 
#     Share — copy and redistribute the material in any medium or format
#     Adapt — remix, transform, and build upon the material
#     The licensor cannot revoke these freedoms as long as you follow the license terms.
# 
# Under the following terms:
# 
#     Attribution — You must give appropriate credit, provide a link to
#     the license, and indicate if changes were made. You may do so in
#     any reasonable manner, but not in any way that suggests the
#     licensor endorses you or your use.
# 
#     NonCommercial — You may not use the material for commercial purposes.
# 
#     ShareAlike — If you remix, transform, or build upon the material,
#     you must distribute your contributions under the same license as
#     the original.
# 
#     No additional restrictions — You may not apply legal terms or
#     technological measures that legally restrict others from doing
#     anything the license permits.
# 
# Notices:
# 
#     You do not have to comply with the license for elements of the
#     material in the public domain or where your use is permitted by an
#     applicable exception or limitation.  No warranties are given. The
#     license may not give you all of the permissions necessary for your
#     intended use. For example, other rights such as publicity,
#     privacy, or moral rights may limit how you use the material.

read.CIRI<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  ciri<-read_tsv(file_name, comment = "", col_names = T,col_types = cols())
  select(ciri,(chr:`#junction_reads`),strand)%>%unite(circ_id,chr,circRNA_start,circRNA_end,strand)%>%rename(!!sample_name:=`#junction_reads`)
}
read.CE<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  ce<-read_tsv(file_name, comment = "", col_names = F,col_types = cols())
  select(ce,(X1:X3),X6,X13)%>%mutate(X2=X2+1)%>%unite(circ_id,X1,X2,X3,X6)%>%rename(!!sample_name:=X13)
}
read.tsv.file<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  df<-read_tsv(file_name, comment = "", col_names = F,col_types = cols())
  unite(df,circ_id,X1,X2,X3,X4)%>%rename(!!sample_name:=X5)
}

parse.files<-function(file_list, file_path, input_type){
  if(input_type=="CIRI"){
    func=read.CIRI
    }else if(input_type=="CE"){
      func=read.CE
      }else if(input_type=="tsv"){
        func=read.tsv.file
        }else stop("input type is exclusively: CIRI, CE or tsv")
  Reduce(function(x,y) {full_join(x,y, by="circ_id")},map(paste0(file_path,file_list),func))
 
}
process.BSJs<-function(cdf,sample_table){
  #stands for circular data frame
  cdf[is.na(cdf)] <- 0
  #cdf<-column_to_rownames(cdf,var = "circ_id")
  cdf[,sample_table$sample_name]<-(cdf[,sample_table$sample_name]/sample_table$lib_size)*10e6
  #cdf<-rownames_to_column(cdf,var = "circ_id")
  cdf$meanc<-rowMeans(cdf[,c(sample_table$sample_name[sample_table$treatment=="control"])])
  cdf$meanRR<-rowMeans(cdf[,c(sample_table$sample_name[sample_table$treatment=="enriched"])])
 # rownames(cdf[cdf$meanc<cdf$meanr,])
  cdf[cdf$meanc<cdf$meanRR,]
}
combine.two.BSJ.tables<-function(ce_bsjs,ciri_bsjs){
  #colnames(ciri_bsjs)<-colnames(ce_bsjs)
  #ciri_bsjs$circ_id<-rownames(ciri_bsjs)
  #ce_bsjs$circ_id<-rownames(ce_bsjs)
  full_table<-bind_rows(ciri_bsjs,ce_bsjs)
  full_table[is.na(full_table)]<-0
  table_circ<-full_table%>%group_by(circ_id)%>%summarize_all(max)
  #table_circ<-table_circ%>%separate(circ_id, into = paste("V", 1:4, sep = ""),sep = "_")
  table_circ<-table_circ%>%separate(circ_id, into = c("chr","start","end","strand"),sep = "_")
  table_circ
}
make.BSJ.gr<-function(BSJ_set){
  BSJ_BED<-as.data.frame(str_split(BSJ_set,"_",simplify = T),stringsAsFactors=F)
  BSJ_BED$V2<-as.numeric(BSJ_BED$V2)
  BSJ_BED$V3<-as.numeric(BSJ_BED$V3)
  BSJ_gr <- GRanges(BSJ_BED$V1, IRanges(start=BSJ_BED$V2,end=BSJ_BED$V3),
                    strand=BSJ_BED$V4, seqlengths=NULL)
  BSJ_gr@elementMetadata@listData[["gene_id"]]<-paste0("circ",seq_along(1:length(BSJ_gr)))
  BSJ_gr
}
plotRanges2 <- function(...) {
  arg_v <- c(as.list(environment()), list(...)) #argument values 
  arg_n <- as.list(match.call()) #argument names 
  t<-paste(arg_n)
  t[1]="c"
  x<-eval(as.call( parse(text=t))) # it needs to be so complicated to save the order of the ranges 
  x_l<-lengths(arg_v)
  col=rep(t[-1],x_l)
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  dat <- cbind(as.data.frame(x), bin = bins, col=col)
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(dat) + 
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = bin, ymax = bin + 0.9, fill=col)) +
    theme_bw()+scale_fill_manual(values=cbPalette)+theme(axis.text.y=element_blank(),
                                                         axis.ticks.y=element_blank())
}
overlap.SG.BSJ<-function(sgfc_pred,BSJ_gr){
  sg_gr<-rowRanges(sgfc_pred) #stands for splice graph granges
  sg_gr_e<-sg_gr[sg_gr@type=="E"]

  BSJ_sg<-SGFeatures(BSJ_gr,type=rep("E",length(BSJ_gr)),splice5p =rep(F,length(BSJ_gr)),
                     splice3p = rep(F,length(BSJ_gr)), featureID =rep(1L,length(BSJ_gr)),
                     geneID = rep(1L,length(BSJ_gr)), txName = mcols(BSJ_gr)$geneID,
                     geneName = mcols(BSJ_gr)$geneID)
  
  sg_gr_e2<-sg_gr_e[sg_gr_e%over%BSJ_sg]  #this is to remove some mapping artifacts
  sg_gr_e<-sg_gr_e[sg_gr_e@geneID%in%sg_gr_e2@geneID] #and keep only the genes that overlap with BSJ regions
  #BSJ_sg2<-GenomicRanges::reduce(BSJ_sg)
  combined<-c(BSJ_sg,sg_gr_e)
  #combined2<-c(BSJ_sg2,sg_gr_e)
  disjoint_ranges<-GenomicRanges::disjoin(combined, with.revmap=T)
  #disjoint_ranges2<-GenomicRanges::disjoin(combined2, with.revmap=T)
  disjoint_ranges<-disjoint_ranges[disjoint_ranges%over%sg_gr_e] #this is needed to remove the intronic artefacts from the split
  disjoined_sg<-SGFeatures(disjoint_ranges,type=rep("E",length(disjoint_ranges)),splice5p =rep(F,length(disjoint_ranges)),
                           splice3p = rep(F,length(disjoint_ranges)), featureID = rep(1L,length(disjoint_ranges)),
                           geneID = rep(1L,length(disjoint_ranges)), txName = mcols(disjoint_ranges)$geneID,
                           geneName = mcols(disjoint_ranges)$geneID)
  disjoined_sg@geneID<-sg_gr_e@geneID[queryHits(findOverlaps(sg_gr_e,disjoined_sg))] #this makes sure we keep geneID info in the new bins
  disjoined_sg@geneName<-sg_gr_e@geneName[queryHits(findOverlaps(sg_gr_e,disjoined_sg))] 
  disjoined_sg@txName<-sg_gr_e@txName[queryHits(findOverlaps(sg_gr_e,disjoined_sg))] 
  #disjoined_sg@featureID<-sg_gr_e@featureID[queryHits(findOverlaps(sg_gr_e,disjoined_sg))] 
  #disjoined_sg@featureID<-seq(1:length(disjoined_sg))
  disjoined_sg@featureID<-seq(max(sg_gr_e@featureID),(max(sg_gr_e@featureID)+length(disjoined_sg)-1))
  disjoined_sg
}
make.BSJ.sg<-function(circ_sg,BSJ_gr){
  BSJ_sg<-SGFeatures(BSJ_gr,type=rep("E",length(BSJ_gr)),splice5p =rep(F,length(BSJ_gr)),
                     splice3p = rep(F,length(BSJ_gr)), featureID =rep(1L,length(BSJ_gr)),
                     geneID = rep(1L,length(BSJ_gr)), txName = mcols(BSJ_gr)$geneID,
                     geneName = mcols(BSJ_gr)$geneID)
  BSJ_sg@geneID<-circ_sg@geneID[queryHits(findOverlaps(circ_sg,BSJ_sg,type = c("start")))] 
  BSJ_sg
}
filter.bam<-function(BSJ_gr,sample_table,samtools_prefix){
  #gn_circ<-gene_ranges[gene_ranges%over%BSJ_gr]
  #unannot<-BSJ_gr[BSJ_gr%outside%gene_ranges]
  #filter_ranges<-c(gn_circ,unannot)
  bed_file<-paste0(bam_file_prefix,"BSJ.bed")
  #export.bed(con = bed_file,object =  filter_ranges)
  export.bed(con = bed_file,object =  BSJ_gr)
  #now we need to use the bed file to filter the bam files; i suggest doing it directly in the command line 
  list_file<-paste0(bam_file_prefix,"list")
  write(sample_table$sample_name,list_file)
  print("Files are trimming")
  trim_command<-paste(paste0(samtools_prefix,"samtools view -b -L ",bed_file," ", sample_table$file_bam," > ",bam_file_prefix,sample_table$sample_name,"_sorted_trimmed.bam"), collapse = "; ")
  system(trim_command)
  print("Files are indexing")
  index_command<-paste(paste0(samtools_prefix,"samtools index ",bam_file_prefix,sample_table$sample_name,"_sorted_trimmed.bam"), collapse = "; ")
  system(index_command)
  trimmed_bams<-BamFileList(gsub("sorted", "sorted_trimmed", sample_table$file_bam), asMates = F)
  trimmed_bams
}
recount.features<-function(full_sg,sample_table){
  feature_info<-as.data.frame(full_sg,stringsAsFactors = F)
  # we make a feature reference for counting
  saf<-feature_info[feature_info$type=="E",c("featureID","seqnames","start","end","strand")]
  colnames(saf)<- c("GeneID",		"Chr",	"Start",	"End",	"Strand")
  #saf$GeneID<-c(1,2,3,4)
  saf.counts<-Rsubread::featureCounts(gsub("sorted", "sorted_trimmed", sample_table$file_bam), annot.ext=saf, isPairedEnd=F,  countChimericFragments= T, 
                                      countMultiMappingReads =T ,juncCounts=T,allowMultiOverlap=T)
  saf.fc<-saf.counts$counts
  saf.fc
}
get.seqs<-function(full_sg, bs_genome=Dmelanogaster){
  seqs<-unname(BSgenome::getSeq(bs_genome, paste0("chr",full_sg@seqnames), start = start(full_sg), end = end(full_sg), strand=strand(full_sg),as.character=T))
  seqs
}
get.seqs.corrected<-function(full_sg, bs_genome=Dmelanogaster){
  seqs<-get.seqs(full_sg)
  short_ranges<-full_sg[width(full_sg)<35]
  long_ranges<-full_sg[!width(full_sg)<35]
  pair_ranges<-long_ranges[nearest(short_ranges,long_ranges)]
  corrected.seqs<-paste0(get.seqs(short_ranges),get.seqs(pair_ranges))
  seqs[width(full_sg)<35]<-corrected.seqs
  seqs
}
get.seqs.positive<-function(full_sg, bs_genome=Dmelanogaster){
  strand(full_sg)<-"+"
  get.seqs(full_sg)
}
#this one needs work to add the true/false if-s and test it; also figure out how to save the model 
RPKM.calc<-function(count_matrix, sg ,bsj_granges, bs_genome, sample_table,feature_type, fsj_overhang=3, bsj_overhang=15, eff_length_correction=T, gc_correction=T){
  bsj_overhang_adj=bsj_overhang-fsj_overhang
  eff_lengths<- width(sg)
  #the scaling is capped based on the  read length, because in case if something maps it is at least a read size 
  eff_length_limit<-max(sample_table$read_len[sample_table$treatment=="enriched"])
  eff_lengths[eff_lengths<eff_length_limit]<-eff_length_limit
  if(feature_type=="j") {eff_lengths<-rep(1,length(sg))}
  if(feature_type=="e") {
    if(eff_length_correction){
      #adjustment
      #min required eff. length adjustment is based on the overhang and then we build up upon that 
      #FSJ fix 
      eff_lengths<-eff_lengths-(2*fsj_overhang)
      #BSJ fix
      eff_lengths[start(sg)%in%start(bsj_granges)]<-eff_lengths[start(sg)%in%start(bsj_granges)]-bsj_overhang_adj
      eff_lengths[end(sg)%in%end(bsj_granges)]<-eff_lengths[end(sg)%in%end(bsj_granges)]-bsj_overhang_adj
      
    }
  }  
  ############################# GC-content adhustment 
  #calculating CG content
  #seqs<-unname(BSgenome::getSeq(bs_genome, paste0("chr",full_sg@seqnames), start = start(full_sg), end = end(full_sg), strand=strand(full_sg),as.character=T))
  #fix that later
  #full_gr<-GRanges(full_sg)
  #seqlevelsStyle(full_gr) <- "UCSC"
  #seqs<-BSgenome::getSeq(bs_genome, full_gr)
  if(gc_correction){
    seqs<-get.seqs.corrected(sg)
    GC_content <- (str_count(seqs, "G") + str_count(seqs, "C")) / str_length(seqs) * 100 
    #again a rough fix; i am capping it within those ranges they need to be tested 
    GC_content[GC_content<35]<-35
    GC_content[GC_content>75]<-75
    # adjusting scaled counts for CG content
    #this section needs to be reworked; may be save the model as rda; in the packege folder 
    data("loessfit7")
    fit.df<-as.data.frame(loessfit7)
    loess_mod <- loess(y~x, data = fit.df, control=loess.control(surface="direct"), span=0.3)
    shifts<-predict(loess_mod, GC_content/100)
    shifts<-2^shifts
    eff_lengths_gc<-eff_lengths*shifts
  }
  else{eff_lengths_gc<-eff_lengths}
  #here is with shifts based on the log2
  #shifts <- predict(loess_mod,count_per_feature_exons_circ$GC_content/100)
  colnames(count_matrix)<-sample_table$sample_name
  full_fc_adj<-sweep(count_matrix,2,sample_table$lib_size/10^9, FUN="/")
  full_fc_adj<-sweep(full_fc_adj,1,eff_lengths_gc, FUN="/")
  #full_fc_adj<-sweep(full_fc_adj,1,shifts, FUN="/")
  full_fc_adj<-apply (full_fc_adj, c (1, 2), function (x) {
    (as.integer(x))
  })
  full_fc_adj
}
find.depleted.features<-function(circ_fc_adj,sample_table,circ_sg, test="DEX"){#may need to add feature type
  if(length(sample_table$sample_name)<4) {test="comparison"}
  if(test=="comparison"){
    cdf<-as.data.frame(circ_fc_adj)
    cdf$meanc<-rowMeans(as.data.frame(cdf[,c(sample_table$sample_name[sample_table$treatment=="control"])]))
    cdf$meanRR<-rowMeans(as.data.frame(cdf[,c(sample_table$sample_name[sample_table$treatment=="enriched"])]))
    cdf[cdf$meanc>cdf$meanRR,]
    rownames(cdf)
  }
  else{
    feature_dex<- DEXSeqDataSet(circ_fc_adj, sample_table, 
                                design= ~ sample + exon + treatment:exon, 
                                as.character(circ_sg@featureID),  as.character(circ_sg@geneID), 
                                featureRanges=GRanges(circ_sg), 
                                transcripts=circ_sg@txName, 
                                alternativeCountData=NULL)
    sizeFactors(feature_dex)<-1
    #if(feature_type=="e") {sizeFactors(feature_dex)<-1}
    #if(feature_type=="j") {feature_dex = estimateSizeFactors(feature_dex)}
    feature_dex = estimateDispersions( feature_dex, fitType="local") # due to the small number of features; we use the local fit type 
    feature_dex = testForDEU( feature_dex)
    feature_dex = estimateExonFoldChanges( feature_dex,fitExpToVar="treatment")
    feature_dex_res = as.data.frame(DEXSeqResults( feature_dex ))
    #filtering of the enriched features
    feature_dex_res_depleted<-feature_dex_res$featureID[feature_dex_res$log2fold_enriched_control<(-1)&feature_dex_res$padj<=0.005]
    #feature_dex_res_depleted_extra<-feature_dex_res$featureID[feature_dex_res$log2fold_enriched_control<(-1)]
    #feature_dex<-union(feature_dex_res_depleted,feature_dex_res_depleted_extra)
    feature_dex_res_depleted<-feature_dex_res_depleted[!is.na(feature_dex_res_depleted)]
    feature_dex_res_depleted #exports the feature id; not to be confused with the index 
  }
}
prep.output<-function(transcript_features,circ_exons){
  transcript_features<-rev(transcript_features)
  transcript_features<-transcript_features[!transcript_features==""]
  csp.results.col<-c("circID","chr","start","end","strand","gene","exon_starts","exon_ends","seq")
  csp.results<-data.frame(matrix(0,nrow=length(transcript_features),ncol=length(csp.results.col)))
  colnames(csp.results)<-csp.results.col
  seqs<-get.seqs.positive(circ_exons)
  for(i in 1:length(transcript_features)){
    #str_split(transcript_features[i],"_",simplify = T)
    temp.features<-unlist(strsplit(transcript_features[i],"_"))
    temp.feature.table<-as.data.frame(circ_exons[featureID(circ_exons)%in%temp.features])
    #the reasonable output part
    csp.results[i,2]<-temp.feature.table$seqnames[1]
    csp.results[i,3]<-temp.feature.table$start[1]
    csp.results[i,4]<-tail(temp.feature.table$end, n=1)
    csp.results[i,5]<-as.character(temp.feature.table$strand[1])
    csp.results[i,6]<-ifelse(length(unlist(temp.feature.table$geneName[1]))==0,temp.feature.table$geneID[1],unlist(temp.feature.table$geneName[1]))
    csp.results[i,7]<-paste(temp.feature.table$start,collapse = ",")
    csp.results[i,8]<-paste(temp.feature.table$end,collapse = ",")
    csp.results[i,1]<-paste(c(i,csp.results[i,2],csp.results[i,3],csp.results[i,4]),collapse ="_")
    csp.results[i,9]<-paste(seqs[featureID(circ_exons)%in%temp.features],collapse = "")
    if(temp.feature.table$strand[1]=="-"){ csp.results[i,9]<-as.character(reverseComplement(DNAString(csp.results[i,9])))}
  }
  csp.results
}
prep.output.gtf<-function(qics_out,circ_exons,annot_list=NULL){
  #check if annotation is given 
  if(is.null(annot_list)){
    mapped_ids<-unique(qics_out$gene)
    names(mapped_ids)<-mapped_ids
  }else{mapped_ids<-annot_list}
  csp.results<-qics_out
  #the unreasonable, weird gtf part
  temp.gene.id<-""
  gtf.results.col<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
  gtf.results<-data.frame(matrix(ncol=length(gtf.results.col)))
  colnames(gtf.results)<-gtf.results.col
  for(i in 1:length(csp.results$circID)){
    #a<-lapply(1:length(csp.results$circID),function(i){
    #gene
    geneid<-geneID(circ_exons)[start(circ_exons)==csp.results[i,3]]
    temp.gene.name<-ifelse(is.na(mapped_ids[csp.results[i,6]]),paste0("NA_",geneid),mapped_ids[csp.results[i,6]])
    if(temp.gene.id!=geneid){
      temp.gene.id<-csp.results[i,1]
      temp.gtf.gene<-data.frame(matrix(0,nrow=1,ncol=length(gtf.results.col)))
      colnames(temp.gtf.gene)<-gtf.results.col
      temp.gtf.gene$seqname<-csp.results[i,2]
      temp.gtf.gene$feature<-"gene"
      temp.gtf.gene$start<-min(start(circ_exons[geneID(circ_exons)==geneid]))
      temp.gtf.gene$end<-max(end(circ_exons[geneID(circ_exons)==geneid]))
      temp.gtf.gene$strand<-csp.results[i,5]
      temp.gtf.gene$attribute<-paste0('gene_id "circ',csp.results[i,6],'"; gene_name "',temp.gene.name,'"; gene_source "CSP"; gene_biotype "circRNA";')
      gtf.results<-rbind(gtf.results,temp.gtf.gene)   
    } 
    #transcript
    temp.gtf.trascript<-data.frame(matrix(0,nrow=1,ncol=length(gtf.results.col)))
    colnames(temp.gtf.trascript)<-gtf.results.col
    temp.gtf.trascript$feature<-"transcript"
    temp.gtf.trascript$seqname<-csp.results[i,2]
    temp.gtf.trascript$start<-min(str_split(csp.results[i,7],",",simplify =T))
    temp.gtf.trascript$end<-max(str_split(csp.results[i,8],",",simplify =T))
    temp.gtf.trascript$strand<-csp.results[i,5]
    temp.gtf.trascript$attribute<-paste0('gene_id "circ',csp.results[i,6],'"; transcript_id "',csp.results[i,1],'"; gene_name "',temp.gene.name,'"; gene_source "CSP"; gene_biotype "circRNA"; transcript_name "',csp.results[i,1],'"; transcript_source "CSP"; transcript_biotype "circRNA";')
    gtf.results<-rbind(gtf.results,temp.gtf.trascript)
    #exons
    temp.gtf.exons<-data.frame(matrix(0,nrow=length(str_split(csp.results[i,7],",",simplify =F)[[1]]),ncol=length(gtf.results.col)))
    colnames(temp.gtf.exons)<-gtf.results.col
    temp.gtf.exons$feature<-"exon"
    temp.gtf.exons$seqname<-csp.results[i,2]
    temp.gtf.exons$start<-str_split(csp.results[i,7],",",simplify =F)[[1]]
    temp.gtf.exons$end<-str_split(csp.results[i,8],",",simplify =F)[[1]]
    temp.gtf.exons$strand<-csp.results[i,5]
    if(temp.gtf.exons$strand[1]=="-") {temp.gtf.exons<-arrange(temp.gtf.exons, rev(rownames(temp.gtf.exons)))}
    temp.gtf.exons$attribute<-paste0('gene_id "circ',csp.results[i,6],'"; transcript_id "',csp.results[i,1],'"; exon_number "',seq(1:length(str_split(csp.results[i,7],",",simplify =F)[[1]])),'"; gene_name "',temp.gene.name,'"; gene_source "CSP"; gene_biotype "circRNA"; transcript_name "',csp.results[i,1],'"; transcript_source "CSP"; transcript_biotype "circRNA"; exon_id "',csp.results[i,1],'-E',seq(1:length(str_split(csp.results[i,7],",",simplify =F)[[1]])),'"')
    gtf.results<-rbind(gtf.results,temp.gtf.exons)
  }
  gtf.results<-gtf.results[-1,]
  gtf.results$source<-"QICS"
  gtf.results$score<-"."
  gtf.results$frame<-"."  
  gtf.results
}
assemble.transcripts.per.sample<-function(i){
  #the junction counts do not scale well to the exon counts; so the rpkm of the junctions is scaled to the rpkm of the exons; we use median ratios scaling 
  #first i am calulating a pseudo count for a junction out of the 2 surrounding exons
  pseudo_ref<-(circ_exons_counts[match(start(circ_junc),end(circ_exons)),i]+circ_exons_counts[match(end(circ_junc),start(circ_exons)),i])/2
  #there are NAs introduced from junctions connecting deplted exons; we can remove them or set them to 0; settign them to 0 would shift the median 
  #pseudo_ref<-pseudo_ref[complete.cases(pseudo_ref)]
  ratios<-circ_junc_counts[,i]/pseudo_ref
  ratios<-ratios[complete.cases(ratios)]
  ratios<-ratios[ratios>1]
  #ratios[is.infinite(ratios)]<-0
  #now we need to calculate the size factors 
  size_factor<-median(ratios)
  #we apply the size factor 
  circ_junc_counts[,i]<-round(circ_junc_counts[,i]/size_factor)
  
  #########################################
  #start recontruction here
  #we need to create the splice graph 
  transcript_features<-c()
  #it is OK to use for loops here since there should not be more than 5000 iterations 
  for(j in unique(circ_exons@geneID)){
    #j=11 #for now
    #print(j)
    sg_exons_df<-cbind(as.data.frame(circ_exons[circ_exons@geneID==j]), circ_exons_counts[circ_exons@geneID==j,i])
    colnames(sg_exons_df)[13]<-"counts"
    sg_junc_df<-cbind(as.data.frame(circ_junc[circ_junc@geneID==j]), circ_junc_counts[circ_junc@geneID==j,i])
    colnames(sg_junc_df)[13]<-"counts"
    #this is needed for the removal of the linear residual counts
    ###avg_weight_exons_lin<-round(mean(lin_rpkm[lin_sg@geneID==j,i]))
    #for now avg_weight_exons_lin=0 #for now 
    #ths is needed for the calcualtion of the threshold to cut an exon toghether with another 
    avg_weight_exons_circ<-mean(sg_exons_df$counts)
    ###min_weight_exon_circ<-min(sg_exons_df$counts)
    #gettign the retained introns 
    #sg_exons_ri<-sg_exons_df[sg_exons_df$start%in%(sg_exons_df$end+1),]
    #sg_exons_ri<-sg_exons_ri[sg_exons_ri$end%in%(sg_exons_ri$start-1),]
    #feature_id_ri<-sg_exons_ri$featureID
    #feature_id_ri<-""
    #sg_exons_df$counts[!sg_exons_df$featureID%in%feature_id_ri]<-(sg_exons_df$counts-avg_weight_exons_lin)[!sg_exons_df$featureID%in%feature_id_ri]
    #creation of the graph
    sg_exons_df<-sg_exons_df[,c("featureID","start","end","counts")]
    sg_junc_df<-sg_junc_df[,c("featureID","start","end","counts")]
    sg_df<-rbind(sg_exons_df,sg_junc_df)
    #adjusting the counts; avoiding retained introns 
    ###sg_df$counts[!sg_df$featureID%in%feature_id_ri]<-(sg_df$counts-avg_weight_exons_lin)[!sg_df$featureID%in%feature_id_ri]
    ###sg_df$counts[sg_df$counts<min_weight_exon_circ]<-min_weight_exon_circ
    sg_df2 <- sg_df[order(sg_df$start),]
    
    #----------------------------------
    # fix of Retain Introns and 5/3 prime splicing
    # we make a data frame of junctions to connect the hanging pieces in the graph 
    fix_df<-NULL
    fix_ends<-sg_exons_df[sg_exons_df$end%in%sg_exons_df$start-1,"start"]
    if(length(fix_ends)>0){
      fix_df<-data.frame(fix_ends-1,fix_ends,rep(max(sg_df2$counts),length(fix_ends)))
      colnames(fix_df)<-c("start","end","weight")
      rownames(fix_df)<-paste0(seq(max(sg_df2$featureID),max(sg_df2$featureID)+length(fix_ends)-1),".1")
    }
    #-----------------------------
    
    rownames(sg_df2)<-sg_df2$featureID
    
    sg_df2<-sg_df2[,c("start","end","counts")]
    colnames(sg_df2)<-c("start","end","weight")
    sg_df2<-rbind(sg_df2,fix_df)
    sg_df2 <- sg_df2[order(sg_df2$start),] 
    #sg_df2<-sg_df2[sg_df2$weight>0,]
    sg<-graph.data.frame(sg_df2)
    #initiation of the graph reconstruction suporting parameters
    curr_edge_exons<-intersect(edge_features,rownames(sg_df2))
    #curr_circ<-as.data.frame(BSJ_sg[BSJ_sg@geneID==j])
    #curr_circ<-unname(curr_circ)
    curr_circ<-as.data.frame(table_circ[table_circ$start%in%sg_exons_df$start&table_circ$end%in%sg_exons_df$end,])
    if(length(curr_circ$start)==0){next}
    curr_circ$used<-0
    tracking_circ<-curr_circ
    #transcript_features<-c()
    min.feature<-min(sg_df2$weight)
    #starts here 
    repeat{
      sg<-graph.data.frame(sg_df2)
      #print(sg_df2)
      #plot(sg, layout=layout.circle, edge.label=paste0(rownames(sg_df2),": ",sg_df2$weight ),edge.label.cex=0.8)
      #check is the circles can be reconstructed, based on avalability of edge exons; save the circles that can use the exons as current 
      curr_circ<-curr_circ[curr_circ$end%in%sg_df2$end,]
      curr_circ<-curr_circ[curr_circ$start%in%sg_df2$start,]
      #the lowest coverage exon
      curr_sg_exons<-sg_df2[as.character(intersect(rownames(sg_df2),sg_exons_df$featureID)),]
      #sg_df2_plot<-sg_df2[rownames(curr_sg_exons),]
      #sg_df2_plot[rownames(sg_df2_plot)%in%rownames(sg_df2),3]<-sg_df2[rownames(sg_df2)%in%rownames(sg_df2_plot),3]
      #sg_df2_plot[!rownames(sg_df2_plot)%in%rownames(sg_df2),3]<-0
      #barplot(sg_df2_plot$weight, ylim=c(0,300000),col="darkblue", names.arg = rownames(sg_df2_plot))
      min_exon<-curr_sg_exons[curr_sg_exons$weight==min(curr_sg_exons$weight),][1,]
      min_exon_weight<-min_exon$weight
      #the circles containing that exon 
      work_circ<-curr_circ[curr_circ$start<=min_exon$start&curr_circ$end>=min_exon$end,]
      #quick check if the circle hasnt been depleted 
      if (length(work_circ[,1])==0) {
        sg_df2<-sg_df2[rownames(sg_df2)!=rownames(min_exon),]
      }else{
        #the minimum quantitiy circle containing that exon 
        work_circ<-work_circ[work_circ$count==min(work_circ$count),]
        work_circ<-work_circ[1,] #just on case there are the same number of BSJ
        #in case of very low levels of circle, it can get depleted before the actual use of the BSJ; this is a rough fix for that case; continues at the end of the loop 
        #increasing a tracker counter when a circle is used 
        tracking_circ[tracking_circ$chr==work_circ$chr&tracking_circ$start==work_circ$start&tracking_circ$end==work_circ$end,"used"]<-tracking_circ[tracking_circ$chr==work_circ$chr&tracking_circ$start==work_circ$start&tracking_circ$end==work_circ$end,"used"]+1
        curr_exons<-sg_exons_df[sg_exons_df$start>=work_circ$start&sg_exons_df$end<=work_circ$end,"featureID"]
        if (length(curr_exons)<3) {
          transcript_features<-c(paste(curr_exons,collapse = "_"),transcript_features)
          sg_df2$weight[rownames(sg_df2)%in%curr_exons]<-(sg_df2$weight-as.numeric(min_exon_weight))[rownames(sg_df2)%in%curr_exons]
        }else{  
          #max flow trough the min exon ; seperated into 2 max flows 
          a<-max_flow(sg,as.character(work_circ[1,2]),as.character(min_exon[1,1]),capacity = sg_df2$weight)
          b<-max_flow(sg,as.character(min_exon[1,2]),as.character(work_circ[1,3]),capacity = sg_df2$weight)
          #if the path is imposible reverse the increase in the tracker 
          if(all(b$flow==0)|all(a$flow==0)) {
            tracking_circ[tracking_circ$chr==work_circ$chr&tracking_circ$start==work_circ$start&tracking_circ$end==work_circ$end,"used"]<-tracking_circ[tracking_circ$chr==work_circ$chr&tracking_circ$start==work_circ$start&tracking_circ$end==work_circ$end,"used"]-1
            sg_df2<-sg_df2[rownames(sg_df2)!=rownames(min_exon),]
            next
          }
          curr_flow<-c(rownames(sg_df2)[as.logical(a$flow)],rownames(min_exon),rownames(sg_df2)[as.logical(b$flow)])
          
          # the list of reaconstructed transcripts as feature ids
          #the if is sanity check for lack of edge exons; if the path is imposible the length will be less then 3 
          if(length(setdiff(curr_flow,rownames(fix_df)))>=3) {transcript_features<-c(paste(intersect(curr_flow,sg_exons_df$featureID),collapse = "_"),transcript_features)}
          
          #depletion of the coverage 
          sg_df2<-sg_df2[rownames(sg_df2)!=rownames(min_exon),]
          #sg_df2$weight[rownames(sg_df2)%in%curr_exons]<-(sg_df2$weight-as.numeric(min_exon_weight))[rownames(sg_df2)%in%curr_exons]
          #this switches on and off the depletion of the junctions
          #curr_flow<-intersect(curr_flow,curr_exons)
          sg_df2$weight[rownames(sg_df2)%in%curr_flow]<-(sg_df2$weight-as.numeric(min_exon_weight))[rownames(sg_df2)%in%curr_flow]
          #sg_exons_df<-sg_exons_df[sg_exons_df$weight>0.2*avg_weight_exons_circ]
          #sg_df2<-sg_df2[sg_df2$weight>0.05*avg_weight_exons_circ,]
        }
      }
      #print(transcript_features)
      #sg_df2<-sg_df2[sg_df2$weight>0.04*avg_weight_exons_circ,]
      sg_df2<-sg_df2[sg_df2$weight>0.05*min.feature,]
      #sg_df2$weight[sg_df2$weight>0.1*avg_weight_exons_circ&rownames(sg_df2)%in%curr_flow]<-0
      #sg_df2<-sg_df2[sg_df2$weight>0,]
      #print(sg_df2)
      if (length(intersect(rownames(sg_df2),curr_edge_exons))==0) {break}
    }
    # here is the continuation of the fix from above  for the prematurely depleted circles 
    tracking_circ<-tracking_circ[tracking_circ$used==0,]
    if (length(tracking_circ$used)>0) {  
      #print("enter tracker")
      for(k in 1:length(tracking_circ$used)){
        work_circ<-tracking_circ[k,]
        curr_exons<-sg_exons_df[sg_exons_df$start>=work_circ$start&sg_exons_df$end<=work_circ$end,]
        #curr_exons<-setdiff(curr_exons,feature_id_ri)
        curr_exons<-curr_exons[curr_exons$counts>0.2*avg_weight_exons_circ,"featureID"]
        transcript_features<-c(paste(curr_exons,collapse = "_"),transcript_features)
      }
    }
  }
  # just in case 
  transcript_features<-unique(transcript_features)
  transcript_features<-transcript_features[transcript_features!=""]
  transcript_features
}
transcripts.per.sample<-function(i){
  transcript_features<-assemble.transcripts.per.sample(i)
  qics_out<-prep.output(transcript_features,circ_exons)
  qics_out
}
merge.qics<-function(qics1,qics2){
  qics_out_merged<-rbind(qics1,qics2)
  qics_out_merged$circID<-1
  qics_out_merged<-unique(qics_out_merged)
  b<-qics_out_merged[qics_out_merged$seq%in%qics_out_merged$seq[duplicated(qics_out_merged$seq)],]
  a<-qics_out_merged[!qics_out_merged$seq%in%qics_out_merged$seq[duplicated(qics_out_merged$seq)],]
  b<-b[order(b$exon_starts),]
  c<-b[!duplicated(b$seq),]
  qics_out_merged<-rbind(a,c)
  qics_out_merged<-qics_out_merged[order(qics_out_merged$gene),]
  qics_out_merged
  #qics_out_sim_75$chr<-sgfc_pred@rowRanges@seqnames@values[qics_out_sim_75$chr]
  #qics_out_sim_75$circID<-paste(1:length(qics_out_sim_75$circID),qics_out_sim_75$chr,qics_out_sim_75$start,qics_out_sim_75$end, sep = "_")
  
}
merge.fasta<-function(qics_fa,known_fa){
  work_fa<-setdiff(qics_fa,known_fa)
  work_fa<-union(known_fa,qics_fa)
}