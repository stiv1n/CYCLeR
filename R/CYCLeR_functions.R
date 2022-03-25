## Copyright (C) 2021 Stefan Stefanov and Irmtraud M. Meyer (www.e-rna.org)
## Contact information: irmtraud.meyer@cantab.net

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License,or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.




##' Parse BSJ files from CIRI, CIRCexplorer2 or a TSV file
##'
##' This processes BSJ prediction files and prepares them for the next step of the pipeline.
##' \code{input_type} is essential for the correct parsing of the files. 
##' 
##' @title Parse BSJ input
##' @param file_list list with file names
##' @param file_path string object with file path, clould be an empty string
##' @param input_type \code{CIRI} for CIRI2 input, \code{CE} for CIRCexplorer2 input and \code{tsv} 
##'  for TSV formatted input
##' @return \code{Tibble} object with combined BSJ coordinate and number of junction spanning 
##'  reads across sample
##' @keywords parse
##' @author Stefan Stefanov
##' @export

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

read.CIRI<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  ciri<-readr::read_tsv(file_name, comment = "", col_names = T, col_types = cols(.default = "c" ))
  dplyr::select(ciri,(chr:`#junction_reads`),strand)%>%unite(circ_id,chr,circRNA_start,circRNA_end,strand)%>%
    dplyr::rename(!!sample_name:=`#junction_reads`)%>%mutate_at(2,as.double)
}
read.CE<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  ce<-read_tsv(file_name, comment = "", col_names = F,col_types = cols(.default = "c" ))
  dplyr::select(ce,(X1:X3),X6,X13)%>%mutate(X2=as.numeric(X2)+1)%>%unite(circ_id,X1,X2,X3,X6)%>%
    dplyr::rename(!!sample_name:=X13)%>%mutate_at(2,as.double)
}
read.tsv.file<-function(file_name){
  sample_name<-tail(str_split(file_name,"_", simplify = T)[1,], n=1)
  df<-read_tsv(file_name, comment = "", col_names = F,col_types = cols(.default = "c" ))
  dplyr::unite(df,circ_id,X1,X2,X3,X4)%>%dplyr::rename(!!sample_name:=X5)%>%mutate_at(2,as.double)
}

##' process the BSJ table and select high confidence BSJs
##'
##' Filters BSJ based on comparison of the average CPM values of BSJs 
##' 
##' @title Process BSJs
##' @param cdf tibble produced by \code{parse.files}
##' @param file_path string object with file path, clould be an empty string
##' @param sample_table sample table formatted according to the manual,
##'   Must contain \dQuote{sample_name} \dQuote{treatment} \dQuote{file_bam} \dQuote{lib_size} 
##'   \dQuote{read_len}; NB the values in column \dQuote{treatment} can only be \dQuote{control} and 
##'   \dQuote{enriched}
##' @return \code{Tibble} object with combined filtered BSJ coordinate and number of junction spanning 
##'  reads across sample. 
##' @keywords filter BSJ
##' @author Stefan Stefanov
##' @export
process_BSJs<-function(cdf,sample_table){
  #stands for circular data frame
  cdf[is.na(cdf)] <- 0
  #cdf<-column_to_rownames(cdf,var = "circ_id")
  sample_index<-c(sample_table$sample_name[sample_table$treatment=="enriched"])
  cdf2<-cdf[,sample_index]
  cdf<-cdf[rowSums(cdf2>1)>0,]
  cdf[,sample_table$sample_name]<-(cdf[,sample_table$sample_name]/sample_table$lib_size)*10e6
  #cdf<-rownames_to_column(cdf,var = "circ_id")
  cdf$meanc<-rowMeans(cdf[,c(sample_table$sample_name[sample_table$treatment=="control"])])
  cdf$meanRR<-rowMeans(cdf[,c(sample_table$sample_name[sample_table$treatment=="enriched"])])
 # rownames(cdf[cdf$meanc<cdf$meanr,])
  cdf[cdf$meanc<cdf$meanRR,]
}


##' Combine 2 BSJ tables
##'
##' Just a combination of BSJ tables to make sure we have a complete set of BSJs. 
##'  The variable names do not actually matter since the all tables have the same formatting.
##' 
##' @title combine BSJs
##' @param ce_bsjs BSJ table 1
##' @param ciri_bsjs BSJ table 2
##' @return \code{Tibble} object with combined filtered BSJ coordinate and number of junction spanning 
##'  reads across sample. 
##' @keywords combine BSJ
##' @author Stefan Stefanov
##' @export
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

##'Convert BSJ string to GRanges obejct
##'
##' Convert BSJ string to \code{GRanges} obejct
##' 
##' @title Convert BSJ string to GRanges obejct
##' @param BSJ_set a list of BSJ ID records procudes by \code{process.BSJs} or \code{combine.two.BSJ.tables}
##' @return \code{GRanges} object indicating BSJ loci 
##' @keywords GRanges BSJ
##' @author Stefan Stefanov
##' @export
make.BSJ.gr<-function(BSJ_set){
  BSJ_BED<-as.data.frame(str_split(BSJ_set,"_",simplify = T),stringsAsFactors=F)
  BSJ_BED$V2<-as.numeric(BSJ_BED$V2)
  BSJ_BED$V3<-as.numeric(BSJ_BED$V3)
  BSJ_gr <- GRanges(BSJ_BED$V1, IRanges(start=BSJ_BED$V2,end=BSJ_BED$V3),
                    strand=BSJ_BED$V4, seqlengths=NULL)
  BSJ_gr@elementMetadata@listData[["gene_id"]]<-paste0("circ",seq_along(1:length(BSJ_gr)))
  BSJ_gr
}

##'Plots GRanges objects 
##'
##' ggplot of multiple GRanges object. Every object is auto assigned a colour from colorblind friendly scheme
##' 
##' @title Plot ranges
##' @return ggplot of multiple GRanges objects 
##' @keywords GRanges plot
##' @author Stefan Stefanov
##' @export
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

##'Creates a disjointed set of exons based on a \code{SGSeq} obejct and a BSJ \code{GRanges} object
##'
##' Creates a disjointed set of exons based on a \code{SGSeq} obejct and a BSJ \code{GRanges} object. 
##'  The function keeps the \code{SGSeq} metadata
##' 
##' @title Overlap of BSJ and a splice graph 
##' @param BSJ_gr a GRange of BSJ cooredinates
##' @param sgfc_pred \code{SGSeq} prediction object
##' @return \code{SGSeq} with disjoint exon bins
##' @keywords GRanges overlap
##' @author Stefan Stefanov
##' 
##' @export
overlap.SG.BSJ<-function(sgfc_pred,BSJ_gr,sg_annot){
  sg_gr<-rowRanges(sgfc_pred) #first we process the predicted features
  sg_gr_e<-sg_gr[sg_gr@type=="E"]

  BSJ_sg<-SGFeatures(BSJ_gr,type=rep("E",length(BSJ_gr)),splice5p =rep(F,length(BSJ_gr)),
                     splice3p = rep(F,length(BSJ_gr)), featureID =rep(1L,length(BSJ_gr)),
                     geneID = rep(1L,length(BSJ_gr)), txName = mcols(BSJ_gr)$geneID,
                     geneName = mcols(BSJ_gr)$geneID)
  
  sg_gr_e2<-sg_gr_e[sg_gr_e%over%BSJ_sg]  #this is to remove some mapping artifacts
  sg_gr_e<-sg_gr_e[sg_gr_e@geneID%in%sg_gr_e2@geneID] #and keep only the genes that overlap with BSJ regions
  #BSJ_sg2<-GenomicRanges::reduce(BSJ_sg)
  sg_an<-sg_annot #the same for annotated features 
  sg_an_e<-sg_an[sg_an@type=="E"]
  sg_an_e2<-sg_an_e[sg_an_e%over%BSJ_sg]  #this is to remove some mapping artifacts
  sg_an_e<-sg_an_e[sg_an_e@geneID%in%sg_an_e2@geneID] #and keep only the genes that overlap with BSJ region
  combined<-c(BSJ_sg,sg_gr_e,sg_an_e)
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
##'Selection of the exons based on BSJ set 
##'
##' Selection of the exons based on BSJ set
##' 
##' @title Preparation of the BSJ-specific splice graphs  
##' @param BSJ_gr a GRange of BSJ coordinates
##' @param circ_sg \code{SGSeq} prediction object
##' @return \code{SGSeq} containing exons belonging to BSJ loci
##' @keywords GRanges overlap
##' @author Stefan Stefanov
##' 
##' @export
make.BSJ.sg<-function(circ_sg,BSJ_gr){
  BSJ_sg<-SGFeatures(BSJ_gr,type=rep("E",length(BSJ_gr)),splice5p =rep(F,length(BSJ_gr)),
                     splice3p = rep(F,length(BSJ_gr)), featureID =rep(1L,length(BSJ_gr)),
                     geneID = rep(1L,length(BSJ_gr)), txName = mcols(BSJ_gr)$geneID,
                     geneName = mcols(BSJ_gr)$geneID)
  BSJ_sg@geneID<-circ_sg@geneID[queryHits(findOverlaps(circ_sg,BSJ_sg,type = c("start")))] 
  BSJ_sg
}
##' A wrapper function for samtools use to trim the files
##'
##' This function removes the BAM file reads that do not overlap with the BSJ loci.
##'  This significantly speeds us the feature detection and lowers the virtual memory requirements
##' 
##' @title BAM file filter  
##' @param BSJ_gr a GRange of BSJ cooredinates
##' @param sample_table sample table formatted according to the manual,
##'   Must contain \dQuote{sample_name} \dQuote{treatment} \dQuote{file_bam} \dQuote{lib_size} 
##'   \dQuote{read_len}; NB the values in column \dQuote{treatment} can only be \dQuote{control} and 
##'   \dQuote{enriched}
##' @param samtools_prefix a string that corresponds to user's samtools run prefix 
##' @return \code{BAMFileList} object with info on the trimmed files 
##' @keywords Bam Filter
##' @author Stefan Stefanov
##' 
##' @export
filter_bam<-function(BSJ_gr,sample_table,samtools_prefix){
  #gn_circ<-gene_ranges[gene_ranges%over%BSJ_gr]
  #unannot<-BSJ_gr[BSJ_gr%outside%gene_ranges]
  #filter_ranges<-c(gn_circ,unannot)
  bed_file<-paste0(bam_file_prefix,"/BSJ.bed")
  #export.bed(con = bed_file,object =  filter_ranges)
  rtracklayer::export.bed(con = bed_file,object =  BSJ_gr)
  #now we need to use the bed file to filter the bam files; i suggest doing it directly in the command line 
  list_file<-paste0(bam_file_prefix,"list")
  write(sample_table$sample_name,list_file)
  print("Files are trimming")
  trim_command<-paste(paste0(samtools_prefix,"samtools view -b -L ",bed_file," ", sample_table$file_bam," > ",bam_file_prefix,"/",sample_table$sample_name,"_sorted_trimmed.bam"), collapse = "; ")
  system(trim_command)
  print("Files are indexing")
  index_command<-paste(paste0(samtools_prefix,"samtools index ",bam_file_prefix,"/",sample_table$sample_name,"_sorted_trimmed.bam"), collapse = "; ")
  system(index_command)
  trimmed_bams<-BamFileList(paste0(bam_file_prefix,"/",sample_table$sample_name,"_sorted_trimmed.bam"), asMates = F)
  trimmed_bams
}
##' A wrapper function for Rsubread
##'
##' This function performs requantification of the exon bins with specifically selected parameters
##' 
##' @title Re-count of the reads per exon bin
##' @param full_sg a \code{SGSeq} object of exon bins
##' @param sample_table sample table formatted according to the manual,
##'   Must contain \dQuote{sample_name} \dQuote{treatment} \dQuote{file_bam} \dQuote{lib_size} 
##'   \dQuote{read_len}; NB the values in column \dQuote{treatment} can only be \dQuote{control} and 
##'   \dQuote{enriched}
##' @param paired_end a binary for pair-end info
##' @return \code{BAMFileList} object with info on the trimmed files 
##' @keywords Bam Filter
##' @author Stefan Stefanov
##' 
##' @export
recount.features<-function(full_sg,sample_table,paired_end=T){
  feature_info<-as.data.frame(full_sg,stringsAsFactors = F)
  # we make a feature reference for counting
  saf<-feature_info[feature_info$type=="E",c("featureID","seqnames","start","end","strand")]
  colnames(saf)<- c("GeneID",		"Chr",	"Start",	"End",	"Strand")
  #saf$GeneID<-c(1,2,3,4)
  saf.counts<-Rsubread::featureCounts(sample_table$file_bam, annot.ext=saf, isPairedEnd=paired_end,  countChimericFragments= T, 
                                      countMultiMappingReads =T ,juncCounts=T,countReadPairs=F,requireBothEndsMapped=F,allowMultiOverlap=T)
  saf.fc<-saf.counts$counts
  saf.fc
}
##' A wrapper function for BSgenome subsequencing
##'
##' Extracts sequence based on \code{SGSeq} object and \dQuote{BSgenome} name
##' 
##' @title Extract sequence per exon bin
##' @param full_sg a \code{SGSeq} object of exon bins
##' @param bs_genome \dQuote{BSgenome} name
##' @return sequence list 
##' @keywords seqs
##' @author Stefan Stefanov
##' 
##' @export
get.seqs<-function(full_sg, bs_genome=Dmelanogaster){
  seqs<-unname(BSgenome::getSeq(bs_genome, paste0("chr",full_sg@seqnames), start = start(full_sg), end = end(full_sg), strand=strand(full_sg),as.character=T))
  seqs
}
get.seqs.corrected<-function(full_sg, bs_genome=Dmelanogaster){
  seqs<-get.seqs(full_sg,bs_genome)
  short_ranges<-full_sg[width(full_sg)<35]
  long_ranges<-full_sg[!width(full_sg)<35]
  pair_ranges<-long_ranges[nearest(short_ranges,long_ranges)]
  corrected.seqs<-paste0(get.seqs(short_ranges),get.seqs(pair_ranges))
  seqs[width(full_sg)<35]<-corrected.seqs
  seqs
}
get.seqs.positive<-function(full_sg, bs_genome=Dmelanogaster){
  strand(full_sg)<-"+"
  get.seqs(full_sg,bs_genome)
}
#this one needs work to add the true/false if-s and test it; also figure out how to save the model 
##' RPKM calculation for the genomic features
##'
##' This function performs RPKM calculations for the exonic features. The RPKM calculation is
##'  performed based on the exact sequences for the exons. For junctions, the sequences are selected based 
##'  on the exons, flanking the junction. The function takes into account the needed effective 
##'  length corrections.
##' 
##' @title RPKM calculation for the genomic features 
##' @param count_matrix count matrix corresponding to the features 
##' @param sg \code{SGSeq} object supplying feature info
##' @param bsj_granges  GRange of BSJ cooredinates
##' @param bs_genome a \code{BSGenome} object used for extracting the sequences 
##' @param sample_table sample table formatted according to the manual,
##'   Must contain \dQuote{sample_name} \dQuote{treatment} \dQuote{file_bam} \dQuote{lib_size} 
##'   \dQuote{read_len}; NB the values in column \dQuote{treatment} can only be \dQuote{control} and 
##'   \dQuote{enriched}
##' @param feature_type either \dQuote{e} for exons ot \dQuote{j} for junctions
##' @param fsj_overhang the FJS overhand used in the mapping a.k.a. anchor
##' @param bsj_overhang the BSJ overhand used in the chimeric detection
##' @param eff_length_correction whether or not to apply effective length correction
##' @param gc_correction whether or not to apply GC-content correction; requires further testing
##' @return \code{BAMFileList} object with info on the trimmed files
##' @keywords RPKM
##' @author Stefan Stefanov
##' 
##' @export
RPKM.calc<-function(count_matrix, sg ,bsj_granges, bs_genome, sample_table,feature_type, fsj_overhang=3, 
                    bsj_overhang=15, eff_length_correction=T, gc_correction=F){
  bsj_overhang_adj=bsj_overhang-fsj_overhang
  eff_lengths<- width(sg)
  #the scaling is capped based on the  read length, because in case if something maps it is at least a 
  #read size 
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
##' CircRNA feature selection 
##'
##' This function works in 2 ways: direct comparison of average quantities or as a wrapper of DEXSeq.
##' In case of dataset with replicates, the suggested approach is the use of DEXSeq statistical test.
##' 
##' @title CircRNA feature selection 
##' @param circ_fc_adj count matrix corresponding to the circRNA features 
##' @param circ_sg \code{SGSeq} object supplying feature info
##' @param sample_table sample table formatted according to the manual,
##'   Must contain \dQuote{sample_name} \dQuote{treatment} \dQuote{file_bam} \dQuote{lib_size} 
##'   \dQuote{read_len}; NB the values in column \dQuote{treatment} can only be \dQuote{control} and 
##'   \dQuote{enriched}
##' @param test either \dQuote{DEX} for DEXSeq based feature selection or \dQuote{comparison} simple average
##'   comaparison
##' @return vector of featureID
##' @keywords depleted
##' @author Stefan Stefanov
##' 
##' @export
find.depleted.features<-function(circ_fc_adj,sample_table,circ_sg, test="DEX"){#may need to add feature type
  if(length(sample_table$sample_name)<4) {test="comparison" 
                                          print("test strategy has been switched to comparison")}
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
    #feature_dex_res_depleted<-feature_dex_res$featureID[feature_dex_res$log2fold_enriched_control<(-0.5)&feature_dex_res$padj<=0.05]
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
  seqs<-get.seqs.positive(circ_exons,bs_genome = bs_genome)
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
##' Creation of GTF based on CYCLeR results 
##'
##' This function takes the 
##' 
##' @title Creation of GTF based on CYCLeR results 
##' @param qics CYCLeR table  of intermediate results
##' @param circ_exons \code{SGSeq} object supplying feature info
##' @param annot_list ORG package; soon to be expanded 
##' @return GTF-like table
##' @keywords GTF
##' @author Stefan Stefanov
##' 
##' @export
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
  #print(circ_junc_counts[1,i])
  #########################################
  #start recontruction here
  #we need to create the splice graph 
  transcript_features<-c()
  #it is OK to use for loops here since there should not be more than 10000 iterations 
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
    #sg_junc_df$counts[sg_junc_df$counts<min(sg_exons_df$counts)]<-round(min(sg_exons_df$counts)+0.1*avg_weight_exons_circ)
    #sg_exons_df$counts[sg_exons_df$featureID%in%edge_features]<-round(min(sg_exons_df$counts)+0.1*avg_weight_exons_circ)
    
    #sg_junc_df$counts<-rep(23031,length(sg_junc_df$counts))
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
        edge_exons<-c(head(curr_exons,n=1),tail(curr_exons,n=1))
        if (length(curr_exons)<3) {
          transcript_features<-c(paste(curr_exons,collapse = "_"),transcript_features)
          sg_df2$weight[rownames(sg_df2)%in%curr_exons]<-(sg_df2$weight-as.numeric(min_exon_weight))[rownames(sg_df2)%in%curr_exons]
        }else{  
          #max flow trough the min exon ; seperated into 2 max flows 
          a<-max_flow(sg,as.character(work_circ[1,2]),as.character(min_exon[1,1]),capacity = sg_df2$weight)
          b<-max_flow(sg,as.character(min_exon[1,2]),as.character(work_circ[1,3]),capacity = sg_df2$weight)
          #print(rownames(sg_df2)[as.logical(a$flow)])
          #print(rownames(min_exon))
          #print(rownames(sg_df2)[as.logical(b$flow)])
          #if the path is imposible reverse the increase in the tracker 
          if((all(b$flow==0)|all(a$flow==0))&!is.element(rownames(min_exon),edge_exons)) {
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
      sg_df2<-sg_df2[sg_df2$weight>0.001*avg_weight_exons_circ,]
      #sg_df2<-sg_df2[sg_df2$weight>0.05*min.feature,]
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
        #curr_exons<-curr_exons[curr_exons$counts>0.002*avg_weight_exons_circ,"featureID"]
        curr_exons<-curr_exons[curr_exons$counts>0.05*avg_weight_exons_circ,"featureID"]
        transcript_features<-c(paste(curr_exons,collapse = "_"),transcript_features)
      }
    }
  }
  # just in case 
  transcript_features<-unique(transcript_features)
  transcript_features<-transcript_features[transcript_features!=""]
  transcript_features
}
##' Transcript assembly per sample 
##' Transcript assembly per sample based on sample name in the \dQuote{sample_table}
##' @title Transcript assembly
##' @param i name of the sample 
##' @return \code{data.frame} of transcript information in flat format 
##' @keywords assembly
##' @author Stefan Stefanov
##' 
##' @export
transcripts.per.sample<-function(i){
  transcript_features<-assemble.transcripts.per.sample(i)
  qics_out<-prep.output(transcript_features,circ_exons)
  qics_out
}
##' Pair-wise merging 2 assemblies
##' Pair-wise merging 2 assemblies
##' @title Merging 2 assemblies 
##' @param qics1 assembly 1
##' @param qics2 assembly 2
##' @return \code{data.frame} of transcript information in flat format 
##' @keywords assembly
##' @author Stefan Stefanov
##' 
##' @export
merge_qics<-function(qics1,qics2){
  qics_out_merged<-rbind(qics1,qics2)
  qics_out_merged$circID<-1
  qics_out_merged<-unique(qics_out_merged)
  b<-qics_out_merged[qics_out_merged$seq%in%qics_out_merged$seq[duplicated(qics_out_merged$seq)],]
  a<-qics_out_merged[!qics_out_merged$seq%in%qics_out_merged$seq[duplicated(qics_out_merged$seq)],]
  b<-b[order(b$exon_starts),]
  c<-b[!duplicated(b$seq),]
  qics_out_merged<-rbind(a,c)
  qics_out_merged<-qics_out_merged[order(qics_out_merged$gene),]
  qics_out_merged$chr<-sgfc_pred@rowRanges@seqnames@values[qics_out_merged$chr]
  qics_out_merged$circID<-paste(1:length(qics_out_merged$circID),qics_out_merged$chr,qics_out_merged$start,qics_out_merged$end, sep = "_")
  qics_out_merged
}
merge_fasta<-function(qics_fa,known_fa){
  work_fa<-setdiff(qics_fa,known_fa)
  work_fa<-union(known_fa,qics_fa)
}
