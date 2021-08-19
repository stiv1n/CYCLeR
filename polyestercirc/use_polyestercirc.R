#polyestercirc follows the established polyester workflow
#the only difference being generating a circular molecule with a random nick in the sequence
#note that the extra code significantly decreases the speed of the algorithm 

library("scales")
library("Biostrings")
library("polyester")


countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

# recalculate a data frame for frag_GC_bias
{data("loessfit7")
fit.df<-as.data.frame(loessfit7)
loess_mod <- loess(y~x, data = fit.df, control=loess.control(surface="direct"), span=0.3)
pdf <- data.frame(x = seq(0, 1, by = 0.01))
loess_df <- data.frame(Percent = pdf, count_adjustment = predict(loess_mod, newdata = pdf))
gc_prob<-rescale(loess_df$count_adjustment, to=c(0,1))
loess_df$prob<-gc_prob
frag_GC_bias<-as.matrix(cbind(gc_prob,gc_prob,gc_prob,gc_prob))
}

# FASTA annotation
fasta_lin<-readDNAStringSet("linear.fa")
fasta_lin_nascend<-readDNAStringSet("nascend.fa")
fasta_circ_high<-readDNAStringSet("/home/sstefan/data/circ_high.fa")
fasta_circ_low<-readDNAStringSet("/home/sstefan/data/circ_high.fa")

#calculating number of reads per sample; 
#the sample values act as a factor of what the ratio of linear and circular transcripts will be
readspertx_lin<-sample(100:300, length(fasta_lin_selected), replace=TRUE)
readspertx_circ_high<-sample(3:20, length(fasta_circ_high), replace=TRUE)
readspertx_circ_low<-sample(3:4, length(fasta_circ_low), replace=TRUE)
readspertx_lin_nascend<-sample(10:15, length(fasta_lin_nascend), replace=TRUE)
#number of reads to per transcript
readspertx_lin<-round((width(fasta_lin_selected) / 500)*readspertx_lin)
readspertx_circ_high<-round((width(fasta_circ_high) / 50)*readspertx_circ_high)
readspertx_circ_low<-round((width(fasta_circ_low) / 50)*readspertx_circ_low)
readspertx_lin_nascend<-round((width(fasta_lin_nascend) / 500)*readspertx_lin_nascend)
fasta_lin_full<-c(fasta_lin_selected,fasta_lin_nascend)
readspertx_lin_full<-c(readspertx_lin,readspertx_lin_nascend)
fasta_circ_full<-c(fasta_circ_high,fasta_circ_low)
readspertx_circ_full<-c(readspertx_circ_high,readspertx_circ_low)


#good statistics to keep track of 
readspertx_circ_full[width(fasta_circ_full)<600]
readspertx_lin_full[width(fasta_lin_full)<2000]
sum(readspertx_circ_full)
sum(readspertx_lin_full)

#fold_changes_lin<-matrix(c(rep(1,length(readspertx_lin)),rep(0.2,length(readspertx_lin))),nrow = length(readspertx_lin))
#fold_changes_nascend<-matrix(c(rep(1,length(readspertx_lin_nascend)),rep(0.2,length(readspertx_lin_nascend))),nrow = length(readspertx_lin_nascend))

fold_changes_lin<-matrix(c(rep(1,length(readspertx_lin_full)),rep(0.2,length(readspertx_lin_full))),nrow = length(readspertx_lin_full))
fold_changes_circ<-matrix(c(rep(1,length(readspertx_circ_full)),rep(2.5,length(readspertx_circ_full))),nrow = length(readspertx_circ_full))


###########final prep################
writeXStringSet(fasta_circ_full,'circles_seq_test.fa')
writeXStringSet(fasta_lin_full,'circles_seq_lin_selected.fa')



# simulation with the modified frag_GC_bias

Sys.time()
polyestercirc::simulate_experiment('circles_seq_test.fa', reads_per_transcript=readspertx_circ_full,   strand_specific=TRUE,
                                   error_model="illumina5",
                                   readlen=75, #distr= "empirical" ,
                                   fraglen=280, fragsd=50,  
                                   frag_GC_bias=frag_GC_bias,
                                   #bias="rnaf",
                                   num_reps=c(2,2), fold_changes=fold_changes_circ, outdir='.../simulation/circ') 
Sys.time()

polyester::simulate_experiment('circles_seq_lin_selected.fa', reads_per_transcript=readspertx_lin_full, 
                               bias="rnaf", 
                               frag_GC_bias=frag_GC_bias,
                               error_model="illumina5",
                               readlen=75, strand_specific=TRUE, #distr= "empirical" ,
                               fraglen=300, fragsd=50,  
                               num_reps=c(2,2), fold_changes=fold_changes_lin, outdir='.../simulation/lin') 

Sys.time()
#the files need to be merged subsequently
#suggested marking the circular reads with "c" makes for easier troubleshooting, when benchmarking 
# sed 's/[>]/>c/' sample_01_1.fasta > csample_01_1.fasta
