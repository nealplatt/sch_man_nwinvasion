#conda install -c r r r-dplyr  r-chron r-maps
#conda install -c conda-forge r-ade4 r-ape r-httpuv r-spdep r-spam r-vegan r-miniui r-devtools
#conda install -c bioconda r-seqinr
#install.packages("ape", dependencies = TRUE)
#install.packages("spdep", dependencies = TRUE)
#install.packages("chron", dependencies = TRUE)
#install.packages("rrcovNA", dependencies = TRUE)
#install.packages("fitdistrplus", dependencies = TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("qvalue")

library("dplyr")
library("rrcovNA")
library(MINOTAUR)
library("MASS")
library("fitdistrplus")
library("qvalue")

setwd("/master/nplatt/sch_man_nwinvasion/results/dcms")

for (pop in c("senegal", "niger", "tanzania", "brazil"))
{
    print(pop)
    #---------------------------------------------------------------------------
    #read in sel stats
    csv_infile = paste(pop, "_raw_data.csv", sep="")
    sel_table<-read.table(csv_infile, header=TRUE, sep=",")

    # h (hscan)
    h_rank<-percent_rank(sel_table$h)
    h_norm<-qnorm((rank(h_rank,na.last="keep")-0.5)/sum(!is.na(h_rank)))
    h_p<-pnorm(h_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

    # fst
    fst_rank<-percent_rank(sel_table$fst)
    fst_norm<-qnorm((rank(fst_rank,na.last="keep")-0.5)/sum(!is.na(fst_rank)))
    fst_p<-pnorm(fst_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

    # lr (sweepfinder2)
    lr_rank<-percent_rank(sel_table$lr)
    lr_norm<-qnorm((rank(lr_rank,na.last="keep")-0.5)/sum(!is.na(lr_rank)))
    lr_p<-pnorm(lr_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


    dfv<-data.frame(sel_table$h,  sel_table$fst,  sel_table$lr)
    dfp<-stat_to_pvalue( dfv, 
                         column.nums  = 1:ncol(dfv), 
                         subset       = 1:nrow(dfv), 
                         two.tailed   = rep(FALSE, 3), 
                         right.tailed = rep(TRUE, 3)
                       )


    dfv<-na.omit(dfv)
    dfp<-na.omit(dfp)

    #calculate a cov matrix on the pvalues
    cov_matrix<-CovNAMcd(dfp, alpha = 0.75, nsamp = 100000)

    #decorelatced composite test
    dcms<-DCMS( dfv, 
                subset = 1:nrow(dfv),
                S=cov_matrix$cov,
                column.nums=1:3,
                dfp,
                column.nums.p =1:3
               )           

    #convert DCMS to p
    dcms_ps<-pnorm(dcms, mean = 0 , sd = 2, lower.tail = FALSE, log.p = FALSE)

    #convert DCMS to q
    dcms_qs<-qvalue(dcms_ps, fdr=0.01)
    dcms_logq=-log10(dcms_qs$qvalues)

    ##rough plot
    png_outfile = paste(pop, "_dcms.png", sep="")
    png(file=png_outfile)
        plot(dcms_logq,
             pch=19,
             cex=0.5,
             col="blue"
             )
    dev.off()

    out_df<-data.frame(sel_table$chrom, sel_table$pos, dfv, dfp, dcms, dcms_ps, dcms_qs$qvalues, dcms_logq, sel_table$fig_x_pos_x)
    colnames(out_df)<-c("chrom", "pos", "h", "fst", "lr", "h_p", "fst_p", "lr_p", "dcms", "dcms_p", "dcms_q", "dcms_logq", "fig_x_pos")   

    csv_outfile = paste(pop, "_dcms.csv", sep="")
    write.csv(out_df, csv_outfile, row.names=FALSE, quote=FALSE)
}

# #make bed
# sig_quantile<-quantile(out_data$dcms, 0.995)
# sig_snps<-subset(out_data, dcms >= sig_quantile)
# cols<-c("chrom", "pos", "pos")
# sig_snps<-sig_snps[cols]
# write.table(sig_snps, "nw_vs_ea_dcms_sig.bed", sep="\t", row.names=FALSE, quote=FALSE)

