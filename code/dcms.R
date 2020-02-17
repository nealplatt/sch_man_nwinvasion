library("rrcovNA")
library(MINOTAUR)
library("MASS")
library("fitdistrplus")

setwd("C:/Users/rplatt/Dropbox (TX Biomed)/projects/sch_man_nwinvasion/results/dcms")

#----------------------------------------------------------------------------------------------------------------------
#   NEW WORLD VS EAST AFRICA
#----------------------------------------------------------------------------------------------------------------------
#read in sel stats
sel_table<-read.table("nw_vs_ea.csv", header=TRUE, sep=",")

#for now only keep certain statistics
keeps <- c("chrom", "pos", "smoothed_pi", "smoothed_td", "smoothed_ihs", "smoothed_h12", "smoothed_sxpehh")

pi_rank<-percent_rank(sel_table$smoothed_pi)
pi_norm<-qnorm((rank(pi_rank,na.last="keep")-0.5)/sum(!is.na(pi_rank)))
pi_p<-pnorm(pi_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


td_rank<-percent_rank(sel_table$smoothed_td)
td_norm<-qnorm((rank(td_rank,na.last="keep")-0.5)/sum(!is.na(td_rank)))
td_p<-pnorm(td_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

abs_ihs<-abs(sel_table$smoothed_ihs)
ihs_rank<-percent_rank(abs_ihs)
ihs_norm<-qnorm((rank(ihs_rank,na.last="keep")-0.5)/sum(!is.na(ihs_rank)))
ihs_p<-pnorm(ihs_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


h12_rank<-percent_rank(sel_table$smoothed_h12)
h12_norm<-qnorm((rank(h12_rank,na.last="keep")-0.5)/sum(!is.na(h12_rank)))
h12_p<-pnorm(h12_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

abs_sxpehh<-abs(sel_table$smoothed_sxpehh)
xpehh_rank<-percent_rank(abs_sxpehh)
xpehh_norm<-qnorm((rank(xpehh_rank,na.last="keep")-0.5)/sum(!is.na(xpehh_rank)))
xpehh_p<-pnorm(xpehh_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


dfp<-data.frame(pi_p, td_p, ihs_p, h12_p, xpehh_p)

#for now only keep certain statistics
dfv<-sel_table[keeps]
dfv<-na.omit(dfv)
dfp<-na.omit(dfp)

#calculate a cov matrix on the pvalues
cov_matrix<-CovNAMcd(dfp, alpha = 0.75, nsamp = 100000)

#decorelatced composite test
dcms<-DCMS(dfv, 
           subset = 1:nrow(dfv),
           S=cov_matrix$cov,
           column.nums=3:7,
           dfp,
           column.nums.p =1:5
)

plot(dcms,
     pch=19,
     cex=0.5,
     col="blue"
     )

out_data<-data.frame(dfv, dfp, dcms)

#write.csv(out_data, "nw_vs_ea_dcms.csv", row.names=FALSE, quote=FALSE)

#make bed
sig_quantile<-quantile(out_data$dcms, 0.995)
sig_snps<-subset(out_data, dcms >= sig_quantile)
cols<-c("chrom", "pos", "pos")
sig_snps<-sig_snps[cols]
write.table(sig_snps, "nw_vs_ea_dcms_sig.bed", sep="\t", row.names=FALSE, quote=FALSE)


#----------------------------------------------------------------------------------------------------------------------
#   WEST AFRICA VS EAST AFRICA
#----------------------------------------------------------------------------------------------------------------------
#read in sel stats
sel_table<-read.table("wa_vs_ea.csv", header=TRUE, sep=",")

#for now only keep certain statistics
keeps <- c("chrom", "pos", "smoothed_pi", "smoothed_td", "smoothed_ihs", "smoothed_h12", "smoothed_sxpehh")

pi_rank<-percent_rank(sel_table$smoothed_pi)
pi_norm<-qnorm((rank(pi_rank,na.last="keep")-0.5)/sum(!is.na(pi_rank)))
pi_p<-pnorm(pi_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


td_rank<-percent_rank(sel_table$smoothed_td)
td_norm<-qnorm((rank(td_rank,na.last="keep")-0.5)/sum(!is.na(td_rank)))
td_p<-pnorm(td_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

abs_ihs<-abs(sel_table$smoothed_ihs)
ihs_rank<-percent_rank(abs_ihs)
ihs_norm<-qnorm((rank(ihs_rank,na.last="keep")-0.5)/sum(!is.na(ihs_rank)))
ihs_p<-pnorm(ihs_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


h12_rank<-percent_rank(sel_table$smoothed_h12)
h12_norm<-qnorm((rank(h12_rank,na.last="keep")-0.5)/sum(!is.na(h12_rank)))
h12_p<-pnorm(h12_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

##### VERY IMPROTANT TAIL CHANGE
abs_sxpehh<-abs(sel_table$smoothed_sxpehh)
xpehh_rank<-percent_rank(abs_sxpehh)
xpehh_norm<-qnorm((rank(xpehh_rank,na.last="keep")-0.5)/sum(!is.na(xpehh_rank)))
xpehh_p<-pnorm(xpehh_norm, mean = 0, sd = 2, lower.tail = TRUE, log.p = FALSE)


dfp<-data.frame(pi_p, td_p, ihs_p, h12_p, xpehh_p)

#for now only keep certain statistics
dfv<-sel_table[keeps]
dfv<-na.omit(dfv)
dfp<-na.omit(dfp)

#calculate a cov matrix on the pvalues
cov_matrix<-CovNAMcd(dfp, alpha = 0.75, nsamp = 100000)

#decorelatced composite test
dcms<-DCMS(dfv, 
           subset = 1:nrow(dfv),
           S=cov_matrix$cov,
           column.nums=3:7,
           dfp,
           column.nums.p =1:5
)

plot(dcms,
     pch=19,
     cex=0.5,
     col="green"
)

out_data<-data.frame(dfv, dfp, dcms)

#write.csv(out_data, "wa_vs_ea_dcms.csv", row.names=FALSE, quote=FALSE)

#make bed
sig_quantile<-quantile(out_data$dcms, 0.995)
sig_snps<-subset(out_data, dcms >= sig_quantile)
cols<-c("chrom", "pos", "pos")
sig_snps<-sig_snps[cols]
write.table(sig_snps, "wa_vs_ea_dcms_sig.bed", sep="\t", row.names=FALSE, quote=FALSE)


#----------------------------------------------------------------------------------------------------------------------
#   NEW WORLD VS WEST AFRICA
#----------------------------------------------------------------------------------------------------------------------
#read in sel stats
sel_table<-read.table("nw_vs_wa.csv", header=TRUE, sep=",")

#for now only keep certain statistics
keeps <- c("chrom", "pos", "smoothed_pi", "smoothed_td", "smoothed_ihs", "smoothed_h12", "smoothed_sxpehh")

pi_rank<-percent_rank(sel_table$smoothed_pi)
pi_norm<-qnorm((rank(pi_rank,na.last="keep")-0.5)/sum(!is.na(pi_rank)))
pi_p<-pnorm(pi_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


td_rank<-percent_rank(sel_table$smoothed_td)
td_norm<-qnorm((rank(td_rank,na.last="keep")-0.5)/sum(!is.na(td_rank)))
td_p<-pnorm(td_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

abs_ihs<-abs(sel_table$smoothed_ihs)
ihs_rank<-percent_rank(abs_ihs)
ihs_norm<-qnorm((rank(ihs_rank,na.last="keep")-0.5)/sum(!is.na(ihs_rank)))
ihs_p<-pnorm(ihs_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


h12_rank<-percent_rank(sel_table$smoothed_h12)
h12_norm<-qnorm((rank(h12_rank,na.last="keep")-0.5)/sum(!is.na(h12_rank)))
h12_p<-pnorm(h12_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)

abs_sxpehh<-abs(sel_table$smoothed_sxpehh)
xpehh_rank<-percent_rank(abs_sxpehh)
xpehh_norm<-qnorm((rank(xpehh_rank,na.last="keep")-0.5)/sum(!is.na(xpehh_rank)))
xpehh_p<-pnorm(xpehh_norm, mean = 0, sd = 2, lower.tail = FALSE, log.p = FALSE)


dfp<-data.frame(pi_p, td_p, ihs_p, h12_p, xpehh_p)

#for now only keep certain statistics
dfv<-sel_table[keeps]
dfv<-na.omit(dfv)
dfp<-na.omit(dfp)

#calculate a cov matrix on the pvalues
cov_matrix<-CovNAMcd(dfp, alpha = 0.75, nsamp = 100000)

#decorelatced composite test
dcms<-DCMS(dfv, 
           subset = 1:nrow(dfv),
           S=cov_matrix$cov,
           column.nums=3:7,
           dfp,
           column.nums.p =1:5
)

plot(dcms,
     pch=19,
     cex=0.5,
     col="blue"
)

out_data<-data.frame(dfv, dfp, dcms)

#write.csv(out_data, "nw_vs_wa_dcms.csv", row.names=FALSE, quote=FALSE)

#make bed
sig_quantile<-quantile(out_data$dcms, 0.995)
sig_snps<-subset(out_data, dcms >= sig_quantile)
cols<-c("chrom", "pos", "pos")
sig_snps<-sig_snps[cols]
write.table(sig_snps, "nw_vs_wa_dcms_sig.bed", sep="\t", row.names=FALSE, quote=FALSE)


