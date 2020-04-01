import os
import shutil
import allel
import math
import yaml
import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
#from scipy import stats

os.chdir("/master/nplatt/sch_man_nwinvasion")

with open('data/pop_assign.yml') as yaml_file:
    pop_assign = yaml.load(yaml_file, Loader=yaml.FullLoader)

#-------------------------------------------------------------------------------
# get genotype info per population

#read in vcf
phased_callset=allel.read_vcf('results/phasing/auto_beagle.vcf')

#now get an index for each sample/population
samples = phased_callset["samples"]

i=0 
for sample in samples:  
pop_idxs = defaultdict(list)   
     pop_idxs[pop_assign[sample]].append(i) 
     i=i+1 

pops= list(pop_idxs.keys()) 

#get genotypes
gt=allel.GenotypeArray(phased_callset['calldata/GT'])

#now get allele count per population
ac=gt.count_alleles()

pop_ac={}
for pop in pops:
    pop_ac[pop] = gt.count_alleles(subpop=pop_idxs[pop])

#-------------------------------------------------------------------------------
################################################################################
#-------------------------------------------------------------------------------
# get accessible bases from bedfile

#initialize an list the lenght of each contig to fale
accessible_bases = {}
chrom_length = {}

with open('/master/nplatt/sch_man_nwinvasion/data/genomes/Smansoni_v7.fa.fai', 'r') as fai:
    for entry in fai:
        chrom, length, *offset = entry.rstrip().split("\t")
        chrom_length[chrom] = int(length)
        accessible_bases[chrom]=[False] * int(length)


#now read the bed
with open('data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed', 'r') as in_bed_file:
    for bed_entry in in_bed_file:
        chrom, start, stop = bed_entry.rstrip().split("\t")
        for base in range(int(start) - 1, int(stop)):
             accessible_bases[chrom][base]=True
#-------------------------------------------------------------------------------
# pop calcs

for pop in pops:
    accessible_genome_size = 0
    pi=0
    td=0
    theta=0
    mu=8.1e-9

    #now loop through each chromosome
    for chrom in list(set(phased_callset['variants/CHROM'])) :
        target_sites = phased_callset['variants/CHROM'] == chrom

     
        chr_poss = phased_callset['variants/POS'][target_sites]
        chr_acs  = pop_ac[pop][target_sites]
        chr_len  = len(accessible_bases[chrom])
        

        chr_pi = allel.sequence_diversity(chr_poss, chr_acs, start=1, stop=chr_len, is_accessible=accessible_bases[chrom])
        chr_theta = allel.watterson_theta(chr_poss, chr_acs, is_accessible=accessible_bases[chrom])
        chr_td = allel.tajima_d(chr_acs, pos=chr_poss, start=1, stop=chr_len, min_sites=3)


        theta = theta + (chr_theta * sum(accessible_bases[chrom]))
        pi = pi + (chr_pi * sum(accessible_bases[chrom]))
        td = td + (chr_td * sum(accessible_bases[chrom]))


        accessible_genome_size = accessible_genome_size + sum(accessible_bases[chrom])

    pi = pi/accessible_genome_size
    td = td/accessible_genome_size
    theta = theta/accessible_genome_size
    
    ne = theta/(4 * mu)

    print("{}: {} {} {} {}".format(pop, pi, theta, td, ne))
