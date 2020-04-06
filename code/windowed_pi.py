import os
import shutil
import allel
import math
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy import stats


# To Do
# - make the vcf read in in a module
# - graph fst output

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
pop_idxs = defaultdict(list)   
for sample in samples:  
     pop_idxs[pop_assign[sample]].append(i) 
     i=i+1 

pops= list(pop_idxs.keys()) 

#get genotypes
gt=allel.GenotypeArray(phased_callset['calldata/GT'])

#now get allele count per population
ac=gt.count_alleles()

#for simplicity add maf info to callset data
maf=ac[:, :2].min(axis=1)/ac[:, :2].sum(axis=1)
phased_callset['maf']=maf 

pop_ac={}
for pop in pops:
    pop_ac[pop] = gt.count_alleles(subpop=pop_idxs[pop])
    
#-------------------------------------------------------------------------------
#generate windows
window=100_000

#define an array of window start and stops
window_starts = [int(x - (window/2)) for x in phased_callset['variants/POS']]
window_stops  = [int(x + (window/2)) for x in phased_callset['variants/POS']]

#make sure that window starts are all gt 1
window_starts = [1 if i < 1 else i for i in window_starts]


#make sure that all stops are not gt chrom length
chr_length = {}
#genome_size = 0
with open('/master/nplatt/sch_man_nwinvasion/data/genomes/Smansoni_v7.fa.fai', 'r') as fai:
    for entry in fai:
        chrom, length, *offset = entry.rstrip().split("\t")
        chr_length[chrom]=int(length)
        #genome_size = genome_size + chr_length[chrom]
    
i=0
for stop in window_stops:
    chrom = phased_callset['variants/CHROM'][i]
    
    if stop > chr_length[chrom]:
        window_stops[i]=chr_length[chrom]
    i=i+1
    
windows = np.column_stack((np.array(window_starts), 
                           np.array(window_stops)))

phased_callset['windows']=windows

#-------------------------------------------------------------------------------
# fst calculations

pops=["new_world", "east_africa", "west_africa"]

#make comparisons between population
for pop in pops:

    chrom_s          = []
    pos_s            = []
    pi_s             = []
    smoothed_pi_s    = []
    pi_calc_window_s = []
    pi_count_s       = []
    z_s              = []
    p_s              = []
    reject_s         = [] 
    q_s              = [] 
    sidak            = []
    bonf             = []
    logq_s           = []

    #create empty dataframe to store data    
    headers = ["chrom", "pos", "pis", "smoothed_pi", "window", "num_snps", "zscore", "pvalue", "q", "-log10(q)"]
    df=pd.DataFrame(columns=headers) 

    #now loop through each chromosome
    for chrom in list(set(phased_callset['variants/CHROM'])) :

        #get the snps per chromsome (and relevnat other info)
        target_sites = [phased_callset['variants/CHROM'] == chrom]  

        chr_acs  = pop_ac[pop][tuple(target_sites)]
        chr_poss = phased_callset['variants/POS'][tuple(target_sites)]
        chr_wins = phased_callset['windows'][tuple(target_sites)]

        #calculate pi
        chr_pis, chr_windows, chr_n_bases, chr_n_snps = allel.windowed_diversity(chr_poss, chr_acs, windows=chr_wins, is_accessible=accessible_bases[chrom])
  
        #smooth
        chr_smoothed_pis=scipy.signal.medfilt(chr_pis, kernel_size = 101)

        #add data to lists (of data from other chroms)
        chrom_s.extend([chrom]*len(chr_poss))
        pos_s.extend(chr_poss)
        pi_s.extend(chr_pis)
        smoothed_pi_s.extend(chr_smoothed_pis)
        pi_calc_window_s.extend(chr_windows)
        pi_count_s.extend(chr_n_snps)
        

    #calculate z, p, and q stats
    z_s=scipy.stats.zscore(pi_s)
    p_s=scipy.special.ndtr(-z_s)     
    reject_s, q_s, sidak, bonf =sm.stats.multitest.multipletests(p_s, alpha=0.01, method="fdr_bh")
    logq_s=-np.log10(q_s)

    #add all info to dataframe/table 
    df = pd.DataFrame(data = [chrom_s, pos_s, pi_s, smoothed_pi_s, pi_calc_window_s, pi_count_s, z_s, p_s, q_s, logq_s]).T
    df.columns=headers

    #save to a csv file
    csv_file = "./results/genomewide-pi/{}_windowed-pi.csv".format(pop)
    df = df.sort_values(["chrom", "pos"], ascending = (True, True))
    df.to_csv(csv_file, index=False, header=True, mode='w')

    #-------------------------------------------------------------------------------
    # generate figure
   
    #get cumul positions
    cumul_start={}
    cumul_start['SM_V7_1']=0
    cumul_start['SM_V7_2']= cumul_start['SM_V7_1'] + chr_length['SM_V7_1']
    cumul_start['SM_V7_3']= cumul_start['SM_V7_2'] + chr_length['SM_V7_2']
    cumul_start['SM_V7_4']= cumul_start['SM_V7_3'] + chr_length['SM_V7_3']
    cumul_start['SM_V7_5']= cumul_start['SM_V7_4'] + chr_length['SM_V7_4']
    cumul_start['SM_V7_6']= cumul_start['SM_V7_5'] + chr_length['SM_V7_5']
    cumul_start['SM_V7_7']= cumul_start['SM_V7_6'] + chr_length['SM_V7_6']
    scanned_size = cumul_start['SM_V7_7'] + chr_length['SM_V7_7']

    #iterate over the df and get the chrom, cumul_pos, and log10pvalue
    chrom_colors={}
    chrom_colors["SM_V7_1"] = "black"
    chrom_colors["SM_V7_2"] = "darkgray"
    chrom_colors["SM_V7_3"] = "black"
    chrom_colors["SM_V7_4"] = "darkgray"
    chrom_colors["SM_V7_5"] = "black"
    chrom_colors["SM_V7_6"] = "darkgray"
    chrom_colors["SM_V7_7"] = "black"
    sig_color = "blue"

    #remove snp windows with lt 10 snps (in 100kb) and create
    # modified positions along the x axis and add color based on chrom
    fig_x_pos_s=[]
    fig_chrom_colors =[]
    for index, row in df.iterrows(): 
        fig_x_pos_s.append(row["pos"] + cumul_start[row['chrom']])
        fig_chrom_colors.append(chrom_colors[row['chrom']])
 
    ##change color of points if q value is significant
    #for index, row in df.iterrows(): 
    #    if row.q<sidak:
    #        fig_chrom_colors[index]=sig_color

    #set tick marks in the figure (where chrom labels will be)
    ticks = [ (cumul_start['SM_V7_1'] + cumul_start['SM_V7_2'] )/2,
              (cumul_start['SM_V7_2'] + cumul_start['SM_V7_3'] )/2,
              (cumul_start['SM_V7_3'] + cumul_start['SM_V7_4'] )/2,
              (cumul_start['SM_V7_4'] + cumul_start['SM_V7_5'] )/2,
              (cumul_start['SM_V7_5'] + cumul_start['SM_V7_6'] )/2,
              (cumul_start['SM_V7_6'] + cumul_start['SM_V7_7'] )/2,
              (cumul_start['SM_V7_7'] + scanned_size )/2 ]

    tick_lbls = [ "1", "2", "3", "4", "5", "6" ,"7"]

    #build the plot
    plt.scatter(fig_x_pos_s, df['smoothed_pi'], marker =".", s = 1, c = fig_chrom_colors) 
    plt.ylabel("Pi") 
    plt.xticks(ticks, tick_lbls) 
    #plt.axhline(y=np.log10(sidak)*-1, color="black", linestyle=":", linewidth=1)
    plt.title("{}".format(pop))
    plt.savefig('results/genomewide-pi/{}_windowed-pi.png'.format(pop), dpi=300) 
    plt.close()

