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

comps = [ [ "new_world",   "west_africa" ],
          [ "east_africa", "west_africa" ],
          [ "new_world",   "east_africa" ] ]

#make comparisons between population
for comp in comps:
    pop1 = comp[0]
    pop2 = comp[1]

    fst_s             = []
    fst_calc_window_s = []
    fst_count_s       = []

    #create empty dataframe to store data    
    headers = ["chrom", "pos", "fst", "smoothed_fst", "window", "num_snps", "zscore", "pvalue", "-log10(p)"]
    df=pd.DataFrame(columns=headers) 

    #now loop through each chromosome
    for chrom in list(set(phased_callset['variants/CHROM'])) :
        target_sites = np.logical_and( phased_callset['maf'] < 0.05, 
                                       phased_callset['variants/CHROM'] == chrom )  

        chr_gts  = gt[target_sites]
        chr_poss = phased_callset['variants/POS'][target_sites]
        chr_wins = phased_callset['windows'][target_sites]
        pop1_idx = pop_idxs[pop1]
        pop2_idx = pop_idxs[pop2]
        
        fsts, fst_calc_windows, fst_counts =allel.windowed_weir_cockerham_fst(chr_poss, chr_gts, subpops=[pop1_idx, pop2_idx], windows=chr_wins )

        #get rid of nan values
        useful_values = np.logical_and( np.isfinite(fsts), fst_counts>=10) 

        fsts = fsts[useful_values]
        fst_calc_windows = fst_calc_windows[useful_values]
        fst_counts = fst_counts[useful_values]
        chr_poss = chr_poss[useful_values]

        #set negative fst values to 0
        i=0
        for fst in fsts:
            if fst <0:
                fsts[i]=0
            i=i+1        
        
        #smooth
        smoothed_fsts=scipy.signal.medfilt(fsts, kernel_size = 101)

        #calculate pvalue
        zs = scipy.stats.zscore(fsts)
        ps = scipy.stats.norm.cdf(zs)
        logps = [-1*math.log10(p) for p in ps] 

        #add data to dataframe/table
        data = list(zip([chrom]*len(fsts), chr_poss, fsts, smoothed_fsts, fst_calc_windows, fst_counts, zs, ps, logps))
        chr_df=pd.DataFrame(data, columns=headers)
        df = df.append(chr_df)
        

    #save data to csv file
    csv_file = "./results/fst_per_window/{}_vs_{}_windowed_fst.csv".format(pop1, pop2)
    df = df.sort_values(["chrom", "pos"], ascending = (True, True))
    csv_file = "./results/fst_per_window/{}_vs_{}_windowed_fst.csv".format(pop1, pop2)
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
    fig_x_pos=[]
    fig_chrom_colors =[]

    chrom_colors={}
    chrom_colors["SM_V7_1"] = "red"
    chrom_colors["SM_V7_2"] = "blue"
    chrom_colors["SM_V7_3"] = "red"
    chrom_colors["SM_V7_4"] = "blue"
    chrom_colors["SM_V7_5"] = "red"
    chrom_colors["SM_V7_6"] = "blue"
    chrom_colors["SM_V7_7"] = "red"

    for index, row in df.iterrows(): 
        fig_x_pos.append(row["pos"] + cumul_start[row['chrom']])
        fig_chrom_colors.append(chrom_colors[row['chrom']])
 
    
    ticks = [ (cumul_start['SM_V7_1'] + cumul_start['SM_V7_2'] )/2,
              (cumul_start['SM_V7_2'] + cumul_start['SM_V7_3'] )/2,
              (cumul_start['SM_V7_3'] + cumul_start['SM_V7_4'] )/2,
              (cumul_start['SM_V7_4'] + cumul_start['SM_V7_5'] )/2,
              (cumul_start['SM_V7_5'] + cumul_start['SM_V7_6'] )/2,
              (cumul_start['SM_V7_6'] + cumul_start['SM_V7_7'] )/2,
              (cumul_start['SM_V7_7'] + scanned_size )/2 ]

    tick_lbls = [ "1", "2", "3", "4", "5", "6" ,"7"]

    plt.scatter(fig_x_pos, df['smoothed_fst'], marker =".", s = 1, c = fig_chrom_colors) 
    plt.ylabel("-log10(p)") 
    plt.xticks(ticks, tick_lbls) 
    plt.title("Fst ({} vs {})".format(pop1, pop2))
    plt.savefig('results/fst_per_window/{}_v_{}.png'.format(pop1, pop2), dpi=300) 
    plt.close()

