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
from scipy import signal


# To Do
# - make the vcf read in in a module
# - graph fst output

os.chdir("/master/nplatt/sch_man_nwinvasion")

with open('data/pop_assign.yml') as yaml_file:
    pop_assign = yaml.load(yaml_file, Loader=yaml.FullLoader)

#-----------------------------------
# get lengths from cumul positions
#make sure that all stops are not gt chrom length
chr_length = {}
#genome_size = 0
with open('/master/nplatt/sch_man_nwinvasion/data/genomes/Smansoni_v7.fa.fai', 'r') as fai:
    for entry in fai:
        chrom, length, *offset = entry.rstrip().split("\t")
        chr_length[chrom]=int(length)

    cumul_start={}
    cumul_start['SM_V7_1']=0
    cumul_start['SM_V7_2']= cumul_start['SM_V7_1'] + chr_length['SM_V7_1']
    cumul_start['SM_V7_3']= cumul_start['SM_V7_2'] + chr_length['SM_V7_2']
    cumul_start['SM_V7_4']= cumul_start['SM_V7_3'] + chr_length['SM_V7_3']
    cumul_start['SM_V7_5']= cumul_start['SM_V7_4'] + chr_length['SM_V7_4']
    cumul_start['SM_V7_6']= cumul_start['SM_V7_5'] + chr_length['SM_V7_5']
    cumul_start['SM_V7_7']= cumul_start['SM_V7_6'] + chr_length['SM_V7_6']
    scanned_size = cumul_start['SM_V7_7'] + chr_length['SM_V7_7']


#-------------------------------------------------------------------------------
# get genotype info per population

#read in vcf
callset=allel.read_vcf('results/variant_filtration/smv7_ex_autosomes.vcf')

#now get an index for each sample/population
samples = callset["samples"]

i=0 
pop_idxs = defaultdict(list)   
for sample in samples:  
     pop_idxs[pop_assign[sample]].append(i) 
     i=i+1 

pops= list(pop_idxs.keys()) 

#get genotypes
gt=allel.GenotypeArray(callset['calldata/GT'])

#now get allele count per population
ac=gt.count_alleles()

#for simplicity add maf info to callset data
maf=ac[:, :2].min(axis=1)/ac[:, :2].sum(axis=1)
callset['maf']=maf 

pop_ac={}
for pop in pops:
    pop_ac[pop] = gt.count_alleles(subpop=pop_idxs[pop])
    
#-------------------------------------------------------------------------------
#generate windows
window=100_000

#define an array of window start and stops
window_starts = [int(x - (window/2)) for x in callset['variants/POS']]
window_stops  = [int(x + (window/2)) for x in callset['variants/POS']]

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
    chrom = callset['variants/CHROM'][i]
    
    if stop > chr_length[chrom]:
        window_stops[i]=chr_length[chrom]
    i=i+1
    
windows = np.column_stack((np.array(window_starts), 
                           np.array(window_stops)))

callset['windows']=windows

#-------------------------------------------------------------------------------
# fst calculations

pops = ["brazil", "tanzania", "niger", "senegal" ]

idx_comps = {"brazil":   [pop_idxs["brazil"],   pop_idxs["tanzania"] + pop_idxs["niger"]  + pop_idxs["senegal"] ],
             "tanzania": [pop_idxs["tanzania"], pop_idxs["brazil"]   + pop_idxs["niger"]  + pop_idxs["senegal"] ],
             "niger":    [pop_idxs["niger"],    pop_idxs["tanzania"] + pop_idxs["brazil"] + pop_idxs["senegal"] ],
             "senegal":  [pop_idxs["senegal"],  pop_idxs["tanzania"] + pop_idxs["niger"]  + pop_idxs["brazil"] ]}

#make comparisons between population
for pop in idx_comps.keys():
    pop1_idx = idx_comps[pop][0]
    pop2_idx = idx_comps[pop][1]

    fst_s             = []
    fst_calc_window_s = []
    fst_count_s       = []

    #create empty dataframe to store data    
    headers = ["chrom", "pos", "fst", "smoothed_fst", "window", "num_snps", "zscore", "pvalue", "-log10(p)"]
    df=pd.DataFrame(columns=headers) 

    #now loop through each chromosome
    for chrom in list(set(callset['variants/CHROM'])) :
        target_sites = np.logical_and( callset['maf'] < 0.05, 
                                       callset['variants/CHROM'] == chrom )  

        chr_gts  = gt[target_sites]
        chr_poss = callset['variants/POS'][target_sites]
        chr_wins = callset['windows'][target_sites]

        
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
        smoothed_fsts=signal.medfilt(fsts, kernel_size = 101)

        #calculate pvalue
        zs = stats.zscore(fsts)
        ps = stats.norm.cdf(zs)
        logps = [-1*np.log10(p) for p in ps] 

        #add data to dataframe/table
        data = list(zip([chrom]*len(fsts), chr_poss, fsts, smoothed_fsts, fst_calc_windows, fst_counts, zs, ps, logps))
        chr_df=pd.DataFrame(data, columns=headers)
        df = df.append(chr_df)

    #add cumul positions
    fig_x_pos_s=[]
    for index, row in df.iterrows(): 
        fig_x_pos_s.append(int(row["pos"]) + int(cumul_start[row['chrom']]))

    df['fig_x_pos']=fig_x_pos_s

    #save data to csv file
    csv_file = "./results/fst_per_window/{}_vs_all_windowed_fst.csv".format(pop)
    df = df.sort_values(["fig_x_pos"], ascending = True)
    df.to_csv(csv_file, index=False, header=True, mode='w')


#-------------------------------------------------------------------------------
# generate figure
   
for pop in pops:
    csv_file = "./results/fst_per_window/{}_vs_all_windowed_fst.csv".format(pop)
    df = pd.read_csv(csv_file,sep=",")

    #remove sites lt 10 snps?

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
    chrom_colors["SM_V7_1"] = "black"
    chrom_colors["SM_V7_2"] = "darkgray"
    chrom_colors["SM_V7_3"] = "black"
    chrom_colors["SM_V7_4"] = "darkgray"
    chrom_colors["SM_V7_5"] = "black"
    chrom_colors["SM_V7_6"] = "darkgray"
    chrom_colors["SM_V7_7"] = "black"
    sig_color = "blue"

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
    plt.ylabel("smoothed_fst") 
    plt.xticks(ticks, tick_lbls) 
    plt.title("Fst ({})".format(pop))
    plt.savefig('results/fst_per_window/{}_v_all.png'.format(pop), dpi=300) 
    plt.close()

