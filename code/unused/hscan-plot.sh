################################################################################
#make table and plot 

#calculate p-values
import scipy.stats
import statsmodels as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

for pop in ["new_world", "east_africa", "west_africa"]:

    chrom_s=[]
    pos_s=[]
    h_s=[]

    for i in range(1,8):
        chrom = "SM_V7_{}".format(i)
        
        with open("results/hscan/{}/{}_{}.hscan-out".format(pop, chrom, pop), 'r') as hscan_file:
            next(hscan_file)
            for calc in hscan_file:
                
                pos, h = calc.split("\t")
                
                chrom_s.append(chrom)
                h_s.append(h)
                pos_s.append(pos)

    chrom_s=np.array(chrom_s)
    pos_s=np.array(pos_s)
    h_s=np.array(h_s).astype(np.float)
    
    #calculate qs    
    z_s=scipy.stats.zscore(h_s)
    p_s=scipy.special.ndtr(-z_s)     
    reject_s, q_s, sidak, bonf =sm.stats.multitest.multipletests(p_s, alpha=0.01, method="fdr_bh")
    logq_s=-np.log10(q_s)

    #now add all info to pop data table
    columns = ["chrom", "pos", "h", "z", "p", "q", "-log10(q)"]
    df = pd.DataFrame(data = [chrom_s, pos_s, h_s, z_s, p_s, q_s, logq_s]).T
    df.columns=columns

    #get cumul positions
    fig_x_pos_s=[]
    for index, row in df.iterrows(): 
        fig_x_pos_s.append(int(row["pos"]) + int(cumul_start[row['chrom']]))

    df['fig_x_pos']=fig_x_pos_s

    #save data to csv file
    csv_file ="results/hscan/{}/{}_hscan.csv".format(pop, pop)
    df = df.sort_values(["fig_x_pos"], ascending = True)
    df.to_csv(csv_file, index=False, header=True, mode='w')

    #-------------------------------------------------------------------------------
    # generate figure
   
    #iterate over the df and get the chrom, cumul_pos, and log10pvalue
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
        fig_chrom_colors.append(chrom_colors[row['chrom']])
 
    #update colors of significant peaks
    for index, row in df.iterrows(): 
        if row.q<sidak:
            fig_chrom_colors[index]=sig_color
    
    #plot log data for H
    ticks = [ (cumul_start['SM_V7_1'] + cumul_start['SM_V7_2'] )/2,
              (cumul_start['SM_V7_2'] + cumul_start['SM_V7_3'] )/2,
              (cumul_start['SM_V7_3'] + cumul_start['SM_V7_4'] )/2,
              (cumul_start['SM_V7_4'] + cumul_start['SM_V7_5'] )/2,
              (cumul_start['SM_V7_5'] + cumul_start['SM_V7_6'] )/2,
              (cumul_start['SM_V7_6'] + cumul_start['SM_V7_7'] )/2,
              (cumul_start['SM_V7_7'] + scanned_size )/2 ]

    tick_lbls = [ "1", "2", "3", "4", "5", "6" ,"7"]

    plt.scatter(df['fig_x_pos'], df['-log10(q)'], marker =".", s = 1, c = fig_chrom_colors) 
    plt.ylabel("-log10(q)")
    plt.axhline(y=np.log10(sidak)*-1, color="black", linestyle=":", linewidth=1)
    plt.xticks(ticks, tick_lbls) 
    plt.title("{} H".format(pop))
    plt.savefig("results/hscan/{}/{}_hscan.png".format(pop, pop), dpi=300) 
    plt.close()

