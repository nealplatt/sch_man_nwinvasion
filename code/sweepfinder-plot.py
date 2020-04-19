################################################################################
#make table and plot 

#calculate p-values
import scipy.stats
import statsmodels.api as sm
import statsmodels as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

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

for pop in ["niger", "senegal", "brazil", "tanzania"]:
    print(pop)
    chrom_s=[]
    pos_s=[]
    lr_s=[]

    for i in range(1,8):
        chrom = "SM_V7_{}".format(i)
        
        with open("results/sweepfinder/{}/{}_{}.sw2out".format(pop, chrom, pop), 'r') as sf2_file:
            next(sf2_file)
            for calc in sf2_file:
                
                pos, lr, alpha = calc.rstrip().split("\t")
                
                chrom_s.append(chrom)
                lr_s.append(lr)
                pos_s.append(round(float(pos)))

    chrom_s=np.array(chrom_s)
    pos_s=np.array(pos_s)
    lr_s=np.array(lr_s).astype(np.float)
    
    ##calculate qs    
    #z_s=scipy.stats.zscore(h_s)
    #p_s=scipy.special.ndtr(-z_s)     
    #reject_s, q_s, sidak, bonf =sm.stats.multitest.multipletests(p_s, alpha=0.01, method="fdr_bh")
    #logq_s=-np.log10(q_s)

    #calculate max neutral h score
    lr_max = 0.0
    sim_files = glob.glob("results/sweepfinder/{}-sim/*sw2out".format(pop))     
    for sim_file in sim_files:
        with open(sim_file, 'r') as in_sim_file:
            next(in_sim_file)
            for sim_entry in in_sim_file:
                 x, sim_lr, sim_alpha = sim_entry.rstrip().split("\t")
                 if (float(sim_lr)>float(lr_max)):
                    lr_max=float(sim_lr)

    #calculate qs (with simulated hmax
    z_s=scipy.stats.zscore(np.insert(lr_s, 0, lr_max))
    p_s=scipy.special.ndtr(-z_s)   
    reject_s, q_s, sidak, bonf =sm.stats.multitest.multipletests(p_s, alpha=0.01, method="fdr_bh")
    logq_s=-np.log10(q_s)

    #now add all info to pop data table
    columns = ["chrom", "pos", "lr", "z", "p", "q", "-log10(q)"]
    df = pd.DataFrame(data = [chrom_s, pos_s, lr_s, z_s[1:], p_s[1:], q_s[1:], logq_s[1:]]).T
    df.columns=columns

    #get cumul positions
    fig_x_pos_s=[]
    for index, row in df.iterrows(): 
        fig_x_pos_s.append(int(row["pos"]) + int(cumul_start[row['chrom']]))

    df['fig_x_pos']=fig_x_pos_s

    #save data to csv file
    csv_file ="results/sweepfinder/{}/{}_sw2.csv".format(pop, pop)
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
        if row.q<sidak and row.lr>lr_max:
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
    #plot FDR line
    plt.axhline(y=np.log10(sidak)*-1, color="black", linestyle=":", linewidth=1)
    #plt.text(df['fig_x_pos'][-1:], np.log10(sidak)*-1, 'FDR=0.01', horizontalalignment='right')
    #plot max_h line)
    plt.axhline(y=logq_s[0], color="black", linestyle="--", linewidth=1)
    #plt.text(df['fig_x_pos'][-1:], logq_s[0], 'max(h) neutral', horizontalalignment='right')
    plt.xticks(ticks, tick_lbls) 
    plt.title("{} LR".format(pop))
    plt.savefig("results/sweepfinder/{}_sweepfinder_q.png".format(pop, pop), dpi=300) 
    plt.close()


    plt.scatter(df['fig_x_pos'], df['lr'], marker =".", s = 1, c = fig_chrom_colors) 
    plt.ylabel("lr")
    #plt.axhline(y=min(h_s[reject_s[1:]]), color="black", linestyle=":", linewidth=1)
    plt.axhline(y=lr_max, color="black", linestyle="--", linewidth=1)
    plt.xticks(ticks, tick_lbls) 
    plt.title("{} LR".format(pop))
    plt.savefig("results/sweepfinder/{}_sweepfinder_lr.png".format(pop, pop), dpi=300) 
    plt.close()



