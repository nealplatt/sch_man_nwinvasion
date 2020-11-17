import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#make sure that all stops are not gt chrom length
chr_length = {}
#genome_size = 0
with open('data/genomes/Smansoni_v7.fa.fai', 'r') as fai:
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

#set plot colors
chrom_colors={}
chrom_colors["SM_V7_1"] = "black"
chrom_colors["SM_V7_2"] = "darkgray"
chrom_colors["SM_V7_3"] = "black"
chrom_colors["SM_V7_4"] = "darkgray"
chrom_colors["SM_V7_5"] = "black"
chrom_colors["SM_V7_6"] = "darkgray"
chrom_colors["SM_V7_7"] = "black"
sig_color = "blue"

for pop in ["senegal", "niger", "tanzania", "brazil"]:
    print(pop)
    
    #open dcms csv
    df = pd.read_csv("results/dcms/{}_dcms.csv".format(pop))
    
    #plot
    fig_chrom_colors =[]

    for index, row in df.iterrows(): 
        fig_chrom_colors.append(chrom_colors[row['chrom']])
 
    #update colors of significant peaks
    for index, row in df.iterrows(): 
        if row.dcms_p<2.056684e-05:
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

    plt.scatter(df['fig_x_pos'], -np.log10(df['dcms_p']), marker =".", s = 1, c = fig_chrom_colors) 
    plt.ylabel("-log10(p)")
    #plot FDR line
    plt.axhline(y=-np.log10(2.056684e-05), color="black", linestyle=":", linewidth=1)
    #plt.text(df['fig_x_pos'][-1:], np.log10(0.1)*-1, 'FDR=0.05', horizontalalignment='right')
    plt.xticks(ticks, tick_lbls) 
    plt.title("{} DCMS".format(pop))
    plt.savefig("results/dcms/{}_dcms.png".format(pop), dpi=300) 
    plt.close()


