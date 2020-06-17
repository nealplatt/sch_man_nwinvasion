from collections import defaultdict
from intervaltree import IntervalTree
import pandas as pd

for pop in ["senegal", "niger", "brazil", "tanzania"]:

    print(pop)
    #designate files
    hscan_file = "results/hscan/{}/{}_hscan.csv".format(pop, pop)
    sweep_file = "results/sweepfinder/{}/{}_sw2.csv".format(pop, pop)
    fst_file = "results/fst_per_window/{}_vs_all_windowed_fst.csv".format(pop)

    #read in to dfs    
    hscan_df = pd.read_csv(hscan_file, sep=",")
    sweep_df = pd.read_csv(sweep_file, sep=",")
    fst_df   = pd.read_csv(fst_file,   sep=",")

    #merge data frames
    df = hscan_df
    df = df.merge(fst_df, how="outer", left_on=['chrom', 'pos'], right_on=['chrom', 'pos'])

    for chrom in range(1,8):
        chrom="SM_V7_{}".format(chrom)

        chrom_df       = df.loc[df['chrom'] == chrom]     
        chrom_sweep_df = sweep_df.loc[sweep_df['chrom'] == chrom]     
       
    #build interval tree of all sweepfinder values
    sweep_trees=defaultdict(lambda: IntervalTree()) 
    for index, row in sweep_df.iterrows():
        sweep_trees[row.chrom][row.pos:row.pos+1001] = row.lr

    lrs=[]
    for index, row in df.iterrows():
        lr = sorted(sweep_trees[row.chrom].at(row.pos))[0].data   
        lrs.append(lr)
    
    #add sweepfinder values to dataframe
    df["lr"]=lrs

    #subsample dataframe for desired values
    sub_df = df[["chrom", "pos", "h", "fst", "lr", "fig_x_pos_x"]]

    #get complete records
    sub_df = sub_df.dropna()

    #save as new csv for dcms
    csv_file = "./results/dcms/{}_raw_data.csv".format(pop)
    sub_df = sub_df.sort_values(["fig_x_pos_x"], ascending = True)
    sub_df.to_csv(csv_file, index=False, header=True, mode='w')













