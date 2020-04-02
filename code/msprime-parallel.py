import msprime
from tqdm import tqdm                                                                                                                                 
import os
import subprocess
import multiprocessing

#################################
def sim_tree(iteration, pop, length, recomb_rate, mut_rate, seed):
    tree_sequence = msprime.simulate( sample_size=n_samples[pop], 
        Ne=ne[pop], 
        length=88.9e6, 
        recombination_rate=3.4e-8,
        mutation_rate=8.1e-9,
        random_seed=12345)

    tree = tree_sequence.first()
    print("iter {}: {} completed".format(iteration, pop))
    #print(tree.draw(format="unicode"))

    with open("results/sch_man_nwinvasion/msprime/{}/chr1_{}_rep_{}.vcf".format(pop, pop, iteration), "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=2)
########################################


#calcualted by watersons_theta.py
ne = { "new_world"   : 18292,
       "west_africa" : 29644,
       "east_africa" : 52234 }

n_samples = { "new_world"   : 46,
              "west_africa" : 36,
              "east_africa" : 56  }

#set up a tuple for all parameters
params=[]
for pop in ["new_world", "east_africa", "west_africa"]:
    out_dir = "results/msprime/{}".format(pop)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i in range(0, 11): 
        params.append((i, pop, 88.9e6, 3.4e-8, 8.1e-9, 12345)) 

#run sims in parallel
with multiprocessing.Pool(processes=11) as pool:
    pool.starmap(sim_tree, params)
      

#-------------------------------------------------------------------------------
# extract only "probed" regions from the simuluated chromosomes/vcfs

#print new genome index file 
#with open("results/sch_man_nwinvasion/msprime/sim_chr1.genome", 'w') as bed_genome:
#    bed_genome.write("1\t1\t8890000")

#now create bed for the new "sim" chr1
#with open('data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed', 'r') as in_bed:
#    with open('results/sch_man_nwinvasion/msprime/sim_probes.bed', 'w') as out_bed:
#        for bed_entry in in_bed:
#            chrom, start, stop = bed_entry.rstrip().split("\t")
#            if chrom == "SM_V7_1":
#                out_bed.write("1\t{}\t{}\n".format(start, stop))

#now loop through all of the sim vcf files to get snps at probed regions
#for pop in ["new_world", "east_africa", "west_africa"]:
#    out_dir = "results/msprime/{}".format(pop)
#    for i in tqdm(range(0, 11), desc="{} intersect".format(pop):
#
#        sim_vcf = "results/msprime/new_world/chr1_{}_rep_{}.vcf".format(pop, i)
#        int_vcf = "results/msprime/new_world/chr1_probed_{}_rep_{}.vcf".format(pop, i)
#        bed = 'results/sch_man_nwinvasion/msprime/sim_probes.bed'
#
#        vcf_cmd = "vcftools --vcf {} --bed {} --recode --recode-INFO-all --stdout >{}".format(sim_vcf, bed, int_vcf)
#
#        #run vcf cmd
#        process = subprocess.Popen(vcf_cmd.split(""),
#                             stdout=subprocess.PIPE, 
#                             stderr=subprocess.PIPE)


