import msprime
import os
import subprocess
import sys

pop = sys.argv[1]
iteration = sys.argv[2]

#################################
def sim_tree(iteration, pop, length, recomb_rate, mut_rate, seed):
    tree_sequence = msprime.simulate( sample_size=n_samples[pop], 
        Ne=ne[pop], 
        length=88.9e6, 
        recombination_rate=3.4e-8,
        mutation_rate=8.1e-9,
        random_seed=seed)

    print("iter {}: {} completed".format(iteration, pop))
    print(tree.draw(format="unicode"))

    with open("results/sch_man_nwinvasion/msprime/{}/chr1_{}_rep_{}_seed_{}.vcf".format(pop, pop, iteration, seed), "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=2)
########################################


#calcualted by watersons_theta.py
ne = { "new_world"   : 18292,
       "west_africa" : 29644,
       "east_africa" : 52234 }

n_samples = { "new_world"   : 46,
              "west_africa" : 36,
              "east_africa" : 56  }

sim_tree(iteration, pop, 88.9e6, 3.4e-8, 8.1e-9, 12345)



