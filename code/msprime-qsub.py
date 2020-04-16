import msprime
import os
import subprocess
import sys
from random import randint

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
    #print(tree.draw(format="unicode"))

    with open("results/msprime/{}/chr1_{}_rep_{}_seed_{}.vcf".format(pop, pop, iteration, seed), "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=2)
########################################

#calcualted by est_pop_params.py
ne = { "brazil"      : 18293,
       "niger"       : 19971,
       "senegal"     : 22421,
       "tanzania"    : 52234 }

n_samples = { "brazil"   : 48*2,
              "senegal"  : 25*2,
              "niger"    : 10*2,
              "tanzania" : 58*2  }

seed = randint(0,1e6)
sim_tree(iteration, pop, 88e-9, 3.4e-8, 8.1e-9, seed)



