import msprime
import os
import subprocess
import sys
import random

pop = sys.argv[1]
iteration = sys.argv[2]
random = random.randrange(1e6)

#calcualted by watersons_theta.py
ne = { "new_world"   : 18292,
       "west_africa" : 29644,
       "east_africa" : 52234 }

n_samples = { "new_world"   : 46,
              "west_africa" : 36,
              "east_africa" : 56  }

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

    with open("results/msprime/{}/chr1_{}_rep_{}_seed_{}.vcf".format(pop, pop, iteration, seed), "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=2)
########################################

sim_tree(iteration, pop, 88.9e6, 3.4e-8, 8.1e-9, random)





#CONDA="conda activate sch_man_nwinvasion-msprime;"
#QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 3 "    
#
#for I in $(seq 1 342); do
#    for POP in new_world east_africa west_africa; do
#        MSPRIME_CMD="code/msprime_qsub.py $POP $I"       
#        echo "$CONDA $MSPRIME_CMD" | $QSUB -N $POP"_"$I -o results/msprime/logs/$POP"_"$I.log
#    done
#done
