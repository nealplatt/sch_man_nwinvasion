import msprime
from tqdm import tqdm                                                                                                                                 
import os

#calcualted by watersons_theta.py
ne = { "new_world"   : 18292,
       "west_africa" : 29644,
       "east_africa" : 52234 }

n_samples = { "new_world"   : 46,
              "west_africa" : 36,
              "east_africa" : 56  }

#for pop in ["new_world", "east_africa", "west_africa"]:
#    out_dir = "results/msprime/{}".format(pop)
#    if not os.path.exists(out_dir):
#        os.makedirs(out_dir)
#
#    tree_sequence = msprime.simulate( sample_size=n_samples[pop], 
#                                      Ne=ne[pop], 
#                                      length=1e6, 
#                                      recombination_rate=3.4e-8,
#                                      mutation_rate=8.1e-9,
#                                      random_seed=12345)
#    tree = tree_sequence.first()
#    print(tree.draw(format="unicode"))
#
#    with open("{}/{}_rep_{}.vcf".format(out_dir, pop, i), "w") as vcf_file:
#        tree_sequence.write_vcf(vcf_file, ploidy=2)

for pop in ["new_world", "east_africa", "west_africa"]:
    out_dir = "results/msprime/{}".format(pop)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #run the equivilant of 100 genomes worth of simulations
    for i in tqdm(range(0, 343), desc="{} sim".format(pop)):
        for pop in ["new_world", "east_africa", "west_africa"]:
            tree_sequence = msprime.simulate( sample_size=n_samples[pop], 
                                              Ne=ne[pop], 
                                              length=88.9e6, 
                                              recombination_rate=3.4e-8,
                                              mutation_rate=8.1e-9,
                                              random_seed=12345)
            tree = tree_sequence.first()
            print(tree.draw(format="unicode"))

            with open("results/sch_man_nwinvasion/msprime/{}/chr1_{}_rep_{}.vcf".format(pop, pop, i), "w") as vcf_file:
                tree_sequence.write_vcf(vcf_file, ploidy=2)


