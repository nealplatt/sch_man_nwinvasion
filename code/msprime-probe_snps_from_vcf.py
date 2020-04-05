import msprime
from tqdm import tqdm
import os
import subprocess
import glob


#now create bed for the new "sim" chr1
with open('data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed', 'r') as in_bed:
    with open('results/msprime/sim_probes.bed', 'w') as out_bed:
        for bed_entry in in_bed:
            chrom, start, stop = bed_entry.rstrip().split("\t")
            if chrom == "SM_V7_1":
                out_bed.write("1\t{}\t{}\n".format(start, stop))

bed = 'results/msprime/sim_probes.bed'
#now loop through all of the sim vcf files to get snps at probed regions
for pop in ["new_world", "east_africa", "west_africa"]:
    out_dir = "results/msprime/{}".format(pop)
    
    sim_vcfs = glob.glob("{}/chr1_*_rep_*.vcf".format(out_dir))
    for sim_vcf in sim_vcfs:

        probed_vcf = sim_vcf.replace(".vcf", "_probed.vcf")         
        jid = "probe_{}".format(probed_vcf.split("/")[-1])
        log = "{}/logs".format(out_dir)

        vcf_cmd = "vcftools --vcf {} --bed {} --recode --recode-INFO-all --stdout >{}".format(sim_vcf, bed, probed_vcf)
        qsub_cmd =  "qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 3 -N {} -o {}".format(jid, log)
        conda_cmd = "conda activate sch_man_nwinvasion-msprime"

        cmd ="echo \"{}; {}\" | {}".format(conda_cmd, vcf_cmd, qsub_cmd)

        #run vcf cmd
        #process = subprocess.Popen(cmd.split(""),
        #                     stdout=subprocess.PIPE, 
        #                     stderr=subprocess.PIPE)
        !{cmd}
