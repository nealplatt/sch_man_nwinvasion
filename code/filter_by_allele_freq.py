import re
import os
import numpy as np

#open in/out files
vcf_infile  = open("results/variant_filtration/smv7_ex_autosomes.vcf", 'r')
snp_outfile = open("results/variant_filtration/rodhaini_variants.list", 'w')



for entry in vcf_infile:
    
    #ignore the headers
    if entry.find("#") != 0:
        
        #declare list to store genotypes in (and empty it)        
        gts=[]
        #get the genotype info for all samples
        chr, pos, id, ref, alt, qual, filter, site_info, fields, *gt_data = entry.split("\t")
    
        #extract the genotypes        
        for gt_datum in gt_data:
            gt=re.split('\/|\|', gt_datum.split(":")[0])
            gts.append(gt)
        
        #for each species
        rod_gts = np.array(gts[1:10])
        man_gts = np.array(gts[10:156])
       
        #count the ref and alt alleles as well as missing data
        rod={}
        for allele in "0", "1", ".":
            rod[allele] = np.count_nonzero(rod_gts == allele)
           
        man={}
        for allele in "0", "1", ".":
            man[allele] = np.count_nonzero(man_gts == allele)

        #get allele frequencies
        rod_gtd=rod["0"] + rod["1"]
        if(rod_gtd > 0):
            rod_ref_freq = rod["0"] / rod_gtd
            rod_alt_freq = rod["1"] / rod_gtd

            man_gtd=man["0"] + man["1"]
            man_ref_freq = man["0"] / man_gtd
            man_alt_freq = man["1"] / man_gtd

            #snp_outfile.write(id + "," + rod_ref_freq + "," + rod_alt_freq + "," + man_ref_freq + "," + man_alt_freq + "\n")
            #here we capture ANY variation in rodent allele freqs
            if rod_alt_freq < 1 and rod_alt_freq > 0:
                snp_outfile.write(id + "\n")
            #here we look for private alleles in man v rod (private)
            elif rod_alt_freq == 1 and man_alt_freq != 1:
                snp_outfile.write(id + "\n")
            elif rod_ref_freq == 1 and man_ref_freq != 0:
                snp_outfile.write(id + "\n")

vcf_infile.close()
snp_outfile.close()

