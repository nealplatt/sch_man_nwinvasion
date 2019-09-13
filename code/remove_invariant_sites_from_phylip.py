### RNPlatt
### 05 Sept 2019
### USAGE: remove_invariant_sites_from_phylip.py <list of invariants> <phylip> <outfile>
### 
### This will take a list of invariant sites predetermined by raxml in a single
###     site per line file and the original phylip file used in raxml and return
###     a phylip file with all the invariant sites removed.
###
### REQUIRES: numpy
###
import sys
import numpy as np
#-------------------- GET CMD LINE OPTIONS--------------------------------------
invariant_file=sys.argv[1]
phylip_file=sys.argv[2]
out_file=sys.argv[3]


#-------------------- GET INVARIANT SITES --------------------------------------
#read in the file of invariant sites generated by raxml
inv_sites_infile=open(invariant_file, "r")
inv_sites=inv_sites_infile.readlines()

i=0
#update each site to 0-based (python) index from 1-based (raxml)
for entry in inv_sites:
    inv_sites[i]=int(entry.rstrip())-1
    i=i+1


#-------------------- GET SEQUENCE DATA AND TRIM -------------------------------
#read in the dat from the untrimmed phylip file
phylip_infile=open(phylip_file, "r")
phylip_data=phylip_infile.readlines()

#get the num samples and sites from the header line
num_samples, num_sites=phylip_data[0].rstrip('\n').split(" ")

#cycle through the seqeunce data into two seperate lists
# sample ids are in a list
# sequences is a 2d list with each base a seperate position
sample_ids=[]
sequences=[]
for entry in phylip_data[1:]:
    sample_id, sequence = entry.rstrip('\n').split()
    sample_ids.append(sample_id)
    sequences.append(list(sequence))

#convert to 2d array
sequences=np.array(sequences)

#trim invariant sites
trimmed_seqs=np.delete(sequences, inv_sites, 1)

#now turn into strings
seq_strings=[]
for trimmed_seq in trimmed_seqs:

    #convert trimmed array/list to a string
    as_string=''.join(list(trimmed_seq))
    
    #add to new list of trimmed sequences
    seq_strings.append(as_string)


#-------------------- CREATING THE OUTPUT FILE ---------------------------------
#create an output file
trimmed_phylip_outfile=open(out_file, "w")
num_sites_after_trimming=len(seq_strings[0])

#print header line
trimmed_phylip_outfile.write(str(num_samples) + " " + str(num_sites_after_trimming) + "\n")

#print trimmed info to outfile
i=0
for sample_id in sample_ids:
    trimmed_phylip_outfile.write(sample_id + " " + seq_strings[i] + "\n")
    i=i+1


#-------------------- CLOSING FILES  -------------------------------------------
inv_sites_infile.close()
phylip_infile.close()
trimmed_phylip_outfile.close()

