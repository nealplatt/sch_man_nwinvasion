#vcf2hscan.py                                                                                                                                         
import vcf                                                                                                                                            
import sys                                                                                                                                            
import os                                                                                                                                             
from tqdm import tqdm                                                                                                                                          

in_vcf    = sys.argv[1]                                                                                                                               
out_hscan = sys.argv[2]                                                                                                                               
                                                                                                                                                      
#delete outfile if it exists                                                                                                                          
if os.path.exists(out_hscan):                                                                                                                         
    os.remove(out_hscan)                                                                                                                              
                                                                                                                                                      
#read in the vcf                                                                                                                                      
vcf_reader = vcf.Reader(open(in_vcf, 'r'))     

outlines=[]
#start converting
for record in tqdm(vcf_reader, desc=out_hscan):
    
    pos   = record.POS
    ref   = record.REF
    alt   = str(record.ALT[0])
    
    homo_gts=[]
    for sample in record.samples:
        #get unphased genotype
        gt = sample['GT'].replace('|', '/') 
        
        homo_ref = ['0/.', './0', '0/0']
        homo_alt = ['1/.', './1', '1/1']
        hetero   = ['0/1','1/0']
        missing  = ['./.']
        
        if gt in missing:
            homo_gt = 'N'
        elif gt in homo_ref:
            homo_gt = ref 
        elif gt in homo_alt:
            homo_gt = alt
        elif gt in hetero:
            homo_gt = '.'
        else:
            print(gt + " does not appear to be a standard genotype at " + chr + ":" + str(pos))
            sys.exit()
        
        homo_gts.append(homo_gt)

    outline = ','.join([str(pos)] + homo_gts)
    outlines.append(outline)

