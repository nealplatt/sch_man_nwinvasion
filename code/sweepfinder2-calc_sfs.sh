conda activate sch_man_nwinvasion-sweepfinder

cd results/sweepfinder

#convert to sweefinder2 format
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --counts2 \
    --out smv7_ex_autosomes_tmp

tail -n+2 smv7_ex_autosomes_tmp.frq.count \
   | awk -v OFS="\t" '{print $2,$6,$4,"1"}' \
   >smv7_ex_autosomes.in

echo -e 'position\tx\tn\tfolded' \
   | cat - smv7_ex_autosomes.in \
   > temp && mv temp smv7_ex_autosomes.in

#calculate genome-wide sfs
SweepFinder2 -f smv7_ex_autosomes.in smv7_ex_autosomes.SpectFile
