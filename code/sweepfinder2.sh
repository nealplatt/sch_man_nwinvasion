conda activate sch_man_nwinvasion-sweepfinder

cd /master/nplatt/sch_man_nwinvasion

mkdir results/sweepfinder

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

#run sweepfinder on the actual data
code/sweepfinder2-run_with_real_data.sh

#get the msprime data

#run sweepfinder on simulated data
code/sweepfinder2-run_with_simulated_data.sh

#coallate and plot



