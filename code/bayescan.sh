conda activate sch_man_nwinvasion-bayescan

cd /master/nplatt/sch_man_nwinvasion

mkdir results/bayescan
cd results/bayescan

#get bayescan2.1.0
#http://cmpg.unibe.ch/software/BayeScan/download.html

# the VCF_to_BAYESCAN spid file is saved in data
cp ../../VCF_to_BAYESCAN.spid .

#create a pop file
#<indiv>\t<pop>


#get cleaned up vcf file
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --keep pops \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >smv7_ex_autosomes_maf05.vcf


for POP in brazil niger senegal tanzania; do

    ##create spid file
    sed "s/pops/$POP.pops/" VCF_to_BAYESCAN.spid >"$POP"_VCF_to_BAYESCAN.spid

    ##convert vcf to bayescan fmt
    java -jar ~/sch_man_nwinvasion/bin/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
        -inputfile smv7_ex_autosomes_maf05.vcf \
        -inputformat vcf \
        -outputfile "$POP"_smv7_ex_autosomes_maf05.bayescan \
        -outputformat GESTE_BAYE_SCAN \
        -spid "$POP"_VCF_to_BAYESCAN.spid

    #run bayescan
    CMD="conda activate sch_man_nwinvasion-bayescan; \
            ~/sch_man_nwinvasion/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
                "$POP"_smv7_ex_autosomes_maf05.bayescan \
                -threads 12 \
                -o "$POP"_smv7_ex_autosomes_maf05_bayescan_pr10_pi10k_bin25k_ngen100k_npb50_thin20 \
                -pr_odds 10 \
                -pilot 10000 \
                -burn 25000 \
                -nbp 50 \
                -n 100000 \
                -thin 20"

    echo $CMD | qsub -V -cwd -S /bin/bash -q all.q -j y -N "$POP"_bayescan -o "$POP"_bayescan.log -pe smp 12
done

