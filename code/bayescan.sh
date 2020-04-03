conda activate sch_man_nwinvasion-bayescan
cd /master/nplatt/sch_man_nwinvasion/results/bayescan

#get bayescan2.1.0

#get cleaned up vcf file
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --maf 0.05 \
    --keep ../lists/new_world.list \
    --keep ../lists/east_africa.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >smv7_ex_autosomes_BR_EA_maf05.vcf

#convert vcf to bayescan fmt
java -jar ~/sch_man_nwinvasion/bin/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
    -inputfile smv7_ex_autosomes_BR_EA_maf05.vcf \
    -inputformat vcf \
    -outputfile smv7_ex_autosomes_BR_EA_maf05.bayescan \
    -outputformat GESTE_BAYE_SCAN \
    -spid VCF_to_BAYESCAN.spid

#run bayescan
CMD="conda activate sch_man_nwinvasion-bayescan; ~/sch_man_nwinvasion/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
        smv7_ex_autosomes_BR_EA_maf05.bayescan \
        -threads 12 \
        -o smv7_ex_autosomes_BR_EA_maf05_bayescan_pr10_pi10k_bin25k_ngen100k_npb50_thin20 \
        -pr_odds 10 \
        -pilot 10000 \
        -burn 25000 \
        -nbp 50 \
        -n 100000 \
        -thin 20"

echo $CMD | qsub -V -cwd -S /bin/bash -q all.q -j y -N bayescan -o bayescan.log -pe smp 12

#-------------------------------------------------------------------------------
#get cleaned up vcf file
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --maf 0.05 \
    --keep ../lists/new_world.list \
    --keep ../lists/west_africa.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >smv7_ex_autosomes_BR_WA_maf05.vcf

#convert vcf to bayescan fmt
java -jar ~/sch_man_nwinvasion/bin/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
    -inputfile smv7_ex_autosomes_BR_WA_maf05.vcf \
    -inputformat vcf \
    -outputfile smv7_ex_autosomes_BR_WA_maf05.bayescan \
    -outputformat GESTE_BAYE_SCAN \
    -spid VCF_to_BAYESCAN.spid

#run bayescan
CMD="conda activate sch_man_nwinvasion-bayescan; ~/sch_man_nwinvasion/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
        smv7_ex_autosomes_BR_WA_maf05.bayescan \
        -threads 12 \
        -o smv7_ex_autosomes_BR_WA_maf05_bayescan_pr10_pi10k_bin25k_ngen100k_npb50_thin20 \
        -pr_odds 10 \
        -pilot 10000 \
        -burn 25000 \
        -nbp 50 \
        -n 100000 \
        -thin 20"

echo $CMD | qsub -V -cwd -S /bin/bash -q all.q -j y -N bayescan -o bayescan.log -pe smp 12

#-------------------------------------------------------------------------------
#get cleaned up vcf file
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --maf 0.05 \
    --keep ../lists/east_africa.list \
    --keep ../lists/west_africa.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >smv7_ex_autosomes_EA_WA_maf05.vcf

#convert vcf to bayescan fmt
java -jar ~/sch_man_nwinvasion/bin/PGDSpider_2.1.1.5/PGDSpider2-cli.jar \
    -inputfile smv7_ex_autosomes_EA_WA_maf05.vcf \
    -inputformat vcf \
    -outputfile smv7_ex_autosomes_EA_WA_maf05.bayescan \
    -outputformat GESTE_BAYE_SCAN \
    -spid VCF_to_BAYESCAN.spid

#run bayescan
CMD="conda activate sch_man_nwinvasion-bayescan; ~/sch_man_nwinvasion/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
        smv7_ex_autosomes_EA_WA_maf05.bayescan \
        -threads 12 \
        -o smv7_ex_autosomes_EA_WA_maf05_bayescan_pr10_pi10k_bin25k_ngen100k_npb50_thin20 \
        -pr_odds 10 \
        -pilot 10000 \
        -burn 25000 \
        -nbp 50 \
        -n 100000 \
        -thin 20"

echo $CMD | qsub -V -cwd -S /bin/bash -q all.q -j y -N bayescan -o bayescan.log -pe smp 12

