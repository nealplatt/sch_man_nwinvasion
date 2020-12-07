conda activate sch_man_nwinvasion-jupyter

mkdir results/mito_network

#get the mito vcf without outgroups
vcftools \
    --vcf results/variant_filtration/smv7_ex_mito.vcf \
    --remove-indv Sro_female_1.1_CCATCCTC \
    --remove-indv Sro_female_1.2_CCGACAAC \
    --remove-indv Sro_female_2.1_CCTAATCC \
    --remove-indv Sro_female_2.2_CCTCTATC \
    --remove-indv Sro_male_1.1_ATCATTCC \
    --remove-indv Sro_male_1.2_ATTGGCTC \
    --remove-indv Sro_male_2.1_CAAGGAGC \
    --remove-indv Sro_male_2.2_CACCTTAC \
    --remove-indv Sm.BR_PdV.1409_rep \
    --remove-indv Sm.BR_PdV.1475_rep \
    --remove-indv Sm.BR_PdV.2406_rep \
    --remove-indv ERX284221 \
    --remove-indv ERR310938 \
    --remove-indv ERR119614 \
    --recode \
    --recode-INFO-all \
    --stdout \
    >results/mito_network/mito.vcf

#convert seq data to nexus
bin/vcf2phylip/vcf2phylip.py -i results/mito_network/mito.vcf -n


Sro_female_1.1_CCATCCTC
Sro_female_1.2_CCGACAAC
Sro_female_2.1_CCTAATCC
Sro_female_2.2_CCTCTATC
Sro_male_1.1_ATCATTCC
Sro_male_1.2_ATTGGCTC
Sro_male_2.1_CAAGGAGC
Sro_male_2.2_CACCTTAC
ERX284221 - margrebowiei
ERR310938 - sro-burundi
ERR046038 -sm.pr
ERR539847 -sm.guad
ERR539848 -sm.guad
ERR103050 -sm.cam
ERR103049 -sm.sengal
ERR119615 -sm.uganda
ERR997461 -sm.uganda
