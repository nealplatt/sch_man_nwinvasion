#treemix

mkdir results/treemix

#remove rodhaini and margrebowiei
vcftools \
    --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
    --mac 2 \
    --remove-indv ERX284221 \
    --remove-indv ERR310938 \
    --remove-indv Sro_female_1.1_CCATCCTC \
    --remove-indv Sro_female_1.2_CCGACAAC \
    --remove-indv Sro_female_2.1_CCTAATCC \
    --remove-indv Sro_female_2.2_CCTCTATC \
    --remove-indv Sro_male_1.1_ATCATTCC \
    --remove-indv Sro_male_1.2_ATTGGCTC \
    --remove-indv Sro_male_2.1_CAAGGAGC \
    --remove-indv Sro_male_2.2_CACCTTAC \
    --recode \
    --recode-INFO-all \
    --stdout \
    >results/treemix/smv7_ex_autosomes_mansoni.vcf

#ld filter
plink \
    --vcf results/treemix/smv7_ex_autosomes_mansoni.vcf \
    --double-id \
    --allow-extra-chr \
    --indep-pairwise 250kb 1 0.20 \
    --out results/treemix/smv7_ex_autosomes_mansoni_ld

vcftools \
    --vcf results/treemix/smv7_ex_autosomes_mansoni.vcf \
    --exclude results/treemix/smv7_ex_autosomes_mansoni_ld.prune.out \
    --recode \
    --recode-INFO-all \
    --stdout \
    >results/treemix/smv7_ex_autosomes_mansoni_ld.vcf


cp results/skyline/brazil/brazil.list .
cp results/skyline/tanzania/tanzania.list .
cp results/skyline/niger/niger.list .
cp results/skyline/senegal/senegal.list .


for POP in brazil niger senegal tanzania; do
    vcftools \
        --vcf results/treemix/smv7_ex_autosomes_mansoni_ld.vcf \
        --keep results/treemix/$POP.list \
        --counts \
        --out results/treemix/smv7_ex_autosomes_mansoni_ld_$POP
done

#paste them all together
#clean up into treemix format

#By default, TreeMix assumes biallelic sites. The input le is a gzipped le that consists of a header
#with a space-delimited list of the names of populations, followed by lines containing the allele counts
#at each SNP. It is assumed that the order of the SNPs in the le is the order of the SNPs in the
#genome. The line is space delimited between populations, and the two allele within the population
#are comma-delimited. For example:

#pop1 pop2 pop3 pop4
#5,1 1,1 4,0 0,4
#3,3 0,2 2,2 0,4
#1,5 0,2 2,2 1,3
