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

for POP in brazil niger senegal tanzania; do
    FRQ_FILE=results/treemix/smv7_ex_autosomes_mansoni_ld_"$POP".frq.count
    TMX_FILE=results/treemix/smv7_ex_autosomes_mansoni_ld_"$POP".frq.tmx
    cat $FRQ_FILE | cut -f5,6 | sed 's/[A|T|C|G]://g' | awk '{print $1","$2}' | sed 1d | sed "1 i\\$POP" >$TMX_FILE

done

paste results/treemix/smv7_ex_autosomes_mansoni_ld_*.frq.tmx \
    | gzip \
    >results/treemix/smv7_ex_autosomes_mansoni_ld.treemix.gz

#looking for format like:
#pop1 pop2 pop3 pop4
#5,1 1,1 4,0 0,4
#3,3 0,2 2,2 0,4
#1,5 0,2 2,2 1,3

treemix -i $FILE.treemix.frq.gz \
    -m $i \
    -o $FILE.$i \
    -root GoldenJackal \
    -bootstrap \
    -k 500 \
    -noss > treemix_${i}_log &

treemix -i results/treemix/smv7_ex_autosomes_mansoni_ld.treemix.gz -k 500 -noss -o out -root tanzania -m 1

wget -P bin/ https://bitbucket.org/nygcresearch/treemix/downloads/treemix-1.13.tar.gz
tar -xvzf bin/treemix-1.13.tar.gz -C bin/treemix

source("bin/treemix-1.13/src/plotting_funcs.R")
plot tree("out")



