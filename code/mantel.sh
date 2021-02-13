#https://github.com/BGI-shenzhen/VCF2Dis


conda activate sch_man_nwinvasion-mantel

cd ~/sch_man_nwinvasion/results/mantel

#get target sample list
cut -f1 locations.tsv >exome_samples.list


#only keep vcf wtih target exome samples
vcftools \
    --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
    --keep exome_samples.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >smv7_ex_autosomes_exomes.vcf


#update the header (max name length = 12 char).
cat smv7_ex_autosomes_exomes.vcf \
    | head -n 5000 \
    | grep "#" \
    | tail -n1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >smv7_ex_autosomes_exomes.list


seq 1 135 >reheader.list

paste reheader.list smv7_ex_autosomes_exomes.list >sample_name_key.list

bcftools reheader \
    --samples reheader.list \
    -o smv7_ex_autosomes_exomes_renamed.vcf \
    smv7_ex_autosomes_exomes.vcf

#calculate genetic distance (p-distance)
~/sch_man_nwinvasion/code/VCF2Dis/bin/VCF2Dis \
    -InPut smv7_ex_autosomes_exomes_renamed.vcf \
    -OutPut smv7_ex_autosomes_exomes_renamed_genetic_distance_matrix.tsv

sed -i 1d smv7_ex_autosomes_exomes_renamed_genetic_distance_matrix.tsv


