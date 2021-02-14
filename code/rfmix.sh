conda activate sch_man_nwinvasion-rfmix

cd ~/sch_man_nwinvasion

mkdir results/rfmix

cd results/rfmix

#get a vcf of references and one of query haplotypes

#first remove unwanted samples (margrebowie and kenya):
echo -e "ERX284221\nERR119614" >remove.list

vcftools \
    --vcf ../phasing/auto_beagle.vcf \
    --remove remove.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >auto_beagle_rfmix.vcf

#choose 5 random samples from each population
cat auto_beagle_rfmix.vcf | head -n 5000 | grep "#" | tail -n1 | cut -f10- | sed 's/\t/\n/g'  | shuf

#manually cleaned to show in ref_samples.key and samples in ref_samples.list
#Sm.BR_PdV.2556.1 brazil
#Sm.BR_PdV.1039.1 brazil
#Sm.BR_PdV.2368.1 brazil
#Sm.BR_PdV.2450.1 brazil
#Sm.BR_PdV.2039.1    brazil
#Sm.SN_Te3.1 senegal
#Sm.SN_Te26.1    senegal
#Sm.SN_Te55.1    senegal
#Sm.SN_Nd56.1    senegal
#Sm.SN_Te68.1    senegal
#Sm.TZ_134.2.2   tanzania
#Sm.TZ_134.1.1   tanzania
#Sm.TZ_077.7.3   tanzania
#Sm.TZ_055.7.1   tanzania
#Sm.TZ_141.6.1   tanzania
#Sm.NE_Di297.1   niger
#Sm.NE_Na381.1   niger
#Sm.NE_Na40.1    niger
#Sm.NE_Na376.2   niger
#Sm.NE_Di186.1   niger
#Sro_female_1.1_CCATCCTC rodhaini
#Sro_female_1.2_CCGACAAC rodhaini
#Sro_female_2.1_CCTAATCC rodhaini
#Sro_female_2.2_CCTCTATC rodhaini
#Sro_male_1.1_ATCATTCC   rodhaini
#Sro_male_1.2_ATTGGCTC   rodhaini
#Sro_male_2.1_CAAGGAGC   rodhaini
#Sro_male_2.2_CACCTTAC   rodhaini
#ERR310938   rodhaini
#ERR046038   puerto_rico
#ERR539847   guadeloupe2
#ERR539848   guadeloupe3
#ERR103050   cameroon
#ERR103049   senegal
#ERR119615   uganda
#ERR997461   uganda

vcftools \
    --vcf auto_beagle_rfmix.vcf \
    --remove ref_samples.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >auto_beagle_rfmix_query.vcf

vcftools \
    --vcf auto_beagle_rfmix.vcf \
    --keep ref_samples.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >auto_beagle_rfmix_ref.vcf


#all rodhaini

rfmix \
    --query-file=auto_beagle_rfmix_query.vcf \
    --reference-file=auto_beagle_rfmix_ref.vcf \
    --sample-map= \
    --genetic-map= \
    --output-basename= \
    --n-threads=12 \
    --random-seed=98765321 \
    --crf-spacing= \    
    --rf-minimum-snps= \
    
