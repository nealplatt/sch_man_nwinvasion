mkdir skylines

cp ../../results/variant_filtration/smv7_ex_autosomes.vcf .

vcftools --vcf smv7_ex_autosomes.vcf --thin 100000 --recode --recode-INFO-all --stdout >smv7_ex_autosomes_100k-thnned.vcf


python ../../bin/vcf2phylip/vcf2phylip.py -i smv7_ex_autosomes_100k-thinned.vcf

sed 1d smv7_ex_autosomes_100k-thinned.min4.phy \
    | awk '{print ">"$1"\n"$2}' \
    >smv7_ex_autosomes_100k-thinned.fasta

#merge baits w/in 500 bp of one another
bedtools merge \
    -d 500 \
    -i data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed \
    >results/skyline/merged_baits_500bp.bed


for POP in brazil niger senegal tanzania; do

    ##get biallelic snps for each population
    vcftools \
        --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
        --keep results/skyline/$POP/$POP.list \
        --mac 2 \
        --recode \
        --recode-INFO-all \
        --stdout \
        >results/skyline/$POP/smv7_ex_autosomes_$POP.vcf

    #remove uninformative loci from each population
    vcftools \
        --vcf results/skyline/$POP/smv7_ex_autosomes_$POP.vcf \
        --singletons \
        --stdout \
        >results/skyline/$POP/smv7_ex_autosomes_singletons_$POP.tbl
        
    cut -f1,2 results/skyline/$POP/smv7_ex_autosomes_singletons_$POP.tbl \
        >results/skyline/$POP/smv7_ex_autosomes_singletons_$POP.list

    vcftools \
        --vcf results/skyline/$POP/smv7_ex_autosomes_$POP.vcf \
        --exclude-positions results/skyline/$POP/smv7_ex_autosomes_singletons_$POP.list \
        --recode \
        --recode-INFO-all \
        --stdout \
        >results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf
   

    #find snps in 500bp regions (from baits)
    bedtools intersect \
        -c \
        -a results/skyline/merged_baits_500bp.bed \
        -b results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
        > results/skyline/$POP/"$POP"_intersect_merged_baits_500bp_counts.bed
    
    ##count the number of genotyped and filter for too may and too few
    cat results/skyline/$POP/"$POP"_intersect_merged_baits_500bp_counts.bed \
        | awk '{ if ($4>3 && $4<50) print $0 }' \
            > results/skyline/$POP/"$POP"_loci_gt3_lt50.bed

    #extract snps from target loci
    vcftools \
        --vcf results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
        --bed results/skyline/$POP/"$POP"_loci_gt3_lt50.bed \
        --recode \
        --recode-INFO-all \
        --stdout \
        >results/skyline/$POP/"$POP"_loci_gt3_lt50.vcf &
    
done

CONDA="conda activate sch_man_nwinvasion-nbanalyses; "
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12 "

for POP in brazil niger senegal tanzania; do

    mkdir results/skyline/$POP/loci_vcfs
    mkdir results/skyline/$POP/logs/

    BED_FILE="results/skyline/$POP/"$POP"_loci_gt3_lt50.bed"
    while read -r BED_ENTRY
    do

        CHROM=$(echo $BED_ENTRY | cut -f1 -d" ")
        START=$(echo $BED_ENTRY | cut -f2 -d" ")
        STOP=$(echo $BED_ENTRY  | cut -f3 -d" ")

        VCF_CMD="vcftools \
            --vcf results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
            --chr $CHROM \
            --from-bp $START \
            --to-bp $STOP\
            --recode \
            --recode-INFO-all \
            --stdout \
            >results/skyline/$POP/loci_vcfs/$POP"_"$CHROM"_"$START"_"$STOP.vcf"

        echo $CONDA $VCF_CMD | $QSUB -N $POP"_"$CHROM"_"$START"_"$STOP -o results/skyline/$POP/logs/$POP"_"$CHROM"_"$START"_"$STOP.vcf.log

    done < $BED_FILE

done


for POP in brazil niger senegal tanzania; do

    mkdir results/skyline/$POP/loci_vcfs
    mkdir results/skyline/$POP/logs/

    BED_FILE="results/skyline/$POP/"$POP"_loci_gt3_lt50.bed"
    while read -r BED_ENTRY
    do

        CHROM=$(echo $BED_ENTRY | cut -f1 -d" ")
        START=$(echo $BED_ENTRY | cut -f2 -d" ")
        STOP=$(echo $BED_ENTRY  | cut -f3 -d" ")

        VCF_CMD="vcftools \
            --vcf results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
            --chr $CHROM \
            --from-bp $START \
            --to-bp $STOP\
            --recode \
            --recode-INFO-all \
            --stdout \
            >results/skyline/$POP/loci_vcfs/$POP"_"$CHROM"_"$START"_"$STOP.vcf"

        echo $CONDA $VCF_CMD | $QSUB -N $POP"_"$CHROM"_"$START"_"$STOP -o results/skyline/$POP/logs/$POP"_"$CHROM"_"$START"_"$STOP.vcf.log

    done < $BED_FILE

done


#physical phasing via whatshap
#phase brazil
BRZ_CMD="whatshap phase \
    -o results/skyline/$POP/"$POP"_loci_gt3_lt50_phased.vcf \
    results/skyline/$POP/"$POP"_loci_gt3_lt50.vcf \
    $(ls results/mapped_reads/Sm.BR_PdV.*processed.bam)"

#phase niger
NGR_CMD="whatshap phase \
    -o results/skyline/$POP/"$POP"_loci_gt3_lt50_phased.vcf \
    results/skyline/$POP/"$POP"_loci_gt3_lt50.vcf \
    $(ls results/mapped_reads/Sm.NE*processed.bam)"

#phase tanzania
TNZ_CMD="whatshap phase \
    -o results/skyline/$POP/"$POP"_loci_gt3_lt50_phased.vcf \
    results/skyline/$POP/"$POP"_loci_gt3_lt50.vcf \
    $(ls results/mapped_reads/Sm.NE*processed.bam)"

#phase senegal
SNG_CMD="whatshap phase \
    -o results/skyline/$POP/"$POP"_loci_gt3_lt50_phased.vcf \
    results/skyline/$POP/"$POP"_loci_gt3_lt50.vcf \
    $(ls results/mapped_reads/Sm.SN*processed.bam)"


#submit jobs to scheduler
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12 "

echo $BRZ_CMD | $QSUB -N BRZ.whatshap -o results/skyline/brazil/whatshap.log
echo $SNG_CMD | $QSUB -N SNG.whatshap -o results/skyline/senegal/whatshap.log
echo $NGR_CMD | $QSUB -N NGR.whatshap -o results/skyline/niger/whatshap.log
echo $TNZ_CMD | $QSUB -N TNZ.whatshap -o results/skyline/tanzania/whatshap.log












#bgzip phased.vcf
#tabix phased.vcf.gz
#bcftools consensus -H 1 -f reference.fasta phased.vcf.gz > haplotype1.fasta
#bcftools consensus -H 2 -f reference.fasta phased.vcf.gz > haplotype2.fasta

