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

        #echo $CONDA $VCF_CMD | $QSUB -N $POP"_"$CHROM"_"$START"_"$STOP -o results/skyline/$POP/logs/$POP"_"$CHROM"_"$START"_"$STOP.vcf.log

    done < $BED_FILE

done

#create pop specific lists of bam file
ls results/mapped_reads/Sm.BR*processed.bam >results/skyline/brazil/bams.list
ls results/mapped_reads/Sm.NE*processed.bam >results/skyline/niger/bams.list
ls results/mapped_reads/Sm.SN*processed.bam >results/skyline/senegal/bams.list
ls results/mapped_reads/Sm.TZ*processed.bam >results/skyline/tanzania/bams.list

QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1 "
#physically phase each locus
for POP in brazil niger senegal tanzania; do

    mkdir results/skyline/$POP/phased_vcfs
    BAM_LIST=$(cat results/skyline/$POP/bams.list)

    for VCF in $(ls results/skyline/$POP/loci_vcfs/*.vcf); do

        LOCUS=$(basename $VCF .vcf | cut -f2- -d"_")
        
        WHATSHAP_CMD="$CONDA whatshap phase \
            -o results/skyline/$POP/phased_vcfs/"$POP"_"$LOCUS"_phased.vcf \
            $VCF \
            $BAM_LIST"

        echo $WHATSHAP_CMD | $QSUB -N $POP"_"$LOCUS"_whatshap" -o results/skyline/$POP/logs/$POP"_"$LOCUS"_whatshap.vcf.log"

    done 
done

#merge all population vcf files

#create pop specific lists of bam file
for POP in brazil niger senegal tanzania; do
    ls results/skyline/$POP/phased_vcfs/*.vcf \
        >results/skyline/$POP/phased_vcfs/phased_vcfs.list

    bin/gatk-4.1.2.0/gatk --java-options "-Xmx4g" \
        MergeVcfs \
            --MAX_RECORDS_IN_RAM 500000 \
            -I results/skyline/$POP/phased_vcfs/phased_vcfs.list \
            -O results/skyline/$POP/phased.vcf
done


#for POP in tanzania; do
#    for CHR in 1 2 3 4 5 6 7; do
#        ls tanzania_SM_V7_$CHR*.vcf \
#            >$CHR.list

#        ../../../../bin/gatk-4.1.2.0/gatk --java-options "-Xmx4g" \
#            MergeVcfs \
#                --MAX_RECORDS_IN_RAM 500000 \
#                -I $CHR.list \
#                -O $CHR.vcf
#    done
#done

#ls ?.vcf >chrs.list
#        ../../../../bin/gatk-4.1.2.0/gatk --java-options "-Xmx4g" \
#            MergeVcfs \
#                --MAX_RECORDS_IN_RAM 500000 \
#                -I chrs.list \
#                -O phased.vcf


#bgzip and taibx pop vcf
#consensus for each individual
for POP in brazil niger tanzania senegal; do
    
    #bgzip results/skyline/$POP/phased.vcf
    #tabix results/skyline/$POP/phased.vcf.gz

    #mkdir results/skyline/$POP/phased_genomes
    rm -r results/skyline/$POP/phased_loci_sequences
    mkdir results/skyline/$POP/phased_loci_sequences

    while read -r INDIV; do

        for HAP in 1 2; do
            #consensus
            bcftools consensus \
                -H $HAP \
                -M "?" \
                --sample $INDIV \
                -f data/genomes/Smansoni_v7.fa \
                results/skyline/$POP/phased.vcf.gz \
                >results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta

            #getfasta (tab) from only selected loci per population
             bedtools getfasta \
                -tab \
                -fi results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta \
                -bed results/skyline/$POP/$POP"_loci_gt3_lt50.bed" \
                -fo results/skyline/$POP/phased_loci_sequences/$INDIV"_H"$HAP.phased.tab

            ##print to loci files
            awk -v header="$INDIV"_H"$HAP" -v outdir="results/skyline/$POP/phased_loci_sequences/" \
                '{print ">"header"#"$1"\n"$2 >>outdir$1".fas"}' \
                results/skyline/$POP/phased_loci_sequences/$INDIV"_H"$HAP.phased.tab
        done

    done < results/skyline/$POP/$POP.list

done


#now set up replicate runs per population
for POP in brazil niger tanzania senegal; do
    for REP in 1 2 3; do
        
        rm -r results/skyline/$POP/replicate_$REP
        mkdir results/skyline/$POP/replicate_$REP

        ##randomize list of loci 
        for FILE in $(ls results/skyline/$POP/phased_loci_sequences/SM*fas | shuf | head -n 50); do
            echo $FILE >>results/skyline/$POP/replicate_$REP/locus.list

            FILE_NAME=$(basename $FILE)
            cp $FILE results/skyline/$POP/replicate_$REP/$FILE_NAME
        
        done
    done
done

#prep beauti files for skyline runs

