mkdir skylines

bedtools merge \
    -d 0 \
    -i data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed \
    >results/skyline/merged_baits.bed


for POP in brazil niger senegal tanzania; do

    #mkdir results/skyline/$POP

    ##get biallelic snps for each population
    #vcftools \
    #    --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
    #    --keep results/skyline/$POP/$POP.list \
    #    --mac 2 \
    #    --recode \
    #    --recode-INFO-all \
    #    --stdout \
    #    >results/skyline/$POP/smv7_ex_autosomes_$POP.vcf &

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
        -a results/skyline/merged_baits.bed \
        -b results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
        > results/skyline/$POP/"$POP"_snps_per_bait.bed
    
    ##count the number of genotyped and filter for too may and too few
    cat results/skyline/$POP/"$POP"_snps_per_bait.bed \
        | awk '{ if ($4>2 && $4<150) print $0 }' \
            > results/skyline/$POP/"$POP"_loci_gt2_lt150.bed

    #extract snps from target loci
    vcftools \
        --vcf results/skyline/$POP/smv7_ex_autosomes_informative_$POP.vcf \
        --bed results/skyline/$POP/"$POP"_loci_gt2_lt150.bed \
        --recode \
        --recode-INFO-all \
        --stdout \
        >results/skyline/$POP/"$POP"_loci_gt2_lt150.vcf &
    
done

CONDA="conda activate sch_man_nwinvasion-nbanalyses; "
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1"

for POP in brazil niger senegal tanzania; do

    mkdir results/skyline/$POP/loci_vcfs
    mkdir results/skyline/$POP/logs/

    BED_FILE="results/skyline/$POP/"$POP"_loci_gt2_lt150.bed"
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

        #dont overload the cluster        
        NUM_JOBS=$(qstat | wc -l)

        while [ $NUM_JOBS -gt 7000 ]; do
            sleep 10s
            NUM_JOBS=$(qstat | wc -l)
        done    


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

        #dont overload the cluster        
        NUM_JOBS=$(qstat | wc -l)

        while [ $NUM_JOBS -gt 7000 ]; do
            sleep 10s
            NUM_JOBS=$(qstat | wc -l)
        done    

    done 
done

#merge all population vcf files

#create pop specific lists of bam file
for POP in brazil niger senegal tanzania; do
    ls results/skyline/$POP/phased_vcfs/*.vcf \
        >results/skyline/$POP/phased_vcfs/phased_vcfs.list

    bin/gatk-4.1.2.0/gatk --java-options "-Xmx40g" \
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



#generate a phased genome and extract individual loci from them (in tab)
for POP in brazil niger tanzania senegal; do
    
    QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12"
    #bgzip results/skyline/$POP/phased.vcf
    #tabix -f results/skyline/$POP/phased.vcf.gz

    rm -r results/skyline/$POP/phased_genomes
    mkdir results/skyline/$POP/phased_genomes
    rm -r results/skyline/$POP/phased_loci_sequences
    mkdir results/skyline/$POP/phased_loci_sequences
    rm -r results/skyline/$POP/logs
    mkdir results/skyline/$POP/logs
   
    #select  random indivuals per population
    #shuf results/skyline/$POP/$POP.list | head -n 10 >results/skyline/$POP/"$POP"_random.list

    while read -r INDIV; do
        for HAP in 1 2; do
            #consensus
            CONS_CMD="$CONDA bcftools consensus \
                -H $HAP \
                -M \"?\" \
                --sample $INDIV \
                -f data/genomes/Smansoni_v7.fa \
                results/skyline/$POP/phased.vcf.gz \
                >results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta"

            BED_CMD="$CONDA bedtools getfasta \
                -tab \
                -fi results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta \
                -bed results/skyline/$POP/"$POP"_loci_gt2_lt150.bed \
                -fo results/skyline/$POP/phased_loci_sequences/$INDIV"_H"$HAP.phased.tab"

            echo $CONS_CMD | $QSUB -N $INDIV"_H"$HAP.phased -o results/skyline/$POP/logs/$INDIV"_H"$HAP.phased.log
            echo $BED_CMD  | $QSUB -N $INDIV"_H"$HAP.bed    -o results/skyline/$POP/logs/$INDIV"_H"$HAP.bed.log -hold_jid $INDIV"_H"$HAP.phased     
        done

    done <results/skyline/$POP/$POP.list
done


for POP in brazil niger tanzania senegal; do
    
    QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12"
    #bgzip results/skyline/$POP/phased.vcf
    #tabix -f results/skyline/$POP/phased.vcf.gz

    rm -r results/skyline/$POP/phased_genomes
    mkdir results/skyline/$POP/phased_genomes
    rm -r results/skyline/$POP/phased_loci_sequences
    mkdir results/skyline/$POP/phased_loci_sequences
    rm -r results/skyline/$POP/logs
    mkdir results/skyline/$POP/logs
   
    #select  random indivuals per population
    #shuf results/skyline/$POP/$POP.list | head -n 10 >results/skyline/$POP/"$POP"_random.list

    while read -r INDIV; do
        for HAP in 1 2; do
            #consensus
            CONS_CMD="$CONDA bcftools consensus \
                -H $HAP \
                -M \"?\" \
                --sample $INDIV \
                -f data/genomes/Smansoni_v7.fa \
                results/skyline/$POP/phased.vcf.gz \
                >results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta"

            BED_CMD="bedtools getfasta \
                -tab \
                -fi results/skyline/$POP/phased_genomes/$INDIV"_H"$HAP.phased.fasta \
                -bed results/skyline/$POP/"$POP"_loci_gt2_lt150.bed \
                -fo results/skyline/$POP/phased_loci_sequences/$INDIV"_H"$HAP.phased.tab"

            echo $CONS_CMD | $QSUB -N $INDIV"_H"$HAP.phased -o results/skyline/$POP/logs/$INDIV"_H"$HAP.phased.log
            echo $BED_CMD  | $QSUB -N $INDIV"_H"$HAP.bed    -o results/skyline/$POP/logs/$INDIV"_H"$HAP.bed.log -hold_jid $INDIV"_H"$HAP.phased     
        done

    done <results/skyline/$POP/$POP.list
done


#create a fasta file for each locus
for POP in brazil niger tanzania senegal; do
    for HAP in 1 2; do
        while read -r INDIV; do
            ##print to loci files
            awk -v header="$INDIV"_H"$HAP" -v outdir="results/skyline/$POP/phased_loci_sequences/" \
                '{print ">"header"#"$1"\n"$2 >>outdir$1".fas"}' \
                results/skyline/$POP/phased_loci_sequences/$INDIV"_H"$HAP.phased.tab
        done <results/skyline/$POP/$POP.list
    done
done


#find phylogenetically informative loci
for POP in brazil niger tanzania senegal; do

    #summarize each alignment
    ./bin/amas/amas/AMAS.py summary \
        -i results/skyline/$POP/phased_loci_sequences/*.fas \
        -f fasta \
        -d dna \
        -o results/skyline/$POP/phased_loci_summary.tsv

    mkdir results/skyline/$POP/informative_phased_loci
    rm results/skyline/$POP/informative_phased_loci/*


    for PARS_INF_LOCUS in $(sed 1d results/skyline/$POP/phased_loci_summary.tsv | awk '{if ($9 > 5 && $3 < 500) print $1}'); do
        cp results/skyline/$POP/phased_loci_sequences/$PARS_INF_LOCUS \
            results/skyline/$POP/informative_phased_loci/ 
    done 
    
done

#git clone https://github.com/marekborowiec/AMAS.git bin/

#now set up replicate runs per population
for POP in brazil niger tanzania senegal; do
    for REP in 1 2 3; do
        
        rm -r results/skyline/$POP/replicate_$REP
        mkdir results/skyline/$POP/replicate_$REP

        ##randomize list of loci 
        for FILE in $(ls results/skyline/$POP/informative_phased_loci/SM*fas | shuf | head -n 50); do
            echo $FILE >>results/skyline/$POP/replicate_$REP/locus.list

            FILE_NAME=$(basename $FILE)
            cp $FILE results/skyline/$POP/replicate_$REP/$FILE_NAME
        done
        
        #modify headers in each fasta file so they are for the indiviudal and conssitent across loci
        sed -i 's/\#.*//' $FILE results/skyline/$POP/replicate_$REP/*.fas

        #modify file names so that trees can be stored properly (in beast2)
        rename : - results/skyline/$POP/replicate_$REP/*.fas

        #make concatenated nexus alignment, partition, and finished files
        ./bin/amas/amas/AMAS.py concat  \
            -i results/skyline/$POP/replicate_$REP/*.fas \
            -f fasta \
            -d dna \
            --out-format nexus \
            --part-format nexus \
            --concat-part results/skyline/$POP/replicate_$REP/$POP"_rep"$REP"_partition.nexus" \
            --concat-out results/skyline/$POP/replicate_$REP/$POP"_rep"$REP"_alignment.nexus"

        cat results/skyline/$POP/replicate_$REP/$POP"_rep"$REP"_alignment.nexus" \
            results/skyline/$POP/replicate_$REP/$POP"_rep"$REP"_partition.nexus" \
            >results/skyline/$POP/replicate_$REP/$POP"_rep"$REP".nexus"
    done
done

#wget -P bin/ http://hudson.cs.auckland.ac.nz/job/BEAST2/lastSuccessfulBuild/artifact/build/dist/beast.jar

#make the xml files 


#run beast
#java -jar beast.jar -working -seed 12345 -threads 1 -overwrite results/skyline/brazil/replicate_1/brazil_rep1.xml


#java -cp ./bin/beast.jar beast.app.tools.EBSPAnalyser \
#    -i results/skyline/niger/replicate_1/EBSP.log \
#    -burnin 20 \
#    -type linear \
#    -o results/skyline/niger/replicate_1/EBSP_linear_20p.tsv


