mkdir results/astral
mkdir results/astral/logs/

bedtools merge \
    -d 0 \
    -i data/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed \
    >results/astral/merged_baits.bed


#get genomes form vcfs for each individual
#generate a genome with IUPAC
#generate a phased genome and extract individual loci from them (in tab)

#prep to generate indivudal genomes
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 12"
mkdir results/astral/individual_genomes
mkdir results/astral/locus_sequences

#process the vcf file
cp results/variant_filtration/smv7_ex_autosomes.vcf results/astral/

#get list of samples
grep "#" results/astral/smv7_ex_autosomes.vcf \
    | tail -n 1 \
    | cut -f10- \
    | sed 's/\t/\n/g' \
    >results/astral/samples.list

#process for bcftools
bgzip results/astral/smv7_ex_autosomes.vcf
tabix -f results/astral/smv7_ex_autosomes.vcf.gz

#for each individual generate a genome
while read -r INDIV; do
    #consensus
    CONS_CMD="$CONDA bcftools consensus \
        -M \"?\" \
        -I \
        --sample $INDIV \
        -f data/genomes/Smansoni_v7.fa \
        results/astral/smv7_ex_autosomes.vcf.gz \
        >results/astral/individual_genomes/"$INDIV".fasta"

    BED_CMD="$CONDA bedtools getfasta \
        -tab \
        -fi results/astral/individual_genomes/"$INDIV".fasta \
        -bed results/astral/merged_baits.bed \
        -fo results/astral/locus_sequences/"$INDIV".tab"

        echo $CONS_CMD | $QSUB -N "$INDIV"_phased -o results/astral/logs/"$INDIV"_phased.log
        echo $BED_CMD  | $QSUB -N "$INDIV"_bed    -o results/astral/logs/"$INDIV"_bed.log -hold_jid "$INDIV"_phased     
done <results/astral/samples.list


#create a fasta file for each locus
while read -r INDIV; do

    ##print to loci files
    awk -v header=$INDIV"_" -v outdir="results/astral/locus_sequences/" \
        '{print ">"header"#"$1"\n"$2 >>outdir$1".fas"}' \
        results/astral/locus_sequences/"$INDIV".tab

done <results/astral/samples.list

#summarize each alignment to find phylogenetically informative loci
mkdir results/astral/informative_loci


for CHR in 1 2 3 4 5 6 7; do

./bin/amas/amas/AMAS.py summary \
    -i results/astral/locus_sequences/SM_V7_$CHR*.fas \
    -f fasta \
    -d dna \
    -c 12 \
    -o results/astral/SM_V7_"$CHR"_loci_summary.tsv

done 

#create a clean loci_summary file (remove excess headers)
cat results/astral/SM_V7_1_loci_summary.tsv >results/astral/loci_summary.tsv
for CHR in 2 3 4 5 6 7; do
    sed 1d results/astral/SM_V7_"$CHR"_loci_summary.tsv \
        >>results/astral/loci_summary.tsv
done 

rm results/astral/SM_V7_?_loci_summary.tsv

#copy to new dir of informative loci
for PARS_INF_LOCUS in $(sed 1d results/astral/loci_summary.tsv | awk '{if ($9 > 1) print $1}'); do
    cp results/astral/loci_sequences/$PARS_INF_LOCUS \
        results/astral/informative_loci/ 
done 


#now do a raxml run for each locus (fast with 100k bootstrap replicates)

#collapse unsupported branches


