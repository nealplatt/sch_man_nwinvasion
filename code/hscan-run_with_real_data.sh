conda activate sch_man_nwinvasion-hscan

CONDA="conda activate sch_man_nwinvasion-hscan"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1 "    


cd ~/sch_man_nwinvasion

mkdir -p results/hscan/logs

#REAL DATA
for POP in brazil tanzania senegal niger; do
    mkdir results/hscan/$POP
        for i in $(seq 1 7 ); do
            CHROM="SM_V7_"$i
    
            #get chrom vcf
            VCF_CMD="vcftools \
                --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
                --keep results/lists/$POP.list \
                --chr $CHROM \
                --recode \
                --recode-INFO-all \
                --stdout \
                >results/hscan/$POP/"$CHROM"_"$POP".vcf"

            #convert vcf to hscan
            CONV_CMD="python code/vcf2hscan.py \
                results/hscan/$POP/"$CHROM"_"$POP".vcf \
                results/hscan/$POP/"$CHROM"_"$POP".hscan-in"

            #convert to hscan
            HSCAN_CMD="bin/H-scan \
                -i results/hscan/$POP/"$CHROM"_"$POP".hscan-in \
                -g 10000 \
                >results/hscan/$POP/"$CHROM"_"$POP".hscan-out"

        echo "$CONDA; $VCF_CMD; $CONV_CMD; $HSCAN_CMD" | $QSUB -N $POP"_"$CHROM -o results/hscan/logs/$POP"_"$CHROM.log

    done
done



