conda activate sch_man_nwinvasion-msprime
cd ~/sch_man_nwinvasion

mkdir -p results/msprime/logs

CONDA="conda activate sch_man_nwinvasion-msprime;"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1 "    

for POP in tanzania senegal niger brazil; do
    mkdir results/msprime/$POP
    for I in $(seq 1 342); do
        MSPRIME_CMD="python code/msprime-qsub.py $POP $I"       
        echo "$CONDA $MSPRIME_CMD" | $QSUB -N $POP"_"$I -o results/msprime/logs/$POP"_"$I.log
    done
done

#now get the vcfs
python code/msprime-probe_snps_from_vcf.py

