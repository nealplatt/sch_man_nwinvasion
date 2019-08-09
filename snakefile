#using snakemake v5.5.4
#using conda v4.7.11

#set main project dir and work from there
LOGS='logs'
RESULTS='results'
DATA='data'
BIN='bin'
GENOME_FILE = DATA + "/genomes/Smansoni_v7.fa"

configfile: "config/config.yml"

localrules: 
    all, 
    create_list_of_contigs,
    prep_gdbimport,

rule all:
    input:
        RESULTS + "/genotype/cohort_raw.vcf",
        #RESULTS + "/variant_filtration/recal_snps_culled_indv.vcf",
        #RESULTS + "/variant_filtration/schMan_v7_exome_snps_filtered.vcf",
        #RESULTS + "/phasing/schMan_v7_exome_phased.vcf"

rule bwa_index_genome:
    input:
        GENOME = GENOME_FILE
    output:
        expand( GENOME_FILE + ".{ext}", ext=["pac", "ann", "amb", "bwt", "sa" ])
    log:
        LOGS + "/bwa_index_genome.log"
    conda:
        "config/env.yml"
    shell:
        """
        bwa index {input.GENOME}
        """

rule seqdict_genome:
    input: 
        GENOME_FILE
    output:
        DATA + "/genomes/Smansoni_v7.dict"
    threads:
        1
    log:
        LOGS + "/seqdict_genome"
    conda:
        "config/env.yml"
    shell:
        """
        bin/gatk-4.1.2.0/gatk CreateSequenceDictionary -R {input}
        """

rule faidx_genome:
    input:
        rules.seqdict_genome.output,
        GENOME = GENOME_FILE
    output:
        GENOME_FILE + ".fai"
    log:
        LOGS + "/faidx_genome"
    conda:
        "config/env.yml"
    shell:
        """
        samtools faidx {input.GENOME}
        """

rule filter_exome_reads:
    input:
        R1 = DATA + "/exomes/{id}_R1.fastq.gz",
        R2 = DATA + "/exomes/{id}_R2.fastq.gz"
    output:
        R1_PE = temp( RESULTS + "/filtered_reads/{id}_filtered_R1.fastq.gz"    ),
        R1_SE = temp( RESULTS + "/filtered_reads/{id}_filtered_SE_R1.fastq.gz" ),
        R2_PE = temp( RESULTS + "/filtered_reads/{id}_filtered_R2.fastq.gz"    ),
        R2_SE = temp( RESULTS + "/filtered_reads/{id}_filtered_SE_R2.fastq.gz" ),
        RX    = temp( RESULTS + "/filtered_reads/{id}_filtered_RX.fastq.gz"    )
    threads:
        12
    log:
        LOGS + "/filter_exome_reads#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        trimmomatic \
            PE \
            -threads {threads} \
            -phred33 \
            {input.R1} \
            {input.R2} \
            {output.R1_PE} \
            {output.R1_SE} \
            {output.R2_PE} \
            {output.R2_SE} \
            LEADING:10 \
            TRAILING:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36

        zcat {output.R1_SE} {output.R2_SE} | gzip >{output.RX} 
        """

rule filter_sra_reads:
    input:
        R1 = DATA + "/sra/{id}_1.fastq.gz",
        R2 = DATA + "/sra/{id}_2.fastq.gz"
    output:
        R1_PE = temp( RESULTS + "/filtered_reads/{id}_filtered_R1.fastq.gz"    ),
        R1_SE = temp( RESULTS + "/filtered_reads/{id}_filtered_SE_R1.fastq.gz" ),
        R2_PE = temp( RESULTS + "/filtered_reads/{id}_filtered_R2.fastq.gz"    ),
        R2_SE = temp( RESULTS + "/filtered_reads/{id}_filtered_SE_R2.fastq.gz" ),
        RX    = temp( RESULTS + "/filtered_reads/{id}_filtered_RX.fastq.gz"    )
    threads:
        12
    log:
        LOGS + "/filter_sra_reads#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        trimmomatic \
            PE \
            -threads {threads} \
            -phred33 \
            {input.R1} \
            {input.R2} \
            {output.R1_PE} \
            {output.R1_SE} \
            {output.R2_PE} \
            {output.R2_SE} \
            LEADING:10 \
            TRAILING:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36

        zcat {output.R1_SE} {output.R2_SE} | gzip >{output.RX} 
        """

rule bwa_map:
    input:
        rules.bwa_index_genome.output,
        FASTQ_FILE = RESULTS + "/filtered_reads/{id}_filtered_{read}.fastq.gz",
        REFERENCE  = GENOME_FILE
    output:
        temp( RESULTS + "/mapped_reads/{id}_{read}.sai" )
    threads:
        12
    log:
        LOGS + "/bwa_map_{read}#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        bwa aln \
            -n 0.15 \
            -t {threads} \
            -f {output} \
            {input.REFERENCE} \
            {input.FASTQ_FILE}
        """

rule bwa_sampe_samse:
    input:
        REFERENCE = GENOME_FILE,
        SAI_R1 = RESULTS + "/mapped_reads/{id}_R1.sai",
        SAI_R2 = RESULTS + "/mapped_reads/{id}_R2.sai",
        SAI_RX = RESULTS + "/mapped_reads/{id}_RX.sai",
        PE_R1  = RESULTS + "/filtered_reads/{id}_filtered_R1.fastq.gz",
        PE_R2  = RESULTS + "/filtered_reads/{id}_filtered_R2.fastq.gz",
        SE_RX  = RESULTS + "/filtered_reads/{id}_filtered_RX.fastq.gz"
    output:
        PE_BAM = RESULTS + "/mapped_reads/{id}_PE.bam",
        SE_BAM = RESULTS + "/mapped_reads/{id}_SE.bam",
        FILTERED_PE_BAM = temp( RESULTS + "/mapped_reads/{id}_filtered_PE.bam" ),
        FILTERED_SE_BAM = temp( RESULTS + "/mapped_reads/{id}_filtered_SE.bam" )
    threads:
        12
    log:
        LOGS + "/bwa_sampe_samse#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        bwa sampe \
            {input.REFERENCE} \
            {input.SAI_R1} \
            {input.SAI_R2} \
            {input.PE_R1} \
            {input.PE_R2} \
            |samtools view \
                -Sb \
                - \
                >{output.PE_BAM}

        samtools view \
                -Sb \
                -F 4 \
                {output.PE_BAM} \
                >{output.FILTERED_PE_BAM}

        bwa samse \
            {input.REFERENCE} \
            {input.SAI_RX} \
            {input.SE_RX} \
            | samtools view \
                -Sb \
                -F 4 \
                - \
                >{output.SE_BAM}

         samtools view \
            -Sb \
            -F 4 \
            {output.SE_BAM} \
            >{output.FILTERED_SE_BAM}
        """

rule sort_merge_bam:
    input:
        FILTERED_PE_BAM = RESULTS + "/mapped_reads/{id}_filtered_PE.bam",
        FILTERED_SE_BAM = RESULTS + "/mapped_reads/{id}_filtered_SE.bam"
    output:
        PE_BAM_SORTED = temp( RESULTS + "/mapped_reads/{id}_sorted_PE.bam" ),
        SE_BAM_SORTED = temp( RESULTS + "/mapped_reads/{id}_sorted_SE.bam" ),
        MERGED_BAM    = temp( RESULTS + "/mapped_reads/{id}_merged.bam"    )
    threads:
        12
    log:
        LOGS + "/sort_merge_bam#{id}"
    conda:
        "config/env.yml"
    shell:
      """
        samtools sort -o {output.SE_BAM_SORTED} {input.FILTERED_SE_BAM}
        samtools sort -o {output.PE_BAM_SORTED} {input.FILTERED_PE_BAM}
            
        samtools merge \
            {output.MERGED_BAM} \
            {output.PE_BAM_SORTED} \
            {output.SE_BAM_SORTED}            
        """

rule add_readgroups_to_bam:
    input:
        MERGED_BAM = RESULTS + "/mapped_reads/{id}_merged.bam"
    output:
        RG_BAM = temp(RESULTS + "/mapped_reads/{id}_rg.bam"),
    threads:
        12
    log:
        LOGS + "/add_readgroups_to_bam#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        bin/gatk-4.1.2.0/gatk --java-options \"-Xmx4g -Xms4g\" AddOrReplaceReadGroups \
            --INPUT={input.MERGED_BAM} \
            --OUTPUT={output.RG_BAM} \
            --RGPU=unk \
            --RGLB=library1 \
            --RGPL=illumina \
            --RGSM={wildcards.id} \
            --RGID={wildcards.id}
        """

rule sort_bam_post_readgroup:
    input:
        RG_BAM = RESULTS + "/mapped_reads/{id}_rg.bam",
    output:
        RG_SORTED_BAM = temp( RESULTS + "/mapped_reads/{id}_rg_sorted.bam" )
    threads:
        12
    log:
        LOGS + "/sort_bam_post_readgroup#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        samtools sort -o {output.RG_SORTED_BAM} {input.RG_BAM}
        """

rule mark_duplicates_in_bam:
    input:
        RG_SORTED_BAM = RESULTS + "/mapped_reads/{id}_rg_sorted.bam",
    output:
        METRICS = temp( RESULTS + "/mapped_reads/{id}_dupmetrics.log" ),
        BAM     = RESULTS + "/mapped_reads/{id}_processed.bam"
    threads:
        12
    log:
        LOGS + "/mark_duplicates_in_bam#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        bin/gatk-4.1.2.0/gatk --java-options \"-Xmx4g -Xms4g\" MarkDuplicates \
            --INPUT {input.RG_SORTED_BAM} \
            --OUTPUT {output.BAM} \
            --METRICS_FILE {output.METRICS} \
            --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900 \
            --ASSUME_SORT_ORDER coordinate
        """

rule index_bam:
    input:
        BAM = RESULTS + "/mapped_reads/{id}_processed.bam"
    output:
        INDEX = RESULTS + "/mapped_reads/{id}_processed.bam.bai"
    threads:
        1
    log:
        LOGS + "/index_bam#{id}"
    conda:
        "config/env.yml"
    shell:
      """
        samtools index {input.BAM}
      """

#################################################################################
###call snps with gatk4
##################################################################################
# Create list of contigs for Sman genome. this will be used frequently for gatk
rule create_list_of_contigs:
    input: 
        GENOME_FILE + ".fai"
    output:
        temp( RESULTS + "/lists/contigs.list" ),
    threads:
        1
    log:
        LOGS + "/create_list_of_contigs.log"
    shell:
        """
        awk '{{print "results/genotype/"$0".vcf"}}' {input} >{output}
        """

#### initial SNPs with haplotype caller
rule haplotype_caller:
    input: 
        DATA + "/genomes/Smansoni_v7.dict",
        GENOME_FILE + ".fai",
        RESULTS + "/mapped_reads/{id}_processed.bam.bai",
        BAM            = RESULTS + "/mapped_reads/{id}_processed.bam",
        REFERENCE      = GENOME_FILE,
        TARGET_REGIONS = DATA + "/renamed-sma_agilent_baits.v7.0.chr_reorderd.bed",
    output:
        VCF     = temp( RESULTS + "/haplotype_caller/{id}.hc.vcf"     ),
        VCF_IDX = temp( RESULTS + "/haplotype_caller/{id}.hc.vcf.idx" )
    threads:
        12
    log:
        LOGS + "/haplotype_caller#{id}"
    conda:
        "config/env.yml"
    shell:
        """
        bin/gatk-4.1.2.0/gatk HaplotypeCaller \
            --input {input.BAM} \
            --output {output.VCF} \
            -reference {input.REFERENCE} \
            -L {input.TARGET_REGIONS} \
            --emit-ref-confidence GVCF
        """

####make a list of contigs and the gdb parent directory
rule prep_gdbimport:
    input:
        VCF_FILES = expand( RESULTS + "/haplotype_caller/{id}.hc.vcf",     id=config["SAMPLE_IDS"] ),
        VCF_IDXS  = expand( RESULTS + "/haplotype_caller/{id}.hc.vcf.idx", id=config["SAMPLE_IDS"] ) 
    output:
        LIST = temp( RESULTS + "/lists/hc.list" ),
        DB   = temp( RESULTS + "/gdbimport/tmp" )
    threads:
        1
    conda:
        "config/env.yml"
    log:
        LOGS + "/prep_gdbimport_r1.log"
    shell:
        """
        ls {input.VCF_FILES} >{output.LIST}
        touch {output.DB}
        """

###gdbimport
rule gdbimport:
    input:
        HC_VCF_LIST = RESULTS + "/lists/hc.list",
        VCF_FILES   = expand(RESULTS + "/haplotype_caller/{id}.hc.vcf", id=config["SAMPLE_IDS"]),
        VCF_IDXS    = expand(RESULTS + "/haplotype_caller/{id}.hc.vcf.idx", id=config["SAMPLE_IDS"]) 
    output:
        DB = temp( directory( RESULTS + "/gdbimport/{contigs}" ))
    threads:
        12
    conda:
        "config/env.yml"
    log:
        LOGS + "/gdimport#{contigs}"
    shell:
        """
        bin/gatk-4.1.2.0/gatk --java-options \"-Xmx4g -Xms4g\" GenomicsDBImport \
                -V {input.HC_VCF_LIST} \
                --genomicsdb-workspace-path {output.DB} \
                -L {wildcards.contigs} \
                --reader-threads {threads} \
                --batch-size 12
        """

###genotype each contg
rule genotype:
    input:
        REFERENCE = GENOME_FILE,
        DB = RESULTS + "/gdbimport/{contigs}"
    output:
        VCF = temp( RESULTS + "/genotype/{contigs}.vcf"     ),
        IDX = temp( RESULTS + "/genotype/{contigs}.vcf.idx" )
    threads:
        4
    conda:
        "config/env.yml"
    log:
        LOGS + "/genotype#{contigs}"
    shell:
        """
        bin/gatk-4.1.2.0/gatk GenotypeGVCFs \
                -R {input.REFERENCE} \
                -V gendb://{input.DB} \
                -new-qual \
                -O {output.VCF}
        """

rule merge_genotyped_vcfs_pre_recal:
    input:
        LIST = RESULTS + "/lists/contigs.list",
        VCFS = expand( RESULTS + "/genotype/{contig}.vcf",     contig=config["CONTIGS"]),
        IDXS = expand( RESULTS + "/genotype/{contig}.vcf.idx", contig=config["CONTIGS"]),
    output:
        MERGED_VCF = RESULTS + "/genotype/cohort_raw.vcf",
        MERGED_IDX = RESULTS + "/genotype/cohort_raw.vcf.idx",
    threads:
        12
    conda:
        "config/env.yml"
    log:
        LOGS + "/merge_genotyped_vcfs_pre_recal"
    shell:
        """
        bin/gatk-4.1.2.0/gatk --java-options "-Xmx4g" \
            MergeVcfs \
                --MAX_RECORDS_IN_RAM 500000 \
                -I {input.LIST} \
                -O {output.MERGED_VCF}
        """

rule variant_recalibration:
    input:
        CROSS_1_TRUTH_VCF   = DATA + "/sch_man_snp_training_panel/Cross_1.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_1_UNTRUTH_VCF = DATA + "/sch_man_snp_training_panel/Cross_1.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_2_TRUTH_VCF   = DATA + "/sch_man_snp_training_panel/Cross_2.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_2_UNTRUTH_VCF = DATA + "/sch_man_snp_training_panel/Cross_2.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_3_TRUTH_VCF   = DATA + "/sch_man_snp_training_panel/Cross_3.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_3_UNTRUTH_VCF = DATA + "/sch_man_snp_training_panel/Cross_3.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_4_TRUTH_VCF   = DATA + "/sch_man_snp_training_panel/Cross_4.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_4_UNTRUTH_VCF = DATA + "/sch_man_snp_training_panel/Cross_4.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_1_TRUTH_IDX   = DATA + "/sch_man_snp_training_panel/Cross_1.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_1_UNTRUTH_IDX = DATA + "/sch_man_snp_training_panel/Cross_1.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_2_TRUTH_IDX   = DATA + "/sch_man_snp_training_panel/Cross_2.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_2_UNTRUTH_IDX = DATA + "/sch_man_snp_training_panel/Cross_2.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_3_TRUTH_IDX   = DATA + "/sch_man_snp_training_panel/Cross_3.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_3_UNTRUTH_IDX = DATA + "/sch_man_snp_training_panel/Cross_3.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        CROSS_4_TRUTH_IDX   = DATA + "/sch_man_snp_training_panel/Cross_4.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_truth.vcf",
        CROSS_4_UNTRUTH_IDX = DATA + "/sch_man_snp_training_panel/Cross_4.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls-dp20.MI_var_untruth.vcf",
        REFERENCE           = GENOME_FILE,
        TARGETED_VCF        = RESULTS + "/genotype/cohort_target_regions.vcf",
        TARGETED_VCF_IDX    = RESULTS + "/genotype/cohort_target_regions.vcf.idx"
    output:
        RECAL_VCF = temp( RESULTS + "/variant_filtration/snp_recal.vcf" ),
        TRANCHES  = RESULTS + "/variant_filtration/snp_recal_tranches.csv",
        RSCRIPT   = RESULTS + "/variant_filtration/snp_recal_plots.R"
    threads:
        12
    conda:
        "config/env.yml"
    log:
        LOGS + "/final_merge_of_genotyped_vcfs"
    shell:
        """
        bin/gatk-4.1.2.0/gatk VariantRecalibrator \
            -R {input.REFERENCE} \
            -V {input.TARGETED_VCF} \
            --resource:cross1_truth,known=false,training=true,truth=true,prior=10.0 {input.CROSS_1_TRUTH_VCF} \
            --resource:cross1_untruth,known=false,training=true,truth=false,prior=10.0 {input.CROSS_1_UNTRUTH_VCF}\
            --resource:cross2_truth,known=false,training=true,truth=true,prior=10.0 {input.CROSS_2_TRUTH_VCF}\
            --resource:cross2_untruth,known=false,training=true,truth=false,prior=10.0 {input.CROSS_2_UNTRUTH_VCF}\
            --resource:cross3_truth,known=false,training=true,truth=true,prior=10.0 {input.CROSS_3_TRUTH_VCF}\
            --resource:cross3_untruth,known=false,training=true,truth=false,prior=10.0 {input.CROSS_3_UNTRUTH_VCF}\
            --resource:cross4_truth,known=false,training=true,truth=true,prior=10.0 {input.CROSS_4_TRUTH_VCF}\
            --resource:cross4_untruth,known=false,training=true,truth=false,prior=10.0 {input.CROSS_4_UNTRUTH_VCF}\
            -an SOR \
            -an MQ \
            -an MQRankSum \
            -an ReadPosRankSum \
            -mode SNP \
            -O {output.RECAL_VCF} \
            --tranches-file {output.TRANCHES} \
            --rscript-file {output.RSCRIPT} \
            --truth-sensitivity-tranche 100.0 \
            --truth-sensitivity-tranche 99.5 \
            --truth-sensitivity-tranche 99.0 \
            --truth-sensitivity-tranche 97.5 \
            --truth-sensitivity-tranche 95.0 \
            --truth-sensitivity-tranche 90.0
           """

rule apply_varient_recal_and_filter:
    input:
        RECAL_VCF      = RESULTS + "/variant_filtration/snp_recal.vcf",
        TRANCHES       = RESULTS + "/variant_filtration/snp_recal_tranches.csv",
        REFERENCE      = GENOME_FILE,
        TARGET_VCF     = RESULTS + "/genotype/cohort_target_regions.vcf",
        TARGET_VCF_IDX = RESULTS + "/genotype/cohort_target_regions.vcf.idx"
    output:
        SOFT_RECAL_VCF    = RESULTS + "/variant_filtration/snp_recal_soft.vcf",
        HARD_FILTERED_VCF = temp( RESULTS + "/variant_filtration/snp_recal_hard.vcf" ),
    params:
        TRANCHE_LEVEL = "97.5"    
    threads:
        12
    conda:
        "config/env.yml"
    log:
        LOGS + "/apply_varient_recal_and_filter"
    shell:
        """
        bin/gatk-4.1.2.0/gatk ApplyVQSR \
            -R {input.REFERENCE} \
            -V {input.TARGET_VCF} \
            -O {output.SOFT_RECAL_VCF} \
            --truth-sensitivity-filter-level {params.TRANCHE_LEVEL} \
            --tranches-file {input.TRANCHES} \
            --recal-file {input.RECAL_VCF} \
            -mode SNP

        bin/gatk-4.1.2.0/gatk SelectVariants \
            -V {output.SOFT_RECAL_VCF} \
            -select-type SNP \
            --exclude-filtered \
            -O {output.HARD_FILTERED_VCF} \
            -R {input.REFERENCE} \
        """



