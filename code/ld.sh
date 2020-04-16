conda activate sch_man_nwinvasion-ld

cd /master/nplatt/sch_man_nwinvasion

mkdir results/ld

#get pop specific VCF files
for POP in senegal niger tanzania brazil; do
    if [ ! -f results/ld/smv7_ex_"$POP"_ld_dist.tbl ]; then
        rm results/ld/smv7_ex_"$POP"_ld_dist.tbl
    fi
    
    for CHR in "1" "2" "3" "4" "5" "6" "7" ; do
    
        CHR="SM_V7_$CHR"
    
        #get pop specific chr vcf
        vcftools \
            --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
            --chr $CHR \
            --maf 0.05 \
            --keep results/lists/$POP.list \
            --recode \
            --recode-INFO-all \
            --stdout \
            >results/ld/smv7_ex_maf05_"$POP"_"$CHR".vcf
    
        plink \
            --threads 6 \
            --vcf results/ld/smv7_ex_maf05_"$POP"_"$CHR".vcf \
            --out results/ld/smv7_ex_maf05_"$POP"_"$CHR" \
            --double-id\
            --recode12 \
            --allow-extra-chr

        #calculate R2 between all snps on a chr
        plink \
            --threads 6 \
            --r2 \
            --file results/ld/smv7_ex_maf05_"$POP"_"$CHR" \
            --out results/ld/smv7_ex_maf05_"$POP"_"$CHR" \
            --double-id \
            --allow-extra-chr \
            --ld-window-r2 0.0 \
            --ld-window 1000000 \
            --ld-window-kb 1000
            
        #calc distance and generate a single table
        awk '{print $0"\t"$5-$2}' \
            results/ld/smv7_ex_maf05_"$POP"_"$CHR".ld \
            >>results/ld/smv7_ex_maf05_"$POP"_ld_dist.tbl
                
        #clean up
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".ld
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".map
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".ped
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".log
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".nosex
        rm results/ld/smv7_ex_maf05_"$POP"_"$CHR".vcf

    done
    
        #remove the header line from the table (duplicated from each chrom)
        sed -i -e '1p' -e '/CHR/d' results/ld/smv7_ex_maf05_"$POP"_ld_dist.tbl 
        sed -i -e '1s/0/BP_DISTANCE/' results/ld/smv7_ex_maf05_"$POP"_ld_dist.tbl
        
        #create a file with distances lt 500kb
        awk '{if ($8 <= 500000) print $0}' results/ld/smv7_ex_maf05_"$POP"_ld_dist.tbl \
            >results/ld/smv7_ex_maf05_"$POP"_ld_dist_lt500kb.tbl
done






%%R
#in R generate summary stats for bins and plot

#500bp bins for 500kb
breaks <- seq( 0, 5e5, 500)

centers   <- vector()
means     <- vector()
pops      <- vector()
smootheds <- vector()

for (pop in c("oman", "eafrica", "wafrica", "caribbean", "brazil_x")) {
    
    #read in lt 500kb r2 table from vcftools
    ld_table <- read.table(paste("results/ld/smv7_ex_maf05_", pop, "_ld_dist_lt500kb.tbl", sep=""), 
                           header=FALSE)
    
    #bin r2 values and calculate stats
    ld_binned <- stats.bin(ld_table$V8, ld_table$V7, breaks = breaks)
    
    #created regression line
    loessMod  <- loess(ld_binned$stats["mean",] ~ ld_binned$centers, span=0.50)
    smoothed  <- predict(loessMod)
    
    #append all data to vectors
    centers   <- append(centers, ld_binned$centers)
    means     <- append(means, ld_binned$stats["mean",])
    pops      <- append(pops, rep(pop, length(ld_binned$centers)))
    smootheds <- append(smootheds, smoothed)
}    

#build the dataframe and save to csv
r2_df <- data.frame(centers, means, pops, smootheds)
write.csv(r2_df, file = paste("results/ld/ld_dist_lt500kb_binned_smoothed.csv", sep="") ,row.names=FALSE)

#subset desired populations
major_groups <- subset(r2_df, pops == "oman"  | 
                              pops == "eafrica" | 
                              pops == "wafrica"   |
                              pops == "caribbean" |
                              pops == "brazil_x" ) 






%%R

pop_colors <- c("eafrica"   = "green",
                "oman"      = "yellow", 
                "wafrica"   = "red",
                "brazil_x"  = "purple",
                "caribbean" = "blue")

#start plotting
p <- ggplot(major_groups, aes(x     = centers, 
                              y     = means, 
                              color = pops))

#adjust colors
p <- p + scale_colour_manual(values = pop_colors)

#plot data
p <- p + geom_point(alpha = 0.0)

#smooth data with loess
p <- p + geom_smooth(span   = 0.75, 
                     method = "loess", 
                     lwd    = 1, 
                     se     = FALSE)

#modify x and y axis
p <- p + scale_x_continuous(name   = "Distance between SNPs (Kb)", 
                            breaks = seq(0,500000,100000), 
                            labels = c("0", "100", "200", "300", "400", "500"), 
                            expand = c(0,0),
                            )
p <- p + scale_y_continuous(name   = "Mean R2", 
                            expand = c(0,0), 
                            limits = c(0,1))

#change theme
p <- p + theme_bw()

#change fonts on axis elements and titles
p <- p + theme(axis.text  = element_text(size = 12),
               axis.title = element_text(size = 14,
                                         face = "bold"))

#removing gridlines
p <- p + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())

#set up plot title etc.
p <- p + ggtitle("LD Decay")
p <- p + theme(plot.title = element_text(hjust = 0.5, 
                                         vjust = 0.5, 
                                         face  = 'bold', 
                                          size  = 18))

#modify legend
p <- p + theme(legend.title         = element_text(size = 14, 
                                                   face = "bold"),
               legend.text           = element_text(size = 12),
               legend.position       = c(0.85, 0.85),
               legend.box.background = element_rect(colour = "black"))
p <- p + labs(col = "Population")
p <- p + scale_color_manual(labels = c("Brazil",
                                       "Caribbean",
                                       "E. Africa",
                                       "Oman",
                                       "W. Africa" ), 
                            values = c(pop_colors["brazil_x"], 
                                       pop_colors["caribbean"], 
                                       pop_colors["eafrica"],
                                       pop_colors["oman"],
                                       pop_colors["wafrica"] ))

#save the figure
svg_img <- "results/ld/ld_decay.svg"
png_img <- "results/ld/ld_decay.png"
ggsave(png_img, plot = p, dpi = 600)
ggsave(svg_img, plot = p)


#display in notebook
print(p)
