#Purpose: Estimate GAP1 copy numbers and CNV breakpoints from read depth of WGS data of isolated clones from the #EE_GAP1_ArchMuts_2021. These clones were isolated from generation 79 (g79).

#Upstream of this, I ran the mapcheck.sh bash script
#which outputted `{NAME}_RD.txt` file. Column headers are not there but are chromosome number, coordinate, and read depth.

#takes output from samtools depth
#`samtools depth -a ${bam} > ${ID}_RD.txt`

#Set up working directory and load packages
#setwd("~/Documents/WGS_clones_g79")
setwd("C:/Users/cassa/OneDrive/Documents/GitHub/Sequencing/Analysis/RdataFiles/Cov_JC/")

library(tidyverse)
library(docstring) #like docstring in python for R

#Note, the ref sequence used for alignment and mapping was
# /scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna which contains the GFP KanMX GAP1 CNV reporter (JC manually searched for the sequences to verify)
# Coordinates in this R script use this ref sequence (i think)

#### functions ###

get_rd_data = function(x) {
  #' This function reads in a file, names columns, removes the mitochondrial chromosome, and changes chromosome names from RefSeq chromosomes to numeric
  #' @param x {name}_RD.txt file outputted from samtools depth -a ${bam} > {NAME}_RD.txt
  #' @return a tibble with three columns - chromosome number, coordinate, and read depth
  #' @example DGY2071_read_data = get_rd_data("DGY2017_RD.txt")
  rd <- read.table(x)
  names(rd) <- c("chromo", "coordinate", "depth")
  rd <- rd[rd$chromo !="NC_001224.1",] #mitochondrial
  #change chromosome names
  rd <- rd %>% mutate(chromosome = case_when(chromo == 'NC_001133.9' ~ 1,
                                             chromo == 'NC_001134.8' ~ 2,
                                             chromo == 'NC_001135.5' ~ 3,
                                             chromo == 'NC_001136.10' ~ 4,
                                             chromo == 'NC_001137.3' ~ 5,
                                             chromo == 'NC_001138.5' ~ 6,
                                             chromo == 'NC_001139.9' ~ 7,
                                             chromo == 'NC_001140.6' ~ 8,
                                             chromo == 'NC_001141.2' ~ 9,
                                             chromo == 'NC_001142.9' ~ 10,
                                             chromo == 'NC_001143.9' ~ 11,
                                             chromo == 'NC_001144.5' ~ 12,
                                             chromo == 'NC_001145.3' ~ 13,
                                             chromo == 'NC_001146.8' ~ 14,
                                             chromo == 'NC_001147.6' ~ 15,
                                             chromo == 'NC_001148.4' ~ 16)) %>%
    select(chromosome, coordinate, depth)
  return(rd)
}

#' The next set of functions calculate the relative copy number at GAP1 and other specified genes
#' Genes are specified by giving the coordinates from the reference sequence
#' @param x the tibble returned from get_rd_data()
#'
#' @return numeric; copy number at specified gene relative to mean whole genome coverage
#' @export
#'
#' @examples get_gap1_reporter_copies(DGY2071_read_data)
get_gap1_reporter_copies = function(x) {
  #coordinates are adjusted to take into account the PrACT1-GFP-KanMX
  mean((x$depth[x$chromosome=="11" & x$coordinate>=513946 & x$coordinate<=520246])/mean(x$depth))
}
#this function gets relative GFP copy number
get_GFP_copies = function(x){
  #coordinates of GFP (717bp) located upstream of GAP1
  mean((x$depth[x$chromosome=="11" & x$coordinate>=514934 & x$coordinate<=515650])/mean(x$depth))

}

#this function gets the relative copy number at GAP1
get_gap1_copies = function(x) {
  #coordinates are GAP1 CDS location our reference sequence GCF_000146045.2_R64_genomic_GAP1.fna
  mean((x$depth[x$chromosome=="11" & x$coordinate>=518438 & x$coordinate<=520246])/mean(x$depth))
}

#this function gets the relative copy number at rDNA
get_rDNA_copies = function(x) {
  mean((x$depth[x$chromosome=="12" & x$coordinate>=451000 & x$coordinate<=470000])/mean(x$depth))
}

#this function gets the relative copy number at DUR3
get_dur3_copies = function(x) {
  mean((x$depth[x$chromosome=="8" & x$coordinate>=72037 & x$coordinate<=74244])/mean(x$depth))
}

#this function gets the relative copy number at HXT 6/7
get_hxt_copies = function(x) {
  mean((x$depth[x$chromosome=="4" & x$coordinate>=1154216 & x$coordinate<=1161320])/mean(x$depth))
}


# A function that generates a table with all read depth and copy numbers for genes of interest.
# Genes of interest can be added in the mutate() argument.

make_summary_table = function(x, name) {
  #' this function generates a table with all read depth and copy numbers for genes of interest
  #' genes of interest can be added in the mutate() argument
  #' @param x should be tibble outputted from get_rd_data()
  #' @param name should be a string stating the sample name
  #'
  #' @return Tibble with read depth and relative copy numbers of each gene of interest for the sample
  #'
  #' @examples DGY2071_summary_table = make_summary_table(DGY2071_read_data, "DGY2071"))
  summary = tibble(get_gap1_copies(x))
  names(summary) = 'gap1_copies'
  summary = summary %>%
    mutate(sample = name,
           GFP_copies = get_GFP_copies(x),
           gap1_reporter_copies = get_gap1_reporter_copies(x),
           rDNA_copies = get_rDNA_copies(x),
           dur3_copies = get_dur3_copies(x),
           hxt_copies = get_hxt_copies(x),
           mean_rd = mean(x$depth), median_rd = median(x$depth), sd_rd = sd(x$depth),
           min_rd = min(x$depth), max_rd = max(x$depth)) %>%
    relocate(sample)
}


#This function makes a plots of read depth normalized to chromosome 11.
#This lets us see GAP1 CNVs.

make_rd_plot_chr11 = function(df, output, dir) {
  #' this function makes plots read depth normalized to chromosome 11 read depth
  #'
  #' @param df is a dataframe of read data,
  #' @param output name of outputted plot file
  #' @param dir is the directory in which to save it
  #'
  #' @return save a .tiff plot read depth normalized to chromosome 11
  #'
  #' @examples
  #check to see if the directory exists and if it doesn't, make it
  if(dir.exists(dir) == FALSE) {
    dir.create(dir)
  }
  #get chromosome 11 only
  chr11 <- df[which(df$chromosome == 11),]
  avg11 <- mean(chr11$depth)
  linez <- runmed(chr11$depth/avg11,25001)
  newdata <- cbind(chr11, linez)
  newdata <- data.frame(newdata)
  #downsample for ease of plotting
  chr11plot <- newdata[sample(nrow(newdata), 20000),]
  #plot
  plot <- ggplot(chr11plot, aes(x=chr11plot[,2],y=chr11plot[,3]/avg11)) +
    geom_point(size = .25) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), 'Relative Depth (compared to mean chr11)', limits=c(-0.5, 6)) +
    scale_x_continuous(expand = c(0, 0), limits=c(5000,655000), breaks = scales::pretty_breaks(n = 5), "Position on Chromosome 11 (kb)") +
    geom_line(aes(x= chr11plot[,2], y=chr11plot[,4]), color = "red")
  ggsave(filename = paste0(output, "_rd_chr11_plot.tiff"), plot = plot, path = dir)
}


#This function makes read depth plots for chromosome 11 normalized to the entire genome.
#I think this lets us detect aneuploidies (gains/losses) of chromosome 11.

make_rd_plot_genome = function(df, output, dir) {

  #this function makes read depth plots for chromosome 11 normalized to the entire genome
  #' Title
  #'
  #' @param df
  #' @param output
  #' @param dir
  #'
  #' @return
  #' @export
  #'
  #' @examples
  #' #df is a dataframe of rd, dir is the directory in which to save it
  #check to see if the directory exists and if it doesn't, make it
  if(dir.exists(dir) == FALSE) {
    dir.create(dir)
  }
  #get mean rd for the genome
  avg <- mean(df$depth)
  #get chromosome 11 only
  chr11 <- df[which(df$chromosome == 11),]
  linez <- runmed(chr11$depth/avg,25001)
  newdata <- cbind(chr11, linez)
  newdata <- data.frame(newdata)
  #downsample for ease of plotting
  chr11plot <- newdata[sample(nrow(newdata), 20000),]
  #plot
  plot <- ggplot(chr11plot, aes(x=chr11plot[,2],y=chr11plot[,3]/avg)) +
    geom_point(size = .25) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), 'Relative Depth (compared to mean genome)', limits=c(-0.5, 6)) +
    scale_x_continuous(expand = c(0, 0), limits=c(5000,655000), breaks = scales::pretty_breaks(n = 5), "Position on Chromosome 11 (kb)") +
    geom_line(aes(x= chr11plot[,2], y=chr11plot[,4]), color = "red")
  ggsave(filename = paste0(output, "_rd_genome_plot.tiff"), plot = plot, path = dir)
}


#This function makes read depth plots of Chromosome 4 normalized to the entire genome
#Let's us detect CNVs of the hxt6/7 gene.
make_rd_plot_genome_hxt = function(df, output, dir) {
  #df is a dataframe of rd, dir is the directory in which to save it
  #check to see if the directory exists and if it doesn't, make it
if(dir.exists(dir) == FALSE) {
  dir.create(dir)
}
#get mean rd for the genome
avg <- mean(df$depth)
#get chromosome 4 only
chr4 <- df[which(df$chromosome == 4),]
linez <- runmed(chr4$depth/avg,25001)
newdata <- cbind(chr4, linez)
newdata <- data.frame(newdata)
#downsample for ease of plotting
chr4plot <- newdata#[sample(nrow(newdata), 20000),]
#plot
plot <- ggplot(chr4plot, aes(x=chr4plot[,2],y=chr4plot[,3]/avg)) +
  geom_point(size = .25) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), 'Relative Depth (compared to mean genome)', limits=c(-0.5, 6)) +
  scale_x_continuous(expand = c(0, 0), limits=c(1155000,1164000), breaks = scales::pretty_breaks(n = 5), "Position on Chromosome 4 (kb)") +
  geom_line(aes(x= chr4plot[,2], y=chr4plot[,4]), color = "red")
ggsave(filename = paste0(output, "_hxt67_rd_genome_plot.tiff"), plot = plot, path = dir)
}


##### run the functions #####

#get all files
files = list.files(path = './', pattern = '*RD.txt')

#perform get_rd_data() in all files
clones_rd_data=map(paste0('./', files), get_rd_data)

#make one big summary table
summary_data_list=map2(clones_rd_data, files, ~make_summary_table(.x, str_sub(.y, 1, -8)))
summary_data_table = bind_rows(summary_data_list)

#write file out, change name of outfile in needed
write_csv(summary_data_table, file = 'cnv_g79clones_rd_summary.csv') #Write a name for outfile.

#make read depth plots relative on chromosome 11 and whole genome
walk2(clones_rd_data, files, ~make_rd_plot_chr11(.x, str_sub(.y, 1, -8), 'chr11_rd_plots'))
walk2(clones_rd_data, files, ~make_rd_plot_genome(.x, str_sub(.y, 1, -8), 'genome_rd_plots'))

####################################
#CB Testing
####################################

read.table(files[1])
