# Sequencing
UPDATE: pipeline now uses Nextflow, READ.ME for which is in Nextflow_Pipeline/

## Nextflow Start from Scratch

1. Prepare reference: we used R64 as a reference as file ```GCF_000146045.2_R64_genomic.fna``` in folder Reference. Alignment requires that the reference have an index, which can be completed using ```samtools faidx Reference/*fna```. The BQSR step requires a dictionary, which can be completed using picard tools. Ensure that the files at the top of the Nextflow pipeline are all being referenced; just ```ll <copy path> | wc -l``` to check

```
java -jar /share/apps/picard/2.17.11/picard.jar CreateSequenceDictionary
R=GCF_000146045.2_R64_genomic.fna O=GCF_000146045.2_R64_genomic.dict
``` 

2. Prepare parent sequences; Oak is ```SRR5331805```, Wine is ```SRR5331804```
```
module load sra-tools/2.10.9
prefetch SRR5331804
fastq-dump --split-files --fasta 60 SRR5331804
```
3. Run the following on the parent sequences, except for the end where you just make a VCF of these locations to call.
   
4. Change the directory that you're creating by making a new .config file (```cp x.config y.config``` and then ```vi y.config```)
5. Make a new .q executable file for running this new config file; run this on the HPC using slurm:
   ```
   sbatch newexecutablefile.q
   ```
   and then monitor using ```watch sbatch --me```
   
6. Once the new folder has the trimmed, aligned, sorted, and bqsr folders, with the correct number of files in each, copy the files in Nextflow_Pipeline/ into that new folder and ```cd``` into it
7. Merge and split them all: ```sbatch CB_2.1_merge.split.q```
8. Index the files: ```for i in *bam; do sbatch CB_3.0_Index.q $i; done```
9. Call variants in a loop: ```for i in *bam; do sbatch CB_4.0_CallVariants_T.q $i; done```
10. Sort and merge all of the files, with last argument being the final name. This will also run the gatk vcftotable so that the final output can be loaded directly into R: ```sbatch CB_5.0_zip.concat.sort.q HVYTYDRX2``` 

## Pipelines
1. NZ: Adjusted from Naomi Ziv's 2017 published pipeline using samtools
2. GATK (g): Adjusted from Mohammed Khalfan's GATK variant calling pipeline
3. **CB_Pipeline**: using a variety of gatk tools, but adjusted for specific sequences
4. **Analysis**: uses QTLSeqR for analysis

Sequencing pipelines for Chromosome Substitution BSA analysis. Reads must be aligned to reference and then variants called through GATK ```haplotypecaller``` before analyzing variants using QTLSeqR. Both CB_Pipeline and Analysis folders are used here. ```NZ_CB_Pipeline``` is obsolete.

### CB_Pipeline

***

#### Trim, align, and sort files
_Use any of these based on which data you're using_
```
CB_1.0_trim.align.sort
CB_1.1_trim.align.sort_all
CB_1.2_trim.align.sort_all
```
#### Trim, align, and sort using BQSR adjustment 
_GOLD STANDARD; use both of these scripts instead of any of the previous_
```
CB_1.3_trim.align.sort_all_BQSR
CB_1.5_ApplyBQSR
```

#### Reformatting .bam files
_Merge bam files into one, then split by reference (so by chromosome) to array variant calling_
```
CB_2.0_merge.split
CB_2.1_.split
CB_3.0_Index
```

#### Variant Calling
_Either of the two scripts; can use excluded regions for telomeres (T)_
```
CB_4.0_CallVariants
CB_4.0_CallVariants_T
```

#### Zip, concatenate split outputs, and sort
_Input: VCF for each chromosome; Output: sorted VCF to input into gatk table_
```
CB_5.0_zip.concat.sort
```

#### Make GATK table for R
```

```

## Analysis
***

#### Use summary table for all runs to keep uniform
```
MQCRuns <- read.csv("C:\\Users\\cassa\\OneDrive\\Documents\\SiegalLab\\Sequencing_Tuboweb\\AllMultiQCRuns.csv")
MQCRuns %>% select(Pool, ShortName, VCF_Table) %>% distinct() -> RawFiles

for(i in 1:length(RawFiles$VCF_Table)){
  < convert gatk vcf tables to alleles >
  }

```
#### Convert gatk vcf tables to alleles
```
cybrInputGATKTable(RawFiles$VCF_Table[i]) %>% mutate(Coverage = as.numeric(AD.REF) + as.numeric(AD.ALT)) %>%
  select(POS, CHROM, Dataset, GQ, AD.REF, AD.ALT, Coverage) %>% mutate(Pool = RawFiles$Pool[i])-> rawdata
```
#### Calling specific SNPs based on parent sequences
Parent (Oak and Wine) sequences were aligned to the same reference genome, and SNPs called for each. Those which are alternate in Oak and not Wine are called "Oak" alleles, those in Wine but not Oak are "Wine" alleles, and those which are shared by both strains are not called as they will not be different in the final bulks. The reference is then used as the opposite strain, and the same script run so that the directions of each are consistent:
```
rawdata %>% #head(3000) %>%
  merge(parentSNPids) %>% mutate(REFW = as.numeric(Type == "Wine"), REFO = as.numeric(Type == "Oak")) %>%
  group_by(Dataset, CHROM, POS) %>%
  mutate(Wine = max(REFW * as.numeric(AD.REF), REFO * as.numeric(AD.ALT)),
         Oak = max(REFW * as.numeric(AD.ALT), REFO * as.numeric(AD.REF))) %>%
  select(Pool, Dataset, POS, CHROM, Coverage, Wine, Oak) %>%
  pivot_longer(c(Oak, Wine), names_to = "Allele", values_to = "Reads") -> rawdata_called
```
#### Smooth the data by finding rolling median, centered
```
rawdata_called %>% group_by(Pool, Dataset, CHROM, Allele) %>% arrange(POS) %>%
  reframe(POS = POS, SmoothCount = ceiling(frollapply(Reads, n = 200, FUN = median, align = "center"))) -> rawdata_smoothed
```

### GLM

#### glm_cb Function on its own
This will take in the data columns, weights (W), formula for glm, and output length and return the coefficients specified (return = "Z")
```
glm_cb2_short <- function(..., W, formula, numgroups = FALSE, outputlength = 4, return = c("Z")) {
  glm_fit <- glm(glm_formula, data = as.data.frame(data), weights = W, family = binomial)

  #Return specific ones like z-score only
  if(return %in% "Z"){
    output <- summary(glm_fit)$coefficients[((length(summary(glm_fit)$coefficients)*0.5)+1):((length(summary(glm_fit)$coefficients)*0.75))]
  }
  
  if(length(output) == outputlength){
    return(output)
  }else{
    return(rep(NA, outputlength))
  }
  
}
```
#### Test function for each experiment
Because each experiment has a different number of replicates, this script should be adjusted for each so that it has the correct output and labels. The formula can be changed for adding a replicate or not, and if using glmer, the function itself must be changed to reflect that.
```
testdata <- data %>% na.omit() %>% filter(CHROM == "I", POS == 34991) %>% arrange(Bulk, Parent)
formula = "Allele ~ Bulk * Parent"
glm(formula = formula, family = "binomial", data = testdata, weights = AvgCount)
```

#### Run parallelized across positions and chromosomes
```
data %>% na.omit() %>% 
  #REMOVE CHROMOSOME OF FIXATION (don't waste resources on that)
  filter(CHROM %in% c("M", "I") == FALSE) %>%
  group_by(CHROM, POS) %>%
  mutate_if(is.character, as.factor) %>%

  summarize(Summary = glm_cb2_short(Allele = Allele,
                             Bulk = Bulk,
                             Parent = Parent,
                             W = AvgCount,
                             formula = "Allele ~ Bulk * Parent",
                            numgroups = 1, outputlength = 4),
            Factor = (c("Intercept", "Bulk", "Parent", "Interaction"))) -> dataglm
```

#### Notes

* nextflow processes each FILE on its own, one file for each flask
* for 8 samples in a sp run, this is ~8 hours
* using CB_2.q, each file is merged together, indexed as one, and THEN SEPARATED by chromosome. This took 10 hours for a 1s Novaseq run
* CB_3.q takes a couple of minutes for the by-chromosome samples
* CB_4.q will take > 8 hours for all except the smallest samples if running 1s (I swear it was shorter for the sp run, maybe 7 max for the large ones)
* CB_4.q finishes for all Chr at: [>9hr]
