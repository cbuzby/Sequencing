# Sequencing
 Sequencing pipelines for Chromosome Substitution BSA analysis. Reads must be aligned to reference and then variants called through GATK ```haplotypecaller``` before analyzing variants using QTLSeqR. Both CB_Pipeline and Analysis folders are used here. ```NZ_CB_Pipeline``` is obsolete.

## Pipelines
1. NZ: Adjusted from Naomi Ziv's 2017 published pipeline using samtools
2. GATK (g): Adjusted from Mohammed Khalfan's GATK variant calling pipeline
3. **CB_Pipeline**: using a variety of gatk tools, but adjusted for specific sequences
4. **Analysis**: uses QTLSeqR for analysis

### CB_Pipeline

Trim, align, and sort files

_Use any of these based on which data you're using_
```
CB_1.0_trim.align.sort
CB_1.1_trim.align.sort_all
CB_1.2_trim.align.sort_all
```
Trim, align, and sort using BQSR adjustment 

_GOLD STANDARD; use both of these scripts instead of any of the previous_
```
CB_1.3_trim.align.sort_all_BQSR
CB_1.5_ApplyBQSR
```

Reformatting .bam files

_Merge bam files into one, then split by reference (so by chromosome) to array variant calling_
```
CB_2.0_merge.split
CB_2.1_.split
CB_3.0_Index
```

Variant Calling

_Either of the two scripts; can use excluded regions for telomeres (T)_
```
CB_4.0_CallVariants
CB_4.0_CallVariants_T
```

Zip, concatenate split outputs, and sort

_Input: VCF for each chromosome; Output: sorted VCF to input into gatk table_
```
CB_5.0_zip.concat.sort
```

Make GATK table for R
```

```

### Analysis

Set parameters and create function
```
mindepth = 50
maxdepth = 1500
minsampledepth = 50
mingq = 0.98

PipelineFunc <- function(HighBulk, LowBulk, rawData = "mergedCuSO4.REF.SortedCat.vcf.output.table",
                         mindepth = mindepth,
                         maxdepth = maxdepth,
                         minsampledepth = minsampledepth,
                         mingq = mingq,
                         windowSize = 2e4,
                         Chroms = c("NC_001134.8", "NC_001135.5", "NC_001136.10", 
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6", 
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5", 
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1")
                         )
  ```
  
**Within Function**

Import Data
  ```
  {
  mytitle <- paste(HighBulk, LowBulk, sep = " vs ")
  HNGLCDRXY <- read.table(rawData, header = TRUE)

  df <- importFromGATK(
        file = rawData,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
        )

df %>% merge(.,ChromKey) %>% 
  group_by(CHROM) %>% mutate(Start = min(POS) + 350, End = max(POS) - 350) %>% 
  as.data.frame() %>% na.omit() -> df

colnames(df)[1] <- "NC_Chrom"
colnames(df)[which(colnames(df) == "chromosomes")] <- "CHROM"

df$CHROM <- factor(df$CHROM, levels = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"))
```

Run analysis
```
#Filter SNPs based on some criteria
df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.10, #0.20
        minTotalDepth = mindepth, #100
        maxTotalDepth = maxdepth, #400
        minSampleDepth = minsampledepth, #40
        minGQ = mingq #99
    )

  #df %>% merge(.,ChromKey)
  #df_filt <- df
  
  #Run G' analysi
  df_filt <- runGprimeAnalysis(
      SNPset = df_filt,
      windowSize = windowSize, #1e6
      outlierFilter = "deltaSNP")
  
  #Run QTLseq analysis
  df_filt <- runQTLseqAnalysis(
      SNPset = df_filt,
      windowSize = windowSize,
      popStruc = "F2",
      bulkSize = 1000, #c(25, 25)
      replications = 10000,
      intervals = c(95, 99)
  )
```

Run stats

```
  df_filt$idu <- row.names(df_filt)
  q <- 0.01
  fdrT <- getFDRThreshold(df_filt$pvalue, alpha = q)
  GprimeT <- df_filt[which(df_filt$pvalue == fdrT), "Gprime"]
```

Return dataframe
```
    
  return(as.data.frame(df_filt))
}
```
