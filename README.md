# Sequencing
 Sequencing pipelines for Chromosome Substitution BSA analysis. Reads must be aligned to reference and then variants called through GATK ```haplotypecaller``` before analyzing variants using QTLSeqR. Both CB_Pipeline and Analysis folders are used here. ```NZ_CB_Pipeline``` is obsolete.

## Pipelines
1. NZ: Adjusted from Naomi Ziv's 2017 published pipeline using samtools
2. GATK (g): Adjusted from Mohammed Khalfan's GATK variant calling pipeline
3. **CB_Pipeline**: using a variety of gatk tools, but adjusted for specific sequences
4. **Analysis**: uses QTLSeqR for analysis

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

### Analysis
***
#### Combine oak and wine parental alleles, and define the bulks, parents, and replicates
```
setwd("../../../GitHub/Sequencing/Analysis/")

parentSNPids <- cybrConvertParentalAlleles(ParentFiles = c("Wine_VCF.txt", "Oak_VCF.txt"), Truncate = TRUE)

parentSNPids %>% group_by(CHROM, POS) %>% summarize(Type = Type, Unique = length(unique(Type))) %>% filter(Unique == 1) -> pSNPs

#CHANGE THIS
mydatatotest = "Data/HKTFTDRX2.SortedCat.vcf.output.table"

FilteredData <- cybrInputGATKTable(mydatatotest) %>% 
  cybrQualityFilter() %>% 
  cybrIDAlleles(BSAdfstart = ., Parentdf = pSNPs, yeast = TRUE) %>% 
  na.omit()

#Using Gsub for this
gsub(FilteredData$Dataset, "HKTFTDRX2_n01_", "") #CHANGE THIS

FilteredData %>% mutate(DShort = gsub("HKTFTDRX2_n01_", "", Dataset),
                       DS = gsub(".fastq", "", DShort)) %>% select(-Dataset, -DShort) -> tempFilteredData

tempFilteredData$Bulk <- NA
tempFilteredData$Parent <- NA
tempFilteredData$Rep <- NA

tempFilteredData$Bulk[grep("C", tempFilteredData$DS)] <- "CuSO4" #CHANGE THIS
tempFilteredData$Bulk[grep("D", tempFilteredData$DS)] <- "Dilute"

tempFilteredData$Rep[grep("a", tempFilteredData$DS)] <- "A"
tempFilteredData$Rep[grep("b", tempFilteredData$DS)] <- "B"

tempFilteredData$Parent[grep("O", tempFilteredData$DS)] <- "Oak"
tempFilteredData$Parent[grep("W", tempFilteredData$DS)] <- "Wine"

tempFilteredData$ReadCount <- as.numeric(tempFilteredData$ReadCount)
  
# #THIS IGNORES REPLICATES
tempFilteredData %>% select(CHROM, POS, PAllele, ReadCount, Bulk, Parent, Rep) %>% distinct %>% 
 pivot_wider(names_from = c(Bulk, Parent, Rep, PAllele), values_from = ReadCount) -> cybr2Data
 ```
 
#### Check log alleles per flask and coverage across sequencing run
```
cybr2Data %>% pivot_longer(c(-CHROM, -POS), names_to = c("Bulk", "Parent", "Rep", "Allele"), names_sep = "_") %>%
  pivot_wider(names_from = Allele, values_from = value) %>% mutate(Coverage = Wine + Oak, logWineOak = log(Wine/Oak)) -> RawCountSummary
```
 
#### Smooth data by rolling mean or median

```
#Use rolling average of 100 SNPs, finding the mean
cybr2Data %>% cybr2_rollmean() -> rollmeanData

#Find the rolling median or change n instead
cybr2Data %>% pivot_longer(c(-CHROM, -POS), names_to = "label") %>% 
  group_by(CHROM, label) %>% arrange(POS) %>% 
  summarize(POS = POS, 
            CHROM = CHROM, 
            SmoothCount = ceiling(frollapply(value, n = 100, FUN = median))) %>% #CHANGE THIS LINE
  na.omit() %>% 
  pivot_wider(names_from = label,values_from = SmoothCount) -> rollData
```

#### Caluclate GLM of rolling data

```
#Change for different datasets
mydata <- rollData

#Run GLM - options are glmfixed_rep(), glmfixed_rep3(), and glmfixed()
#IF NOT USING REPS, change 1:5 to 1:4, and remove "Rep" from labels
mydata %>% group_by(CHROM, POS) %>% 
  summarise(summary = glmfixed_rep(HOOa = CuSO4_Oak_A_Oak, 
                           HOWa = CuSO4_Oak_A_Wine, 
                           HWOa = CuSO4_Wine_A_Oak,
                           HWWa = CuSO4_Wine_A_Wine,
                           LOOa = Dilute_Oak_A_Oak,
                           LOWa = Dilute_Oak_A_Wine,
                           LWOa = Dilute_Wine_A_Oak, 
                           LWWa = Dilute_Wine_A_Wine,
                           
                           HOOb = CuSO4_Oak_B_Oak, 
                           HOWb = CuSO4_Oak_B_Wine, 
                           HWOb = CuSO4_Wine_B_Oak,
                           HWWb = CuSO4_Wine_B_Wine,
                           LOOb = Dilute_Oak_B_Oak,
                           LOWb = Dilute_Oak_B_Wine,
                           LWOb = Dilute_Wine_B_Oak, 
                           LWWb = Dilute_Wine_B_Wine)[1:5],
                                                   label = c("intercept", "Bulk", "Parent", "Rep", "Interaction")) -> GLMdata
```

### Visualizing
***
#### Single GLM Plot for this Data for reference

```
GLMdata %>% 
  filter(label != "intercept", CHROM != "I", CHROM != "M") %>% 
  
  ggplot(aes(x = POS, y = summary, color = label)) + geom_line() + 
  
  facet_grid(~CHROM, scales = "free", space = "free") +
  theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) + ggtitle("GLM of Data")
```

#### Log Odds of Alleles Plot for reference

```
RawCountSummary %>% 
  filter(CHROM != "I", CHROM != "M") %>% 
  
  ggplot(aes(x = POS, y = logWineOak, shape = paste(Bulk, Parent, Rep, sep = "_"), color = Bulk)) + 
  geom_point(alpha = 0.3) + 
  
  facet_grid(~CHROM, scales = "free", space = "free") +
  theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
  scale_color_manual(values = c("Violet", "Black"))
```

#### Notes

* nextflow processes each FILE on its own, one file for each flask
* for 8 samples in a sp run, this is ~8 hours
* using CB_2.q, each file is merged together, indexed as one, and THEN SEPARATED by chromosome. This took 10 hours for a 1s Novaseq run
* CB_3.q takes a couple of minutes for the by-chromosome samples
* CB_4.q will take > 6 hours for even the smallest samples if running 1s (I swear it was shorter for the sp run, maybe 7 max for the large ones)
