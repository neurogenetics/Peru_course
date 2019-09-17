# Lima, Peru Course: ISIA

**Date:** Tuesday 17 September 2019
**Creators:** Cornelis Blauwendraat, Sara Bandres-Ciga
**Acknowledgements:** Hail-driven Jupyter Notebooks, very kindly provided by the [Broad Institute]([https://www.broadinstitute.org/](https://www.broadinstitute.org/)) (Thank you very much!!!)

## Table of Contents
### [0. Getting Started](#0)
 1. Software that will be used
 2. Login information 
### [1. Cleaning Data for GWAS](#1)
 1. Brief overview of what the scripts do
 2. Going to the proper directory 
 3. Edit the STEP1 files to clean using PLINK and GCTA 
 4. Running the cleaning script (.sh)
 5. *(Optional)* Visualization in R
 6. Make principle components 
### [2. Running a GWAS](#2)
 1. Going to the proper directory 
 2. Filter imputed data
 3. Using RVtests to run a genome-wide association study 
 4. Checking the results
### [3. Generating a Genetic Risk Score (GRS)](#3)
 1. Going to the proper directory 
 2. Calculate GRS
 3. Using R to calculate p-values
 4. Visualization 
###  [4. Run a Burden Test](#4)
 1. Going to the proper directory 
 2. Using RVtests to run a burden test (One or More Genes)
 3. Checking the results 


<a id="0"></a>
## Part 0: Getting Started 
### Software that will be used: 
 - **plink v1.9** @ https://www.cog-genomics.org/plink/1.9/
 - **GCTA** @ https://cnsgenomics.com/software/gcta/#Overview
 - **RVTests** @ http://zhanxw.github.io/rvtests/
 - **R** @ https://www.r-project.org/

### Login Information:  
 - Login is https://notebook.hail.is 
	 - If there are problems, use https://notebook.hail.is/new
- Password:  Provided during the Course

<a id="1"></a>
## Part 1: Cleaning Data for GWAS

Preparing data for a genome-wide association study (GWAS) is a two step process
1. Clean the data 
2. Harmonize the data with reference panel 
	- No time for this today, but you can check the `STEP2_cleaning.sh` script if interested in how this is done

### Step 1: Brief Overview of What the Scripts Do
1. `sh STEP1_cleaning.sh INPUTFILE` does the following:
	- Call rate filter
	- Heterozygosity filter
	- Sex check
	- Variant filter
	- Variant filter between cases and controls
	- *(Optional)* Ancestry filter
	- *(Optional)* Relatedness filter


2. `sh STEP2_cleaning.sh INPUTFILE` (we don't have time for this one)
does the following:
	- Makes sure the right alleles are used
	- Makes sure the variants between input and reference match
	- Makes zipped `.vcf` files of cleaned input files per chromosome
	- ...and more small things


### Step 2: Going to the Proper Directory
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
ls
```
### Step 3: Edit the STEP1 files to Clean using PLINK and GCTA 
1. Click the `GWAS_course_files` directory and then the `QC_PACKAGE` directory
2. Click and open `STEP1_cleaning.sh` and edit the following: 
	- Change `./plink` to `plink`
	- Change `./gcta` to `gcta64`
	- At the top, `File` and then `Save`
3. Click and open `PCA_in_R.R` and edit all of the following: 
	- `pdf("FILE.pdf")` to `jpeg("FILE.jpeg")`
	- At the top, `File` and then `Save`

### Step 4:  Running the Cleaning Script (.sh)

```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
sh STEP1_cleaning.sh EXAMPLE_DATA
```

### Step 5: *(Optional)* Visualization in R
Check plots with the following script:
```python
from IPython.display import Image
Image(filename="GWAS_course_files/QC_PACKAGE/raw_hapmap_plot.jpeg")
```
Check files with **one** of the following scripts:
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
head PCA_filtered_europeans.txt
```
**OR**
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
head GENDER_FAILURES.txt
```
**OR**
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
head CALL_RATE_OUTLIERS.txt
```

### Step 6: Make Principle Components
Use PLINK to generate principle components on cleaned, pruned SNPs using the following script: 
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
plink --bfile FILTERED.EXAMPLE_DATA --maf 0.05 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt \
--make-bed --out pass1

# Filter the variants
plink --bfile pass1 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data

# Extract pruned SNPs and only these variants will be used for PC calculation
plink --bfile pass1 --extract pruned_data.prune.in --make-bed --out pass1_pruned 

# Calculate/Generate PCs based on pruned data set
plink --bfile pass1_pruned --pca --out PCA

# Then look at the .eigenvec file
```
Look at the Eigenvectors using the following script:
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/QC_PACKAGE/
head PCA.eigenvec
```

<a id="2"></a>
## Part 2: Running a GWAS
- The GWAS will be run using the imputed data 
	- We have no time to do the actual imputation during this demo
	- So we will use an example dataset
- ***Warning:*** Imputation imputes a lot of variants, many that might be not high quality

1. First, we create variant list R2>0.3 
2. Run a GWAS on a small piece of the genome

### Step 1: Going to the Proper Directory
Inspect content by doing the following: 
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/GWAS/
ls
```

### Step 2: Filter Imputed Data
- Normally, you would do this for all 22 chromosomes...
	- but for this demo, we are only focusing on chromosome 4

Here just one chromosome part:
```R
%%bash
# Go to the proper directory
cd GWAS_course_files/GWAS/
# Load R
R
library(plyr)
  input <- paste("chr4_short.info.gz")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.01 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  da <- dat[,c("SNP","ALT_Frq","Rsq")]
  write.table(da, paste("maf001rsq03minimums_chr4.info",sep = ""), row.names = F, quote = F, sep = "\t")
  
  input <- paste("chr4_short.info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.01 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
  da <- dat[,c("range")]
  write.table(da, paste("maf001rsq03minimums_chr4.txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
```

### Step 3: Using RVtests to Run a GWAS
- Again, normally we would loop for each chromosome... 
- Analysis takes about 3-4 minutes for 765 cases and 1982 controls and a couple thousand variants

Run a GWAS in RVTests after filtering using the following script:
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/GWAS/
# Run the GWAS 
rvtest --noweb --hide-covar --rangeFile maf001rsq03minimums_chr4.txt \
--out EXAMPLE_DATA_GWAS --single wald \
--inVcf EXAMPLE_DATA.vcf.gz --dosage DS --pheno covariates.txt \
--pheno-name PHENO_PLINK --covar covariates.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5
```
### Step 4: Check the Results 
Have a look at the results file:
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/GWAS/
sort -gk 9 EXAMPLE_DATA_GWAS.SingleWald.assoc | head
```
- The SNP variant `4:90641340` has a P=`4.20E-06` 
	- dbSNP: rs356220 -- One of the top variants for PD
	- This variant has a beta of -0.29, which is on OR of ~1.3
	
<a id="3"></a>
# Part 3: Generating a Genetic Risk Score (GRS)
In this section, you will calculate a GRS 

### Step 1: Going to the Proper Directory
Check overview of files that are in the folder:
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/GRS/
ls
```

### Step 2: Calculate GRS
Calculate a GRS using PLINK by running the following: 
```bash
%%bash
cd GWAS_course_files/GRS/
# Go to the proper directory
plink --bfile NEUROX_GRS_only --score META5_GRS_NEUROX.txt --out NEUROX_GRS
```

### Step 3: Using R to calculate p-values
```R
%%bash
# Go to the proper directory
cd GWAS_course_files/GRS/
R
# Load the necessary packages 
library(dplyr)
library(ggplot2)
# Read in the Data
GRS <- read.table("NEUROX_GRS.profile",header=T)
cov <- read.table("PHENOTYPES.txt",header=T)
cov$PHENO <- NULL
MM2 = merge(GRS,cov,by='IID')
data <- MM2
meanGRS <- mean(data$SCORE)
sdGRS <- sd(data$SCORE)
data$SCOREZ <- (data$SCORE - meanGRS)/sdGRS

# Association Test
data$PHENO <- data$PHENO-1
thisFormula1 <- formula(paste("PHENO ~ SCOREZ + sex  + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
model1 <- glm(thisFormula1, data = data, ,family=binomial)
print(summary(model1))
# Optional if p-value is really low
summary(model1)$coefficients[,4]

## Plotting
# Option 1 (simple)
data$PHENO <- data$PHENO+1
#pdf("MY_PLOT.pdf",width=4)
jpeg(file="MY_PLOT.jpeg")
boxplot(data$SCOREZ~data$PHENO,col=c('grey', 'red'),xlab="1 control, 2 PD-case",ylab="GRS Z-score",main="PD genetic risk score")
grid()
boxplot(data$SCOREZ~data$PHENO, add=TRUE)
dev.off()

# Option 2 (more fancy)
# pdf("MY_PLOT2.pdf",width=6)
jpeg(file="MY_PLOT2.jpeg")
p <- ggplot(data, aes(x=as.factor(PHENO), y=SCOREZ, fill=as.factor(PHENO))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("grey", "red")) + theme_bw() + 
labs(title="PD genetic risk score",x="1 control, 2 PD case", y = "GRS Z-score") + theme(legend.position="none")
dev.off()
```

### Step 4: Visualization 
View the generated plots by running the following:

```python
# View the simple plot
from IPython.display import Image
Image(filename="GWAS_course_files/GRS/MY_PLOT.jpeg")
```

```python
# View the more complex plot
from IPython.display import Image
Image(filename="GWAS_course_files/GRS/MY_PLOT2.jpeg")
```

<a id="4"></a>
## Part 4: Run a Burden Test


In preparation for this demo, we already created some files to make it easier, and this is how it was created: 
```bash
## Subset only the variants of interest in this case 3 GBA variants that have been associated with Parkinson's disease

# plink --bfile NEUROX --extract GBA_BURDEN_variants.txt --make-bed --out NEUROX_GBA

# GBA variants:
	# 1:155206167 => E326K 
	# 1:155206037 => T369M
	# 1:155205634 => N370S

# cat GBA_BURDEN_variants.txt
	# NeuroX_rs2230288
	# exm106220
	# exm106217

## Make from the plink files a .vcf file
# Extract the variants you are interested in (usually coding variants)
# plink --bfile NEUROX_GBA --extract GBA_BURDEN_variants.txt --recode vcf-iid --out BURDEN_INPUT

# Then zip and index the .vcf file (need samtools or htslib for this
# bgzip -c BURDEN_INPUT.vcf > BURDEN_INPUT.vcf.gz
# tabix -p vcf BURDEN_INPUT.vcf.gz
```
### Step 1: Going to the Proper Directory
Check the files that are in the folder:
```bash
%%bash
cd GWAS_course_files/BURDEN/
ls
```
### Step 2: Using RVtests to run a burden test (One or More Genes)
We will be running a burden test using RVtests. This can be done on one gene or on a list of multiple genes by running the following: 
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/BURDEN/

rvtest --noweb --hide-covar --out burden_skat --kernel skato --burden cmc --inVcf BURDEN_INPUT.vcf.gz \
--pheno PHENOTYPES.txt --pheno-name PHENO --covar PHENOTYPES.txt \
--covar-name sex,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--geneFile refFlat_hg19.txt.gz --gene GBA
	# Optional options to add: --freqUpper 0.05
```

### Step 3: Checking the Results 
Have a look at the results by running the following: 
```bash
%%bash
# Go to the proper directory
cd GWAS_course_files/BURDEN/
ls

# Read out the SKAT results 
echo "SKAT"
head burden_skat.SkatO.assoc

# Read out the CMC results
echo "CMC"
head burden_skat.CMC.assoc
```

