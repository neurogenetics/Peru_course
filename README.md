# Peru_course
Tuesday 17 September 2019

Jupyter notebooks, very kindly provided by Broad (Thank you very much!!!)

Login is https://notebook.hail.is if there are problems use https://notebook.hail.is/new

password will be provided during the course.

# Part1: Clean GWAS data....

Two step process
1) cleans data 
2) harmonizes data with reference panel (no time for this today, but you can check STEP2 script)

```
# sh STEP1_cleaning.sh INPUTFILE
does the following:
- call rate filter
- heterozygosity filter
- sex-check
- variant filter
- variant filter between cases and controls

optional
- ancestry filter
- relatedness filter

# sh STEP2_cleaning.sh INPUTFILE (we dont have time for this one)
does the following:
- make sure the right alleles are used
- make sure the variants between input and reference are matches
- makes zipped vcf files of cleaned input files per chromosome
- and more small things
```

## step 1 check to inventory
```
%%bash
cd GWAS_course_files/QC_PACKAGE/
ls
```
## step 2 edit files
```click and open STEP1_cleaning.sh and edit
change ./plink into plink
change ./gcta into gcta64
then file save
```

```click and open PCA_in_R.R and edit
pdf("FILE.pdf") -> jpeg("FILE.jpeg")
then file save
```

## step 3 run it...
```%%bash
cd GWAS_course_files/QC_PACKAGE/
sh STEP1_cleaning.sh EXAMPLE_DATA
```

## step 4 (optional) to vizualize files
``` check images:
from IPython.display import Image
Image(filename="GWAS_course_files/QC_PACKAGE/raw_hapmap_plot.jpeg")
```
``` check files:
%%bash
cd GWAS_course_files/QC_PACKAGE/
head PCA_filtered_europeans.txt
OR
%%bash
cd GWAS_course_files/QC_PACKAGE/
head GENDER_FAILURES.txt
OR
%%bash
cd GWAS_course_files/QC_PACKAGE/
head CALL_RATE_OUTLIERS.txt
```

## Make principal components based on cleaned data....

```
%%bash
cd GWAS_course_files/QC_PACKAGE/
plink --bfile FILTERED.EXAMPLE_DATA --maf 0.05 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt \
--make-bed --out pass1
# filter variants
plink --bfile pass1 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
# Extract pruned SNPs and only these variants will be used for PC calculation
plink --bfile pass1 --extract pruned_data.prune.in --make-bed --out pass1_pruned 
# Calculate/generate PCs based on pruned data set
plink --bfile pass1_pruned --pca --out PCA
# then use the .eigenvec file
```
```
# look at them
%%bash
cd GWAS_course_files/QC_PACKAGE/
head PCA.eigenvec
```

# Part2: Run GWAS....


Part3: Run genetic risk score....


Part4: Run burden test....



