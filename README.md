# HFI
Haplotype Fixation Index (HFI) for crop populations **with homozygous nature**, such as rice.

By Zhuo CHEN, contact: chenomics@163.com or zhuochen@genetics.ac.cn

# Motivation and description:

Considering the homozygous nature of cultivated rice, we designed a haplotype-based estimate HFI for genetic differentiation analysis. This disign of HFI estimate is inspired from fixation index [Weir and Cockerham's FST](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x) but with multiple major changes.

HFI is based on two other haplotype-based estimates, namely hapDiv and hapDist.

* Note: Due to security settings, the formula may not display properly on Microsoft Edge or IE(version 11 and above) browsers. If this happens, please open this page with another browser.

![image](https://latex.codecogs.com/gif.latex?hapDiv=\sum_{i=1}^{n}\sum_{j=i}^{n}x_{i}x_{j}d_{ij})

where *n* is the number of haplotypes in the window; *xi* and *xj* are the allele frequency of the haplotype *i* and *j*; *dij* is the genetic distance between haplotype *i* and *j*. If there are any clear base difference (excluding missing genotype or heterozygous genotype) between the two haplotypes, *dij* will be set as 1; otherwise, it will be set as zero.

![image](https://latex.codecogs.com/gif.latex?hapDist_{ab}=\frac{1}{2}\sum_{i=1}^{n}\sum_{j=1}^{n}x_{ai}x_{bj}d_{ij})

where definition of *n*, *x* and *d* in the equation are the same with those in the equation of hapDiv.

Then:

![image](https://latex.codecogs.com/gif.latex?HFI_{ab}=hapDist_{ab}-\min(hapDiv_{a},hapDiv_{b}))

# Usage:

Typical usage:

`perl HFI.pl --in pop.geno --out out.hfi --list1 pop1.list --list2 pop2.list`

For detailed usage:

`perl HFI.pl`

Input format: tab-deleimited table with header. Each line contains: chr pos geno1 geno2 ... genoX

This program was designed only for populations with homozygous nature. The genotype coding in the input file is:

0 for reference type; 1 for alternative type; - for missing or heterozygous genotype.

The input file is recommended to be created from a VCF format file and pruned with the script SNP_pruning.r2.pl

`perl SNP_pruning.r2.pl --in pop.vcf --out pop.geno`

# Reference:

Weir BS, Cockerham CC. ESTIMATING F-STATISTICS FOR THE ANALYSIS OF POPULATION STRUCTURE. Evolution. 1984;38(6):1358–1370. doi:10.1111/j.1558-5646.1984.tb05657.x

# Citation:

Chen Zhuo et al. Genomic Footprints of Breeding in a Large Chinese Rice Collection. In preparation.
