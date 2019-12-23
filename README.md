# HFI
Haplotype Fixation Index (HFI) for crop populations with homozygous nature, such as rice.

By Zhuo CHEN, contact: chenomics@163.com or zhuochen@genetics.ac.cn

# Motivation and description:

Considering the homozygous nature of cultivated rice, we designed a haplotype-based estimate HFI for genetic differentiation analysis. This disign of HFI estimate is inspired from fixation index [Weir and Cockerham's FST](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x) but with multiple major changes.

HFI is based on two other haplotype-based estimates, namely hapDiv and hapDist.

![](http://latex.codecogs.com/gif.latex?\hapDiv=\sum_{i=1}^{n}\sum_{j=i}^{n}x_{i}x_{j}d_{ij})

# Usage:

Typical usage:

`perl HFI.pl --in pop.geno --out out.hfi --list1 pop1.list --list2 pop2.list`

For detailed usage:

`perl HFI.pl`

Input format: tab-deleimited table with header. Each line contains: chr pos geno1 geno2 ... genoX

This program was designed only for populations with homozygous nature. The genotype coding in the input file is:

0 for reference type; 1 for alternative type; - for missing or heterozygous genotype.

# Reference:

Weir BS, Cockerham CC. ESTIMATING F-STATISTICS FOR THE ANALYSIS OF POPULATION STRUCTURE. Evolution. 1984;38(6):1358–1370. doi:10.1111/j.1558-5646.1984.tb05657.x

# Citation:

Chen Zhuo et al. Genomic Footprints of Breeding in a Large Chinese Rice Collection. In preparation.
