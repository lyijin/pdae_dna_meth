==============================================================================
Scripts and commands used in SNP identification from Bismark deduplicated bams
==============================================================================

Scripts in this folder serve two purposes. 

1. Identifying SNPs from the deduplicated ``*.bam`` files produced by Bismark's ``deduplicate_bismark``, and processing the massive files into something that other scripts can use.

2. Plotting genetic distances based on the SNPs identified in each sample.

For #1, commands that I describe below are mainly for descriptive purposes--as I cannot upload the ``bam`` files used in the input, the commands below would not work without them.

SNP identification using Bis-SNP
--------------------------------
Bis-SNP v1.0.1 requires Java 8, so that has to be installed first.

Bis-SNP also requires the ``*.bam`` files to be sorted and with proper metadata (Bismark outputs unsorted ``*.bam``). Hence, several Picard tools were used to format the bam files into something that Bis-SNP accepted (took me a while to figure that out...). Let's take E8-AD as an example:

``samtools sort E8-AD.bam > E8-AD.sorted.bam && rm -f E8-AD.bam``

``PicardCommandLine AddOrReplaceReadGroups I=E8-AD.sorted.bam O=test.bam LB=lib1 PL=illumina PU=unit1 SM=E8``

``mv -f test.bam E8-AD.sorted.bam``

``samtools index E8-AD.sorted.bam``

Repeat these commands (with a loop) for all the samples involved.

Bis-SNP can simultaneously identify methylated positions (which I didn't need) and SNPs (which I did). Output of the former is controlled by the ``-vfn1`` flag while the latter ``-vfn2``, so the gist of the commands below is that I only retained the VCF files produced by ``vfn2`` .

``for a in dedup_bams/*.sorted.bam; do b=`echo ${a} | sed 's/sorted.bam/vfn/'` && /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Xmx10g -jar BisSNP-1.0.1.jar -T BisulfiteGenotyper -C CG,1 -I ${a} -R ../../raw_data/pdae_genome.v1.fa -vfn1 ${b}1.vcf -vfn2 snp_vcfs/${b}2.vcf & done; wait``

(wait for three days...)

``rm -f *vfn1*``

``cd snp_vcfs/``

``rename 's/.vfn2//' *``

The script ``tabulate_snp_vcfs.py`` makes its debut here.

``python3 tabulate_snp_vcfs.py snp_vcfs/*.vcf > tabulated_genotypes.tsv``

A compressed copy of that output file (``tabulated_genotypes.tsv.gz``) can be found in this folder. This file is one of the input files for ``filter_pos.five_criteria.py``, which blacklists methylated positions that happen to also be a potential SNP in ANY of the 23 files.

As homozygous reference calls ("0") could also mean the inability to call het/hom alt due to low coverage, I wrote a script to filter ``tabulated_genotypes.tsv`` for positions with coverage values of >= 10 and <= 100 in all samples. The table ``tabulated_depths.tsv`` referred to in ``filter_high_cov_pos.py`` is too large to be uploaded here, but it was produced by merging output files produced by ``samtools depth``, which calculated the coverages of SNP positions in every file:

``for a in dedup_bams/*.sorted.bam; do b=`echo $a | sed 's/dedup_bams/depths_at_snp_loci/' | sed 's/sorted.bam/tsv/'` && samtools depth ${a} -b tabulated_genotypes.bed -d 0 > $b & done; wait``

The post-coverage filtering table is available here as well (``tabulated_genotypes.cov10-100.tsv.gz``), and these positions are used to plot graphs, calculate F\ :sub:`ST` values, etc.

Plotting code
-------------
The only plotting script in this folder is ``plot_correlation.snps.py``, which produces Fig. 2a.
