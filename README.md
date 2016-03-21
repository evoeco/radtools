radtools
========

This package provides a set of tools for data conversion and analysis of RAD (restriction site associated DNA markers).

[//]: # (Please refer to our [wiki](https://github.com/evoeco/radtools/wiki) for a full-length manual.)

## Installation

To install/update all of the tools - clone the repository and copy the scripts to the path:

    git clone https://github.com/evoeco/radtools.git
    sudo cp radtools/perl/* /usr/local/bin

## stacks2fasta

This is a utility to convert tsv-files produced by [Stacks](http://creskolab.uoregon.edu/stacks) into a number of other formats, suitable for down-stream analyses.

### Usage

	stacks2fasta INPUT [OPTIONS] > OUTPUT

	   INPUT:
	     tsv file produced by export_sql.pl ('-' for STDIN)
	   OPTIONS:
	     -lim      - minimal number of samples required for a locus to be reported              [Default: 0 (without limit)]
	    -alls      - min,max interval for the number of alleles per locus                       [Default: without limit]
	    -snps      - min,max interval for the number of SNPs per locus                          [Default: without limit]
	    -outg      - take only those loci, which are genotyped in these individuals             [comma-separated list of names, default: no preferences]
	    -excl      - discard these individuals                                                  [comma-separated list of names, default: include everything]
	    -popul     - population assignments                                                     [file with pairs 'sample population']
	    -fill      - substitute data-missing-positions with this symbol                         [Default: '' (the whole locus filled with '-')]
	    -data      - data to output:                                                            [Default: snps]
	       variable positions only ('snps')
	       whole loci (whole)
	       whole loci including fixed positions (whole-all)
	    -explode   - output each nucleotide position as separate locus (for '-data snps' only)  [Default: 0]
	    -form      - output format:                                                             [Default: fa]
	       Arlequin ('arp')
	       Fasta ('fa'/'fas')
	       Mega ('meg')
	       Nexus ('nex')
	       DIYABC GenePop ('pop')
	       DIYABC SNP ('snp')
	       PED ('ped')
	       Stacks TSV ('tsv')
	    -major     - filter out minor alleles with this frequency threshold (absolute number)   [Default: 0 (no filtering)]
	    -nodoubles - filter out SNPs with the same pattern in the same locus                    [Default: 0 (no filtering)]
	    -samplesnp - pick up only one SNP per locus:                                            [Default: 0]
	       first SNP ('first')
	       random SNP ('random') - not implemented yet
	    -ploidy    - expected ploidy level                                                      [Default: 2]

### Citation

If you decide to use this program for your next paper, we recommend you to cite it as follows:

* Macher J et al (2015). Assessing the phylogeographic history of the montane caddisfly *Thremma gallicum* using mitochondrial and restriction-site-associated DNA (RAD) markers. *Ecology and Evolution* 5(3): 648–662, [doi:10.1002/ece3.1366](http://onlinelibrary.wiley.com/doi/10.1002/ece3.1366/abstract)

## arpsampler

This script produces random samples drawn from an alignment in Arlequin (arp) format.

### Usage

    arpsampler.pl <Input> [-min <number>] [-max <number>] [-perm <number>] [-mode <boot|jack>] [-noms <0|1>] [-step <number>] > <Output>

    -min  - min size of the permuted alignments
    -max  - max size of the permuted alignments
    -perm - number of permutations to perform
    -mode - re-sampling mode: boot(strap) or jack(knifing)
    -noms - no 'missing data' allowed
    -outg - focus on these individuals only (names)
    -step - step for permutated alignment sizes

### Citation

If you decide to use this program for your next paper, we recommend you to cite it as follows:

* Macher J et al (2015). Assessing the phylogeographic history of the montane caddisfly *Thremma gallicum* using mitochondrial and restriction-site-associated DNA (RAD) markers. *Ecology and Evolution* 5(3): 648–662, [doi:10.1002/ece3.1366](http://onlinelibrary.wiley.com/doi/10.1002/ece3.1366/abstract)

