
radtools
=======

This package provides a set of tools for data conversion and analysis of RAD (restriction site associated DNA markers).

[//]: # (Please refer to our [wiki](https://github.com/evoeco/radtools/wiki) for a full-length manual.)

**Installation**

To install/update all of the tools - clone the repository and copy the scripts to the path:

        git clone https://github.com/evoeco/radtools.git
        sudo cp radtools/perl/* /usr/local/bin

stacks2fasta
--------------

This is a utility to convert tsv-files produced by [Stacks](http://creskolab.uoregon.edu/stacks) into a number of other formats, suitable for down-stream analyses.

**Usage**

        stacks2fasta.pl INPUT [OPTIONS] > OUTPUT

        INPUT:
         tsv file produced by stacks' export_sql.pl ('-' for STDIN)
        OPTIONS:
           -lim      - minimal number of taxa required for a locus to be reported
          -alls      - min,max interval for the number of alleles per locus
          -outg      - take only those loci, which are genotyped in these individuals
          -excl      - discard these individuals
          -popul     - population assignments
          -fill      - substitute data-missing-positions with this symbol
          -data      - data to output:
              variable positions only (snps), whole loci (whole) or whole loci including fixed positions (whole-all)
          -explode   - output each nucleotide position as separate locus (for -data snps only)
          -form      - output format:
              Arlequin (arp), Fasta (fa/fas), Mega (meg), Nexus (nex), GenePop (pop), GenePop-SNP (snp), PED (ped)
          -major     - filter out minor alleles with this frequency threshold (absolute number)
          -nodoubles - filter out minor SNPs with the same pattern from the same locus
          -dumm      - add 20% of dummy positions in for GenePop-SNP format (DIYABC bug)
          -ploidy    - expected ploidy level

**Citation**

If you decide to use this program for your next paper, we recommend you to cite it as follows:

* Macher J et al (2014). Assessing the phylogeographic history of the montane caddisfly *Thremma gallicum* using mitochondrial and restriction-site associated DNA (RAD) markers. Ecology and Evolution, accepted

arpsampler
--------------

This script produces random samples drawn from an alignment in Arlequin (arp) format.

**Usage**

        arpsampler.pl <Input> [-min <number>] [-max <number>] [-perm <number>] [-mode <boot|jack>] [-noms <0|1>] [-step <number>] > <Output>

        -min  - min size of the permuted alignments
        -max  - max size of the permuted alignments
        -perm - number of permutations to perform
        -mode - re-sampling mode: boot(strap) or jack(knifing)
        -noms - no 'missing data' allowed
        -outg - focus on these individuals only (names)
        -step - step for permutated alignment sizes

**Citation**

If you decide to use this program for your next paper, we recommend you to cite it as follows:

* Macher J et al (2014). Assessing the phylogeographic history of the montane caddisfly *Thremma gallicum* using mitochondrial and restriction-site associated DNA (RAD) markers. Ecology and Evolution, accepted

