
radtools
=======

This package provides a set of tools for data conversion and analysis of RAD (restriction site associated DNA markers).

To clone:

    git clone https://github.com/har-wradim/radtools.git


stacks2fasta
--------------

This is a utility to convert tsv-files produced by Stacks (http://creskolab.uoregon.edu/stacks) into a number of other formats, suitable for down-stream analyses.

**Installation**

The program itself is written in Perl: hence it's cross-platform and requires no installation, but you need perl installed on your computer to run the script. You may optionally make the file executable and place it somewhere in the PATH.

**Usage**

The script could be run as follows (change "stacks2fasta" to "perl stacks2fasta.pl", if the file is not executable).

    stacks2fasta.pl <Input> [-lim dd] [-selz 0/1] [-alls dd,dd] [-outg dd,dd,...] [-fill w] [-snps 0/1] [-form www] > <Output>\n";

     -lim - minimal number of taxa required for a locus to be reported                 [Default: 0 (without limit)]
    -alls - min,max interval for the number of alleles per locus                       [Default: without limit]
    -selz - output zygosity based on filtered loci only                                [Default: 0 (count all loci)]
    -outg - take only those loci, which are genotyped in these individuals             [1-based indices, default: no preferences]
    -excl - discard these individuals                                                  [1-based indices, default: include everything]
    -fill - substitute data-missing-positions with this symbol ('\?' for '?')          [Default: 0 (the whole locus filled with '-')]
    -snps - output variable positions only                                             [Default: 0 (output whole loci)]
    -form - output format: Arlequin (arp), Fasta (fa/fas), Mega (meg), Nexus (nex), GenePop (pop), GenePop-SNP (snp)   [Default: fa]
    -dumm - add 20% of dummy positions in for GenePop-SNP format (DIYABC bug)          [Default: 0 (no dummy positions)]

**Citation**

If you decide to use this program for your next paper, we recommend you to cite it as follows:
Macher J et al (2014). Mol Ecol, under revision

arpsampler
--------------

This script produces random samples drawn from an alignment in Arlequin (arp) format.

**Installation**

Same as for stacks2fasta.

**Usage**

    arpsampler <Input> [-min <number>] [-max <number>] [-perm <number>] [-mode <boot|jack>] [-noms <0|1>] [-step <number>] > <Output>

    -min  - min size of the permuted alignments            [default: 100]
    -max  - max size of the permuted alignments            [default: 100]
    -perm - number of permutations to perform              [default: 10]
    -mode - re-sampling mode: boot(strap) or jack(knifing) [default: boot]
    -noms - no 'missing data' allowed                      [default: 0 (no)]
    -outg - focus on these individuals only (names)        [default: take all of them]
    -step - step for permutated alignment sizes            [default: 1]

**Citation**

If you decide to use this program for your next paper, we recommend you to cite it as follows:
Macher J et al (2014). Mol Ecol, under revision

