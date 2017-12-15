#!/usr/bin/perl

### Written by Andrey Rozenberg (jaera at yandex.com), Ruhr-Universität Bochum

### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.

### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use POSIX;

my $ver = "0.1b";
print STDERR "This is arpsampler script ver. $ver\n";
print STDERR "Written by Andrey Rozenberg, Ruhr-Universitaet Bochum and distributed under GNU GPL v. 3\n";
print STDERR "If you decide to use this program for your next paper, we recommend you to cite it as follows:
Macher J et al (2014). Mol Ecol, under revision\n";

if (!@ARGV) {
	print STDERR "No input file specified\n";
	print STDERR "\narpsampler.pl <Input> [-min <number>] [-max <number>] [-perm <number>] [-mode <boot|jack>] [-noms <0|1>] [-step <number>] > <Output>\n";
	print STDERR "
    -min  - min size of the permuted alignments                     [default: 100]
    -max  - max size of the permuted alignments                     [default: 100]
    -perm - number of permutations to perform                       [default: 10]
    -mode - re-sampling mode: boot(strap) or jack(knifing)          [default: boot]
    -noms - no 'missing data' allowed                               [default: 0 (no)]
    -outg - focus on these individuals only (names)                 [default: take all of them]
    -step - step for permutated alignment sizes                     [default: 1]
    -alig - output alignments: fas(ta), nex(us) or phy(lip)         [default: 0 (no)]
";
	exit();
}

my $min = 100;
my $max = 100;
my $perm = 10;
my $boot = 1;
my @outgs;
my $noms = 0;
my $step = 1;
my $alig = 0;

if ($#ARGV > 0) {
	if ($#ARGV % 2 != 0) {
		print STDERR "Wrong number of parameters\n";
		exit();
	}
	my %args = @ARGV[1..$#ARGV];
	if (defined($args{-min})) {
		if ($args{-min} !~ m/^\d+$/) {
			print STDERR "-min is expected to be a number\n";
			exit();
		}
		$min = $args{-min};
		delete($args{-min});
	}
	if (defined($args{-step})) {
		if ($args{-step} !~ m/^\d+$/) {
			print STDERR "-step is expected to be a number\n";
			exit();
		}
		$step = $args{-step};
		delete($args{-step});
	}
	if (defined($args{-outg})) {
		@outgs = split(/,/, $args{-outg});
		delete($args{-outg});
	}
	if (defined($args{-mode})) {
		if ($args{-mode} ne 'boot' and $args{-mode} ne 'jack') {
			print STDERR "-mode is expected to be either 'boot' or 'jack'\n";
			exit();
		}
		$boot = ($args{-mode} eq 'boot');
		delete($args{-mode});
	}
	if (defined($args{-noms})) {
		$noms = $args{-noms};
		delete($args{-noms});
	}
	if (defined($args{-max})) {
		if ($args{-max} !~ m/^\d+$/) {
			print STDERR "-max is expected to be a number\n";
			exit();
		}
		$max = $args{-max};
		delete($args{-max});
	}
	if (defined($args{-perm})) {
		if ($args{-perm} !~ m/^\d+$/) {
			print STDERR "-perm is expected to be a number\n";
			exit();
		}
		$perm = $args{-perm};
		delete($args{-perm});
	}
	if (defined($args{-alig})) {
		if ($args{-alig} !~ m/^(0|fas|nex|phy)$/) {
			print STDERR "-alig is expected to be '0' or one of ('fas', 'nex', 'phy')\n";
			exit();
		}
		if ($args{-alig} ne '0') { $alig = $args{-alig}; }
		delete($args{-alig});
	}
	foreach my $key (keys %args) {
		print STDERR "$key: unknown parameter\n";
		exit();
	}
}

my %seqs1;
my %seqs2;
my $len;
my $name;
my $sep = "\t";
my $missing_symbol = "N";

while (<>) {
##	print STDERR substr($_, -2, 1)."\n";
	if (substr($_, -2, 1) eq '{') {
		last;
	}
	if (/MissingData=.(.)./) {
		$missing_symbol = $1;
	}
	elsif (/LocusSeparator=(.+)/) {
		if ($1 eq 'TAB' or $1 eq 'NONE') { next; }
		if ($1 eq 'WHITESPACE') { $sep = ""; }
		else { $sep = $1; }
	}
	if (eof()) {
		print STDERR "The input file has incorrect format\n";
		exit();
	}
}


my %heteroz_or_nucl;
my %nodata;
my @prev_row;
my %missing;
my $subset = @outgs;
sub merge {
	my $a = shift;
	my $b = shift;
	my $diff = ($a ne $b);
	if (!$alig) { return $diff; }
	if (!$diff) { return $a; }
	my $seq = "";
	my $len = length($a);
	for (my $i = 0; $i < $len; $i++) {
		my $nuc_a = substr($a, $i, 1);
		my $nuc_b = substr($b, $i, 1);
		if ($nuc_a eq $nuc_b) { $seq .= $nuc_a; }
		elsif (($nuc_a eq 'A' and $nuc_b eq 'T') or ($nuc_b eq 'A' and $nuc_a eq 'T')) { $seq .= 'W'; }
		elsif (($nuc_a eq 'G' and $nuc_b eq 'T') or ($nuc_b eq 'G' and $nuc_a eq 'T')) { $seq .= 'K'; }
		elsif (($nuc_a eq 'A' and $nuc_b eq 'C') or ($nuc_b eq 'A' and $nuc_a eq 'C')) { $seq .= 'M'; }
		elsif (($nuc_a eq 'A' and $nuc_b eq 'G') or ($nuc_b eq 'A' and $nuc_a eq 'G')) { $seq .= 'R'; }
		elsif (($nuc_a eq 'C' and $nuc_b eq 'T') or ($nuc_b eq 'C' and $nuc_a eq 'T')) { $seq .= 'Y'; }
		elsif (($nuc_a eq 'C' and $nuc_b eq 'G') or ($nuc_b eq 'C' and $nuc_a eq 'G')) { $seq .= 'S'; }
		else { $seq .= 'N'; }
	}
	return $seq;
}


while (<>) {
	$_ =~ s/\s+$//;
	my @row = split($sep, $_);
	if ($subset and !@outgs) { last; }
	if (@row > 2) {
		if ($row[0] ne '') {
			$name = $row[0];
			if ($subset and !grep($_ eq $name, @outgs)) { next; }
			@prev_row = @row;
			@{$nodata{$name}} = ();
			@{$heteroz_or_nucl{$name}} = ();
		}
		else {
			if ($subset) {
				if (!grep($_ eq $name, @outgs)) { next; }
				@outgs = grep($_ ne $name, @outgs);
			}
			for (my $i = 2; $i < @row; $i++) {
				if (index($row[$i], $missing_symbol) != -1 or index($prev_row[$i], $missing_symbol) != -1) {
					$nodata{$name}[$i - 2] = 1;
					$missing{$i - 2} = 1;
 					$heteroz_or_nucl{$name}[$i - 2] = ($alig ? $row[$i] : 0);
				}
				else {
					$nodata{$name}[$i - 2] = 0;
					$heteroz_or_nucl{$name}[$i - 2] =  merge($row[$i], $prev_row[$i]);
				}
			}
		}
	}
	if (eof) { last; }
}

my %chosen_nodata;
my %chosen_heteroz_or_nucl;
foreach my $outg (@outgs) {
	print STDERR "Specimen '$outg' not found\n";
	exit();
}

if (%chosen_nodata) {
	%heteroz_or_nucl = %chosen_heteroz_or_nucl;
	%nodata = %chosen_nodata;
}

my @keys = keys(%nodata);

if ($noms) {
	my @missing = sort({$b <=> $a} keys(%missing));
	foreach my $ms (@missing) {
		foreach my $key (@keys) {
			splice(@{$heteroz_or_nucl{$key}}, $ms, 1);
			splice(@{$nodata{$key}}, $ms, 1);
		}
	}
}

$len = @{$nodata{$keys[0]}};
print STDERR "Total number of positions to draw samples from is: $len\n";
if (!$boot and ($max >= $len or $min >= $len)) {
	print STDERR "In the jackknifing mode sampling size should be less than the observed number of loci ($len)\n";
	exit();
}
if (!$alig) { print "Size\tIteration\tName\tHeterozygotes\tLoci\tHeterozygosity"; }

## find max name length
my $max_name_len = 10;
my $keys_count = 0;
my $taxlabels = "";

foreach my $key (@keys) {
	$keys_count++;
	$taxlabels .= " $key";
	if (length($key) > $max_name_len) {
		$max_name_len = length($key);
	}
}
for (my $size = $min; $size <= $max; $size += $step) {
	print STDERR "\rIterating alignmets of size $size (from $min/$max with step $step)";
	my %iter;
	for (my $i = 0; $i < $perm; $i++) {
		my @sample;
		if ($boot) {
			for (my $j = 0; $j < $size; $j++) {
				push(@sample, int(rand($len)));
			}
		}
		else {
			my @numbers = (0..$len-1);
			for (my $j = 0; $j < $size; $j++) {
				my $k = int(rand(@numbers));
				push(@sample, $numbers[$k]);
				splice(@numbers, $k, 1);
			}
		}
		
		if ($alig eq 'nex') {
			print "#NEXUS[arpsampler|$size|$i]\n";
			print "begin taxa;\n";
			print "dimensions ntax=$keys_count;\n";
			print "taxlabels $taxlabels;\n";
			print "end;\n";
			print "begin characters;\n";
			my $nchar = 0;
			foreach my $key (@keys) {
				foreach my $j (@sample) {
					$nchar += length($heteroz_or_nucl{$key}[$j]);
				} last; }
			print "dimensions nchar=$nchar;\n";
			print "format datatype=dna missing=$missing_symbol interleave;\n";
			print "matrix\n";
			foreach my $j (@sample) {
				foreach my $key (@keys) {
					print "$key\t".$heteroz_or_nucl{$key}[$j]."\n";
				}
			}
			print ";\nend;\n";
			next;
		}
		
		my $fst = 1;
		foreach my $key (@keys) {
			if ($alig) {
				my $seq = "";
				my $align_len = 0;
				foreach my $j (@sample) {
					#if ($nodata{$key}[$j]) { next; }
					$seq .= $heteroz_or_nucl{$key}[$j]." ";
					$align_len += length($heteroz_or_nucl{$key}[$j]);
				}
				if ($alig eq 'fas') {
					print ">$key|$size|$i\n$seq\n";
				}
				elsif ($alig eq 'phy') {
					if ($fst) { print " $keys_count $align_len $size.$i\n"; }
					print sprintf("%-${max_name_len}s", $key)." $seq\n";
				}
			}
			else {
				my $heterozygotes = 0;
				my $total = 0;
				foreach my $j (@sample) {
					if ($nodata{$key}[$j]) { next; }
					$heterozygotes += $heteroz_or_nucl{$key}[$j];
					$total++;
				}
				print "\n$size\t$i\t$key\t$heterozygotes\t$total\t".($total > 0 ? $heterozygotes/$total : 'NA');
			}
			$fst = 0;
		}
	}
}

print STDERR "\nDone\n";
