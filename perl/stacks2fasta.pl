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
### along with this program. If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

my $ver = "0.2b";
print STDERR "This is stacks2fasta script ver. $ver\n";
print STDERR "Written by Andrey Rozenberg, Ruhr-Universitaet Bochum and distributed under GNU GPL v. 3\n";
print STDERR "If you decide to use this program for your next paper, we recommend you to cite it as follows:
Macher J et al (2014). Mol Ecol, under revision\n";
my @cols;
my @outgs;
my @excls = ();
my @specimens;
my @chosen_specimens;
my $consensus;
my @snps;
my @snp_cols;
my @genotypes;
my @to_print;
my @to_print2;
my $flag;
my $outg_flag;
my $seq;
my @alleles;
my @output;
my @output2;
my @a;
my @b;
my $selz = 0;
my $fill = '';
my $snps_only = 0;
my $limit = 0;
my $format = 'fa';
my $len;
my $lim_stack;
my $homo;
my $count = 0;
my $totalcount = 0;
my $total_len = 0;
my @include_tax;
my @exclude_tax;
my @heteroz;
my @homoz;
my %actual_alleles;
my $allmin = 0;
my $allmax = 0;
my $add_dummy = 0;
my $minor_thr;
my $arp_minor_thr = 0;
my $arp_no_double = 0;

if (!@ARGV) {
	print STDERR "No input file specified\n";
	print STDERR "\nstacks2fasta.pl <Input> [-lim dd] [-selz 0/1] [-alls dd,dd] [-outg dd,dd,...] [-fill w] [-snps 0/1] [-form www] > <Output>\n";
	print STDERR "
     -lim - minimal number of taxa required for a locus to be reported                 [Default: 0 (without limit)]
    -alls - min,max interval for the number of alleles per locus                       [Default: without limit]
    -selz - output zygosity based on filtered loci only                                [Default: 0 (count all loci)]
    -outg - take only those loci, which are genotyped in these individuals             [1-based indices, default: no preferences]
    -excl - discard these individuals                                                  [1-based indices, default: include everything]
    -fill - substitute data-missing-positions with this symbol ('\\?' for '?')          [Default: 0 (the whole locus filled with '-')]
    -snps - output variable positions only                                             [Default: 0 (output whole loci)]
    -form - output format: Arlequin (arp), Fasta (fa/fas), Mega (meg), Nexus (nex), GenePop (pop), GenePop-SNP (snp)   [Default: fa]
    -dumm - add 20% of dummy positions in for GenePop-SNP format (DIYABC bug)          [Default: 0 (no dummy positions)]\n";
	exit();
}

if ($#ARGV > 0) {
	if ($#ARGV % 2 != 0) {
		print STDERR "Wrong number of parameters\n";
		exit();
	}
	my %args = @ARGV[1..$#ARGV];
	if (defined($args{-fill})) {
		if (length($args{-fill}) > 1) {
			print STDERR "-fill is expected to be a single symbol\n";
			exit();
		}
		$fill = $args{-fill};
		delete($args{-fill});
	}
	if (defined($args{-lim})) {
		if ($args{-lim} !~ m/^\d+$/) {
			print STDERR "-lim is expected to be a number\n";
			exit();
		}
		$limit = $args{-lim};
		delete($args{-lim});
	}
	if (defined($args{-form})) {
		if ($args{-form} eq 'fas') { $args{-form} = 'fa'; }
		if ($args{-form} ne 'fa' and $args{-form} ne 'arp' and $args{-form} ne 'meg' and $args{-form} ne 'nex' and $args{-form} ne 'snp' and $args{-form} ne 'pop') {
			print STDERR "-form: currently supported formats are 'arp', 'fa'/'fas', 'meg', 'nex', 'pop' and 'snp'\n";
			exit();
		}
		$format = $args{-form};
		delete($args{-form});
	}
	if (defined($args{-snps})) {
		$snps_only = $args{-snps};
		delete($args{-snps});
	}
	if (defined($args{-selz})) {
		$selz = $args{-selz};
		delete($args{-selz});
	}
	if (defined($args{-dumm})) {
		$add_dummy = $args{-dumm};
		delete($args{-dumm});
	}
	if (defined($args{-outg})) {
		@outgs = split(/,/, $args{-outg});
		foreach my $outg (@outgs) {
			if ($outg !~ m/^\d+$/)
			{
				print STDERR "-outg is expected to be a sequence of numeric indices\n";
				exit();
			}
		}
		delete($args{-outg});
	}
	if (defined($args{-excl})) {
		@excls = split(/,/, $args{-excl});
		foreach my $excl (@excls) {
			if ($excl !~ m/^\d+$/) {
				print STDERR "-excl is expected to be a sequence of numeric indices\n";
				exit();
			}
		}
		delete($args{-excl});
	}
	if (defined($args{-alls})) {
		my @alls = split(/,/, $args{-alls});
		if (($#alls != 1) or ($alls[0] !~ m/^\d+$/) or ($alls[1] !~ m/^\d+$/) or ($alls[0] >= $alls[1]))
		{
			print STDERR "-alls is expected to be an interval 'min,max' with min < max\n";
			exit();
		}
		$allmin = $alls[0];
		$allmax = $alls[1];
		delete($args{-alls});
	}
	if (defined($args{-arp_minor_thr})) {
		$arp_minor_thr = $args{-arp_minor_thr};
		delete($args{-arp_minor_thr});
	}
	if (defined($args{-arp_no_double})) {
		$arp_no_double = $args{-arp_no_double};
		delete($args{-arp_no_double});
	}
	foreach my $key (keys %args) {
		print STDERR "$key: unknown parameter\n";
		exit();
	}
}

if ($format eq 'arp' and $fill eq '') {
	print STDERR "For output in Arlequin format you have to specify -fill\n";
	exit();
}
if ($format eq 'snp') {
	$allmin = 1;
	$allmax = 2;
	$fill = '-';
	$snps_only = 1;
}

if ($format eq 'pop') {
	$fill = "0";
}

my $phased = ($format eq 'arp' or $format eq 'snp' or $format eq 'pop');
my @all_alleles;
print STDERR "stacks2fasta started\n";
my $j = 0;
while (my $row = <>) {
	if (eof) { last; }
	$totalcount++;
	chomp $row;
	@to_print = ();
	@to_print2 = ();
	%actual_alleles = ();
	@cols = split(/\t/, $row);
	if ($#cols < 12) {
		next;
	}
	if ($cols[0] =~ /Catalog ID/) {
		@specimens = @cols[12..$#cols];
		@include_tax = ();
		@exclude_tax = ();
		for my $i (0..$#specimens) {
			@{$output[$i]} = ();
			$include_tax[$i] = 0;
			$homoz[$i] = 0;
			$heteroz[$i] = 0;
			for my $j (0..$#excls) {
				if ($excls[$j] - 1 == $i) {
					$exclude_tax[$i] = "1";
					last;
				}
			}
			if (!defined($exclude_tax[$i])) {
				push(@chosen_specimens, $i);
			}
		}
		for my $i (0..$#outgs) {
			$include_tax[$outgs[$i] - 1] = 1; }
		next;
	}
	if ($cols[8] eq "") {
		next;
	}
	$consensus = $cols[4];
	@snps = split(/;/, $cols[8]);
	foreach my $snp (@snps) {
		@snp_cols = split(/,/, $snp);
		$snp = $snp_cols[0];
	}
	if ($snps_only) {
		$len = $#snps+1;
	}
	else {
		$len = length($consensus);
	}
	@genotypes = @cols[12..($#specimens + 12)];
	my $g_count = -1;
	$lim_stack = 0;
	$flag = 0;
	$outg_flag = 1;
	my @current_heteroz = ();
	my @current_homoz = ();
	foreach my $genotype (@genotypes) {
		$g_count++;
		if (defined($exclude_tax[$g_count])) {
			next;
		}
		if (!defined($genotype) or $genotype eq '') {
		### if no genotype ###
			if ($include_tax[$g_count]) {
			### this taxon should have been included, but we have no genotype for it ###
				$outg_flag = 0;
				next;
			}
			if ($fill ne "0" and $fill ne "") {
				if ($snps_only) {
					$seq = $fill x $len;
				}
				else {
					$seq = $consensus;
					foreach my $snp (@snps) {
						$seq = substr($seq, 0, $snp).$fill.substr($seq, $snp + 1);
					}
				}
				$to_print[$g_count] = $seq;
				if ($phased) { $to_print2[$g_count] = $seq; }
			}
			else {
				$to_print[$g_count] = "-" x $len;
			}
		}
		else {
		### if genotype is present ###
			$lim_stack++;
			if ($genotype eq 'consensus') {
			### it's like in consensus ###
				$current_homoz[$g_count] = 1;
				if ($snps_only) {
					$seq = "";
					foreach my $snp (@snps) {
						$seq .= substr($consensus, $snp, 1);
					}
				}
				else {
					$seq = $consensus;
				}
				$to_print[$g_count] = $seq;
				if ($phased) { $to_print2[$g_count] = $seq; }
			}
			else {
			### it's specified explicitely ###
				@alleles = split(/\//, $genotype);
				if ($#alleles > 1) {
				### more than 2 alleles - discard the locus ###
					$flag = 0;
					last;
				}
				$flag = 1;
				@a = split('', $alleles[0]);
				$homo = ($#alleles == 0);
				if (!$homo) {
				### if two alleles ###
					$current_heteroz[$g_count] = 1;
					@b = split('', $alleles[1]);
					if (!$phased) {
					### make ambiguities ###
						for my $i (0..$#a) {
							if ($a[$i] eq $b[$i]) { next; }
							   if (($a[$i] eq 'A' and $b[$i] eq 'T') or ($b[$i] eq 'A' and $a[$i] eq 'T'))
							{ $a[$i] = 'W'; }
							elsif (($a[$i] eq 'G' and $b[$i] eq 'T') or ($b[$i] eq 'G' and $a[$i] eq 'T'))
							{ $a[$i] = 'K'; }
							elsif (($a[$i] eq 'A' and $b[$i] eq 'C') or ($b[$i] eq 'A' and $a[$i] eq 'C'))
							{ $a[$i] = 'M'; }
							elsif (($a[$i] eq 'A' and $b[$i] eq 'G') or ($b[$i] eq 'A' and $a[$i] eq 'G'))
							{ $a[$i] = 'R'; }
							elsif (($a[$i] eq 'C' and $b[$i] eq 'T') or ($b[$i] eq 'C' and $a[$i] eq 'T'))
							{ $a[$i] = 'Y'; }
							elsif (($a[$i] eq 'C' and $b[$i] eq 'G') or ($b[$i] eq 'C' and $a[$i] eq 'G'))
							{ $a[$i] = 'S'; }
						}
					}
				}
				else {
					$current_homoz[$g_count] = 1;
				}
				if ($snps_only) {
					$to_print[$g_count] = join('', @a);
					if ($phased) {
						if ($homo) { $to_print2[$g_count] = join('', @a); }
						else { $to_print2[$g_count] = join('', @b); }
					}
				}
				else {
					$seq = $consensus;
					for my $i (0..$#snps) {
						$seq = substr($seq, 0, $snps[$i]).$a[$i].substr($seq, $snps[$i] + 1);
					}
					$to_print[$g_count] = $seq;
					if ($phased) {
						if ($homo) { $to_print2[$g_count] = $seq; }
						else {
							$seq = $consensus;
							for my $i (0..$#snps) {
								$seq = substr($seq, 0, $snps[$i]).$b[$i].substr($seq, $snps[$i] + 1);
							}
							$to_print2[$g_count] = $seq;
						}
					}
				}
			}
			if ($allmax > 0) {
				$actual_alleles{$to_print[$g_count]} = 1;
				if (defined($to_print2[$g_count])) {
					$actual_alleles{$to_print2[$g_count]} = 1;
				}
			}
		}
	}
	if (!$flag) { next; }
	if ($allmax > 0) {
		my $num_alleles = scalar(values %actual_alleles);
		if ($num_alleles < $allmin or $num_alleles > $allmax) {
			next;
		}
	}
	if ($outg_flag and $lim_stack >= $limit or !$selz) {
		foreach my $i (@chosen_specimens) {
			if (defined($current_homoz[$i])) {
				$homoz[$i] += $current_homoz[$i]; }
			elsif (defined($current_heteroz[$i])) {
				$heteroz[$i] += $current_heteroz[$i]; }
		}
	}
	if (!$outg_flag or $lim_stack < $limit) { next; }
	foreach my $i (@chosen_specimens) {
		if ($format eq 'snp') {
			@all_alleles = keys %actual_alleles;
			if (substr($to_print[$i], 0, 1) eq $fill) {
				$to_print[$i] = "9";
			}
			elsif ($to_print[$i] eq $all_alleles[0] and $to_print2[$i] eq $all_alleles[0]) {
				$to_print[$i] = "0";
			}
			elsif ($to_print[$i] eq $all_alleles[1] and $to_print2[$i] eq $all_alleles[1]) {
				$to_print[$i] = "2";
			}
			elsif ($to_print[$i] eq $all_alleles[1] and $to_print2[$i] eq $all_alleles[0] or $to_print[$i] eq $all_alleles[0] and $to_print2[$i] eq $all_alleles[1]) {
				$to_print[$i] = "1";
			}
			else {
				print STDERR join(", ", keys %actual_alleles)."\n";
				die "'$to_print[$i]' eq '$all_alleles[0]' and '$to_print2[$i]' eq '$all_alleles[1]'";
			}
		}
		push(@{$output[$i]}, $to_print[$i]);
		if ($format eq 'arp' or $format eq 'pop') {
			push(@{$output2[$i]}, $to_print2[$i]);
		}
	}
	$count++;
	$total_len += length($to_print[$chosen_specimens[0]]);
}
if ($format eq 'fa') {
	foreach my $i (@chosen_specimens) {
		print ">".$specimens[$i]."\n".join("\n", @{$output[$i]})."\n";
	}
}
elsif ($format eq 'arp')
{
	print "[Profile]
	Title=\"Stacks output\"
	NbSamples=1
	GenotypicData=1
	DataType=DNA
	GameticPhase=0
	LocusSeparator=TAB
	MissingData='".$fill."'
[Data]
[[Samples]]
	SampleName=\"Stacks sample\"
	SampleSize=".($#specimens+1)."
	SampleData={
";
	foreach my $i (@chosen_specimens) {
		print $specimens[$i]."\t1\t".join("\t", @{$output[$i]})."\n";
		print "\t\t".join("\t", @{$output2[$i]})."\n";
	}
}
elsif ($format eq 'snp') {
	print "<NM=1.0NF> SNP data prepared for Stacks output ($count positions)
IND   SEX   POP".(" A"x($#{$output[0]}+1));
	if ($add_dummy) {
		print " A"x($count * 0.2);
	}
	print "\n";
	foreach my $i (@chosen_specimens) {
		print $specimens[$i]."\tF\tP1";
		my @left = @{$output[$i]};
		for my $j (0..$#left) {
			print " ".$left[$j];
		}
		if ($add_dummy) {
			print ((" ".(int(rand(3))))x($count * 0.2));
		}
		print "\n";
	}
}
elsif ($format eq 'pop') {
	print "<NM=1.0NF> SNP data prepared for Stacks output\n";
	for my $i (0..$#{$output[0]}) {
		print "locus stacks_$i\t<A>\n";
	}
	print "POP\n";
	foreach my $i (@chosen_specimens) {
		print $specimens[$i]." ,";
		my @left = @{$output[$i]};
		my @right = @{$output2[$i]};
		for my $j (0..$#left) {
			if (substr($left[$j], 0, 1) eq "-") {
				$left[$j] = "";
				$right[$j] = "";
			}
			print "\t<[".$left[$j]."][".$right[$j]."]>";
		}
		print "\n";
	}
}
elsif ($format eq 'meg') {
	print "#mega
!Title Stacks output;
!Description This file was produced by stacks2fasta with the following setting:
  -lim: $limit
 -fill: $fill
 -snps: $snps_only ;
!Format DataType=DNA;\n";
	foreach my $i (@chosen_specimens) {
		print "#".$specimens[$i]."\n".join("\n", @{$output[$i]})."\n";
	}
}
elsif ($format eq 'nex') {
	print "#NEXUS
begin taxa;
	dimensions ntax=".($#specimens + 1).";
	taxlabels\n";
	foreach my $i (@chosen_specimens) {
		print "\t".$specimens[$i]."\n";
	}
	print ";\nend;\n\nbegin characters;
	dimensions nchar=$total_len;
	format datatype=dna ".($fill ? "missing=$fill " : '')."interleave=yes;
	matrix";
	for my $j (0..@{$output[0]}-1) {
		foreach my $i (@chosen_specimens) {
			print "\n".$specimens[$i]."\t".$output[$i][$j];
		}
		print "\n";
	}
	print "\n;\nend;"
}
print STDERR "Done: $totalcount lines, $count loci\n";
print STDERR "Num\tName\tHomozygotes\tHeterozygotes\n";
foreach my $i (@chosen_specimens) {
	print STDERR ($i + 1)."\t".$specimens[$i]."\t".$homoz[$i]."\t".$heteroz[$i]."\n";
}
