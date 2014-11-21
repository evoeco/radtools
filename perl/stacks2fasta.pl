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

my $ver = "2.1";
print STDERR "This is stacks2fasta script ver. $ver\n";
print STDERR "Written by Andrey Rozenberg, Ruhr-Universitaet Bochum and distributed under GNU GPL v. 3\n";
print STDERR "If you decide to use this program for your next paper, we recommend you to cite it as follows:
Macher J et al (2014). Assessing the phylogeographic history of the montane caddisfly Thremma gallicum using mitochondrial and restriction-site associated DNA (RAD) markers. Ecology and Evolution, accepted\n";
my @cols;
my %outgs;
my %excls;
my @specimens;
my @chosen_specimens;
my $consensus;
my @snps;
my @snp_cols;
my @genotypes;
my @chr1;
my @chr2;
my $flag;
my $outg_flag;
my $seq;
my @alleles;
my @output  = ();
my @output2 = ();
my @loci;
my $fill = '';
my $snps_only = 1;
my $include_fixed = 0;
my $limit = 0;
my $format = 'fa';
my $len;
my $lim_stack;
my $homo;
my $count = 0;
my $totalcount = 0;
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
my $major = 0;
my $nodoubles = 0;
my @current_heteroz;
my @current_homoz;
my $explode = 0;
my $phased = 0;
my %formats = ( arp => 1, fa  => 1, meg => 1, nex => 1, pop => 1, snp => 1, ped => 1);
my @all_alleles;
my %popul;
my $multiallele = 0;
my $locus;
my $ploidy_df = 1;
my $polymorphic = 0;
my $first_specimen;
my $encode_snp = 0;
my $total_len = 0;

my %ambig_map = (
	'AT' => 'W', 'TA' => 'W',
	'GT' => 'K', 'TG' => 'K',
	'AC' => 'M', 'CA' => 'M',
	'AG' => 'R', 'GA' => 'R',
	'CT' => 'Y', 'TC' => 'Y',
	'CG' => 'S', 'GC' => 'S',
);

if (!@ARGV) {
	print STDERR "No input file specified\n";
	print STDERR "\nstacks2fasta.pl INPUT [OPTIONS] > OUTPUT\n";
	print STDERR "
   INPUT:
     tsv file produced by export_sql.pl ('-' for STDIN)
   OPTIONS:
     -lim      - minimal number of taxa required for a locus to be reported                 [Default: 0 (without limit)]
    -alls      - min,max interval for the number of alleles per locus                       [Default: without limit]
    -outg      - take only those loci, which are genotyped in these individuals             [comma-separated list of names, default: no preferences]
    -excl      - discard these individuals                                                  [comma-separated list of names, default: include everything]
    -popul     - population assignments                                                     [comma-separated list name1/pop1,name2/pop2,...]
    -fill      - substitute data-missing-positions with this symbol ('\\?' for '?')          [Default: 0 (the whole locus filled with '-')]
    -data      - data to output:                                                            [Default: snps]
       variable positions only (snps), whole loci (whole) or whole loci including fixed positions (whole-all)
    -explode   - output each nucleotide position as separate locus (for -data snps only)    [Default: 0]
    -form      - output format:                                                             [Default: fa]
       Arlequin (arp), Fasta (fa/fas), Mega (meg), Nexus (nex), GenePop (pop), GenePop-SNP (snp), PED (ped)
    -major     - filter out minor alleles with this frequency threshold (absolute number)   [Default: 0 (no filtering)]
    -nodoubles - filter out minor SNPs with the same pattern from the same locus            [Default: 0 (no filtering)]
    -dumm      - add 20% of dummy positions in for GenePop-SNP format (DIYABC bug)          [Default: 0 (no dummy positions)]
    -ploidy    - expected ploidy level                                                      [Default: 2]
";
	exit;
}

## BEGIN parse arguments

if ($#ARGV > 0) {
	die "Wrong number of parameters\n" if $#ARGV % 2 != 0;
	my %args = @ARGV[1..$#ARGV];
	if (defined($args{-fill})) {
		die "-fill is expected to be a single symbol\n" if length($args{-fill}) > 1;
		$fill = $args{-fill};
		delete($args{-fill});
	}
	if (defined($args{-lim})) {
		die "-lim is expected to be a number\n" if $args{-lim} !~ m/^\d+$/;
		$limit = $args{-lim};
		delete($args{-lim});
	}
	if (defined($args{-ploidy})) {
		die "-ploidy is expected to be a number\n" if $args{-ploidy} !~ m/^\d+$/;
		$ploidy_df = $args{-ploidy} - 1;
		delete($args{-ploidy});
	}
	if (defined($args{-form})) {
		$args{-form} = 'fa' if $args{-form} eq 'fas';
		die "-form: unrecognized format\n" if !defined($formats{$args{-form}});
		$format = $args{-form};
		delete($args{-form});
	}
	if (defined($args{-data})) {
		die "-data is expected to be 'snps', 'whole' or 'whole-all'\n" if $args{-data} ne 'snps' && $args{-data} ne 'whole' && $args{-data} ne 'whole-all';
		$snps_only = $args{-data} eq 'snps';
		$include_fixed = $args{-data} eq 'whole-all';
		delete($args{-data});
	}
	if (defined($args{-dumm})) {
		$add_dummy = $args{-dumm};
		delete($args{-dumm});
	}
	if (defined($args{-outg})) {
		%outgs = map { $_ => 1 } split(/,/, $args{-outg});
		delete($args{-outg});
	}
	if (defined($args{-excl})) {
		%excls = map { $_ => 1 } split(/,/, $args{-excl});
		if (%outgs) {
			foreach my $excl (%excls) {
				die "-excl and -outg intersect\n" if defined($outgs{$excl});
			}
		}
		delete($args{-excl});
	}
	if (defined($args{-alls})) {
		my @alls = split(/,/, $args{-alls});
		die "-alls is expected to be an interval 'min,max' with min < max\n" if $#alls != 1 || $alls[0] !~ /^\d+$/ || $alls[1] !~ /^\d+$/ || $alls[0] >= $alls[1];
		$allmin = $alls[0];
		$allmax = $alls[1];
		delete($args{-alls});
	}
	if (defined($args{-popul})) {
		foreach my $pair (split(/,/, $args{-popul})) {
			my @parts = split(/\//, $pair);
			die "-popul is expected to be a comma-separated list name1/pop1,name2/pop2,...\n" if $#parts != 1;
			$popul{$parts[0]} = $parts[1];
		}
		delete($args{-popul});
	}
	if (defined($args{-major})) {
		$major = $args{-major};
		delete($args{-major});
	}
	if (defined($args{-nodoubles})) {
		$nodoubles = $args{-nodoubles};
		delete($args{-nodoubles});
	}
	if (defined($args{-explode})) {
		$explode = $args{-explode};
		delete($args{-explode});
	}
	foreach my $key (keys %args) {
		die "$key: unknown parameter\n";
	}
}

my %callargs = (
	arp => \&args_arp,
	snp => \&args_snp,
	pop => \&args_pop,
	ped => \&args_ped,
);

my %callprocess = (
	fa  => \&process_fa,
	arp => \&process_arp,
	snp => \&process_snp,
	pop => \&process_pop,
	meg => \&process_meg,
	nex => \&process_nex,
	ped => \&process_ped,
);

my %callbottom = (
	fa  => \&bottom_fa,
	arp => \&bottom_arp,
	snp => \&bottom_snp,
	pop => \&bottom_pop,
	meg => \&bottom_meg,
	nex => \&bottom_nex,
	ped => \&bottom_ped,
);

$callargs{$format}->() if defined($callargs{$format});

die "'-explode 1' is compatible only with '-data snps'\n" if $explode && ($include_fixed || !$snps_only);

if ($nodoubles) {
	die "-nodoubles is incompatible with '-explode 0'\n" if !$explode;
	$include_fixed = 0;
	$major = 1;
}


## END parse arguments

## BEGIN subs

sub args_arp {
	die "For output in Arlequin format you have to specify -fill\n" if $fill eq '';
	$phased = 1;
}
sub args_snp {
	$fill = 9;
	$snps_only = 1;
	$explode = 1;
	$phased = 1;
	$encode_snp = 1;
}
sub args_pop {
	$fill = 0;
	$phased = 1;
}
sub args_ped {
	$explode = 1;
	$fill = 0;
	$phased = 1;
}

sub process_fa  {
	push(@output,  [@chr1]);
}
sub process_arp {
	push(@output,  [@chr1]);
	push(@output2, [@chr2]);
	push(@loci, $locus);
}
sub process_snp {
	push(@output,  [@chr1]);
}
sub process_pop {
	push(@output,  [@chr1]);
	push(@output2, [@chr2]);
	push(@loci, $locus);
}
sub process_meg {
	push(@output,  [@chr1]);
}
sub process_nex {
	push(@output,  [@chr1]);
}
sub process_ped {
	push(@output,  [@chr1]);
	push(@output2, [@chr2]);
	push(@loci, $locus);
}

sub parse_header {
	my $header = shift;
	chomp $header;
	my @cols = split(/\t/, $header);
	return 0 if $cols[0] !~ /Catalog/;
	@specimens = @cols[12..$#cols];
	@include_tax = ();
	@exclude_tax = ();
	SPEC: for my $i (0..$#specimens) {
		my $specimen = $specimens[$i];
		$include_tax[$i] = 0;
		$homoz[$i] = 0;
		$heteroz[$i] = 0;
		if (defined($excls{$specimen})) {
			delete $excls{$specimen};
			$exclude_tax[$i] = 1;
			next SPEC;
		}
		elsif (defined($outgs{$specimen})) {
			$include_tax[$i] = 1;
		}
		push(@chosen_specimens, $i);
		print STDERR "Warning: no population specified for specimen '$specimen'\n" if %popul && !defined($popul{$specimen}); 
	}
	$first_specimen = $chosen_specimens[0];
	return 1;
}

sub no_data {
	my $g_count = shift;
	if ($fill) {
		my $seq;
		if ($snps_only) {
			$seq = $fill x $len;
		}
		else {
			$seq = build_seq(($fill) x ($#snps + 1));
		}
		@{$chr1[$g_count]} = ($seq);
		@{$chr2[$g_count]} = ($seq) if $phased;
	}
	else {
		@{$chr1[$g_count]} = ("-" x $len);
		@{$chr2[$g_count]} = ("-" x $len) if $phased;
	}
}

sub like_in_consensus {
	my $g_count = shift;
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
	@{$chr1[$g_count]} = ($seq);
	@{$chr1[$g_count]} = ($seq) if $phased;
}

sub parse_snps {
	@snps = split(/;/, shift);
	foreach my $snp (@snps) {
		@snp_cols = split(/,/, $snp);
		$snp = $snp_cols[0];
	}
	$len = $snps_only ? $#snps+1 : length($consensus);
}

sub parse_alleles {
	my @alleles = split(/\//, shift);
	return 0 if $#alleles > $ploidy_df;
	my $g_count = shift;
	my @a = split('', $alleles[0]);
	my @b;
	$homo = ($#alleles == 0);
	if (!$homo) {
		$current_heteroz[$g_count] = 1;
		@b = split('', $alleles[1]);
		@a = make_ambig(\@a, \@b) if !$phased;
	}
	else {
		$current_homoz[$g_count] = 1;
	}
	if ($snps_only) {
		$seq = join('', @a);
		@{$chr1[$g_count]} = ($seq);
		@{$chr2[$g_count]} = ($homo ? $seq : join('', @b)) if $phased;
	}
	else {
		$seq = build_seq(@a);
		@{$chr1[$g_count]} = ($seq);
		@{$chr2[$g_count]} = ($homo ? $seq : build_seq(@b)) if $phased;
	}
	return 1;
}

sub build_seq {
	my @genotype = @_;
	my $seq = $consensus;
	for my $i (0..$#snps) {
		$seq = substr($seq, 0, $snps[$i]).$genotype[$i].substr($seq, $snps[$i] + 1);
	}
	return $seq;
}

sub make_ambig {
	my ($_a, $_b) = @_;
	my @a = @$_a;
	my @b = @$_b;
	for my $i (0..$#{a}) {
		my $cmp = $a[$i] cmp $b[$i];
		next if !$cmp;
		$a[$i] = $ambig_map{$a[$i].$b[$i]}
	}
	return @a;
}

sub convert_snps {
	my $a = shift;
	return (9) if substr($a, 0, 1) eq $fill;
	my $b = shift;
	return (0) if $a eq $all_alleles[0] and $b eq $all_alleles[0];
	return (2) if $a eq $all_alleles[1] and $b eq $all_alleles[1];
	return (1) if $a eq $all_alleles[1] and $b eq $all_alleles[0] || $a eq $all_alleles[0] and $b eq $all_alleles[1];
	print  STDERR join(", ", keys %actual_alleles)."\n";
	printf STDERR "'%s' eq '%s' and '%s' eq '%s'", $a, $all_alleles[0], $b, $all_alleles[1];
	exit;
}

sub filter_output {
	return if !$explode && !$major && !$nodoubles && !$encode_snp;
	my @alleles;
	my @patterns;
	my @prev_snps;
	my @out1;
	my @out2;
	my @first_states;
	
	foreach my $i (@chosen_specimens) {
		@{$out1[$i]} = $explode ? split("", $chr1[$i][0]) : ($chr1[$i]);
		@{$out2[$i]} = $explode ? split("", $chr2[$i][0]) : ($chr2[$i]) if $phased;
		next if !$major && !$nodoubles && !$encode_snp;
		for my $k (0..$#{$out1[$i]}) {
			next if $out1[$i][$k] =~ $fill;
			$first_states[$k] = $out1[$i][$k] if !defined($first_states[$k]);
			if ($major) {
				++$alleles[$k]{$out1[$i][$k]};
				++$alleles[$k]{$out2[$i][$k]} if $phased;
			}
			if ($nodoubles) {
				my $prev_snp = defined($prev_snps[$k]) ? $prev_snps[$k] : 1;
				$patterns[$k] .= ($out1[$i][$k] ne $prev_snp);
				$patterns[$k] .= ($out1[$i][$k] ne $out2[$i][$k]) if $phased;;
				$prev_snps[$k] = $out1[$i][$k];
			}
			if ($encode_snp) {
				if ($out1[$i][$k] ne $out2[$i][$k]) {
					$out1[$i][$k] = 1;
				}
				else {
					$out1[$i][$k] = $out1[$i][$k] eq $first_states[$k] ? 0 : 2;
				}
			}
		}
	}
	if ($major) {
		for my $k (0..$#{$out1[$first_specimen]}) {
			my $c = 0;
			foreach my $count (values %{$alleles[$k]}) {
				$c += ($count >= $major);
			}
			next if $c > 1;
			foreach my $i (@chosen_specimens) {
				delete $out1[$i][$k];
				delete $out2[$i][$k] if $phased;
			}
			delete $patterns[$k];
		}
	}
	if ($nodoubles) {
		foreach my $k1 (0..$#patterns - 1) {
			next if !$patterns[$k1];
			foreach my $k2 ($k1 + 1..$#patterns) {
				next if !$patterns[$k2] || $patterns[$k1] ne $patterns[$k2];
				foreach my $i (@chosen_specimens) {
					delete $out1[$i][$k2];
					delete $out2[$i][$k2];
				}
				delete $patterns[$k2];
			}
		}
	}
	@chr1 = @out1;
	@chr2 = @out2 if $phased;
}

sub bottom_fa {
	foreach my $i (@chosen_specimens) {
		printf ">%s\n", $specimens[$i];
		for my $loc (0..$count - 1) {
			printf "%s\n", join(' ', @{$output[$loc][$i]});
		}
	}
}

sub bottom_arp {
	printf "[Profile]
	Title=\"Stacks output\"
	NbSamples=1
	GenotypicData=1
	DataType=DNA
	GameticPhase=0
	LocusSeparator=TAB
	MissingData='%s'
[Data]
[[Samples]]\n", $fill;
	my %samples;
	my %sample_sizes;
	foreach my $i (@chosen_specimens) {
		my $first = sprintf "%s\t1\t", $specimens[$i];
		my $second = "";
		for my $loc (0..$count - 1) {
			for my $k (0..$#{$output[$loc][$i]}) {
				next if !defined($output[$loc][$i][$k]);
				$first  .= " ".(!$explode && $output[$loc][$i][$k]  =~ /$fill/ ? $fill : $output[$loc][$i][$k] );
				$second .= " ".(!$explode && $output2[$loc][$i][$k] =~ /$fill/ ? $fill : $output2[$loc][$i][$k]);
			}
		}
		my $pop = defined($popul{$specimens[$i]}) ? $popul{$specimens[$i]} : "Unassigned";
		$samples{$pop} .= sprintf("%s\n\t\t%s\n", $first, $second);
		$sample_sizes{$pop}++;
	}
	while (my ($pop, $data) = each %samples) {
		printf "	SampleName=\"%s\"
	SampleSize=%d
	SampleData={\n%s}\n", $pop, $sample_sizes{$pop}, $data;
	}
}

sub bottom_snp {
	print "<NM=1.0NF> SNP data prepared from Stacks output\nIND   SEX   POP";
	for my $loc (0..$count - 1) {
		for my $k (0..$#{$output[$loc][$first_specimen]}) {
			printf " A";
		}
	}
	if ($add_dummy) {
		print " A"x($count * 0.3);
	}
	print "\n";
	my %populations;
	my $pop_c = 0;
	foreach my $i (@chosen_specimens) {
		my $pop = defined($popul{$specimens[$i]}) ? $popul{$specimens[$i]} : "_Unassigned";
		my $pop_num = defined($populations{$pop}) ? $populations{$pop} : ++$pop_c;
		$populations{$pop} = $pop_num;
		printf "%s\t%s\tP%d", $specimens[$i], "F", $pop_num;
		for my $loc (0..$count - 1) {
			for my $k (0..$#{$output[$loc][$i]}) {
				print " ".$output[$loc][$i][$k] if defined($output[$loc][$i][$k]);
			}
		}
		if ($add_dummy) {
			print ((" ".(int(rand(3))))x($count * 0.3));
		}
		print "\n";
	}
}

sub bottom_pop {
	print "<NM=1.0NF> SNP data prepared from Stacks output\n";

	for my $loc (0..$count - 1) {
		for my $k (0..$#{$output[$loc][$first_specimen]}) {
			die "$loc\n" if !defined($loci[$loc]);
			printf "locus stacks_%d.%d\t<A>\n", $loci[$loc], $k if defined $output[$loc][$first_specimen][$k];
		}
	}

	my %samples;
	foreach my $i (@chosen_specimens) {
		my $row = sprintf "%s\t1\t", $specimens[$i];
		for my $loc (0..$count - 1) {
			for my $k (0..$#{$output[$loc][$i]}) {
				next if !defined($output[$loc][$i][$k]);
				$row .= sprintf "\t<[%s][%s]>",
					$output[$loc][$i][$k]  =~ /$fill/ ? "" : $output[$loc][$i][$k],
					$output2[$loc][$i][$k] =~ /$fill/ ? "" : $output2[$loc][$i][$k];
			}
		}
		my $pop = defined($popul{$specimens[$i]}) ? $popul{$specimens[$i]} : "_Unassigned";
		$samples{$pop} .= $row."\n";
	}
	while (my ($pop, $data) = each %samples) {
		printf "POP\n%s", $data;
	}
}

sub bottom_ped {
	my @snp_data_l;
	my @snp_data_r;
	print "#Family Individual Paternal Maternal Sex Phenotype";
	for my $loc (0..$count - 1) {
		for my $k (0..$#{$output[$loc][$first_specimen]}) {
			printf " %d.%d.1 2", $loci[$loc], $k if defined $output[$loc][$first_specimen][$k];
		}
	}
	print "\n";
	foreach my $i (@chosen_specimens) {
		my $pop = defined($popul{$specimens[$i]}) ? $popul{$specimens[$i]} : "0";
		printf "%s %s %d %d %d %d", $pop, $specimens[$i], 0, 0, 0, 0;
		for my $loc (0..$count - 1) {
			for my $k (0..$#{$output[$loc][$i]}) {
				printf " %s %s", $output[$loc][$i][$k], $output2[$loc][$i][$k] if defined $output[$loc][$i][$k];
			}
		}
		print "\n";
	}
}

sub bottom_meg {
	print "#mega
!Title Stacks output;
!Description This file was produced by stacks2fasta
!Format DataType=DNA;\n";

	foreach my $i (@chosen_specimens) {
		printf "#%s\n", $specimens[$i];
		for my $loc (0..$count - 1) {
			printf "%s\n", join(' ', @{$output[$loc][$i]});
		}
	}
}

sub bottom_nex {
	printf "#NEXUS
begin taxa;
	dimensions ntax=%d;
	taxlabels\n", $#specimens + 1;
	foreach my $i (@chosen_specimens) {
		printf "\t%s\n", $specimens[$i];
	}
	printf ";\nend;\n\nbegin characters;
	dimensions nchar=%d;
	format datatype=dna %s interleave=yes;
	matrix", $total_len, $fill ? "missing=$fill " : '';

	for my $loc (0..$count-1) {
		foreach my $i (@chosen_specimens) {
			printf "\n%s\t%s", $specimens[$i], join("\t", @{$output[$loc][$i]});
		}
		print "\n";
	}
	print "\n;\nend;";
}

## END subs

print STDERR "stacks2fasta started\n";
my $j = 0;

my $fh;
if ($ARGV[0] eq '-') {
	open $fh, '-' or die "Couldn't read from STDIN";
}
else {
	open $fh, '<', $ARGV[0] or die "Couldn't open $ARGV[0]";
}

my $header = <$fh>;
die("The input seems to have incorrect format\n") if !parse_header($header);

ROW: while (<$fh>) {
	chomp;
	@cols = split(/\t/);
	next if $#cols < 12;
	next if !$include_fixed && $cols[8] eq ""; 
	$consensus = $cols[4];
	parse_snps($cols[8]);
	@genotypes = @cols[12..($#specimens + 12)];

	my $g_count = -1;
	$lim_stack = 0;
	@current_heteroz = ();
	@current_homoz = ();
	foreach my $genotype (@genotypes) {
		$g_count++;
		next if defined($exclude_tax[$g_count]);
		if (!defined($genotype) || $genotype eq '') {
			next ROW if $include_tax[$g_count];
			no_data($g_count);
		}
		else {
			$lim_stack++;
			if ($genotype eq 'consensus') {
				like_in_consensus($g_count);
			}
			else {
				next ROW if !parse_alleles($genotype, $g_count);
			}
			$actual_alleles{$chr1[$g_count]}++;
			$actual_alleles{$chr2[$g_count]}++ if $phased;
		}
	}
	@all_alleles = keys %actual_alleles;
	next if !$include_fixed && !$#all_alleles;
	next if $allmax > 0 && ($#all_alleles <= $allmin || $#all_alleles + 1 > $allmax);
	next if $lim_stack < $limit;

	$polymorphic += ($#all_alleles > 0);
	$locus = $.;
	filter_output();
	$callprocess{$format}->($.) if defined($callprocess{$format});

	foreach my $i (@chosen_specimens) {
		$homoz[$i]   += defined($current_homoz[$i])   ? $current_homoz[$i]   : 0;
		$heteroz[$i] += defined($current_heteroz[$i]) ? $current_heteroz[$i] : 0;
	}
	$count++;
} continue {
	$totalcount++;
	@chr1 = ();
	@chr2 = ();
	%actual_alleles = ();
}

for my $loc (0..$count-1) {
	$total_len += length(join("", @{$output[$loc][$first_specimen]}));
}

$callbottom{$format}->() if defined($callbottom{$format}) && @output;

printf STDERR "Done: %d lines, %d loci, %d total length", $totalcount, $count, $total_len;
printf STDERR ", %d polymorphic loci", $polymorphic if $include_fixed;
printf STDERR "\nNum\tName\tHomozygotes\tHeterozygotes\n", $totalcount, $count;
foreach my $i (@chosen_specimens) {
	printf STDERR "%d\t%s\t%d\t%d\n", $i + 1, $specimens[$i], $homoz[$i], $heteroz[$i];
}
