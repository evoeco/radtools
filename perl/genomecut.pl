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

my %args = @ARGV;
if (!@ARGV or defined($args{-h}) or defined($args{-v})) {
	print STDERR "Usage:
perl genomecut.pl (-r1 <enzyme> | -r1 <enzyme1> -r2 <enzyme2>) -t <task> (-fa <fasta> | -fq <fastq>) \
      (-min <length>) (-max <length>) (-l <length> | -cat <y/n>) (-no <y/n>)

Arguments:
	-r1 - the first (or the sole) enzyme/file for double digest
	-r2 - the second              enzyme/file for double digest
	-fa - input data in (optionally gzipped) fasta format ('-' for STDIN)
	-fq - input data in (optionally gzipped) fastq format ('-' for STDIN)
	-t  - task, one of
		'dist'  - distribution of recognition sites          (single digest)
		'frag'  - sequences        for restriction fragments (double digest)
		'lens'  - sequence lengths for restriction fragments (double digest)
		'ldist' - summary for 'lens'                         (double digest)
	-min - minimal length of restriction fragments
	-max - maximal length of restriction fragments
	-cat - concatenate all input sequences before cutting
	-l   - minimal input sequence length for digest
        -no  - remove N's before cutting the sequences\n";
	exit 0;
}

my ($r1, $r2);
my ($atgc_total, $gc_total, $len_total, $len_total_sq, $seq_n);
my (%frag_total, %frag_total_sq, %frag_n);
my $min = 0;
my $max = 0;
my @seqlens;
my $t = "dist";
my $cat = 0;
my $double = 0;
my $l = 0;
my $file;
my $fq = 0;
my $no = 0;
my %enz;
my (@pat1, @pat2);
my %counts;
my %tsubs = (
	dist  => \&do_dist,
	frag  => \&do_frag,
	lens  => \&do_lens,
	ldist => \&do_ldist,
);

if (defined($args{-r1})) {
	$r1 = $args{-r1};
	delete($args{-r1});
}
else {
	die("-r1 not specified\n");
}
if (defined($args{-r2})) {
	$r2 = $args{-r2};
	delete($args{-r2});
	$double = 1;
}
if (defined($args{-t})) {
	$t = $args{-t};
	(!defined($tsubs{$t})) && die("Unrecognized -t\n");
	delete($args{-t});
}
if (defined($args{-l})) {
	$l = $args{-l};
	($l !~ m/^\d+$/) && die("Incorrect -l\n");
	delete($args{-l});
}
if (defined($args{-min})) {
	$min = $args{-min};
	($min !~ m/^\d+$/) && die("Incorrect -min\n");
	delete($args{-min});
}
if (defined($args{-max})) {
	$max = $args{-max};
	($max !~ m/^\d+$/ || $max < $min) && die("Incorrect -max\n");
	delete($args{-max});
}
if (defined($args{-no})) {
	($args{-no} ne "y" && $args{-no} ne "n") && die("Incorrect -no\n");
	$no = ($args{-no} eq "y");
	delete($args{-no});
}
if (defined($args{-cat})) {
	($args{-cat} ne "y" && $args{-cat} ne "n") && die("Incorrect -cat\n");
	$cat = ($args{-cat} eq "y");
	delete($args{-cat});
}
if (defined($args{-fa})) {
	$file = $args{-fa};
	$fq = 0;
	delete($args{-fa});
}
elsif (defined($args{-fq})) {
	$file = $args{-fq};
	$fq = 1;
	delete($args{-fq});
}
else {
	die("No input sequences specified\n");
}
foreach my $key (keys %args) {
	die("$key: unrecognized parameter\n");
}

sub gc {
	my $dna = shift;
	$dna =~ s/[^ACGTSW]//g;
	my $len = length($dna);
	$dna =~ s/[ATW]//g;
	my $gc  = length($dna);
	return ($len, $gc);
}

sub revcomp {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTMRYKVHDB/TGCAKYRMBDHV/;
	return $revcomp;
}

sub treatseq {
	my $seq = shift;
	$seq = uc($seq);
	$seq =~ tr/XU*/NTN/;
	$seq =~ s/\s//g;
	$seq =~ s/N//g if $no;
	return $seq;
}

sub checkenzyme {
	my $dna = shift;
	$dna = uc($dna);
	$dna =~ tr/XIU^/NHT /;
	$dna =~ s/\s//g;
	($dna =~ m/[^ACGTMRWSYKVHDBN]/) && die("Error! The pattern $dna contains unexpected characters!\n");
	($dna =~ m/[MRWSYKVHDBN]/) && print STDERR "Warning! The pattern $dna contains ambiguities!\n";
	($dna ne revcomp($dna)) && print STDERR "Warning! The pattern $dna is not palindromic!\n";
	return $dna;
}

sub getpattern {
	my $pattern = shift;
	my @in  = ('M',    'R',    'W',    'S',    'Y',    'K',    'V',      'H',      'D',      'B',      'N');
	my @out = ('A|C|M','A|G|R','A|T|W','C|G|S','C|T|Y','G|T|K','A|C|G|V','A|C|T|H','A|G|T|D','C|G|T|B','.');
	for (my $i = 0; $i <= $#in; $i++) {
		$pattern =~ s/\Q$in[$i]/(?:$out[$i])/g;
	}
	return $pattern;
}

sub do_frag {
	my $seq = treatseq(shift);
	my $seq_len = length($seq);
	($seq_len <= $l) && return;
	my ($len, $gc) = gc($seq);
	add_seqlen($len, $gc, $seq_len);
	foreach my $p1 (@pat1) {
		foreach my $p2 (@pat2) {
			while ($seq =~ /($p1(?:(?:(?!$p2|$p1).)+)$p2|$p2(?:(?:(?!$p1|$p2).)+)$p1)/g) {
				printf "%s\t%s\t%s\n", $enz{$p1}, $enz{$p2}, $1;
			}
		}
	}
}

sub do_lens {
	my $seq = treatseq(shift);
	my $seq_len = length($seq);
	($seq_len <= $l) && return;
	my ($len, $gc) = gc($seq);
	my $frag_len;
	add_seqlen($len, $gc, $seq_len);
	foreach my $p1 (@pat1) {
		foreach my $p2 (@pat2) {
			while ($seq =~ /($p1(?:(?:(?!$p2|$p1).)+)$p2|$p2(?:(?:(?!$p1|$p2).)+)$p1)/g) {
				$frag_len = length($1);
				($frag_len < $min || $max && $frag_len > $max) && next;
				printf "%s\t%s\t%s\n", $enz{$p1}, $enz{$p2}, $frag_len;
			}
		}
	}
}

sub do_ldist {
	my $seq = treatseq(shift);
	my $seq_len = length($seq);
	($seq_len <= $l) && return;
	my ($len, $gc) = gc($seq);
	my $frag_len;
	add_seqlen($len, $gc, $seq_len);
	foreach my $p1 (@pat1) {
		foreach my $p2 (@pat2) {
			while ($seq =~ /($p1(?:(?:(?!$p2|$p1).)+)$p2|$p2(?:(?:(?!$p1|$p2).)+)$p1)/g) {
				$frag_len = length($1);
				($frag_len < $min || $max && $frag_len > $max) && next;
				$frag_total{$p1}{$p2}    += $frag_len;
				$frag_total_sq{$p1}{$p2} += $frag_len * $frag_len;
				$frag_n{$p1}{$p2}++;
			}
		}
	}
}

sub add_seqlen {
	$atgc_total += shift;
	$gc_total   += shift;
	my $seq_len = shift;
	$len_total    += $seq_len;
	$len_total_sq += $seq_len * $seq_len;
	$seq_n++;
}

sub do_dist {
	my $seq = treatseq(shift);
	my $seq_len = length($seq);
	($seq_len <= $l) && return;
	my ($len, $gc) = gc($seq);
	add_seqlen($len, $gc, $seq_len);
	push(@seqlens, $seq_len);
	foreach my $p (@pat1, @pat2) {
		my $c = 0;
		$c++ while $seq =~ /$p/g;
		$counts{$p} += $c;
	}
}

sub parse_enz {
	my $r = shift;
	my @pat;
	foreach my $e (split(',', $r)) {
		$e = checkenzyme($e);
		my $p = getpattern($e);
		push(@pat, $p);
		$enz{$p} = $e;
	}
	return @pat;
}

sub stdev {
	my $n   = shift;
	return 0 if $n == 1;
	my $s   = shift;
	my $sq = shift;
	my $stdev = ($sq - $s * $s / $n) / ($n - 1);
	return sqrt($stdev);
}

my $fh;

if ($file eq '-') {
	open($fh, '-') || die "Couldn't read from STDIN\n";
}
else {
	if ($file =~ /\.gz$/) {
		open($fh, "gzip -cd $file |") || die "Couldn't open gzipped file '$file'\n";
	}
	else {
		open($fh, "<",         $file) || die "Couldn't open file '$file'\n";
	}
}

if (-f $r1) {
	die("-r1 is a file: this functionality is not implemented yet\n");
}
else {
	@pat1 = parse_enz($r1);
}

if ($double) {
	if (-f $r2) {
		die("-r2 is a file: this functionality is not implemented yet\n");
	}
	else {
		@pat2 = parse_enz($r2);
	}
}

my $seq = "";
my $cut;
my $name;
my $line;

if ($t eq "lens") {
	printf "%s\t%s\t%s\n", "R1", "R2", "Length";
}

while ($line = <$fh>) {
	($fq && ($.-1)%4 > 1) && next;
	if ($fq && $.%4 == 1 || $line =~ /^>/) {
		$cat && next;
		$tsubs{$t}->($seq);
		$seq = "";
		chomp $line;
		$name = $line;
	}
	else {
		$seq .= $line;
	}
}
$tsubs{$t}->($seq);

(!$len_total) && die("It seems, the input has incorrect format\n");

my $gcf = $gc_total / $atgc_total;
my $atf = 1 - $gcf;

printf STDERR "Input sequences:\n";
printf STDERR "Total number - %d\n", $seq_n;
printf STDERR "Total length - %d bp\n", $len_total;
printf STDERR "Mean length (SD) - %.1f (%.2f)\n", $len_total / $seq_n, stdev($seq_n, $len_total, $len_total_sq);
printf STDERR "GC%% - %.3f\n", $gcf;

if ($t eq "ldist") {
	printf "%-10s\t%-10s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\n", "R1", "R2", "Mean len", "SD", "N", "N/bp", "N/Mb";
	foreach my $p1 (@pat1) {
		foreach my $p2 (@pat2) {
			printf "%-10s\t%-10s\t", $enz{$p1}, $enz{$p2};
			if (!defined($frag_n{$p1}{$p2})) {
				print "0\t0\n0\n";
				next;
			}
			my $n = $frag_n{$p1}{$p2};
			my $s = $frag_total{$p1}{$p2};
			my $f = $n / $len_total;
			printf "%-8.2f\t%-8.2f\t%-8d\t%f\t%-8.1f\n", $s / $n, stdev($n, $s, $frag_total_sq{$p1}{$p2}), $n, $f, $f * 1000000;
		}
	}
}
elsif ($t eq "dist") {
	printf "%-10s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Enzyme", "Obs occ", "Obs occ/bp", "Obs occ/Mb", "Exp occ", "Exp occ/bp", "Exp occ/Mb";
	while (my ($p, $e) = each %enz) {
		my %f = (A => 0, C => 0, G => 0, T => 0,
			 M => 0, R => 0, W => 0, S => 0, Y => 0, K => 0, V => 0, H => 0, D => 0, B => 0);
		$f{$_}++ foreach split('', $e);
		my $exp = (    $gcf/2) ** ($f{C} + $f{G})
			* (    $atf/2) ** ($f{A} + $f{T})
			* (1 - $gcf/2) ** ($f{H} + $f{D})
			* (1 - $atf/2) ** ($f{B} + $f{V})
			* (    $gcf  ) ** ($f{S}        )
			* (    $atf  ) ** ($f{W}        )
			* (     0.5  ) ** ($f{R} + $f{Y} + $f{K} + $f{M});
		my $exp_raw = $exp * $len_total;
		my $obs_raw   = defined($counts{$p}) ? $counts{$p} : 0;
		my $obs       = $obs_raw / $len_total;
		printf "%-10s\t%-7d\t%-10f\t%-10.1f\t%-7.f\t%-10f\t%-10.1f\n", $enz{$p}, $obs_raw, $obs, $obs * 1000000, $exp_raw, $exp, $exp * 1000000;
	}
}

# 'd <- read.table(file("stdin"),TRUE)',
# 'with(d, Hist(Length, scale="frequency", breaks=200, col="darkgray"))'
# "svg(\"$name.svg\",width=3,height=3)",
