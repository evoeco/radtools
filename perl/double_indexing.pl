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
	print STDERR "double_indexing.pl
Use as:
  (perl) double_indexing.pl [arguments]

  Mandatory arguments:
       -l : left fq-files  (comma-separated)
       -r : right fq-files (comma-separated)
 -indices : indices        (comma-separated)
       -i : DBR sequence
              The standard [IUPAC codes](http://www.bioinformatics.org/sms/iupac.html) are expected to denote variable positions

  Optional arguments:
       -x : allowed mismatches in DBR     [1]
       -y : allowed mismatches in barcode [1]
       -t : tolerated barcode shift       [4]
      -gz : gzip output files (y/n)       [y]
";
	exit();
}

my ($orig_dbr, $dbr, $icut);
my (@lfiles, @rfiles, @indices);
my ($indexread, $inlineindex);
my $m = 10;
my $s = 0;
my $x = 1;
my $y = 1;
my $t = 4;
my $lt = 90;
my $rt = 90;
my $gz = 'y';


## START ARGS ##

if (defined($args{-l}) and defined($args{-r})) {
	@lfiles = split(',', $args{-l});
	@rfiles = split(',', $args{-r});
	delete($args{-l});
	delete($args{-r});
}
else {
	die("at least two fastq-files are expected\n");
}
if (defined($args{-indices})) {
	@indices = split(',', $args{-indices});
	delete($args{-indices});
}
else {
	die("Indices not specified\n");
}
if (defined($args{-i})) {
	$orig_dbr = $args{-i};
	delete($args{-i});
}
else {
	die("DBR not specified\n");
}
if (defined($args{-x})) {
	$x = $args{-x};
	delete($args{-x});
}
else {
	print STDERR "Number of mismatches in DBRs not specified: taking the default of $x\n";
}
if (defined($args{-y})) {
	$y = $args{-y};
	delete($args{-y});
}
else {
	print STDERR "Number of mismatches in barcodes not specified: taking the default of $y\n";
}
if (defined($args{-t})) {
	$t = $args{-t};
	delete($args{-t});
}
else {
	print STDERR "Tolerated barcode shift not specified: taking the default of $t\n";
}
if (defined($args{-gz})) {
	$gz = $args{-gz};
	die("Invalid choice '$gz' for -gz\n") if $gz ne 'y' && $gz ne 'n';
	delete($args{-gz});
}

foreach my $key (keys %args) {
	die("$key: unknown parameter\n");
}

if ($#lfiles != $#rfiles) {
	die("Numbers of left and right files are different\n");
}
if ($#lfiles != $#indices) {
	die("Numbers of files and indices are different\n");
}


## END ARGS ##

my $start_time = time();

my %barcode_names;
my $il;
my $bl;
my (%uniq, %pieces12, %pieces13, %pieces23);
my ($lline, $l, $lh, $ls, $lq, $lcut, $bl_chosen, $b_chosen);
my ($rline, $r, $rh, $rs, $rq, $rcut, $i_offset,  $i_chosen);
my ($subseq, $seq, $subseq1, $subseq2, $subseq3);
my (@nodbrs, @i_shifts);
my $nodbr;
my $i_shift;
my %loci;
my ($locus, $mask);
my @i_rescue;

## BEGIN SUBS ##

## make a regex pattern out of DBR ##
sub prepare_dbr {
	my $pattern = shift;
	$pattern =~ s/R/[AG]/g;
	$pattern =~ s/Y/[CT]/g;
	$pattern =~ s/S/[GC]/g;
	$pattern =~ s/W/[AT]/g;
	$pattern =~ s/K/[GT]/g;
	$pattern =~ s/M/[AC]/g;
	$pattern =~ s/B/[CGT]/g;
	$pattern =~ s/D/[AGT]/g;
	$pattern =~ s/H/[ACT]/g;
	$pattern =~ s/V/[ACG]/g;
	$pattern =~ s/N/[ACGTNX]/g;
	return $pattern;
}

## get DBR and create pattern out of it ##
sub read_dbr {
	$il = length($orig_dbr);
	$dbr = prepare_dbr($orig_dbr);
	my $var_dbr = $orig_dbr;
	$var_dbr =~ s/[ATGC]+$//g;
	$icut = length($var_dbr) - length($orig_dbr);
	my @i_bases = split //, $orig_dbr;
	my $pos = -1;
	foreach my $i_base (@i_bases) {
		$pos++;
		next if $i_base eq 'N';
		$i_base =~ tr/ACGTRYMKBDHVYRKM/BDHVYRKMACGTRYMK/;
		push(@i_rescue, prepare_dbr(substr($orig_dbr, 0, $pos).$i_base.substr($orig_dbr, $pos + 1)) );
	}
}

## identify dbr ##
sub find_dbr {
	my $seq = shift;
	foreach my $j (0..$t) {
		my $subseq = substr($seq, $j, $il);
		if ($subseq =~ /^($dbr)/) {
			return $il + $j;
		}
	}
	foreach my $j (0..$t) {
		my $subseq = substr($seq, $j, $il);
		foreach my $rescue (@i_rescue) {
			if ($subseq =~ /^($rescue)/) {
				return $il + $j;
			}
		}
	}
	return -1;
}

## find and check double index ##
sub check_inline {
	my $seq    = shift;
	my $offset = shift;
	foreach my $index (@indices) {
		my $subseq = substr($seq, $offset, length($index));
		$mask = $index ^ $subseq;
		$mask =~ tr/\0//d;
		if (length($mask) <= $y) {
			return $index;
		}
	}
	return 0;
}

sub openin {
	my $fname = shift;
	my $fp;
	if ($fname =~ /\.gz$/) {
		open($fp, "gzip -cd $fname |") or die($!);
	}
	else {
		open($fp, "<",  $fname) or die($!);
	}
	return $fp;
}

sub openout {
	my $fname = shift;
	my $gz    = shift;
	my $fp;
	if ($gz eq 'y') {
		open($fp, "| gzip -c > $fname.gz") or die($!);
	}
	else {
		open($fp, ">", $fname) or die($!);
	}
	return $fp;
}

## END SUBS ##

read_dbr();

## BEGIN FILES ##

my (%LOUT, %ROUT);
my (%total, %correct, %incorrect, %noindex, %nodbr);


for (my $f = 0; $f <= $#indices; $f++) {
	my $index  = $indices[$f];
	$LOUT{$index} = openout($lfiles[$f]."out.fq", $gz);
	$ROUT{$index} = openout($rfiles[$f]."out.fq", $gz);
	$total{$index}   = 0;
	$correct{$index} = 0;
	$incorrect{$index} = 0;
	$noindex{$index} = 0;
	$nodbr{$index}   = 0;
}

my $time;

for (my $f = 0; $f <= $#indices; $f++) {
	my $lin = openin($lfiles[$f]);
	my $rin = openin($rfiles[$f]);
	$indexread  = $indices[$f];
	my $mod;
	my $eof;

	$lline = <$lin>;
	$rline = <$rin>;

	die("Input files have wrong format\n") if (substr($lline, 0, 1) ne '@' or substr($rline, 0, 1) ne '@');
	$lh = substr($lline, 1);
	$rh = substr($rline, 1);
	chomp $lh;
	chomp $rh;

	## fastq parsing: get the lines simultaneously for the left- and the right-read files ##

	while (!$eof) {
		$lline = <$lin>;
		$rline = <$rin>;
		chomp $lline;
		chomp $rline;
		$mod = $. % 4;
		## if it is header ##
		$eof = (eof($lin) or eof($rin) or !$lline or !$rline);
		## if it is the quality line ##
		if ($mod == 0) {
			$lq = $lline;
			$rq = $rline;
			$rq = substr($rq, 0, $i_offset).substr($rq, $i_offset + length($inlineindex)) if $inlineindex;
		}
		if ($mod == 1 || $eof) {
			$l = sprintf("@%s\n%s\n+\n%s\n", $lh, $ls, $lq);
			$r = sprintf("@%s\n%s\n+\n%s\n", $rh, $rs, $rq);
			$total{$indexread}++;
			if ($inlineindex) {
				$correct{$indexread}   += ($inlineindex eq $indexread);
				$incorrect{$indexread} += ($inlineindex ne $indexread);
				print { $LOUT{$inlineindex} } $l;
				print { $ROUT{$inlineindex} } $r;
			}
			else {
				$nodbr{$indexread} += ($i_offset < 0);
				$noindex{$indexread} += ($i_offset >= 0);
				print { $LOUT{$indexread} } $l;
				print { $ROUT{$indexread} } $r;
			}
			$lh = substr($lline, 1);
			$rh = substr($rline, 1);
		}
		## if it is the "+"-line ##
		elsif ($mod == 3) {
			die("Files are not in phase: $lline\n") if ($lline ne "+" or $rline ne "+");
		}
		## if it is the sequence line ##
		elsif ($mod == 2) {
			$ls = $lline;
			$rs = $rline;
			$i_offset = find_dbr($rs);
			$inlineindex = 0;
			next if $i_offset < 0;
			$inlineindex = check_inline($rs, $i_offset);
			next if !$inlineindex;
			$rs = substr($rs, 0, $i_offset).substr($rs, $i_offset + length($inlineindex));
		}
	}
	close($lin);
	close($rin);
}

for (my $f = 0; $f <= $#indices; $f++) {
	$indexread  = $indices[$f];
	close($LOUT{$indexread});
	close($ROUT{$indexread});
}

print STDERR "\nJob finished in ";

$time = time() - $start_time;
if ($time < 3600) {
	printf STDERR "%d sec\n", $time;
}
else {
	printf STDERR "%d min\n", $time / 60;
}

## print out the summaries ##

printf STDERR "%s\t%s\t%s\t%s\t%s\t%s\n", "Index", "Total", "Correct", "Incorrect", "No DBR", "No index";

for (my $f = 0; $f <= $#indices; $f++) {
	my $index  = $indices[$f];
	printf STDERR "%s\t%d\t%d\t%d\t%d\t%d\n", $index, $total{$index}, $correct{$index}, $incorrect{$index}, $nodbr{$index}, $noindex{$index};
}

