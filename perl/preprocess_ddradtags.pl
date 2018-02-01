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

my $verion = "1.6";

if (!@ARGV or defined($args{-h}) or defined($args{-v})) {
	print STDERR "preprocess_ddradtags.pl version $verion
Use as:
  preprocess_ddradtags.pl [arguments]

  Mandatory arguments:
       -l : left fq-file
       -r : right fq-file
       -b : barcode file
              Each line in this file is expected to represent the barcode sequence (with respective Ins concatenated),
              the presence-of-DBR-flag (1/0) and DBR offset (lengths of P7 Ins).
              The first line may optionally contain the actual length of the barcode itself. 
       -i : DBR sequence
              The standard [IUPAC codes](http://www.bioinformatics.org/sms/iupac.html) are expected to denote variable positions

  Optional arguments:
       -m : matching length [10]
       -x : allowed mismatches in DBR [1]
       -y : allowed mismatches in barcode [1]
       -t : tolerated barcode shift [1]
      -lt : trim the left  sequences to this length [90]
      -rt : trim the right sequences to this length [90]
   -lrest : if provided - check for incomplete restriction for the left  restrictase []
   -rrest : if provided - check for incomplete restriction for the right restrictase []
      -gz : gzip output files (y/n) [y]
";
	exit();
}

my ($bfile, $lfile, $rfile, $orig_dbr, $dbr, $icut);
my $m = 10;
my $x = 1;
my $y = 1;
my $t = 1;
my $lt = 90;
my $rt = 90;
my $gz = 'y';

my %skipped;
my %dbrs;

my $lrest = 0;
my $rrest = 0;


## START ARGS ##

if (defined($args{-b})) {
	$bfile = $args{-b};
	delete($args{-b});
}
else {
	die("barcode file not specified\n");
}
if (defined($args{-i})) {
	$orig_dbr = $args{-i};
	die("The DBR sequence '$orig_dbr' seems to be incorrect\n") if $orig_dbr !~ /^[ACGTRYMKBDHVN]+$/;
	delete($args{-i});
}
else {
	die("DBR not specified\n");
}
if (defined($args{-l}) and defined($args{-r})) {
	$lfile = $args{-l};
	$rfile = $args{-r};
	delete($args{-l});
	delete($args{-r});
}
else {
	die("two fastq-files are expected\n");
}
if (defined($args{-m})) {
	$m = $args{-m};
	delete($args{-m});
}
else {
	print STDERR "Matching length not specified: taking the default of $m (x2)\n";
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
if (defined($args{-lt})) {
	$lt = $args{-lt};
	delete($args{-lt});
}
else {
	print STDERR "Left trim size not specified: taking the default of $lt\n";
}
if (defined($args{-rt})) {
	$rt = $args{-rt};
	delete($args{-rt});
}
else {
	print STDERR "Right trim size not specified: taking the default of $rt\n";
}

if (defined($args{-lrest})) {
	$lrest = $args{-lrest};
	die("Invalid restriction site '$lrest'\n") if $lrest !~ /^[ATGC]+$/;
	delete($args{-lrest});
}
if (defined($args{-rrest})) {
	$rrest = $args{-rrest};
	die("Invalid restriction site '$rrest'\n") if $rrest !~ /^[ATGC]+$/;
	delete($args{-rrest});
}
if (defined($args{-gz})) {
	$gz = $args{-gz};
	die("Invalid choice '$gz' for -gz\n") if $gz ne 'y' && $gz ne 'n';
	delete($args{-gz});
}
foreach my $key (keys %args) {
	print STDERR "$key: unknown parameter\n";
	exit();
}

print STDERR "The script was launched as: -l $lfile -r $rfile -b $bfile -i $orig_dbr -m $m -x $x -y $y -t $t -lt $lt -rt $rt -lrest $lrest -lrest $rrest -gz $gz\n";

my $m2 = $m * 2;

## END ARGS ##

my $start_time = time();

my (@barcodes, @barcode_lens);
my %barcode_names;
my $il;
my $bl;
my (%uniq, %pieces12, %pieces13, %pieces23);
my ($lline, $l, $lh, $ls, $lq, $lcut, $bl_chosen, $b_chosen);
my ($rline, $r, $rh, $rs, $rq, $rcut, $i_offset,  $i_chosen);
my ($subseq, $seq, $subseq1, $subseq2, $subseq3);
my (@no_dbrs, @i_shifts);
my ($no_dbr, $is_dup);
my $i_shift;
my %loci;
my ($locus, $mask);
my @i_rescue;

read_barcodes();
read_dbr();

my $lin = openin($lfile);
my $rin = openin($rfile);

my (%LABB, %RABB, %LWBI, %RWBI, %LNOI, %RNOI, $LWOB, $RWOB, %LWOI, %RWOI, %LRUR, %RRUR, %LLUR, %RLUR);

foreach my $k (0..$#barcodes) {
	$b = $barcodes[$k];
	my $lfix = defined($barcode_names{$b}) ? $barcode_names{$b}.".l" : "$lfile.$b";
	my $rfix = defined($barcode_names{$b}) ? $barcode_names{$b}.".r" : "$rfile.$b";
	$LABB{$b} = openout("$lfix.abb");
	$RABB{$b} = openout("$rfix.abb");
	
	if ($no_dbrs[$k]) {
		$LNOI{$b} = openout("$lfix.noi");
		$RNOI{$b} = openout("$rfix.noi");
	}
	else {
		$LWBI{$b} = openout("$lfix.wbi");
		$RWBI{$b} = openout("$rfix.wbi");
		$LWOI{$b} = openout("$lfix.woi");
		$RWOI{$b} = openout("$rfix.woi");
	}
	if ($rrest) {
		$LRUR{$b} = openout("$rfix.rur");
		$RRUR{$b} = openout("$lfix.rur");
	}
	if ($lrest) {
		$LLUR{$b} = openout("$rfix.lur");
		$RLUR{$b} = openout("$lfix.lur");
	}
}

$LWOB = openout("$lfile.wob");
$RWOB = openout("$rfile.wob");

my $wbic = 0;
my $wobc = 0;
my $woic = 0;
my $noic = 0;
my $dupc = 0;
my $abbc = 0;
my $lurc = 0;
my $rurc = 0;
my $total = 0;
my ($lplus, $rplus);

my $time;
while (!eof($lin) && !eof($rin)) {
	report_time() if ++$total % 1000000 == 0;
	$lh = <$lin>;
	$rh = <$rin>;
	$ls = <$lin>;
	$rs = <$rin>;
	$lplus = <$lin>;
	chomp $lplus;
	die "Left file out of phase at line $.: $lplus\n" if $lplus ne "+";
	$rplus = <$rin>;
	chomp $rplus;
	die "Right file out of phase at line $.: $rplus\n" if $rplus ne "+";
	$lq = <$lin>;
	$rq = <$rin>;
	chomp($ls, $rs, $lh, $rh, $lq, $rq);
	$is_dup = 0;
	$no_dbr = 0;
	$bl_chosen = find_similar_barcode($ls);
	if ($bl_chosen >= 0) {
		$ls = substr($ls, $bl_chosen, $lt);
		$lq = substr($lq, $bl_chosen, $lt);
		if ($no_dbr) {
			$rs = substr($rs, $i_shift, $rt);
			$rq = substr($rq, $i_shift, $rt);
			is_duplicate($rs);
		}
		else {
			$i_offset = find_dbr($rs);
			if ($i_offset >= 0) {
				$rs = substr($rs, $i_offset, $rt);
				$rq = substr($rq, $i_offset, $rt);
				$is_dup = is_duplicate($rs);
				$rh .= " ".$i_chosen if !$is_dup;
			}
		}
	}
	if ($is_dup) {
		$skipped{$b_chosen}++;
		$dupc++;
	}
	else {
		$l = sprintf("%s\n%s\n+\n%s\n", $lh, $ls, $lq);
		$r = sprintf("%s\n%s\n+\n%s\n", $rh, $rs, $rq);
		if ($bl_chosen < 0) {
			print { $LWOB } $l;
			print { $RWOB } $r;
			$wobc++;
		}
		elsif ($no_dbr) {
			print { $LNOI{$b_chosen} } $l;
			print { $RNOI{$b_chosen} } $r;
			$noic++;
		}
		elsif ($i_offset < 0) {
			print { $LWOI{$b_chosen} } $l;
			print { $RWOI{$b_chosen} } $r;
			$woic++;
		}
		elsif ($lrest && index($ls, $lrest) > -1) {
			print { $LLUR{$b_chosen} } $l;
			print { $RLUR{$b_chosen} } $r;
			$lurc++;
		}
		elsif ($rrest && index($rs, $rrest) > -1) {
			print { $LRUR{$b_chosen} } $l;
			print { $RRUR{$b_chosen} } $r;
			$rurc++;
		}
		elsif (length($ls) < $lt || length($rs) < $rt) {
			print { $LABB{$b_chosen} } $l;
			print { $RABB{$b_chosen} } $r;
			$abbc++;
		}
		else {
			print { $LWBI{$b_chosen} } $l;
			print { $RWBI{$b_chosen} } $r;
			$wbic++;
		}
	}
}

close($lin);
close($rin);

foreach my $k (0..$#barcodes) {
	$b = $barcodes[$k];
	close($LABB{$b});
	close($RABB{$b});

	if (exists($LWBI{$b})) {
		close($LWBI{$b});
		close($RWBI{$b});
	}
	if (exists($LNOI{$b})) {
		close($LNOI{$b});
		close($RNOI{$b});
	}
	if (exists($LWOI{$b})) {
		close($LWOI{$b});
		close($RWOI{$b});
	}
	if (exists($LLUR{$b})) {
		close($LLUR{$b});
		close($RLUR{$b});
	}
	if (exists($LRUR{$b})) {
		close($LRUR{$b});
		close($RRUR{$b});
	}
}

close($LWOB);
close($RWOB);

print STDERR "\nJob finished in ";

$time = time() - $start_time;
if ($time < 3600) {
	printf STDERR "%d sec\n", $time;
}
else {
	printf STDERR "%d min\n", $time / 60;
}

## print out the summaries ##

printf STDERR "Total number of input reads: %d.\nAmong them: %d abb, %d wbi, %d noi, %d woi, %d wob, %d dup, %d lur, %d rur\n\n", $total, $abbc, $wbic, $noic, $woic, $wobc, $dupc, $lurc, $rurc;

print STDERR "The counts of skipped reads are:\n";

while (my ($b, $c) = each(%skipped)) {
	print STDERR "$b: $c\n";
}

my %RDUP;
foreach my $k (0..$#barcodes) {
 	$b = $barcodes[$k];
	my $rfix = defined($barcode_names{$b}) ? $barcode_names{$b}.".r" : "$rfile.$b";
 	open($RDUP{$b}, ">", "$rfix.dup") or die($!);
}

my $prev_b = "";
my $prev_dup = "";
my %aread_num;
my %uread_num;
my %loci_num;
my $sum = 0;

while (my ($b, $locdbrs) = each (%uniq)) {
	foreach my $locus (0..$#{$locdbrs}) {
		my $udbrs = ${$locdbrs}[$locus];
		$sum = 0;
		while (my ($i, $count) = each (%{$udbrs})) {
			$sum += $count;
			++$uread_num{$b};
		}
		$aread_num{$b} += $sum;
		++$loci_num{$b};
		while (my ($i, $count) = each (%{$udbrs})) {
			printf { $RDUP{$b} } "%s %d %s %s %d %d\n", $b, $locus, $loci{$b}[$locus], $i, $count, $sum;
		}
	}
}
foreach my $k (0..$#barcodes) {
 	$b = $barcodes[$k];
 	close($RDUP{$b});
}

print STDERR "\nReads/loci:\n\n";
print STDERR "barcode all_reads uniq_reads loci\n";
foreach my $b (keys %uniq) {
	printf STDERR "%s %d %d %d\n", $b, defined($aread_num{$b}) ? $aread_num{$b} : 0, defined($uread_num{$b}) ? $uread_num{$b} : 0, defined($loci_num{$b}) ? $loci_num{$b} : 0;
}

open(INT, ">", "$rfile.int") or die($!);
while (my ($i, $count) = each (%dbrs)) {
	printf INT "%s %d\n", $i, $count;
}
print STDERR "\nDone\n";
close(INT);


## read the barcodes from the files ##
sub read_barcodes {
	open(BARCODES, "<", $bfile) or die("Couldn't open barcode file\n");
	my (@all_barcodes, $len, @parts);
	my (%no_dbrs_hash, %i_shifts_hash);
	while (my $line = <BARCODES>) {
		chomp $line;
		if ($line =~ /^\d+$/) {
 			$bl = $line;
			next;
		}
		die("The barcode file has incorrect format. At least three fields expected: the barcode, the presence-of-DBR-flag, P7 offset\n") if ($line !~ /^[ATGC]+\s+(1|0)\s+\d(\s+\w+)?$/);
		@parts = split ' ', $line;
		$b = $parts[0];
		push(@all_barcodes, $b);
		$no_dbrs_hash{$b} = !$parts[1];
		$i_shifts_hash{$b} = $parts[2];
		$len = length($b);
		if (!defined($bl) or $bl > $len) {
			$bl = $len;
		}
		if (defined($parts[3])) {
			$barcode_names{$b} = $parts[3];
		}
	}
	close(BARCODES);
	my $subb;
	foreach my $b (@all_barcodes) {
		$subb = substr($b, 0, $bl);
		if (defined($barcode_names{$b}) && $subb ne $b) {
			$barcode_names{$subb} = $barcode_names{$b};
			delete($barcode_names{$b});
		}
		push(@barcodes, $subb);
		$barcode_lens[$#barcodes] = length($b);
		$no_dbrs[$#barcodes] = $no_dbrs_hash{$b};
		$i_shifts[$#barcodes] = $i_shifts_hash{$b};
		$skipped{$subb} = 0;
		$uniq{$subb} = ();
	}
}

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
	for my $pos (0..$#i_bases) {
		my $i_base = $i_bases[$pos];
		next if $i_base eq 'N';
		$i_base =~ tr/ACGTRYMKBDHV/BDHVYRKMACGT/;
		push(@i_rescue, prepare_dbr(substr($orig_dbr, 0, $pos).$i_base.substr($orig_dbr, $pos + 1)) );
	}
}

## identify barcode in a read ##
sub find_similar_barcode {
	my $seq = shift;
	my $b;
	$b_chosen = "";
	$i_shift = -1;
	foreach my $j (0..$t) {
		$subseq = substr($seq, $j, $bl);
		foreach my $k (0..$#barcodes) {
			$b = $barcodes[$k];
			$mask = $b ^ $subseq;
			$mask =~ tr/\0//d;
			if (length($mask) <= $y) {
				$b_chosen = $b;
				$no_dbr = $no_dbrs[$k];
				$i_shift = $i_shifts[$k];
				return $barcode_lens[$k] + $j;
			}
		}
	}
	return -1;
}

## identify dbr ##
sub find_dbr {
	my $seq = shift;
	$subseq = substr($seq, $i_shift, $il);
	foreach my $j (0..$t) {
		my $subseqj = $j ? substr($seq, $i_shift + $j, $il) : $subseq;
		if ($subseqj =~ /^($dbr)/) {
			$i_chosen = $icut < 0 ? substr($1, 0, $icut) : $1;
			$dbrs{$i_chosen}++;
			return $il + $i_shift + $j;
		}
	}
	foreach my $rescue (@i_rescue) {
		if ($subseq =~ /^($rescue)/) {
			$i_chosen = $icut < 0 ? substr($1, 0, $icut) : $1;
			$dbrs{$i_chosen}++;
			return $il + $i_shift;
		}
	}
	return -1;
}

## detect if this is a duplicate ##
sub is_duplicate {
	$seq = shift;
	$subseq1 = substr($seq, 0,  $m);
	$subseq2 = substr($seq, $m, $m);
	$locus = $pieces12{$b_chosen}{$subseq1}{$subseq2};
	if (!defined($locus)) {
		$subseq3 = substr($seq, $m2, $m);
		$locus = $pieces13{$b_chosen}{$subseq1}{$subseq3};
		if (!defined($locus)) {
			$locus = $pieces23{$b_chosen}{$subseq2}{$subseq3};
			if (!defined($locus)) {
				$locus = push(@{$loci{$b_chosen}}, $subseq1.$subseq2.$subseq3) - 1;
				$pieces23{$b_chosen}{$subseq2}{$subseq3} = $locus;
			}
			$pieces13{$b_chosen}{$subseq1}{$subseq3} = $locus;
		}
		$pieces12{$b_chosen}{$subseq1}{$subseq2} = $locus;
	}
	if ($no_dbr) {
		return (++$uniq{$b_chosen}[$locus]{"no_dbr"} > 1);
	}
	if ($x) {
		foreach my $dbr (keys %{$uniq{$b_chosen}[$locus]}) {
			$mask = $i_chosen ^ $dbr;
			$mask =~ tr/\0//d;
			if (length($mask) <= $x) {
				++$uniq{$b_chosen}[$locus]{$dbr};
				return 1;
			}
		}
	}
	return (++$uniq{$b_chosen}[$locus]{$i_chosen} > 1);
}

sub openin {
	my $fname = shift;
	my $fp;
	my $pipe = ($fname =~ /\.gz$/ ? "gzip -cd $fname |" : "< $fname");
	open($fp, $pipe) or die("$!\n");
	return $fp;
}

sub openout {
	my $fname = shift;
	my $fp;
	my $pipe = ($gz eq 'y' ? "| gzip -c > $fname.gz" : "> $fname");
	open($fp, $pipe) or die("$!\n");
	return $fp;
}

sub report_time {
	$time = time() - $start_time;
	printf STDERR "Read #%d, %d sec, %d reads/sec\n", $total, $time, $total / $time if $time > 0;
}
