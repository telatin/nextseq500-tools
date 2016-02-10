#!/usr/bin/env perl
 
# Andrea Telatin 2016, Computer Laboratory - University of Cambridge

# todo:
# - Allow single end (SE) trimming
# - Implement minimum stretch (trailing G tail) lenght

# issues:

use Time::HiRes qw( time );
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Term::ANSIColor  qw(:constants);
use strict;


# defaults
my $output = "filtered";		# Output basename, will add _R1.fastq / _R2.fastq
my $minInitialLength = 50;		# Discard raw reads shorter than
my $minFinalLength = 50;		# Discard trimmed reads resulting shorter thane 
my $mismatches = 1;				# Mismatches allowed for trailing G-tail.
my $minstretchlen = 10;			# Minimum trailing G length (NOT IMPLEMENTED YET)

sub shortHelp {
	print STDERR "	   
	---------------------------------------------------------------------------------
	   FASTQ QUALITY FILTER FOR ILLUMINA PAIRED END
	   Filters a set of two paired end FASTQ file removing the final G stretch
	   (dark) signal, defined as a [NG] stretch
	---------------------------------------------------------------------------------
	   Usage: [parameters] File1 [File2]
	
		-o, --output       STRING    Output basename (default: $output)
		-m, --mismatches   INT       Maxmimum mismatches (default: $mismatches)		
		-q, --minavgqual   FLOAT     Minimum average quality
		-5, --trimStart    INT       Trim bases at the begin
		-3, --trimEnd      INT       Trim bases at end	
		-n, --nomachine              Remove machine name from Illumina reads
		-f, --finallen     INT       Minimum length (after trimming, default: $minInitialLength)
		--minlen           INT       Minimum initial length (default: $minFinalLength)

";
	exit;
}


my $noname;
my $help;
my $debug;
my $first_file;
my $second_file = '';
my $minavgqual;
my $trimStart = 0;
my $trimEnd = 0;

my $getinfoResult = GetOptions(
   'h|help'          => \$help,
   'd|debug'         => \$debug,
   'o|output=s'      => \$output,
   'm|mismatches=i'  => \$mismatches,
   'q|minavgqual=f'  => \$minavgqual,
   '5|trimStart=i'   => \$trimStart,
   '3|trimEnd=i'     => \$trimEnd,
   'n|nomachine'     => \$noname,
   'minlen=i'        => \$minInitialLength,
   'f|finallength=i'  => \$minFinalLength
);

pod2usage({-exitval => 0, -verbose => 2}) if ($help);

($first_file, $second_file) = @ARGV;
if ($first_file and !$second_file) {
	$second_file = $first_file;
	$second_file  =~s/_R1/_R2/;
}


&shortHelp() if (!$first_file or ! -e "$second_file");

my $R1;
my $R2;
my $O1;
my $O2;
my $global_trimmed = 0;
my $global_input_bases = 0;
my $global_trimmed_bases = 0;
my $global_max_trimmed = 0;
open($R1, "<", "$first_file") || die " FATAL ERROR:\n Unable to read $first_file.\n";
open($R2, "<", "$second_file") || die " FATAL ERROR:\n Unable to read $second_file.\n";

open($O1, ">", "$output\_R1.fastq") || die " FATAL ERROR:\n Unable to write to $output\_R1.fastq\n";
open($O2, ">", "$output\_R2.fastq") || die " FATAL ERROR:\n Unable to write to $output\_R2.fastq\n";

my $countSeq = 0;
my $countPrinted = 0;
my $total_qual1 = 0;
my $total_qual2 = 0;
my $printed_qual1 = 0;
my $printed_qual2 = 0;

sub trimBlankEnd {
	my $seq = shift;
	my $copy = $seq;
	my $tlen = 0;
	for (my $i = 0; $i < $mismatches; $i++) {
		# Detect a [ACT]nnngnggngngngng stretch and change the first non-GN base
		if ($seq=~/[^NG]([NG]+)$/i) {
			my $replace = 'N' . $1;
			$seq=~s/[^NG]([NG]+)$/$replace/;
		}
	}
	
	# Now remove the "altered" stretch
	$seq =~s/([^NG][NG]+)$//;
	$global_trimmed++ if ($seq ne $copy);
	$tlen = length($copy) - length($seq);
	$global_trimmed_bases += $tlen;
	$global_input_bases   += length($copy);
	$global_max_trimmed = $tlen if ($tlen > $global_max_trimmed);
	return ($seq, $tlen);

}

my %trimmed_R1;
my %trimmed_R2;
while (my $name1 = readline($R1)) {
	my $seq1 = readline($R1);
	my $sep1 = readline($R1);
	my $qua1 = readline($R1);
	my $name2= readline($R2);
	my $seq2 = readline($R2);
	my $sep2 = readline($R2);
	my $qua2 = readline($R2);
	$countSeq++;
	
	next if ($minInitialLength and ( length($seq1) < $minInitialLength or length($seq2) < $minInitialLength));
	if ($noname) {
		$name1 = &noName($name1);
		$name2 = &noName($name2);
	}
	
	# Trim G spots
	my $trim1 = 0;
	my $trim2 = 0;
	
	($seq1, $trim1) = trimBlankEnd($seq1);
	$qua1 = substr($qua1, 0, length($seq1)-1);
	
	($seq2, $trim2) = trimBlankEnd($seq2);
	$qua2 = substr($qua2, 0, length($seq2)-1);
	
	$trimmed_R1{$trim1}++;	
	$trimmed_R2{$trim2}++;
	# Trim end
	if ($trimEnd) {
	    $seq1 = substr($seq1, 0, -1*$trimEnd);
    	$qua1 = substr($qua1, 0, -1*$trimEnd);
	    $seq2 = substr($seq2, 0, -1*$trimEnd);
    	$qua2 = substr($qua2, 0, -1*$trimEnd);
	}	
	
	# Trim start
	if ($trimStart) {
		print STDERR " [$countSeq] Warning: messing things up here \n SEQ:$seq1\n Trim: $trimStart\n" if (length($seq1) <= $trimStart);
		$seq1 = substr($seq1,  $trimStart);
        $qua1 = substr($qua1,  $trimStart);
   		$seq2 = substr($seq2,  $trimStart);
        $qua2 = substr($qua2,  $trimStart);   
	}
	
	my $avg1 = string2qual($qua1);
	my $avg2 = string2qual($qua2);
	$total_qual1 += $avg1;
	$total_qual2 += $avg2;
	
	# Quality filter
	if ($avg1 < $minavgqual or $avg2 < $minavgqual) {
		next;
	}
	$countPrinted++;
	$printed_qual1 += $avg1;
	$printed_qual2 += $avg2;	
	if ($countSeq % 10000 == 0) {
		my $r = sprintf("%.2f", 100 * $countPrinted / $countSeq);
		
		my $a1 = sprintf("%.2f",  $total_qual1 / $countSeq);
		my $a2 = sprintf("%.2f",  $total_qual2 / $countSeq);
		my $b1 = sprintf("%.2f",  $printed_qual1 / $countPrinted);
		my $b2 = sprintf("%.2f",  $printed_qual2 / $countPrinted);
		print STDERR " $countSeq parsed, $global_trimmed trimmed. $r% printed. R1 $a1 -> $b1; R2 $a2 -> $b2\r";
	}
	
	
	if ($minFinalLength and (length($seq1) < $minFinalLength or length($seq2) < $minFinalLength)) {
		next;
	}
	print $O1 &fastq($name1, $seq1, $sep1, $qua1);
	print $O2 &fastq($name2, $seq2, $sep2, $qua2);
	
}
my $r = sprintf("%.2f", 100 * $countPrinted / $countSeq);

my $a1 = sprintf("%.2f", 100 * $total_qual1 / $countSeq);
my $a2 = sprintf("%.2f", 100 * $total_qual2 / $countSeq);
my $b1 = sprintf("%.2f", 100 * $printed_qual1 / $countPrinted);
my $b2 = sprintf("%.2f", 100 * $printed_qual2 / $countPrinted);

print STDERR  "

Parsed sequences     $countSeq
Printed sequences    $countPrinted
Trimmed sequences    $global_trimmed (R1 or R2)
Total bases          $global_input_bases
Total trimmed bases  $global_trimmed_bases
Ratio                $r%
R1 Average           $a1
R1 Average filtered  $b1
R2 Average           $a2
R2 Average filtered  $b2

";
for (my $i = 0; $i <= $global_max_trimmed; $i++) {
	print STDERR "Trimmed $i bases\t$trimmed_R1{$i}\t$trimmed_R2{$i}\n";
}


sub noName {
	my $n  = shift;
	#M02007:38:000000000-
	$n =~s/\@(.*?)-/\@/;
	return $n;
}

sub fastq {
	my ($name, $sequence, $separator, $quality) = (@_);
	chomp($name);
	chomp($sequence);
	chomp($separator);
	chomp($quality);
	die "FASTQ: N:$name\nS:$sequence\n+:$separator\nQ:$quality\n" unless ($quality);
	die "FASTQ not starting with \@: N:$name\nS:$sequence\n+:$separator\nQ:$quality\n" unless (substr($name, 0, 1) eq '@');
	return "$name\n$sequence\n$separator\n$quality\n";
}
sub string2qual {
    my $string = shift;
    my $sum;
    for (my $i=0; $i<length($string); $i++) {
        my $q = substr($string, $i, 1);
        $sum += ord($q) - 33;
    }
    return $sum/length($string) if (length($string));
}


__END__
 
=head1 NAME
 
B<filter_NextSeq500_fastq_PE.pl> - Filter paired end Illumina NextSeq500 reads
removing the final stretch of G (dark) in case of short libraries.

=head1 SYNOPSIS
 
filter_NextSeq500_fastq_PE.pl [-m 1 -o <Basename>]  <FirstPE_R1.fastq> [<SecondPE_R2.fastq>]

The name of the second PE file will be inferred replacing _R1 with _R2, if possible.
 
=head1 DESCRIPTION
 
. 
 
=head1 PARAMETERS

=over 12

=item B<-o, --output> BASENAME (default: "filtered")

Basename for the output files (the program will add _R1.fastq and _R2.fastq)

=item B<-m, --mismatches> INT

Number of mismatches allowed in the ending stretch of [N/G]s.

=item B<-3, --trimStart> INT

Remove the first N bases to all sequences

=item B<-5, --trimEnd> INT

Remove the last N bases to all sequences

=item B<-m, --mismatches> INT

Number of mismatches allowed in the ending stretch of [N/G]s.

=item B<-q, --minavgqual> FLOAT

Discard reads having the average quality lower than this number.

=back

=head1 AUTHOR AND BUGS

Written by Andrea Telatin, January 2016. 
Please report bugs to <andrea.telatin@gmail.com>.

=head1 LICENCE

This script is open sourced under GPL3 licence. 
For further details please see <http://www.gnu.org/licenses/gpl-3.0.en.html>.
