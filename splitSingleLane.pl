#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $out_dir = './';
my $optSuccess = GetOptions('o|output-directory=s' => \$out_dir);
our $version = "1.00";
print STDERR "
  NextSeq Split Lanes v. $version 
 --------------------------------------------------------
 This program split a single NextSeq 500 FASTQ file into
 four files (one per virtual lane).
 Usage:
 nextseq_split_lanes.pl [-o OUTPUT_DIR] FASTQ_FILE
 
 Where OUTPUT_DIR is the destination directory for the 
 four files, and FASTQ_FILE is the input file.
 --------------------------------------------------------
";
my $fileName = $ARGV[0]; 
my $baseName = basename($fileName);

exit unless ($fileName);
my @fn = split /_/, $baseName;
my $fnItems = @fn;

#44715-41897_S1_R1_001.fastq.gz
#print STDERR "$fileName: $fnItems\n";
#44715-41897_S1_L001_R1_001.fastq.gz


die "FATAL ERROR: Input file should be in the format CODE_S1_R1_001.fastq.\n" if ($fnItems != 4); 
print STDERR "\n Preparing output files:\n";
open my $f1, ">$out_dir/$fn[0]_$fn[1]_L001_$fn[2].fastq" || die " FATAL ERROR:\n Unable to WRITE to $out_dir/$fn[0]_$fn[1]_L001_$fn[2].fastq\n";
print STDERR " * $out_dir/$fn[0]_$fn[1]_L001_$fn[2].fastq\n";
open my $f2, ">$out_dir/$fn[0]_$fn[1]_L002_$fn[2].fastq" || die " FATAL ERROR:\n Unable to WRITE to $out_dir/$fn[0]_$fn[1]_L002_$fn[2].fastq\n";
print STDERR " * $out_dir/$fn[0]_$fn[1]_L002_$fn[2].fastq\n";
open my $f3, ">$out_dir/$fn[0]_$fn[1]_L003_$fn[2].fastq" || die " FATAL ERROR:\n Unable to WRITE to $out_dir/$fn[0]_$fn[1]_L003_$fn[2].fastq\n";
print STDERR " * $out_dir/$fn[0]_$fn[1]_L003_$fn[2].fastq\n";
open my $f4, ">$out_dir/$fn[0]_$fn[1]_L004_$fn[2].fastq" || die " FATAL ERROR:\n Unable to WRITE to $out_dir/$fn[0]_$fn[1]_L004_$fn[2].fastq\n";
print STDERR " * $out_dir/$fn[0]_$fn[1]_L004_$fn[2].fastq\n";





open my $fh,  $fileName ||
	die "FATAL ERROR:\nUnable to read input file \"$fileName\".\n";


print STDERR " PARSING FILE...";
while (my $name = <$fh>) {
	my $seq  = <$fh>;
	my $space= <$fh>;
	my $qual = <$fh>;
	#@NB501072:102:HLFLWBGXY:1:11101:15667:1046 1:N:0:CGATGT
	my @items = split /:/, $name;
	my $lane = $items[3];
	if ($lane !~/^[1-4]$/) {
		die "Malformed file: $name\nLane: $lane\n";
	} elsif ($lane == 1) {
		print $f1 "$name$seq$space$qual";
	} elsif ($lane == 2) {
		print $f2 "$name$seq$space$qual";
	} elsif ($lane == 3) {
		print $f3 "$name$seq$space$qual";
	} elsif ($lane == 4) {
		print $f4 "$name$seq$space$qual";
	}

}
print STDERR "DONE!\n";
