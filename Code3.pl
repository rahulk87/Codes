## Perl code for calculating variant allele fraction 
#!/usr/bin/perl -w

use strict;
use Getopt::Long;
$| = 1;

my $vcf = "";
my $out = undef;
my $columns = undef;
my @columns=();
my $args = scalar(@ARGV);
my $help;
my $result = GetOptions (
  "vcf=s" => \$vcf,
  "columns=s" => \$columns,
  "out=s" => \$out,
  "help" => \$help,
  );

# print usage message if requested or no args supplied
if(defined($help) || !$args) {
  &usage;
  exit(0);
}

$out = $vcf . '.vaf' unless defined($out);
@columns = split(/,/, $columns) unless !defined($columns);

open VCF, "< $vcf" or die "Can't read vcf input file $vcf: $!\n";
open OUT, "> $out" or die "Can't write to file $out: $!\n";
while(<VCF>){
	next if (/^##/);
  # split VCF line to find fields after FORMAT
	chomp;
	print OUT $_;
	my @data=split(/\t/);
	if ($data[0] =~ m/#CHROM/){
		@columns=(9..$#data) if ($#columns == -1);
		foreach my $col (@columns){
      		print OUT "\t$data[$col].depth\t$data[$col].alt\t$data[$col].vaf";
      	}
		print OUT "\n";
		next;
	}
	foreach my $col (@columns){
		my $freq_string = &get_AD($data[8], $data[$col]);
		print OUT $freq_string;
  	}
	print OUT "\n";
}

close VCF;

close OUT;

sub get_AD{
	my ($label, $input_string) = @_;
	my %infor;
	my $output_string="\tna\tna\tna";
	if($input_string ne './.'){
		my @labels=split(/:/,$label);
		my @fields = split(/:/,$input_string);
		foreach my $i(0..$#labels){
			$infor{$labels[$i]}=$fields[$i];
		}
  		if(defined($infor{'AD'}) && $infor{'AD'} =~ /(\d+),(\d+)/){
  			my $dep=$1+$2;
  			if($dep !=0){
  				my $freq=$2/$dep;
  				$output_string="\t$dep\t$2\t$freq";
  			}
  		}
  	}
  		return $output_string;
  	
}


sub usage() {
        my $usage =<<END;
#------------------------#
#  get_allele_freqs.pl   #
#------------------------#
Usage:
perl get_allele_freqs.pl [options]
Options
--help          Display this message and quit
--vcf           Path to where the vcf file is located           [required]
--columns	Defaults to (zero-based) all samples in the VCF. List of comma-separated integers
--out           Path to output file                             [default: ./inFile.classes]
END

        print $usage;
}
