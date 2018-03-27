#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use feature qw(say);
use Getopt::Long;    #use for adding parameters
use Pod::Usage;      #use for adding parameters

#GLOBALS
my $swissProtFileData;    # I use this to define hla_dat
my $allele;               # This is the allele name from the HLA.DAT database
my @HLAS
  ; # This stores the allele names in a proper formate for the sequence program from the winners file
my $sequence;   # This gets the sequence for the alleles in the HLA.DAT database
my %hla;        # A hash to store alleles and sequences from HLA.DAT database
my $winners_file;    # This is the filename for the winners.txt file
my $hla_dat = '';    # Filename for HLA.DAT
my $HLAS    = '';    # Filename for HLAS file


# This is to define usage of the program
my $usage = "\n\n$0 [options] \n
Options:
This program takes a HLAS file and makes a patient.hlaFasta.fa for polysolver

--winner_in	This is the HLAS file.
--hla_dat_in	This is the hla.dat file.
--help		Show this message.
\n";

#check the flags
GetOptions(
	'winner_in=s'  => \$winners_file,
	'hla_dat_in=s' => \$hla_dat,
	help           => sub { pod2usage($usage); },
) or pod2usage(2);

# Get 'patient' info from filename
my $patient;
if ($winners_file =~ /(.*\/)*(.+)\.w+/){
	$patient = $2;
}


# Get names for output files
my $name      = join ".", $patient, "hlaFasta.fa";
my $hlas_name = join ".", $patient, "hlas";

# Open hla.dat file
unless ( open( hla_dat, "<", $hla_dat ) ) {
	die "Cannot open ", $hla_dat, " for reading: $!";
}

# open winner
unless ( open( WINNER, "<", $winners_file ) ) {
	die "Cannot open ", $winners_file, " for reading: $!";
}

# Open file to write
unless ( open( FILE_OUT, ">", $name ) ) {
	die "Cannot open file for writing: ", $!;
}

# Open file to write
unless ( open( HLAS_OUT, ">", $hlas_name ) ) {
	die "Cannot open file for writing: ", $!;
}

# Change the record seperator
$/ = "\n";

# Get alleles from winner.txt file
while (<WINNER>) {
	getWINNER($_);
}

# close winner file
close WINNER;

# Change the record seperator for HLA.DAT file
$/ = "//\n";

##get the entire file in the variable
while ( $swissProtFileData = <hla_dat> ) {

	# This gets the allele and puts it into a scalar
	$allele = getAllele( $swissProtFileData, 'DE' );

	# This gets the sequence and puts it into a scalar
	$sequence = getSequence($swissProtFileData);

	# This puts the allele and sequence into a hash.
	$hla{$allele} = $sequence;
}

# Change the record seperator for HLAS file
$/ = "\n";

# Get HLAS
getHLAS(@HLAS);

close FILE_OUT;

sub getWINNER {

	# This reformats the winners.hla.txt file to match the hlas format

	# Read file
	my ($Winner_file) = @_;

	# Chomp
	chomp $Winner_file;

	# set variable for an array of alleles
	my @alleles;

	# Split input by whitespace
	@alleles = split /\s+/, $Winner_file;

	# Remove first element because it is a header
	shift @alleles;

	# Add each element to the global array @HLAS
	for my $line (@alleles) {
		push @HLAS, $line;
	}
}

sub getHLAS {

# This matches the winners.hla.txt alleles to the hla.dat alleles and prints them in fasta format
# this also prints the HLAS file with the alleles.

	my (@lines) = @_;
	foreach my $line (@lines) {

		# if the allele exists
		if ( exists $hla{$line} ) {
			print FILE_OUT ">$line\n$hla{$line}\n";
			print HLAS_OUT $line, "\n";
		}

		elsif ( !exists $hla{$line} ) {

			# if the allele has 8 positions then find one with 6 or 4
			if ( $line =~ m/(((.*_.*_\d+_\d+)_\d+)_\d+)/ms ) {

				# Define the 6 position alllele
				my $pos8 = $1;
				my $pos6 = $2;
				my $pos4 = $3;

				# Array with numbers
				my @numbers = ( 01 ... 99 );

				# Go through each number until a match is found
				for my $num (@numbers) {
					my $n = sprintf( "%02d", $num );
					my $line2 = join( "_", $pos6, $n );
					my $line3 = join( "_", $pos4, $n );

					# find 8 position
					if ( exists $hla{$line2} ) {
						print FILE_OUT ">$line2\n$hla{$line2}\n";
						print HLAS_OUT $line2, "\n";
						last;
					}

					# find 6 position
					elsif ( exists $hla{$line3} ) {

						print FILE_OUT ">$line3\n$hla{$line3}\n";
						print HLAS_OUT $line3, "\n";
						last;
					}

					# find 4 position
					elsif ( exists $hla{$pos4} ) {

						print FILE_OUT ">$pos4\n$hla{$pos4}\n";
						print HLAS_OUT $pos4, "\n";
						last;
					}
				}
			}

			# if the allele has 6 positions then find one with 8 or 4
			elsif ( $line =~ m/(((.*_.*_\d+)_\d+)_\d+)/ms ) {

				# Define the 6 position alllele
				my $pos6 = $1;
				my $pos4 = $2;
				my $pos2 = $3;

				# Array with numbers
				my @numbers = ( 01 ... 99 );

				# Go through each number until a match is found
				for my $num (@numbers) {
					my $n = sprintf( "%02d", $num );
					my $line2 = join( "_", $pos6, $n );
					my $line3 = join( "_", $pos4, $n );
					my $line4 = join( "_", $pos2, $n );

					# find 8 position
					if ( exists $hla{$line2} ) {
						print FILE_OUT ">$line2\n$hla{$line2}\n";
						print HLAS_OUT $line2, "\n";
						last;
					}

					# find 6 position
					elsif ( exists $hla{$line3} ) {

						print FILE_OUT ">$line3\n$hla{$line3}\n";
						print HLAS_OUT $line3, "\n";
						last;
					}

					# find 4 position
					elsif ( exists $hla{$pos4} ) {

						print FILE_OUT ">$pos4\n$hla{$pos4}\n";
						print HLAS_OUT $pos4, "\n";
						last;
					}
				}
			}

			# if the allele has 4 positions then find one with 6 or 8
			elsif ( $line =~ m/(.*_.*_\d+_\d+)/ms ) {

				# Define variable
				my @line2;
				my @line3;

				# Define the 6 position alllele
				my $seq = $1;

				# Array with numbers
				my @numbers = ( 01 ... 99 );

				# Go through each number until a match is found
				for my $num (@numbers) {
					my $n = sprintf( "%02d", $num );
					my $line2 = join( "_", $seq, $n );
					push @line2, $line2;

				}

				# Make 6 positions
				for my $lines (@line2) {
					for my $num (@numbers) {
						my $n = sprintf( "%02d", $num );
						my $line3 = join( "_", $lines, $n );
						push @line3, $line3;
					}
				}

				for my $n3 (@line3) {
					if ( exists $hla{$n3} ) {
						print FILE_OUT ">$n3\n$hla{$n3}\n";
						print HLAS_OUT $n3, "\n";
						last;
					}

					else {
						for my $n2 (@line2) {
							if ( exists $hla{$n2} ) {
								print FILE_OUT ">$n2\n$hla{$n2}\n";
								print HLAS_OUT $n2, "\n";
								last;
							}
						}
						last;
					}
				}
			}
		}
		else {
			print FILE_OUT ">$line\nerror\n";
		}
	}
}

sub getSequence {

 # This changes the allele from the fasta format to the HLA.DAT format to search
	my ($GB_file) = @_;

	# Define variable
	my $seq;

# This uses regex and skips the first line and gets everything untill it sees \\
# Then it removes the whitespace and numbers and returns the scalar.
	if ( $GB_file =~ /^SQ.*;$(.*)\/\//ms ) {
		$seq = $1;

		# This removes the whitespaces and numbers
		$seq =~ s/[\s\d]//g;
		return $seq;
	}
	else {
		return "error";
	}
}

sub getAllele {

# This reformats the allele from the hla.dat file to match the winners.hla.txt allele

	# Read file line by line
	my ( $record, $DE ) = @_;

	# set variable for an array of alleles
	my @allele;
	my $allele;

	# Use regex to get the allele
	if ( $record =~ /^$DE\s+(.*?),/ms ) {

		# Set group one to variable allele
		my $allele = $1;

		# Split by -
		@allele = split /-/, $allele;

		# join by _
		$allele = join( '_', @allele );

		# Split by *
		@allele = split /\*/, $allele;

		# join by _
		$allele = join( '_', @allele );

		# Split by :
		my @allele = split /:/, $allele;

		# join by _
		$allele = join( '_', @allele );

		# return lowercase allele
		return lc($allele);
	}
	else { return 'error' }
}

# Close output files
close hla_dat;
close HLAS_OUT;

## For troubleshooting, check the fasta database ##
## Open file to write
#unless ( open( FASTA_OUT, ">", 'hla.dat.fa' ) ) {
#	die "Cannot open file for writing: ", $!;
#}

#i use a while loop to print out the keys and values of the hash
#while ( my ( $key, $value ) = each %hla ) {
#	print FASTA_OUT ">$key\n$value\n";
#}
