The Perl script will use the winners.hla.txt file and format it so it can be used as
an input for LOHHLA. The program with use the HLA alleles in the winners.hla.txt file
and match it to a database of HLA alleles and get the sequence and put it in FASTA
format. If no matches are found, then the program will add or remove positions (from 4
through 8) until a match is found. The script requires two inputs, the winners.hla.txt
file produced from Polysolver and a reference hla.dat file supplied from LOHHLA. The
Perl script takes about 10 seconds to run.

An example of how the script runs is also shown below:

# For help type perl SeqLOHHLA.pl --help
perl SeqLOHHLA.pl \
--winner_in example_winners.hla \
--hla_dat hla.dat


The R script, poly_bcw_compare_best.r, compares two HLA files. One from Blood Center
of Wisconsin (BCW) and one derived from Polysolver. If there are any differences, this
program defaults to the BCW results and prints out the results in Polysolver format.

An example of how the script runs is also shown below:

# For help type Rscript poly_bcw_compare.r --help
Rscript poly_bcw_compare.r \
--winnersFile example_winners.hla \
--bcwFile example_output.csv

