use strict;
use warnings;
use Parallel::ForkManager;


#----------------------------------loading files and finding the genes without known TFs------------------------------------------


my $file_TF = $ARGV[0] or die "Need to get file_TF file on the command line\n";
my $file_seq = $ARGV[1] or die "Need to get martquery_1008224531_649 file on the command line\n";


my @genes_names;
my @seq_ar;
my %check;





open(FILE1,$file_TF) or die "$!\n";
my $r;

while(my $line=<FILE1>){
		chomp $line; 
		if(length($line)>20){
		$r = substr($line, 33, 15);
		#print $r."\n";
		$check{$r}=$r;
		}
	}
close(FILE1);



my $i=0;
my %de_novo_gene;
my @de_novo_ar;
my $c=0;
open(FILE2,$file_seq) or die "$!\n";
open(FILE1,">","de_novo_genes.txt") or die "$!\n";
	while(my $line=<FILE2>){
		chomp $line; 
		#input all genes			
			if(($i%35==0)||($i==0)){

			$r = substr($line, 1, 15);
			if(defined $check{substr($line, 1, 15)}){
				#print substr($line, 1, 15);
					
				
			}else{	
				$de_novo_gene{substr($line, 1, 15)} = substr($line, 1, 15); 
				$de_novo_ar[$c] = $line;
				#print $line;
				print $de_novo_ar[$c]."\n";
				print FILE1 $de_novo_ar[$c]."\n";
				$c++;
			     }
			}
			
			$i++;
			
	}
close(FILE2);
close(FILE1);
