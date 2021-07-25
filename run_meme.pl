use strict;
use warnings;
use Parallel::ForkManager;


#----------------------------------loading files and finding the genes without known TFs------------------------------------------


my $file_genes = $ARGV[0] or die "Need to get related_genes.txt file on the command line\n";
my $file_seq = $ARGV[1] or die "Need to get martquery_1008224531_649.txt file on the command line\n";

my @genes_names;
my @seq_ar;
my %check;





open(FILE1,$file_genes) or die "$!\n";
my $r;
my $i=0;
my @r_genes_str;
while(my $line=<FILE1>){
		chomp $line; 
		push @r_genes_str,$line;
		$i++;
	}
#insert into a 2-Dhash


$i=0;
my $j;
my $sub_seq;
my $primary_gene;
my $related_gene;
my %genes; #genes used for de novo with sequences as keys


while($i<scalar(@r_genes_str)){
		$j=0;
		$primary_gene = substr($r_genes_str[$i], 0, 15);
		$j=16;
		$genes{$primary_gene}{$primary_gene} = "not";
		while($j<length($r_genes_str[$i])){
			$related_gene = substr($r_genes_str[$i], $j, 15);		
			$genes{$primary_gene}{$related_gene} = "not";
			$j=$j+16;
		}
		$i++;
	}



#insert in each key of the hash the sequence


my $l = scalar(@seq_ar);
my @seq_;



#------------extract all genes found to be expressed related to the de novo genes in files

my $seq;
my $flag=0;
system("mkdir files_fasta");
foreach my $primary (keys %genes){
		
		open(FILE1,">","files_fasta/de_novo".$primary.".fasta") or die "$!\n";
		foreach my $related (keys %{$genes{$primary}}){		
			$flag=0;
			$r = $related;
			
			open(FILE2,$file_seq) or die "$!\n";
			$i=0;
			$seq="";
			while(my $line=<FILE2>){
				chomp $line; 
				if($r=~substr($line,1,15)){
					$flag = 1;
				}
					if($flag==1){
						
						if($i<35){
							if($i>0){
							$seq = $seq.$line."\n";
							if($i==34){
							$genes{$primary}{$related} = $seq;
							print FILE1 $seq;
							}
							}else{$seq = substr($line,0,16)."\n";}
							$i++;
								
						}
						else{$flag=0;
						last;
						}
					}
				
			}
			close(FILE2);
			
		}
		close(FILE1);
print "\n".$primary;
}




#run meme
print "\n";
#run parallel
system("mkdir files_meme");
my $MAX_PROCESSES=11;
my $pm = new Parallel::ForkManager($MAX_PROCESSES);

#print "\n";
foreach my $primary (keys %genes){
		#print $primary."\n";
		$r="files_fasta/de_novo".$primary.".fasta";
		$j ="files_meme/de_novo".$primary;
		my $pid = $pm->start and next;
		system("meme '$r' -oc $j -nmotifs 10");
		$pm->finish;
		close(FILE1);
		}
$pm->wait_all_children;


print "\n done with meme\n";
system("mkdir files_mast");
#now we must use mast for de novo motifs
foreach my $primary (keys %genes){
		print $primary."\n";
		$r="files_meme/de_novo".$primary."/meme.txt";
		$j ="files_mast/de_novo".$primary;
		my $pid = $pm->start and next;
		system("mast '$r' '$file_seq' -oc $j -remcorr");
		$pm->finish;
		close(FILE1);
		}
$pm->wait_all_children;





