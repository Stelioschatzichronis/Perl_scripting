use strict;
use warnings;
use Parallel::ForkManager;

#----------------------------------loading files------------------------------------------


my $file_tss = $ARGV[0] or die "Need to get CSV file on the command line\n";
my $file_bed = $ARGV[1] or die "Need to get summit.bed file on the command line\n";
open(FILE1,$file_tss) or die "$!\n";
open(FILE2,$file_bed) or die "$!\n";
my @tss_ar;
my $i=0;
	while(my $line=<FILE1>){
		chomp $line; 
		if($i==1){push @tss_ar,$line;}
		else{$i++;}			
	}

close(FILE1);
$i=0;
my @promoter_region;
	while(my $line=<FILE2>){
		chomp $line; 
		if($i==1){push @promoter_region,$line;}
		else{$i++;}
	}

undef $i;
close(FILE2);


#----------------------------------inserting into hashes------------------------------------------


#must take 2500b threshold for tss and promoter regions

#put data in hash, chromosome, strand and other data 
my %gene;
my $chr;
my $tss;
my @what;
my (%tss_hash,%macs_tss);
	

	#{strand}{chromosome}{position of peak}, value--> geneID
	for( my $i = 0; $i < scalar(@tss_ar); $i++ ){
 		@what = split(' ', $tss_ar[$i]);
 		$tss_hash{$what[5]}{$what[0]}{$what[1]} = $what[4]; 
	}

	#{chromosome}{number} --> value of position of tss
	for( my $i = 0; $i < scalar(@promoter_region); $i++ ){
 		@what = split(' ', $promoter_region[$i]);
		$what[0] =~ s/chr//g; #removing chr from chromosomes' names
 		$macs_tss{$what[0]}{$i} = $what[1]; 
	}

	undef @promoter_region;
	undef @tss_ar;
	

#---------------------find expressed genes ---------------------------
	
	my %expressed_genes;
	print "\n";	
			my $counter=0;	
			my $counter_child=0;
			foreach my $strand (sort keys %tss_hash) {		
			$counter_child++;
			#sockets are not an optimal way of using parallel programming in perl		
				foreach my $chr2 (sort keys %{$tss_hash{$strand}}) {
				foreach my $geneID (sort keys %{$tss_hash{$strand}{$chr2}}) {
				foreach my $number (sort keys %{$macs_tss{$chr2}}) {
				#check if peak and tss are close enough
					if(exists($macs_tss{$chr2}{$number})){
					if(abs($tss_hash{$strand}{$chr2}{$geneID}-$macs_tss{$chr2}{$number}) <=2500){
						$expressed_genes{$strand}{$chr2}{$geneID} = $geneID;
						#using geneID to store only unique results for each gene
						#overwriting existing keys for same geneID	
					}
					}	
				}
				}
				}
			}
				open(FILE,">", "expr_genes.txt") or die "$!\n";	
				foreach my $strand (sort keys %expressed_genes) {
				foreach my $chr2 (sort keys %{$expressed_genes{$strand}}) {
				foreach my $geneID (sort keys %{$expressed_genes{$strand}{$chr2}}) {			
				print FILE $strand."	".$chr2."	".$expressed_genes{$strand}{$chr2}{$geneID}."\n";
				}
				}
				}
					close(FILE);
		
