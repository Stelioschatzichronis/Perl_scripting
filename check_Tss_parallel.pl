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
 		$macs_tss{$what[0]}{$what[1]} = $what[1]; 
	}

	undef @promoter_region;
	undef @tss_ar;
	

#---------------------find expressed genes ---------------------------
	
	my %expressed_pos;
	my %expressed_genes;	
			my $counter=0;	
			my $counter1=0;
			my $counter_child=0;
			my $pm;
			my $MAX_PROCESSES=2;
			$pm = new Parallel::ForkManager($MAX_PROCESSES);
			foreach my $strand (sort keys %tss_hash) {
			
			$counter_child++;
			my $pid = $pm->start and next;
			open(FILE,">", "expr_pos".$counter_child.".txt") or die "$!\n";	
			open(FILE2,">", "expr_genes".$counter_child.".txt") or die "$!\n";	
			#sockets are not an optimal way of using parallel programming in perl		
				foreach my $chr2 (sort keys %{$tss_hash{$strand}}) {
				foreach my $geneID (sort keys %{$tss_hash{$strand}{$chr2}}) {
				foreach my $number (sort keys %{$macs_tss{$chr2}}) {
				#check if peak and tss are close enough
					if(exists($macs_tss{$chr2}{$number})){
					if(abs($tss_hash{$strand}{$chr2}{$geneID}-$macs_tss{$chr2}{$number}) <=2500){
						$expressed_pos{$strand}{$chr2}{$macs_tss{$chr2}{$number}} = $geneID;
						$expressed_genes{$strand}{$chr2}{$geneID} = $geneID;
						#using geneID to store only unique results for each gene
						#overwriting existing keys for same geneID
						#creating file for multiple regions to be inserted in Biomart
		print FILE $chr2.":".$macs_tss{$chr2}{$number}.":".($macs_tss{$chr2}{$number}+2000).":".$strand.",\n";
		print FILE2 $chr2."	".$geneID."	".$macs_tss{$chr2}{$number}."	".$strand.",\n";
						$counter++;
					}
					}	
				}
				}
				}
			
			close(FILE);
			$pm->finish;
			}
	

$pm->wait_all_children;
	

	#-------------------write the sub files in the main expression file
	my @tss_a;
	$i=1;	
	while($counter_child>=$i){
	open(FILE1,"<", "expr_pos".$i.".txt") or die "$!\n";		
	while(my $line=<FILE1>){
		chomp $line; 
		push @tss_a,$line;
		
	}
	close(FILE1);
	unlink "expr_pos".$i.".txt";
	$i++;
	}
	
	undef $i;
	open(FILE2,">", "expr_pos.txt") or die "$!\n";		
	my $string = join "\n", @tss_a;
	print FILE2 $string;  
	close(FILE2);

	#write 2 files into a final file
	open(FILE2,"<", "expr_pos.txt") or die "$!\n";
	$counter=0;		
	while(my $line=<FILE2>){
		chomp $line; 
		push @tss_a,$line;	
		$counter++;
	}
	close(FILE2);
print "\nnumber of tss in genes discovered is ".$counter;
		for( my $i = 0; $i < scalar(@tss_a); $i++ ){
			@what = split(':', $tss_a[$i]);
			$expressed_pos{$what[0]}{$what[1]}{$what[2]} = $what[3];
		}




#-------------------same for gene expression file------------------------------------
my @tss_b;
	$i=1;	
	while($counter_child>=$i){
	open(FILE1,"<", "expr_genes".$i.".txt") or die "$!\n";		
	while(my $line=<FILE1>){
		chomp $line; 
		push @tss_b,$line;
		
	}
	close(FILE1);
	unlink "expr_genes".$i.".txt";
	$i++;
	}
	
	undef $i;
	open(FILE2,">", "expr_genes.txt") or die "$!\n";		
	 $string = join "\n", @tss_b;
	print FILE2 $string;  
	close(FILE2);

	#write 2 files into a final file
	open(FILE2,"<", "expr_genes.txt") or die "$!\n";
	
	while(my $line=<FILE2>){
		chomp $line; 
		push @tss_b,$line;	

	}
	close(FILE2);

		
		my %sh;	
		for( my $i = 0; $i < scalar(@tss_b); $i++ ){
			@what = split('	', $tss_b[$i]);
			$sh{$what[0]}{$what[1]}{$what[3]} = $what[1];
		}

		$counter1=0;
		foreach my $strand (sort keys %sh) {
		foreach my $chr2 (sort keys %{$sh{$strand}}) {
		foreach my $geneID (sort keys %{$sh{$strand}{$chr2}}) {
			$counter1++;	
		}
		}
		}
		print "\nnumber of expressed genes is ".$counter1;
			

		#time perl check_Tss_parallel.pl mart_export.txt file_peak_summits.bed
