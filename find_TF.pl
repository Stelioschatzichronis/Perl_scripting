use strict;
use warnings;
use Parallel::ForkManager;


#----------------------------------loading files------------------------------------------


my $file_JASPAR = $ARGV[0] or die "Need to get CSV file on the command line\n";
my $file_seq = $ARGV[1] or die "Need to get summit.bed file on the command line\n";
open(FILE1,$file_JASPAR) or die "$!\n";
open(FILE2,$file_seq) or die "$!\n";
my @jaspar_ar;
my @TF_names;
my $i=0;
	while(my $line=<FILE1>){
		chomp $line; 
		if(($i%5==0)||($i==0)){ 
			push @TF_names,$line;}
		else{
				push @jaspar_ar,$line;	
		;}
		$i++			
	}
close(FILE1);





$i=0;
my @genes_names;
my @seq_ar;
	while(my $line=<FILE2>){
		chomp $line; 
		if(($i%35==0)||($i==0)){push @genes_names,$line;}
		else{push @seq_ar,$line;}
		$i++;
	}


close(FILE2);













my $l = scalar(@seq_ar);
my @seq_;
#------------merge array----one gene-one string-----------
$i=0;
my $counter;
my $line;
	while($i<$l){
		$counter=$i;
		$line = "";
		while($counter<($i+34)){
			if(($counter%34==0)||($counter==0)){
				$line = $seq_ar[$counter];
				#print $counter."\n";
			}else{
				$line = $line.$seq_ar[$counter];
				#print $counter."\n";	
			}
			#sleep(1);
			$counter++;
		}
		push @seq_,$line;	
		$i = $i+34;				
	}
undef $i;
undef $counter;






#----------------------------------search each gene-----------------------------
#$seq_[0], sequence of gene[0], genes_names[0], name of gene[0]
#$JASPAR_ar, array with PFM, TF_names[0], name of TF[0]
my $maxscore;
my $score;
my $sub_seq;
my $letter;
my @temp1;
my @temp2;
my @temp3;
my @temp4;
my $pm;
my $counter_child=0;
my $MAX_PROCESSES=4;
my $flag=0;
$i=0;
$pm = new Parallel::ForkManager($MAX_PROCESSES);



#one file to copy them all scalar(@genes_names)


print "\nnumber of genes are ".scalar(@genes_names)*scalar(@TF_names);

#checking all TFs 
for( my $j = 0; $j <  scalar(@TF_names)-1; $j++ ){	
	print "\n TF is ".$j."\n";
	
	@temp1 = split(' ', $jaspar_ar[$j*4]);
	@temp2 = split(' ', $jaspar_ar[$j*4+1]);
	@temp3 = split(' ', $jaspar_ar[$j*4+2]);
	@temp4 = split(' ', $jaspar_ar[$j*4+3]);
	$maxscore = 0;
	
	
	#calculating max scores for j TF.
	for( my $o = 0; $o < scalar(@temp1); $o++ ){		
		if($temp1[$o] > $temp2[$o]){
			if($temp1[$o] > $temp3[$o]){
				if($temp1[$o] > $temp4[$o]){$maxscore = $maxscore + $temp1[$o];}	
				else{$maxscore = $maxscore + $temp4[$o];}
			}
			else{	if($temp3[$o] > $temp4[$o]){$maxscore = $maxscore + $temp3[$o];}	
				else{$maxscore = $maxscore + $temp4[$o];}
			}		
		}else{
			if($temp2[$o] > $temp3[$o]){
				if($temp2[$o] > $temp4[$o]){$maxscore = $maxscore + $temp2[$o];}	
				else{$maxscore = $maxscore + $temp4[$o];}
			}else{
				if($temp3[$o] > $temp4[$o]){$maxscore = $maxscore + $temp3[$o];}	
				else{$maxscore = $maxscore + $temp4[$o];}
			}
		}

	}#end of for $o
	
	
	
	

	#checking for TF $j all genes scalar(@genes_names)
	for( my $h = 0; $h < 30; $h++ ){
	$counter_child++;
	my $pid = $pm->start and next;	
	open(FILE,">", "file_TF".$counter_child.".txt") or die "$!\n";	

	


	#checking TF j in sub sequences of gene h
	for( my $k = 0; $k < length($seq_[$h]); $k++ ){
		$sub_seq = substr($seq_[$h], $k, scalar(@temp1));
		$score=0;
		
		#checking sub sequence starting in $k
		for( my $m = 0; $m < length($sub_seq); $m++ ){

			if(substr($sub_seq, $m, 1)=~ /A/){
				$score = $temp1[$m] + $score;
			}
			if(substr($sub_seq, $m, 1)=~ /C/){
				$score = $temp2[$m] + $score;
			}
			if(substr($sub_seq, $m, 1)=~ /G/){
				$score = $temp3[$m] + $score;
			}
			if(substr($sub_seq, $m, 1) =~ /T/){
				$score = $temp4[$m] + $score;
			}
		}
		if($score==$maxscore){
			#print "\n"."gene ".substr($genes_names[$h],1, 13)." is binded by TF ".$TF_names[$j]; 
			print FILE substr($genes_names[$h],1, 13)."	".$TF_names[$j]."\n";			
			$k =length($seq_[$h])+2;	#terminate, tf found in gene		
		}
	
	}
	close(FILE);
	$pm->finish;
	

	
		
	
	}#for ------------- $h check each gene
	$pm->wait_all_children;
	}#----for ----$j check each TF----
	
	
















#-------------------write the sub files in the main file
	my @tss_a;
	$i=1;	
	while($counter_child>=$i){
	open(FILE1,"<", "file_TF".$i.".txt") or die "$!\n";		
	while(my $line1=<FILE1>){
		chomp $line1; 
		push @tss_a,$line1;
		
	}
	close(FILE1);
	unlink "file_TF".$i.".txt";
	$i++;
	}
	
	undef $i;
	open(FILE2,">", "final_TF.txt") or die "$!\n";		
	my $string = join "\n", @tss_a;
	print FILE2 $string;  
	close(FILE2);

	#write into the final file
	open(FILE2,"<", "final_TF.txt") or die "$!\n";
	$counter=0;		
	while(my $line=<FILE2>){
		chomp $line; 
		push @tss_a,$line;	
		$counter++;
	}
	close(FILE2);
	print "\nnumber is ".$counter;		
		









