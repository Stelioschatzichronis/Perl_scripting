use strict;
use warnings;
use Parallel::ForkManager;
use Data::Dumper qw(Dumper);


my $start = time();
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

#remove extra unecessary line
my $nothing = pop @TF_names;



$i=0;
my @genes_names;
my @seq_ar;
	while(my $line=<FILE2>){
		chomp $line; 
		if(($i%35==0)||($i==0)){
			push @genes_names,$line;}
		else{
			push @seq_ar,$line;
		}
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


my $sub_seq;
my @temp1;
my @temp2;
my @temp3;
my @temp4;
my $pm;
my $pid;
my $MAX_PROCESSES=50;
my $flag=0;
my %hash_of_hash;
my $res;
my $m;
my $flag1;
$i=0;

$pm = new Parallel::ForkManager($MAX_PROCESSES);

$pm->run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
    my $q = $data_structure_reference->{input1};
    my $p = $data_structure_reference->{input2};
    $hash_of_hash{$q}{$p} = $data_structure_reference->{result};
});

#one file to copy them all scalar(@TF_names)-1


my $u = (time()-$start)/60;
#checking all TFs 
for( my $j = 0; $j <  scalar(@TF_names); $j++ ){	
	$u = (time()-$start)/60;
	print "TF is ".$j." "."time is ".$u."\n";
	
	@temp1 = split(' ', $jaspar_ar[$j*4]);
	@temp2 = split(' ', $jaspar_ar[$j*4+1]);
	@temp3 = split(' ', $jaspar_ar[$j*4+2]);
	@temp4 = split(' ', $jaspar_ar[$j*4+3]);
	#print "\n TF is ".$j."\n";

	#checking for TF $j all genes scalar(@genes_names)
	for( my $h = 0; $h < scalar(@genes_names); $h++ ){
	$pid = $pm->start and next;
	$res = "non_TF";
	
	#checking TF j in sub sequences of gene h length($seq_[$h])-1
	for( my $k = 0; $k < length($seq_[$h])-scalar(@temp1)+1; $k++ ){
		$sub_seq = substr($seq_[$h], $k, scalar(@temp1));
		$flag1 = 0;
		$m=0;
		#checking sub sequence starting in $k
		while($m < length($sub_seq)){
			#check if letter is different, if yes maxscore cannot reached	
			if(substr($sub_seq, $m, 1)=~ /A/){		
				if(($temp1[$m]>$temp2[$m])&&($temp1[$m]>$temp3[$m])&&($temp1[$m]>$temp4[$m])){
				}else{ $flag1=1;}				
			}else{
				if(substr($sub_seq, $m, 1)=~ /C/){
					if(($temp2[$m]>$temp1[$m])&&($temp2[$m]>$temp3[$m])&&($temp2[$m]>$temp4[$m])){
					}else{ $flag1=1;}
				}
				else{
					if(substr($sub_seq, $m, 1)=~ /G/){
					if(($temp3[$m]>$temp1[$m])&&($temp3[$m]>$temp2[$m])&&($temp3[$m]>$temp4[$m])){
					}else{ $flag1=1;}
					}else{
					if(($temp4[$m]>$temp1[$m])&&($temp4[$m]>$temp2[$m])&&($temp4[$m]>$temp3[$m])){
					}else{ $flag1=1;}
					}
				}
					
			}
			if($flag1==0){
				$m=$m+1;
			}else{$m = length($sub_seq); $flag1=1;} #print "\n".$m;		
		}
		
		if($flag1==0){
			$hash_of_hash{$TF_names[$j]}{$genes_names[$h]} = substr($genes_names[$h],1, 15)."	".$TF_names[$j];		
			$k =length($seq_[$h])+2;	#terminate, tf found in gene		
			$res = $hash_of_hash{$TF_names[$j]}{$genes_names[$h]};
		}
	}
	
	
	
	if($res=~/non_TF/){
	$pm->finish(0, { result => $res , input1 => $flag, input2 => $flag});
	}else{
	$pm->finish(0,{result => $res ,input1 => \$hash_of_hash{$TF_names[$j]},input2 => \$hash_of_hash{$TF_names[$j]}{$genes_names[$h]}});
	}
	
	
	}#for ------------- $h check each gene
	$pm->wait_all_children;
	}#----for ----$j check each TF----
	
	
	#print Dumper \%hash_of_hash;


#foreach my $counter (keys %hash_of_hash){
	#foreach my $name (keys %{$hash_of_hash{$counter}}){
#print $counter." ".$name." ".$hash_of_hash{$counter}{$name}."\n";
#}
#}




	open(FILE1,">", "file_TF.txt") or die "$!\n";		 
		foreach my $counter (keys %hash_of_hash){
		foreach my $name (keys %{$hash_of_hash{$counter}}){
		print FILE1 $counter." ".$name." ".$hash_of_hash{$counter}{$name}."\n";
		}
		}
	close(FILE1);





my $end = time();
my $minutes = ($end-$start)/60;
print "\nTotal Job took $minutes minutes\n";







