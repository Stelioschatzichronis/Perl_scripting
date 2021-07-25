    use strict;
    use warnings;
    use Parallel::ForkManager;
    use Data::Dumper qw(Dumper);
     
    #my $forks = shift or die "Usage: $0 N\n";
     
    my @numbers = (1,2,3,4,5,6,7,3.5);
    my %results;
    my $MAX_PROCESSES=4;
    #print "Forking up to $forks at a time\n";
    #my $pm = Parallel::ForkManager->new($forks);

    my $pm = new Parallel::ForkManager($MAX_PROCESSES);
    $pm->run_on_finish( sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $q = $data_structure_reference->{input};
	my $p = $data_structure_reference->{input2};
        $results{$q}{$p} = $data_structure_reference->{result};
    });

	my $p =0;
    foreach my $q (@numbers) {
        my $pid = $pm->start and next;
        my $res = calc($q);
        $pm->finish(0, { result => $res, input => $q, input2 => $p });
    }
    $pm->wait_all_children;
     
    print Dumper \%results;
     
print "\n";
print "\n";
#my %new_hash = %{$pm};
foreach my $counter (keys %results){
	foreach my $name (keys %{$results{$counter}}){
print "{".$counter."} {".$name."} ".$results{$counter}{$name}."\n";
}
}



    sub calc {
        my ($n) = @_;
        my $sum = 0;
        for (1 .. $n) {
            $sum += 3;
        }
        return $sum;
    }
