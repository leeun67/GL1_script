use strict;

my $file = shift;
open IN, $file;
my %sequences;
my %samples;
while(my $line = <IN>){
	chomp $line;
	next if($line =~ /^##/);
	my @info = split /\t/, $line;
	if($line =~ /^#/){
		for(my $i = 9; $i < @info; $i++){
			$samples{$i} = $info[$i];
		}
	}
	else{
		for(my $i = 9; $i < @info; $i++){
			my ($gt) = $info[$i] =~ /^(.*?):/;
			my $nt = $gt eq "0/0" ? $info[3] : 
				$gt eq "0/1" ? $info[3] : 
				$gt eq "1/1" ? $info[4] : 
				$gt eq "0|0" ? $info[3] :
				$gt eq "0|1" ? $info[3] :
				$gt eq "1|1" ? $info[4] : "-";      
			$sequences{$i} .= $nt;
		}
	}
}

foreach my $id(sort {$a<=>$b} keys %samples){
	print ">$samples{$id}\n$sequences{$id}\n";
}
