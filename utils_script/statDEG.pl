use strict;
use warnings;


my $usage=<<USAGE;
   perl $0 DEGs file
   The original code was utilized to statistically determine the number of cell types \
   demonstrating significant expression differences of certain genes.
USAGE


my $file=shift;
my %hash=();my %count=();
open IN, $file;
<IN>;
while(<IN>){
	chomp;
	my @t=split(/,/);
	next if ($t[-1] eq "No");
	push @{$hash{$t[0]}}, $t[5];
}
foreach my $gene(keys %hash){
	my $n=scalar(@{$hash{$gene}});
	if($n>8){
		$n=">8";
	}
	$count{$n}++;
}
print "ctnumber\tnumber\n";
foreach my $n (keys %count){
	print "$n\t$count{$n}\n";
}
