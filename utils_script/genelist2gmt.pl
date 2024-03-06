use strict;
use warnings;

my $usage=<<USAGE;
   perl $0 genes list 
   This script was used for transform genes list to gmt format.
USAGE
die $usage if @ARGV !=1;


my $Rout=shift;
open IN, $Rout;
<IN>;my %hash=();
my $url="http://www.gsea-msigdb.org/gsea/msigdb/cards";
while(<IN>){
        chomp;
        $_=~s/"//g;
        my @t=split(/\t/);
	#shift @t;
        my $key="$t[0]\t$url";
        push @{$hash{$key}},$t[1];
}
foreach my $key(keys %hash){
        print "$key\t",join("\t",@{$hash{$key}}),"\n";
}

