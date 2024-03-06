use strict; use warnings;

my $dkd=shift;
my $RA=shift;
my $Restype="";
my %hash=();my %flag=();


my $usage=<<USAGE;
   perl $0 dkd_vs_normal(degs)   RA-dkd_vs_dkd(degs)
   This original script was used for obtaining genes exhibiting opposite trends
USAGE
die $usage if @ARGV !=2;


open IN, $dkd;
<IN>;
while(<IN>){
	chomp;
	my @t=split(/,/);
	next if ($t[-1] eq "No");
	my $key="$t[0]\t$t[-2]";
	$hash{$key}=join("\t",@t);
	$flag{$key}=$t[-1];
}
open IN, $RA;
my $t=<IN>;chomp($t);
$t=~s/,/\t/g;
print "$t\tRA\tRescuedType\n";
while(<IN>){
	chomp;
	my @t=split(/,/);
	next if ($t[-1] eq "No");
        my $key="$t[0]\t$t[-2]";
	$_=~s/,/\t/g;
	if(exists($hash{$key}) && $flag{$key} ne $t[-1]){
		if($flag{$key} eq 'Up'){
			$Restype="DKDUp";
		}
		if($flag{$key} eq 'Down'){
			$Restype="DKDDown"
		}
		print "$hash{$key}\t$t[-1]\t$Restype\n";
	}
}

