# By Zhuo CHEN, IGDB, CAS
# Contact: zhuochen@genetics.ac.cn
# Contact: chenomics@163.com

use strict;
use warnings;
use Getopt::Long;

my($genoFile,$outFile,$windowSize,$stepSize,$list1,$list2);

my $usage = "\nCalculating hapDiv, hapDist and HFI(Haplotype Fixation Index) to measure the differentiation genome regions between two homozygous crop (such as rice) populations.\n\n";
$usage .= "USAGE:\nperl $0 --in <geno file> --out <output file> --list1 <sample list 1> --list2 <sample list 2>\n";
$usage .= "--geno <geno file>: a tab-deleimited file with header line. Each line consists of:\n";
$usage .= "\t[CHROM] [POS] [REF] [ALT] [GENO 1] [GENO 2] ... [GENO n]\n";
$usage .= "\tOnly bi-allelic SNPs are allowed, coded as 0 for REF, 1 for ALT and - for HET or missing.\n";
$usage .= "\t[CHROM] names shall not contain ':'\n";
$usage .= "--window <window size>. Default=10000\n";
$usage .= "--step <step size>. Default=<window size>\n";
$usage .= "--keep <keep sample list>\n";
$usage .= "\nBy Zhuo CHEN, contact: chenomics\@163.com or zhuochen\@genetics.ac.cn\n\n";

GetOptions(
	"in=s" => \$genoFile,
	"out=s" => \$outFile,
	"window=s" => \$windowSize,
	"step=s" => \$stepSize,
	"list1=s" => \$list1,
	"list2=s" => \$list2,
) or die $usage;

die $usage unless(defined $genoFile and defined $outFile and defined $list1 and defined $list2);

unless(defined $windowSize){
	$windowSize = 10 * 1000;
}
unless(defined $stepSize){
	$stepSize = $windowSize;
}

my %hash_list1;
my %hash_list2;
& readList($list1,\%hash_list1);
& readList($list2,\%hash_list2);

open(IN,"<$genoFile") or die $!;
open(OUT,">$outFile");
print OUT "#CHROM\tSTART\tEND\thapPi1\thapPi2\thapDist\tHFI0\tHFI\n";

my %hash_bin = ();
my $bin_tmp = "null:-1";
my @keepRanks = ();
my @allSamples = ();
# define recycle variables
my @binpool;
my %hash_bin1;
my %hash_bin2;
my $sum1 = 0;
my $sum2 = 0;

while(<IN>){
	chomp;
	
	my($chr,$pos,$ref,$alt,@datas) = split/\t/;
	if($_ =~ /^#/){
		for(my $i = 0; $i < @datas; $i++){
			next unless(defined $hash_list1{$datas[$i]} or defined $hash_list2{$datas[$i]});
			push @keepRanks, $i;
		}
		@allSamples = @datas;
		print "# Reamined sample: ".@keepRanks."\n";
		next;
	}

	# CHROM names shall not contain ':'
	$chr =~ s/:/_/g;

	my $binRank = int(($pos-1)/$stepSize);
	my $bin = "$chr:$binRank";
	
	# output data when entering a new bin
	if($bin ne $bin_tmp){
		# initiate new bins and missing bins
		my($chr_tmp,$rank_tmp) = split/:/,$bin_tmp;
		if($chr_tmp eq $chr or $chr_tmp eq "null"){
			$rank_tmp++;
			for(my $i = $rank_tmp; $i <= $binRank; $i++){
				push @binpool, "$chr:$i";
			}
		}else{
			push @binpool, $bin;
		}
		# check the ended bins and output
		for(my $i = 0; $i < @binpool; $i++){
			my $bin_act = $binpool[$i];
			my($chr_,$rank_) = $bin_act =~ /(\S+):(\d+)/;
			my $start_ = $rank_ * $stepSize + 1;
			my $end_ = $rank_ * $stepSize + $windowSize;
			last if($chr eq $chr_ and $pos <= $end_);
			#print "out: $bin_act\n";
			my @values = binPipe(\%hash_list1,\%hash_list2,\%{$hash_bin{$bin_act}});
			my $line = "$chr_\t$start_\t$end_\t".join("\t",@values)."\n";
			print OUT $line;
			splice(@binpool,$i,1);
			$i--;
			delete($hash_bin{$bin_act});
		}
		RESET:
		$bin_tmp = $bin;
	}
	
	# data input
	foreach my $bin_act(@binpool){
		foreach my $rank(@keepRanks){
			push @{$hash_bin{$bin_act}{$allSamples[$rank]}{genos}}, $datas[$rank];
		}
	}
}
close IN;

if(@binpool > 0){
	for(my $i = 0; $i < @binpool; $i++){
		my $bin_act = $binpool[$i];
		my($chr_,$rank_) = $bin_act =~ /(\S+):(\d+)/;
		my $start_ = $rank_ * $stepSize + 1;
		my $end_ = $rank_ * $stepSize + $windowSize;
		my @values = binPipe(\%hash_list1,\%hash_list2,\%{$hash_bin{$bin_act}});
		my $line = "$chr_\t$start_\t$end_\t".join("\t",@values)."\n";
		print OUT $line;
		splice(@binpool,$i,1);
		$i--;
		delete($hash_bin{$bin_act});
	}	
}
close OUT;

sub binPipe{
	my($hashList1,$hashList2,$hashBin) = @_;
	my %hash_list1 = %{$hashList1};
	my %hash_list2 = %{$hashList2};
	my %hash_bin = %{$hashBin};
	my %hash_bin1 = ();
	my %hash_bin2 = ();
	my $sum1 = 0;
	my $sum2 = 0;
	foreach my $sample(sort keys %hash_bin){
		my $haplotype = join("|",@{$hash_bin{$sample}{genos}});
		my $sampleName = $sample;
		if(exists $hash_list1{$sampleName}){
			$sum1++;
			$hash_bin1{$haplotype}{freq}++;
		}
		if(exists $hash_list2{$sampleName}){
			$sum2++;
			$hash_bin2{$haplotype}{freq}++;
		}
	}
	
	# calculate Pi and Dist
	my $Dist = hapDist($sum1,$sum2,\%hash_bin1,\%hash_bin2);
	my $Pi1 = hapPi($sum1,\%hash_bin1);
	my $Pi2 = hapPi($sum2,\%hash_bin2);
	my $fst1 = $Dist - ($Pi1+$Pi2);
	my $PiMin = $Pi1;
	if($PiMin > $Pi2){
		$PiMin = $Pi2;
	}
	my $fst2 = $Dist - $PiMin;
	return($Pi1,$Pi2,$Dist,$fst1,$fst2);
}

sub hapDist{
	my($sum1,$sum2,$hash1,$hash2) = @_;
	my $Dist = 0;
	my @haps1 = sort {${$hash1}{$b}{freq} <=> ${$hash1}{$a}{freq}} keys %{$hash1};
	my @haps2 = sort {${$hash2}{$b}{freq} <=> ${$hash2}{$a}{freq}} keys %{$hash2};
	for(my $i = 0; $i < @haps1; $i++){
		my $hap1 = $haps1[$i];
		my $freq1 = ${$hash1}{$hap1}{freq}/$sum1;
		for(my $j = 0; $j < @haps2; $j++){
			my $hap2 = $haps2[$j];
			my $freq2 = ${$hash2}{$hap2}{freq}/$sum2;
			next if($hap1 eq $hap2);
			my($both,$diff) = hapCompare($hap1,$hap2);
			next unless($diff > 0);
			$Dist += $freq1 * $freq2;
		}
	}
	return($Dist/2);
}

sub hapPi{
	my($sum,$in) = @_;
	my %hash = %{$in};
	my @haps = sort {$hash{$b}{freq} <=> $hash{$a}{freq}} keys %hash;
	my $Pi = 0;
	for(my $i = 0; $i < @haps; $i++){
		my $hapi = $haps[$i];
		my $freqi = $hash{$hapi}{freq}/$sum;
		for(my $j = $i + 1; $j < @haps; $j++){
			my $hapj = $haps[$j];
			my $freqj = $hash{$hapj}{freq}/$sum;
			my($both,$diff) = hapCompare($hapi,$hapj);
			next unless($diff > 0);
			$Pi += $freqi * $freqj;
		}
	}
	return($Pi);
}

sub hapCompare{
	my($hap1,$hap2) = @_;
	my @arr1 = split/\|/, $hap1;
	my @arr2 = split/\|/, $hap2;
	my $both = 0;
	my $diff = 0;
	if(@arr1 != @arr2){
		die "#ERROR: two haplotypes in hapCompare must be of equal length!\n";
	}
	for(my $i = 0; $i < @arr1; $i++){
		my $base1 = $arr1[$i];
		my $base2 = $arr2[$i];
		next unless($base1 =~ /[0-9]/ and $base2 =~ /[0-9]/);
		$both++;
		next if($base1 eq $base2);
		$diff++;
	}
	return($both,$diff);
}

sub readList{
	my($list,$hash) = @_;
	open(IN,"<$list") or die $!;
	while(<IN>){
		chomp;
		next if($_ =~ /^#/);
		my($sample,$others) = split/\t/, $_, 2;
		${$hash}{$sample}{input} = 0;
	}
	close IN;
}
