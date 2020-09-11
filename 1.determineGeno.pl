#!/usr/bin/perl -w

die "Usage: perl $0 input.vcf.gz\n" if(!defined $ARGV[0]);
my %genodb;
open(TMP, "> format.txt") or die"";
open(IN, "gunzip -dc $ARGV[0]|grep -v '##' |") or die"";
my $line = <IN>;
my @hdb  = split(/\s+/,$line); 
print TMP "SNPid	FhF3	FhF4	FhF6	FhHM1	FhHM2	FhHM3	Xgenotype	Ygenotype\n";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	next if($data[0] ne "Chr10");
	next if($data[3]=~/,/);
	next if($data[4]=~/,/);
	my $fm3  = $data[9]; $fm3 =~ s/:.*//g;
	my $fm4  = $data[10];$fm4 =~ s/:.*//g;
	my $fm6  = $data[11];$fm6 =~ s/:.*//g;
	my $hm1  = $data[12];$hm1 =~ s/:.*//g;
	my $hm2  = $data[13];$hm2 =~ s/:.*//g;
	my $hm3  = $data[14];$hm3 =~ s/:.*//g;
	
	my %maledb = ();
	$maledb{$hm1}++; $maledb{$hm2}++; $maledb{$hm3}++;
	next if(exists($maledb{"0/0"}) or exists($maledb{"1/1"}) or exists($maledb{"2/2"}));
	my %femaledb = ();
	$femaledb{$fm3}++; $femaledb{$fm4}++; $femaledb{$fm6}++;
	next if(exists($femaledb{"0/1"}));
	next if(exists($femaledb{"0/0"}) and exists($femaledb{"1/1"}));
	next if(exists($femaledb{"./."}) and ($femaledb{"./."}==3));
	my $Xgeno = ""; my $Ygeno = "";
	if(exists($femaledb{"0/0"})){
		$Xgeno  = $data[3];
		$Ygeno  = $data[4];
	}elsif(exists($femaledb{"1/1"})){
		$Xgeno  = $data[4];
		$Ygeno  = $data[3];		
		}
  my $key   = $data[0]."_".$data[1];
  print TMP	"$key	$fm3	$fm4	$fm6	$hm1	$hm2	$hm3	$Xgeno	$Ygeno\n"; 
  $genodb{$key}->{'X'} = $Xgeno;
  $genodb{$key}->{'Y'} = $Ygeno;
	}
close IN;
close TMP;


my %blockdb;
open(OUT, "> pb.variants.genotype.txt") or die"";
open(IN, "phased.info") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $key  = $data[0]."_".$data[1];
	my $Bid  = (split/:/,$data[4])[-1];
	print OUT "$_	";
	print OUT "$genodb{$key}->{'X'}	$genodb{$key}->{'Y'}" if(exists($genodb{$key}));
	print OUT "\n";
	$blockdb{$Bid}->{'tSNPs'}++;
	$blockdb{$Bid}->{'pSNPs'}++  if(exists($genodb{$key}));
	}
close IN;
close OUT;

print "BLOCKid	totalSNPs	phasedSNPs\n";
foreach my $i (sort keys %blockdb){
	my $num_tSNPs = $blockdb{$i}->{'tSNPs'};
	my $num_pSNPs = exists($blockdb{$i}->{'pSNPs'})?$blockdb{$i}->{'pSNPs'}:0;
	print "$i	$num_tSNPs	$num_pSNPs\n";
	}


