#!/usr/bin/perl -w

my %readb;
open(IN, "result.bwa.filtered.csv") or die"";
while(<IN>){
	chomp;
	my @data = split(/,/,$_);
	my $key  = $data[0].",".$data[1];
	$readb{$key}->{'0'} = $data[3];
	$readb{$key}->{'1'} = $data[4];
	}
close IN;


my %posidb;
my %P1XYdb;
my $RBid = 1;
my %BlOCKdb;
open(OUT, "> tmp.txt") or die"";
print OUT "ID	CHRN	POSI	REF	ALT	HA	HB	Xbase	Ybase	Xgeno	Ygeno\n";
open(IN, "pb.variants.genotype.txt") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $chrn = $data[0];
	my $posi = $data[1];
	my $refN = $data[2];
	my $altN = $data[3];
	my $geno = (split/:/,$data[4])[0];
	my ($P1,$P2) = split(/\|/,$geno);
	my $Bid  = (split/:/,$data[4])[-1];
	my $X    = (defined $data[5])?$data[5]:"NA";
	my $Y    = (defined $data[6])?$data[6]:"NA";
	my %sexdb;
	if($X eq $refN and $Y eq $altN){
		$x     = 0;
		$y     = 1;
		$sexdb{'0'} = "X";
		$sexdb{'1'} = "Y";
	}elsif($X eq $altN and $Y eq $refN){
		$x     = 1;
		$y     = 0;
	  $sexdb{'1'} = "X";
		$sexdb{'0'} = "Y";
	}else{
		$x     = "NA";
		$y     = "NA";
	  $sexdb{'0'} = "-";
		$sexdb{'1'} = "-";
		}
	print OUT "$Bid	$chrn	$posi	$refN	$altN	$P1|$sexdb{$P1}	$P2|$sexdb{$P2}	$X	$Y	$x	$y\n";
	if(exists($BLOCKdb{$Bid})){
		$posidb{$RBid}->{'ONAME'} = $Bid;
		$posidb{$RBid}->{'POSI'} .= $posi." ";
		$P1XYdb{$RBid}->{'X'}++  if($sexdb{$P1} eq "X");
		$P1XYdb{$RBid}->{'Y'}++  if($sexdb{$P1} eq "Y");		
	}else{
		$RBid++;
		$posidb{$RBid}->{'ONAME'} = $Bid;
		$posidb{$RBid}->{'POSI'} .= $posi." ";
		$P1XYdb{$RBid}->{'X'}++  if($sexdb{$P1} eq "X");
		$P1XYdb{$RBid}->{'Y'}++  if($sexdb{$P1} eq "Y");	
		$BLOCKdb{$Bid}++;	
		}

	}
close IN;
close OUT;

open(OUT, "> BLOCK.info") or die"";
foreach my $id (sort {$a<=>$b} keys %posidb){
	my @posidb = split(/\s+/,$posidb{$id}->{'POSI'});
	my $oname  = $posidb{$id}->{'ONAME'};
	my $pa     = $posidb[0];
	my $pb     = $posidb[-1];
	my $numX   = (exists($P1XYdb{$id}->{'X'}))?$P1XYdb{$id}->{'X'}:0;
	my $numY   = (exists($P1XYdb{$id}->{'Y'}))?$P1XYdb{$id}->{'Y'}:0;
	my $tNum   = $numX+$numY;
	   $tNum   = 0.1 if($tNum==0);
	my $Xratio = sprintf("%.2f",$numX/$tNum);
	my $Yratio = sprintf("%.2f",$numY/$tNum);
  next if($Xratio<0.7 and $Yratio<0.7);
	print OUT ">$id	$oname	X=HA	Y=HB\n" if($Xratio>=0.7);
	print OUT ">$id	$oname	X=HB	Y=HA\n" if($Yratio>=0.7);
	print OUT "Block position:	$pa	$pb\n";
	print OUT "Number of SNPs that support P1 as X:	$numX\n";
	print OUT "Number of SNPs that support P1 as Y:	$numY\n";
	
	}
close OUT;

%sexdb = ();
open(IN, "grep '>' BLOCK.info|") or die"";
while(<IN>){
	chomp;
	my ($Bid,$X, $Y) = (split/\s+/,$_)[1,2,3];
	$sexdb{$Bid}->{'X'} = $X;
	$sexdb{$Bid}->{'Y'} = $Y;
	}
close IN;

my $Xlist = ();
my $Ylist = ();
open(OUT, "> sex-phasedBLOCKs.csv") or die"";
open(IN, "tmp.txt") or die"";
<IN>;
print OUT "ID,CHRN,POSI,REF,ALT,X,Y,X-reads,Y-reads\n";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $Bid  = $data[0];
	my $CHRN = $data[1];
	my $POSI = $data[2];
	my $REFN = $data[3];
	my $ALTN = $data[4];
	my $HA   = $data[5];
	my $HB   = $data[6];
	my $X    = (exists($sexdb{$Bid}->{'X'}))?$sexdb{$Bid}->{'X'}:"NA";
	my $Y    = (exists($sexdb{$Bid}->{'Y'}))?$sexdb{$Bid}->{'Y'}:"NA";
	if($X=~/HA/ and $Y=~/HB/){
		$HA    = $HA.":".$X;
		$HB    = $HB.":".$Y;
	}elsif($X=~/HB/ and $Y=~/HA/){
		$HA    = $HA.":".$Y;
		$HB    = $HB.":".$X;
		}
	my @hadb = split("",$HA);
	my @hbdb = split("",$HB);
	next if(@hadb!=8);
	my $ha   = $hadb[0]."|".$hadb[4];
	my $hb   = $hbdb[0]."|".$hbdb[4];
	my $x    = ($ha=~/X/)?$ha:$hb;
	my $y    = ($hb=~/Y/)?$hb:$ha;
	my $key  = $CHRN.",".$POSI;
	my $xgeno = (split/\|/,$x)[0];
	my $ygeno = (split/\|/,$y)[0];
	my $xreads = $readb{$key}->{$xgeno};
	my $yreads = $readb{$key}->{$ygeno};
	   $xreads = "NA" if(!defined($xreads));
	   $yreads = "NA" if(!defined($yreads));
	print OUT "$Bid,$CHRN,$POSI,$REFN,$ALTN,$x,$y,$xreads,$yreads\n";
	$xreads =~ s/NA//g;
	$yreads =~ s/NA//g;
	$Xlist .= $xreads."|";
	$Ylist .= $yreads."|";
  }
close IN;
close OUT;

open(OUT, "|sort -u > x-reads.list") or die"";
$Xlist =~ s/\|/\n/g;
print OUT "$Xlist\n";
close OUT;
open(OUT, "|sort -u > y-reads.list") or die"";
$Ylist =~ s/\|/\n/g;
print OUT "$Ylist\n";
close OUT;

