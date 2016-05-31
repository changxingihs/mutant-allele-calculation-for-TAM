#!/usr/bin/perl
use warnings; use strict;

my($list,$bam,$small,$large) = @ARGV;
die " Usage:
	$0 <snp_list> <bam_file> <start_site> <end_site>
	out put file will be :
		<bam_file>.types.count.txt
		<bam_file>.read.count.txt\n" unless @ARGV;

my %list;
open LI, $list or die $!;
while(<LI>){
	chomp;
	my($chr,$pos,$ref) = split /\t/,$_;
	$list{$pos} = 1;
}

my ($s,$e) = (sort {$a <=> $b} (keys %list))[0,-1];
close LI;

$small = $small?$small:$s-5;
$large = $large?$large:$e+5;

open BAM, "samtools view $bam | " or die $!;
open TYPES, ">$bam.reads.class.txt" or die $!;

my %hash;
my $total;
while (<BAM>){
	chomp;
	my $rd = $_;
	my @ar = split /\t/, $_;
	my($id,$pos,$cig,$seq,$md) = @ar[0,3,5,9,12];
	$pos --;
	my $pos_e;
	while ($cig =~ /(\d+)[MD]/g){
		$pos_e += $1;
	}
	next unless  ($pos <= $small  and $pos_e >= $large);	
	$total ++;
# transfer read sequence to reference sequence 	
	my $adds = 0;
	my %mus;
	while ($cig =~ /(\d+)([MDI])/g){
		my $sl = $1;
		my $tag = $2;
		if($tag =~ /D/){
			my $ma = "N" x $sl;
			substr($seq,$adds,0,$ma);
		}elsif($tag =~ /I/){
			my $ins = substr($seq,$adds,$sl,"");
			$mus{$adds} = "^$ins";
		}
		$adds += $sl;
	}
	my $lps;
	while($md =~ /(\d+)([ATCG\^]+)/g){
		$lps += $1;
		my $ref = $2;
		
		if($ref !~ /\^/){
			my $alt = substr($seq,$lps,length($ref));
			$lps += length($ref);
			$mus{$lps}= "$ref->$alt"
		}else{
			#print "jiang\n";
			$ref =~ s/\^/\-/;
			$mus{$lps} = $ref;
			$lps += (length($ref) - 1);
			
		}
	}

	my $bk;
	foreach my $k (sort {$a<=>$b} keys %mus){
		my $v = $mus{$k};
		my $p = $pos + $k;
		next unless ($list{$p});
		$bk .= "$p=$v;";
	}
	if($bk){
		$hash{$bk}{num} ++;
		print TYPES "$bk\t$rd\n";
	}else{
		print TYPES "Zero\t$rd\n";
	}
}
open COUNT, ">$bam.types.count.txt" or die $!;
print COUNT "Zero mutation\t0\t$total\n";
for my $bk (sort keys %hash){
	my $c = ($bk =~ tr/=/=/);
	my $num =  $hash{$bk}{num};
	print COUNT "$bk\t$c\t$num\n";
}
