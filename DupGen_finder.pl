#!/usr/bin/perl

use Getopt::Std;

%options=();
getopts("i:t:o:c:d:k:g:s:e:m:w:a:x:", \%options);
if(!exists $options{i} || !exists $options{t} || !exists $options{o} || !exists $options{c})
{
print "Usage: perl DupGen_finder.pl -i data_directory -t target_species -c outgroup_species -o output_directory\n";
print "#####################\nOptional:\n";
print "-a 1 or 0(are segmental duplicates ancestral loci or not? default: 1, yes)\n";
print "-d number_of_genes(maximum distance to call proximal, default: 10)\n";
print "#####################\nThe following are optional MCScanX parameters:\n";
print "-k match_score(cutoff score of collinear blocks for MCScanX, default: 50)\n";
print "-g gap_penalty(gap penalty for MCScanX, default: -1)\n";
print "-s match_size(number of genes required to call a collinear block for MCScanX, default: 5)\n";
print "-e e_value(alignment significance for MCScanX, default: 1e-05)\n";
print "-m max_gaps(maximum gaps allowed for MCScanX, default: 25)\n";
print "-w overlap_window(maximum distance in terms of gene number, to collapse BLAST matches for MCScanX, default: 5)\n";
exit;
}
$options{i}=~s/\/$//;
$options{o}=~s/\/$//;
$good=1;
$t_gff="$options{i}\/$options{t}\.gff";
$t_bla="$options{i}\/$options{t}\.blast";
unless(-e $t_gff)
{
$good=0;
print "Error: Cannot find the gff file for the target species at $t_gff\n";
}
unless(-e $t_bla)
{
$good=0;
print "Error: Cannot find the blast file for the target species at $t_bla\n";
}
@ogrp=split("\,",$options{c});
############################
$n_epoch=1;
###########################
for($i=0;$i<=$#ogrp;$i++)
{
$o_gff="$options{i}\/$options{t}_$ogrp[$i]\.gff";
$o_bla="$options{i}\/$options{t}_$ogrp[$i]\.blast";
unless(-e $o_gff)
{
$good=0;
print "Error: Cannot find the gff file between the target species and the outgroup species \"$ogrp[$i]\" at $o_gff\n";
}
unless(-e $o_bla)
{
$good=0;
print "Error: Cannot find the blast file between the target species and the outgroup species \"$ogrp[$i]\" at $o_bla\n";
}
}
if($good==0)
{
print "##############################################################\n";
print "!!!The execution is terminated due to incorrect input files!!!\n";
print "##############################################################\n";
exit;
}
unless(-d $options{o})
{
mkdir $options{o} or die "Error: Cannot create the output directory\n";
}
@para=("k","g","s","e","m","w");
$mcscanx_para="";
for($i=1;$i<6;$i++)
{
if(exists $options{$para[$i]})
{
$mcscanx_para=$mcscanx_para." -$para[$i] $options{$para[$i]}";
}
}
########################################################################
$in="$options{i}\/$options{t}";
$out="$options{o}\/$options{t}";
system("./MCScanX $in $out $mcscanx_para");
for($i=0;$i<=$#ogrp;$i++)
{
$in="$options{i}\/$options{t}_$ogrp[$i]";
$out="$options{o}\/$options{t}_$ogrp[$i]";
system("./MCScanX $in $out $mcscanx_para");
system("rm $out.gff.sorted");
}
########################################################################
$sorted_gff="$options{o}\/$options{t}.gff.sorted";
open(input,$sorted_gff);
%gid=();
%gch=();
%glc=();
%anc=();
%gmd=();
%gmd2=();#qx
$i=0;
while($line=<input>)
{
chomp($line);
@a=split("\t",$line);
$gid{$a[0]}=$i;
$gch{$a[0]}=$a[1];
$glc{$a[0]}=$a[2];
$gmd{$a[0]}=0;
$gmd2{$a[0]}=0;#qx
$anc{$a[0]}=0;
$i++;
}
close(input);
$cut=10;
if(exists $options{d})
{
$cut=$options{d};
}
#############read blastp of target genome#########################
%bla1=();
%bla2=();
%blae=();
%gmd3=();#qx
open(input,$t_bla);
while($line=<input>)
{
@a=split("\t",$line);
if($a[0] ne $a[1] && exists $gid{$a[0]} && exists $gid{$a[1]})
{
$key="$a[0]\t$a[1]";
if(!exists $bla2{$key})
{
#$bla2{$key}=$a[2];#yupeng wang
#####qx####
#Here, the Score or E-value is better than identity value to be used to filter blastp results
#because Score is correlated to identity and seq_length, while identify is not related to seq-length
$bla2{$key}=$a[11];#qx
$gmd3{$key}=0;#qx
}
else
{
#if($a[2]>$bla2{$key})#yupeng wang
if($a[11]>$bla2{$key})#qx
{
#$bla2{$key}=$a[2];#yupeng
$bla2{$key}=$a[11];#qx
}
}
$blae{$key}=$a[10];
if($gid{$a[0]} gt $gid{$a[1]})
{
$key="$a[1]\t$a[0]";
}
$bla1{$key}=0;
}
}
#######qx#######
#open(OUT,">$options{o}\/$options{t}.bla2");
#foreach $key1 (keys %bla2)
#{
#	print OUT "$key1\t$bla2{$key1}\n";
#	#print OUT $bla2{$key1}, "\n";
#}
######qx#######
#######################WGD/segmental#############################################
$in="$options{o}\/$options{t}.collinearity";
$out1=">$options{o}\/$options{t}.wgd.pairs";
$out2=">$options{o}\/$options{t}.wgd.genes";
$incl=1;
$wp=0;
$wg=0;
if(exists $options{a})
{
$incl=$options{a};
}
%h=();
open(input,$in);
open(output1,$out1);
open(output2,$out2);
print output1 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
while($line=<input>)
{
chomp($line);
if($line!~/\#/ && $line ne "")
{
@a=split("\t",$line);
$t_key="$a[1]\t$a[2]";
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
$gmd3{$t_key}=1;#qx
}
else
{
$t_key="$a[2]\t$a[1]";
$temp=$blae{$t_key};
$gmd3{$t_key}=1;#qx
}
print output1 "$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$a[2]\t$gch{$a[2]}\:$glc{$a[2]}\t$temp\n";
$wp++;#qx
$h{$a[1]}=0;
$h{$a[2]}=0;
if($incl==1)
{
$anc{$a[1]}=1;
$anc{$a[2]}=1;
}
$gmd{$a[1]}=1;
$gmd{$a[2]}=1;
$gmd2{$a[1]}=1;#qx
$gmd2{$a[2]}=1;#qx
}
}
foreach $key (sort(keys %h))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
$wg++;#qx
}
####################tandem&proximal#####################
$out1=">$options{o}\/$options{t}.tandem.pairs";
$out2=">$options{o}\/$options{t}.tandem.genes";
$out3=">$options{o}\/$options{t}.proximal.pairs";
$out4=">$options{o}\/$options{t}.proximal.genes";
open(output1,$out1);
open(output2,$out2);
open(output3,$out3);
open(output4,$out4);
print output1 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
print output3 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output4 "Duplicate\tLocation\n";
%tan_p=();
%tan_d=();
$tp=0;
$tg=0;
$pp=0;
$pg=0;
foreach $key (keys %bla1)
{
@a=split("\t",$key);
$dis= abs($gid{$a[1]}-$gid{$a[0]});
if($gch{$a[0]} eq $gch{$a[1]} && $dis<$cut)
{
for($i=0;$i<2;$i++)
{
if(! exists $tan_p{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
$tan_d{$a[$i]}=$dis;
}
else
{
if($dis<$tan_d{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
$tan_d{$a[$i]}=$dis;
}
elsif($dis==$tan_d{$a[$i]})
{
if($gid{$a[1-$i]}>$gid{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
}
}
else
{
}
}
}
}
}
%tan_pairs=();
foreach $key (keys %tan_p)
{
$newkey=$key."\t".$tan_p{$key};
if($gid{$key} gt $gid{$tan_p{$key}})
{
$newkey=$tan_p{$key}."\t".$key;
}
$tan_pairs{$newkey}=$tan_d{$key};
}
%h=();
%k=();
foreach $key (sort(keys %tan_pairs))
{
@a=split("\t",$key);
if($tan_pairs{$key}==1)
{
$t_key="$a[0]\t$a[1]";
$t_key_r="$a[1]\t$a[0]";#qx
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
$gmd3{$t_key}=2;#qx
$gmd3{$t_key_r}=2;#qx
}
else
{
$t_key="$a[1]\t$a[0]";
$t_key_r="$a[0]\t$a[1]";#qx
$temp=$blae{$t_key};
$gmd3{$t_key}=2;#qx
$gmd3{$t_key_r}=2;#qx
}
print output1 "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$temp\n";
$tp++;#qx
for($i=0;$i<2;$i++)
{
if($gmd{$a[$i]}==0)
{
$h{$a[$i]}=2;
$gmd{$a[$i]}=2;
$gmd2{$a[$i]}=2;#qx
}
}
}
else
{
$t_key="$a[0]\t$a[1]";
$t_key_r="$a[1]\t$a[0]";#qx
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
$gmd3{$t_key}=3;#qx
$gmd3{$t_key_r}=3;#qx
}
else
{
$t_key="$a[1]\t$a[0]";
$t_key_r="$a[0]\t$a[1]";#qx
$temp=$blae{$t_key};
$gmd3{$t_key}=3;#qx
$gmd3{$t_key_r}=3;#qx
}
print output3 "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$temp\n";
$pp++;
for($i=0;$i<2;$i++)
{
if($gmd{$a[$i]}==0 || $gmd{$a[$i]}==3)
{
$k{$a[$i]}=3;
$gmd{$a[$i]}=3;
$gmd2{$a[$i]}=3;#qx
}
}
}
}
foreach $key (sort(keys %h))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
$tg++;#qx
}
foreach $key (sort(keys %k))
{
print output4 "$key\t$gch{$key}\:$glc{$key}\n";
$pg++;#qx
}
####################transposed################
%epoch_pair=();
$trp=0;#qx
$trg=0;#qx
for($k=$n_epoch;$k>=1;$k--)
{
for($i=$#ogrp;$i>=$k-1;$i--)
{
$in="$options{o}\/$options{t}_$ogrp[$i].collinearity";
open(input,$in);
while($line=<input>)
{
chomp($line);
if($line!~/\#/ && $line ne "")
{
@a=split("\t",$line);
for($j=1;$j<3;$j++)
{
if(exists $anc{$a[$j]})
{
$anc{$a[$j]}=1;
}
}
}
}
}
%tran_genes=();
%tran_genes2=();#qx
%tran_ident=();
#open(OUT, ">$options{o}\/$options{t}.tran_genes.tmp");#qx
foreach $key (keys %bla2)
{
$good=0;
@a=split("\t",$key);
if($gmd{$a[0]}==0)
{
if($gch{$a[0]} ne $gch{$a[1]})
{
$good=1;
}
else
{
if(abs($gid{$a[0]}-$gid{$a[1]})>=$cut)
{
$good=1;
}
}
}
if($good==1)
{
if($anc{$a[0]}==0 && $anc{$a[1]}==1)
{
if(!exists $tran_genes{$a[0]})
{
$tran_genes{$a[0]}=$a[1];
$tran_ident{$a[0]}=$bla2{$key};
#print OUT "$a[0]\t$tran_genes{$a[0]}\t$tran_ident{$a[0]}\n";#qx
}
else
{
#print OUT "$a[0]\t$a[1]sign\n"; #qx
#if($a[2]>$tran_ident{$a[0]})#yupeng wang
if($bla2{$key}>$tran_ident{$a[0]})#qx modified
{
$tran_genes{$a[0]}=$a[1];
$tran_ident{$a[0]}=$bla2{$key};
#print OUT "AAA\n";#qx
}
}
}
}
}

#############qx###################
###Note the following situation:reserved these trans_pairs
#AT5G37640	AT4G02890	90.43	303	28	1	3	305	1	302	0.0	 535
#AT5G37640	AT5G03240	90.43	303	28	1	3	305	1	302	0.0	 535
#AT5G37640	AT5G20620	90.13	304	29	1	2	305	76	378	0.0	 535
#AT5G37640	AT4G05320	90.13	304	29	1	2	305	76	378	0.0	 535
#However you must avoid the following situation: because the oreders of keys in a hash is random, so the follow second line may be the first input, if then
#the third line as following input, the codes will store the incorrect gene pairs contained in the third line into %tran_genes2
#In actually, we just want the gene pairs contained in the following first line
#AT1G47580	AT2G15690	46.99	183	86	5	66	239	399	579	7e-47	 163
#AT1G47580	AT3G26782	62.16	111	42	0	129	239	549	659	3e-44	 156
#AT1G47580	AT2G22070	54.40	125	57	0	115	239	662	786	7e-44	 156
#AT1G47580	AT1G25360	59.26	108	43	1	133	239	683	790	1e-43	 156
foreach $key (keys %bla2)
{
$good=0;
@a=split("\t",$key);
if($gmd{$a[0]}==0)
{
if($gch{$a[0]} ne $gch{$a[1]})
{
$good=1;
}
else
{
if(abs($gid{$a[0]}-$gid{$a[1]})>=$cut)
{
$good=1;
}
}
}
if($good==1)
{
if($anc{$a[0]}==0 && $anc{$a[1]}==1)
{
if($bla2{$key}==$tran_ident{$a[0]})
{
	if($tran_genes{$a[0]} ne $a[1])
	{
		$tran_genes2{$key}=$bla2{$key};
	}
}
}
}
}
#############qx####################

###yupeng wang
$out1=">$options{o}\/$options{t}.transposed.pairs";
$out2=">$options{o}\/$options{t}.transposed.genes";
open(output1,$out1);
open(output2,$out2);
#open(output3,$out1."bak"); #qx
print output1 "Transposed\tLocation\tParental\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
foreach $key (sort(keys %tran_genes))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
$trg++;#qx
$gmd2{$key}=4;#qx
if($gmd2{$tran_genes{$key}}==0)#qx
{
$gmd2{$tran_genes{$key}}=4;#qx
}
$temp="$key\t$tran_genes{$key}";
#print output3 "$temp\n";#qx
$gmd3{$temp}=4;#qx
$rtemp="$tran_genes{$key}\t$key";#qx
$gmd3{$rtemp}=4;#qx
print output1 "$key\t$gch{$key}\:$glc{$key}\t$tran_genes{$key}\t$gch{$tran_genes{$key}}\:$glc{$tran_genes{$key}}\t$blae{$temp}\n"; #qx modified
$trp++;#qx
$epoch_pair{$temp}=$k;
}
###yupeng wang

#############qx##############
#open(OUT,">tan_genes2_pair");
foreach $key (sort(keys %tran_genes2))#qx
{
	@a=split("\t",$key);
	$temp="$key";#qx
	$rtemp="$a[1]\t$a[0]";
	$gmd3{$temp}=4;
	$gmd3{$rtemp}=4;
	#print OUT "$temp\n";
	print output1 "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$blae{$temp}\n";#qx
	$trp++;
	if($gmd2{$a[1]}==0)
	{
	$gmd2{$a[1]}=4;
	}
	$epoch_pair{$temp}=$k;#qx
}
#############qx##############
}

##############dispersed##############qx
%disp_genes=();
%hash=();
%hash2=();
%ident=();
%disp_genes2=();
%hash4=();
$mdp=0;
$mdg=0;

# Deducting WGD, tandem, proximal, transposed duplications from the BLASTP output ($bla2)
foreach $key (keys %bla2)
{
	@a=split("\t",$key);
	#$j++;
	#if($gmd2{$a[0]}==0 && $gmd2{$a[1]}==0)
	#if($gmd2{$a[0]} == $gmd2{$a[1]})
	if($gmd3{$key} == 1)
	{
		#$w++;
		next;
	}
	elsif($gmd3{$key} == 2)
	{
		#$t++;
		next;
	}
	elsif($gmd3{$key} == 3)
	{
		#$pr++;
		next;
	}
	elsif($gmd3{$key} == 4)
	{
		#$trp++;
		next;
	}
	elsif($gmd3{$key} == 0)
	{
		#$k++;
		$disp_genes{$key}=$bla2{$key};
		$hash{$key}='A';
		if(!exists $hash2{$a[0]}{$a[1]})
		{
			$hash2{$a[0]}{$a[1]}=$bla2{$key};
		}
	}
	else
	{
		#$l++;
		next;
	}
}

#print "blastpairs: $j\ndispersedpairs: $k\nwgd: $w\ntandem: $t\nproximal: $pr\ntransposed: $trp\notherpairs:$l\n";

$out5=">$options{o}\/$options{t}.dispersed.pairs";
$out6=">$options{o}\/$options{t}.dispersed.genes";

open(output5,$out5);
open(output6,$out6);

print output5 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output6 "Duplicate\tLocation\n";

foreach $key (sort(keys %disp_genes))
{
	@a=split("\t",$key);
	if(exists $hash2{$a[0]}{$a[1]}) #remove repetitive lines like ab--ba
	{
		delete $hash2{$a[1]}{$a[0]};
	}
	else
	{
		next;
	}
}

#------------------ Output the dispersed gene pairs in which some genes are ---------------------------------------------
#------------------ overlapped with other modes of duplicated genes but only reserved the best hit for each gene --------

foreach $key1 (sort(keys %hash2))
{
	foreach $key2 (sort(keys %{$hash2{$key1}}))
	{
		#$score=${$hash2{$key1}}{$key2}
		$pair="$key1\t$key2";
		if(!exists $hash4{$key1})
		{
			$hash4{$key1}=$key2;
			$ident2{$key1}=$bla2{$pair};
		}
		else
		{
			if($bla2{$pair}>$ident2{$key1})# Reserve the best hit for each gene
			{
				$hash4{$key1}=$key2;
				$ident2{$key1}=$bla2{$pair};
			}
		}
	}
}


foreach $key (sort(keys %hash4))
{
if($gmd2{$key} == 0 && $gmd2{$hash{$key}} == 0)
{
	print output6 "$key\t$gch{$key}\:$glc{$key}\n$hash4{$key}\t$gch{$hash4{$key}}\:$glc{$hash4{$key}}\n";
	$mdg++;
	$mdg++;
}
elsif($gmd2{$key} == 0 && $gmd2{$hash{$key}} != 0)
{
	print output6 "$key\t$gch{$key}\:$glc{$key}\n";
	$mdg++;
}
elsif($gmd2{$key} != 0 && $gmd2{$hash{$key}} == 0)
{
	print output6 "$hash4{$key}\t$gch{$hash4{$key}}\:$glc{$hash4{$key}}\n";
	$mdg++;
}
$temp="$key\t$hash4{$key}";
print output5 "$key\t$gch{$key}\:$glc{$key}\t$hash4{$key}\t$gch{$hash4{$key}}\:$glc{$hash4{$key}}\t$blae{$temp}\n";
$mdp++;
}

##############Singletons##############qx
$sg=0;
$out=">$options{o}\/$options{t}.singletons";
open(output,$out);
print output "GeneID\tLocation\n";

foreach $key (keys %gmd2)
{
	if($gmd2{$key} == 0)
	{
		print output "$key\t$gch{$key}\:$glc{$key}\n";
		$sg++;
	}
}

#--------------------Output Statistic files for different modes of gene duplications--------------------qx
#Note: for diespersed duplication, the MQD dataset is used in here.
$out1=">$options{o}\/$options{t}.pairs.stats";
#$out2=">$options{o}\/$options{t}.genes.stats";
open(output1,$out1);
#open(output2,$out2);
print output1 "Types\tNumber of gene pairs\n";

print output1 "WGD-pairs\t$wp\nTD-pairs\t$tp\nPD-pairs\t$pp\nTRD-pairs\t$trp\nDSD-pairs\t$mdp\n";
#print output1 "WGD-genes\t$wg\nTD-genes\t$tg\nPD-genes\t$pg\nTRD-genes\t$trg\nDSD-genes\t$mdg\nSingleton-genes\t$sg\n";

#print output2 "Types\tWGD-pairs\tTD-pairs\tPD-pairs\tTRD-pairs\tDD-pairs\n";
#print output2 "$options{t}\t$wp\t$tp\t$pp\t$trp\t$mdp\n";
close output1;

__END__