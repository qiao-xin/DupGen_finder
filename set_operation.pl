if($#ARGV<3)
{
print "Usage: perl set_operation.pl gene_pairs1 gene_pairs2 output operator\n";
}
%h=();
open(input,"$ARGV[0]");
while($line=<input>)
{
$h{$line}=1;
}
open(input,"$ARGV[1]");
open(output,">$ARGV[2]");
if($ARGV[3] eq "+")
{
while($line=<input>)
{
$h{$line}=0;
}
foreach $key (keys %h)
{
print output $key;
}
}
if($ARGV[3] eq "-")
{
while($line=<input>)
{
if(exists $h{$line})
{
$h{$line}=0;
}
}
foreach $key (keys %h)
{
if($h{$key}==1)
{
print output $key;
}
}
}
