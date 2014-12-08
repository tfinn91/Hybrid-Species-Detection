open($fh,'smallerlepus.gb');
open(FH1,">hybridoutput14.xls");
print FH1 "AccessionA";
print FH1 "\t";
print FH1 "SpeciesA";
print FH1 "\t";
print FH1 "\t";
print FH1 "AccessionB";
print FH1 "\t";
print FH1 "SpeciesB";
print FH1 "\t";
print FH1 "\t";
print FH1 "Hyb.Alignment";
print FH1 "\t";
print FH1 "Lowest Intra Alignment";
print FH1 "\t";
print FH1 "Lowest Intra Alignment accession #";
print FH1 "\n";

while ($record=IRS($fh)){
($acc1, $anno1, $species1, $dna1, $len1)= separate($record);
$offset=tell($fh);
push(@annotation, $anno1);
push(@accession,$acc1);
push(@species,$species1);
push(@genome, $dna1);
}
$num=@annotation;
$countb=scalar @genome;


for($k=0;$k<($countb);++$k){
print "countb is $k\n";
$lowest[$k]=1;
$lowacc[$k]="";
for($l=0;$l<($countb);++$l){
## if they are the same species
if ($species[$k] eq $species[$l]){
print "HERES AN INTRA\n";
 ####calculates lowest percent alignment within species

$leng=(length($genome[$k])-2);
print "Genome length is $leng\n\n";

$alignment[$k]=0;
$bpdiffs[$k]=diff($genome[$k],$genome[$l]);
$alignment[$k]=(($leng-$bpdiffs[$k])/$leng);
print "ALIGNMENTK WAS $alignment[$k], LOWESTK WAS $lowest[$k]\n\n";
if ($alignment[$k]<$lowest[$k]){
$lowest[$k]= $alignment[$k];
$lowacc[$k]=$accession[$l];}
print "This is a new lowest: lowest for $species[$k] is: $lowest[$k]\n\n";
## take that individual and compare w. others not same species
}
}

for($m=0;$m<($countb);++$m){
if ($species[$k] ne $species[$m]){
#$alignmentinter[$l]=0;
$bpdiffsinter[$m]=diff($genome[$k], $genome[$m]);
$alignmentinter[$m]=(($leng-$bpdiffsinter[$m])/$leng);
#@print "This is alignmentinter $alignmentinter[$m]\n";
if ($alignmentinter[$m]>$lowest[$k]){
print "EUREKA! alignment: $alignmentinter[$m] comapring $species[$k] and $species[$m] was greater than the lowest intra for $species[$k] which was $lowest[$k]\n";
print FH1 $accession[$k];
print FH1 "\t";
print FH1 $species[$k];
print FH1 "\t";
print FH1 "\t";
print FH1 $accession[$m];
print FH1 "\t";
print FH1 $species[$m];
print FH1 "\t";
print FH1 "\t";
print FH1 $alignmentinter[$m];
print FH1 "\t";
print FH1 "\t";
print FH1 $lowest[$k];
print FH1 "\t";
print FH1 "\t";
print FH1 $lowacc[$k];
print FH1 "\n";;}
}
}
print "HYBRID! $species[$k] had inter $alignmentinter[$1] and intra of $lowest[$1]\n";

#}
print "species [$k] is done. the lowest intra match was $lowest[$k]\n";
}
#}
## if an alignment with another species is higher than the lowest within species, push

print "done";
print @species;

sub IRS{
my($fg)=@_;
$/="//\n\n";
$record=<$fg>;
return $record;
}
sub diff{
($string1,$string2)=@_;
$length=length($string1);
$count=0;
for($position=0;$position<$length;++$position){
if(substr($string1,$position,1) ne substr($string2,$position,1)){
++$count;
}
}
return $count;
}
sub separate{
my($record)=@_;
($anno, $dna)=($record =~/(LOCUS.*ORIGIN\s*\n)(.*\/\/\n)/s); #full annotation and dna
($def)=($record =~ /DEFINITION(.*)ACCESSION.*ORIGIN\s*\n.*\/\/\n/s); #full org def
($cou)=($record=~/country="(.*)"\s*\/collection.*\/\/\n/s); # country recorded
($year)=($record=~/collection_date="(\N*)".*\/\/\n/s); #year recorded
($acc)=($record=~/ACCESSION\s*(\S*)\n.*\/\/\n/s); #accession number
($len)=($record=~/LOCUS\s*\S*\s*(\d*)/s); #sequence length
($spec)=($record=~/ORGANISM\s*(\S*\s\S*)\n.*/s); #Species
$dna=~s/\s//g;
$dna=~s/\d//g;
$def=~s/\n//g;
return($acc, $anno, $spec, $dna, $len);
}
sub search_sequence {
my ($sequence, $term)=@_;
my (@locations)=();
while ($sequence =~ /$term/isg){
push (@locations, pos);
}
return(@locations);
}