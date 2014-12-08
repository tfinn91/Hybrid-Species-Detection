open( $fh, './smallerlepus.gb' );
open( FH1, ">hybridoutput14.xls" );
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

while ( $record = IRS($fh) ) {
    ( $accession_number, $annotation, $species, $dna, $length ) = separate($record);
    $offset = tell($fh);
    push( @annotation, $annotation );
    push( @accession,  $accession_number );
    push( @species,    $species );
    push( @genome,     $dna );
}
#$num    = @annotation;
$countb = scalar @genome; # Count b is the number of dna sequences

$record = IRS($fh);
print "Record: $record\n";



print "annotation list: [@annotation]\n";
 print "accession list: [@accession]\n";
   print "species list: [@species]\n";
    print "genome list: [@genome]\n";


for ( $k = 0 ; $k < ($countb) ; ++$k ) {
    print "countb is $k\n";
    $lowest[$k] = 1;
    $lowacc[$k] = "";
    for ( $l = 0 ; $l < ($countb) ; ++$l ) {
## if they are the same species
        if ( $species[$k] eq $species[$l] ) {
            print "HERES AN INTRA\n";
            ####calculates lowest percent alignment within species

            $leng = ( length( $genome[$k] ) - 2 );
            print "Genome length is $leng\n\n";

            $alignment[$k] = 0;
            $bpdiffs[$k]   = diff( $genome[$k], $genome[$l] );
            $alignment[$k] = ( ( $leng - $bpdiffs[$k] ) / $leng );
            print "ALIGNMENTK WAS $alignment[$k], LOWESTK WAS $lowest[$k]\n\n";
            if ( $alignment[$k] < $lowest[$k] ) {
                $lowest[$k] = $alignment[$k];
                $lowacc[$k] = $accession[$l];
            }
            print
"This is a new lowest: lowest for $species[$k] is: $lowest[$k]\n\n";
## take that individual and compare w. others not same species
        }
    }

    for ( $m = 0 ; $m < ($countb) ; ++$m ) {
        if ( $species[$k] ne $species[$m] ) {

            #$alignmentinter[$l]=0;
            $bpdiffsinter[$m] = diff( $genome[$k], $genome[$m] );
            $alignmentinter[$m] = ( ( $leng - $bpdiffsinter[$m] ) / $leng );

            #@print "This is alignmentinter $alignmentinter[$m]\n";
            if ( $alignmentinter[$m] > $lowest[$k] ) {
                print
"EUREKA! alignment: $alignmentinter[$m] comapring $species[$k] and $species[$m] was greater than the lowest intra for $species[$k] which was $lowest[$k]\n";
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
                print FH1 "\n";
            }
        }
    }
    print
"HYBRID! $species[$k] had inter $alignmentinter[$1] and intra of $lowest[$1]\n";

    #}
    print "species [$k] is done. the lowest intra match was $lowest[$k]\n";
}

#}
## if an alignment with another species is higher than the lowest within species, push

print "done";
print @species;

sub IRS {
    my ($fg) = @_;
    $/      = "//\n\n";
    $record = <$fg>;
    return $record;
}

sub diff {
    ( $string1, $string2 ) = @_;
    $length = length($string1);
    $count  = 0;
    for ( $position = 0 ; $position < $length ; ++$position ) {
        if (
            substr( $string1, $position, 1 ) ne
            substr( $string2, $position, 1 ) )
        {
            ++$count;
        }
    }
    return $count;
}

sub separate {
    my ($record) = @_;
    ( $annotation, $dna ) = ( $record =~ /(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n/s ); #full annotation and dna
    $dna =~ s/\s//g;
    $dna =~ s/\d//g;
    ($accession_number)  = ( $annotation =~ /ACCESSION\s*(\S*)/ );   #accession number
    ($length)  = ( $annotation =~ /LOCUS\s*\S*\s*(\d*)/ );           #sequence length
    ($species) = ( $annotation =~ /ORGANISM\s*(.*)/ );     #Species

    #print "anno: $annotation\n";
    #print "acce: $accession_code\n";
    #print "spec: $species\n";
    #print "dna : $dna\n";
    #print "leng: $length\n";
    return ( $accession_number, # The accession number for the sequence in GenBank
        $annotation, # All text between between LOCUS and ORIGIN\s*\n, inclusive
        $species, # Species name
        $dna, # DNA sequence
        $length ); # Length of the DNA sequence
}

sub search_sequence {
    my ( $sequence, $term ) = @_;
    my (@locations) = ();
    while ( $sequence =~ /$term/isg ) {
        push( @locations, pos );
    }
    return (@locations);
}
