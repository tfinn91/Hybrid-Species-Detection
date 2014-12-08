use strict;
use warnings;

#use Bio::AlignIO;
#use Bio::Align::DNAStatistics;
use Bio::Tree::Tree;
use Bio::Tree::Node;
use lib './lib';
use My::Locus;

# TODO: Come up with a better name for this file.
# TODO: Add the xls output?
# TODO: Use references for efficiency

$| = 1; # Turn on stdout flushing

# --------------------
# Record Processing
# --------------------

# Input Record Separator
sub IRS {
    my ($fg) = @_;
    $/ = "//\n\n";
    my $record = <$fg>;
    return $record;
}

# separate splits each record into the individual parts we want to use
sub separate {
    my ($record) = @_;
    my ( $annotation, $dna ) = ( $record =~ /(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n/s ); # full annotation and dna
    $dna =~ s/\s//g;
    $dna =~ s/\d//g;
    my ($accession_number)  = ( $annotation =~ /ACCESSION\s*(\S*)/ );   # accession number
    my ($species) = ( $annotation =~ /ORGANISM\s*(.*)/ );     # Species name

    return ( $accession_number, # The accession number for the sequence in GenBank
        $species, # Species name
        $dna # DNA sequence
        );
}

sub get_data {
    my ($file_name) = @_;

    my $file_handle;
    open ($file_handle, $file_name);

    if (!$file_handle) {
        die "The file $file_name could not be opened. Make sure it exists.";
    }

    my @loci;
    my %species_index_lists;
    my $current_index = 0;

    while ( my $record = IRS($file_handle) ) {
        my ($accession_number, $species, $dna) = separate( $record );

        # Read data into locus objects
        my $locus = My::Locus->new(
                accession_number => $accession_number,
                species => $species,
                dna => $dna
            );


        # push objects onto an array
        # (these will be their positions in the matrix as well)
        push ( @loci, $locus );

        # Build lists of indices for each species so we can properly
        # pair species when we build the three-way trees
        # (1 list of indices per species)
        if (exists $species_index_lists{$species}) {
            push ( $species_index_lists{$species}, $current_index );
        }
        else {
            $species_index_lists{$species} = [ $current_index ];
        }

        $current_index++;
    }

    my @ks = keys %species_index_lists;

    return (\@loci, \%species_index_lists);
}



# --------------------
# Sequence Processing
# --------------------

# diff counts the number of indices where characters differ between sequences.
# TODO: Are both sequences always the same length?
sub diff {
    my ( $string1, $string2 ) = @_;
    my $length = length($string1);
    my $count  = 0;
    for ( my $position = 0 ; $position < $length ; ++$position ) {
        if (
            substr( $string1, $position, 1 ) ne
            substr( $string2, $position, 1 ) )
        {
            ++$count;
        }
    }
    return $count;
}

# returns a matrix containing the alignment percentage for each sequence pair
sub compute_alignments {
    my ($loci_ref) = @_;
    my @matrix;
    for (my $i = 0; $i < scalar @$loci_ref; ++$i) {
        # j = i is a always a 100% comparison
        $matrix[$i][$i] = 1.0;
        for (my $j = $i + 1; $j < scalar @$loci_ref; ++$j) {
            # for now we just use the length of one of the
            # sequences, since all of the sequences are the
            # same length in our sample data, eventually we
            # should reconsider how we do this.
            my $length = length( $loci_ref->[$i]->dna );
            my $bpDiffs = diff( $loci_ref->[$i]->dna, $loci_ref->[$j]->dna );
            my $alignment_percentage = ($length - $bpDiffs) / $length;
            $matrix[$i][$j] = $alignment_percentage;
            $matrix[$j][$i] = $alignment_percentage;
        }
    }
    return \@matrix;
}

# --------------------
# Tree Analysis
# --------------------

# Returns a three-way tree built from the sequences identified
# by the provided indices in the provided similarity matrix.
# The tree will group nodes by creating edges between the two
# pairs of sequences with the most similarity.
# Each node's id stores the index associated with its sequence
# in the matrix.
sub build_tree {
    my ($indices_ref, $matrix_ref) = @_;
    my $tree = Bio::Tree::Tree->new();

    # Node indices:
    my $A1 = $indices_ref->[0];
    my $A2 = $indices_ref->[1];
    my $B = $indices_ref->[2];

    # Inter-species edges:
    my $A1B = $matrix_ref->[$A1][$B];
    my $A2B = $matrix_ref->[$A2][$B];

    # Cut the lower-valued inter-species edge:
    if ($A1B < $A2B) {
        my $L = Bio::Tree::Node->new( -id => $A1 );
        my $R = Bio::Tree::Node->new( -id => $B );
        my $root = Bio::Tree::Node->new( -id => $A2, -descendents => [$L, $R] );
        $tree->set_root_node($root);
    }
    else { # If A1 and A2 align equally to B, then it is an instance of hybridization no matter which tree we build
        my $L = Bio::Tree::Node->new( -id => $A2 );
        my $R = Bio::Tree::Node->new( -id => $B );
        my $root = Bio::Tree::Node->new( -id => $A1, -descendents => [$L, $R] );
        $tree->set_root_node($root);
    }

    return \$tree;
}

# Returns 0 if the provided tree is not an instance of hybridization,
# and 1 if the provided tree is an instance of hybridization.
# This subroutine flags the tree as an instance of hybridization if
# an edge is incident on sequences from different species and the
# alignment percentage for that edge is greater than the alignment
# percentages for any edges in the tree that are incident on sequences
# from the same species. If both edges are between different species,
# the tree will not be flagged as an instance of hybridization.
sub is_hybridization {
    my ($tree_ref, $loci_ref, $matrix_ref) = @_;

    my $root = $$tree_ref->get_root_node;
    my @descendents = $root->get_all_Descendents;
    # Check if the tree is a three-way (if not, return 0)
    if ( scalar @descendents != 2 ) {
        print "Not a 3-way!\n";
        return 0;
    }

    my $L = $descendents[0];
    my $R = $descendents[1];

    my $speciesRoot = $loci_ref->[$root->id]->species;
    my $speciesL = $loci_ref->[$L->id]->species;
    my $speciesR = $loci_ref->[$R->id]->species;

    if ($speciesRoot ne $speciesL && $speciesRoot eq $speciesR) {
        # If the root and L are closer than or equal to the root and R,
        # this is hybridization (different species closer than same).
        if ($matrix_ref->[$root->id][$L->id] >= $matrix_ref->[$root->id][$R->id]) {
            return 1;
        }
    }
    elsif ($speciesRoot eq $speciesL && $speciesRoot ne $speciesR) {
        # If the root and R are closer than or equal to the root and L,
        # this is hybridization (different species closer than same).
        if ($matrix_ref->[$root->id][$R->id] >= $matrix_ref->[$root->id][$L->id]) {
            return 1;
        }
    }

    return 0;
}

# Returns a hash, where the keys are strings of two comma-separated
# species names each (no spaces between names), with the names
# in alphabetical order.
# The values are trees that were flagged as instances of hybridization.
# You can obtain all of the trees that were flagged as instances of hybridization
# between two species like so:
# $returned_hash{"First species name,Second species name"}
sub find_hybrids {
    my ( $loci_ref, $matrix_ref, $species_index_lists_ref ) = @_;

    my %hybridizations;

    # Sort keys alphabetically so they are in the same order every
    # time we run this program, and so the hash can be accessed
    # in a consistent way.
    my @species_list = sort keys %$species_index_lists_ref;

    # Compare sequences for each pair of species. Since we sorted
    # the species names, species_A always comes before species_B
    # in alphabetical order (j > i in all cases).
    for (my $i = 0; $i < scalar @species_list; ++$i) {

        my $species_A = $species_list[$i];
        my $species_A_indices_ref = $$species_index_lists_ref{ $species_A };

        for (my $j = $i + 1; $j < scalar @species_list; ++$j) {

            my $species_B = $species_list[$j];
            my $species_B_indices_ref = $$species_index_lists_ref{ $species_B };

            my @AB_hybridizations;

            # TODO: Handle edge case where there aren't enough sequences
            #       for a three-way comparison



            # First compare triples with 2 sequences from Species A
            # and 1 sequence from Species B:
            for (my $i = 0; $i < scalar @$species_A_indices_ref; ++$i) {
                my $A1 = $species_A_indices_ref->[$i];

                for (my $j = $i + 1; $j < scalar @$species_A_indices_ref; ++$j) {
                    my $A2 = $species_A_indices_ref->[$j];

                    for (my $k = 0; $k < scalar @$species_B_indices_ref; ++$k) {
                        my $B = $species_B_indices_ref->[$k];

                        # Build tree and check if it is
                        # an instance of hybridization.
                        my $tree_ref = build_tree([$A1, $A2, $B], $matrix_ref);

                        if (is_hybridization( $tree_ref, $loci_ref, $matrix_ref )) {
                            push (@AB_hybridizations, $tree_ref);
                        }
                    }
                }
            }

            # Then 2 From Species B and 1 from Species A
            for (my $i = 0; $i < scalar @$species_B_indices_ref; ++$i) {
                my $B1 = $species_B_indices_ref->[$i];

                for (my $j = $i + 1; $j < scalar @$species_B_indices_ref; ++$j) {
                    my $B2 = $species_B_indices_ref->[$j];

                    for (my $k = 0; $k < scalar @$species_A_indices_ref; ++$k) {
                        my $A = $species_A_indices_ref->[$k];

                        # Build tree and check if it is
                        # an instance of hybridization.
                        my $tree_ref = build_tree([$B1, $B2, $A], $matrix_ref);

                        if (is_hybridization( $tree_ref, $loci_ref, $matrix_ref )) {
                            push (@AB_hybridizations, $tree_ref);
                        }
                    }
                }
            }

            $hybridizations{"$species_A,$species_B"} = \@AB_hybridizations;
        }
    }


    return \%hybridizations;

}

# --------------------
# Output
# --------------------

# Prints the tree, labeling each node with the
# accession number of the associated sequence.
sub print_tree {
    my ($tree_ref, $loci_ref, $matrix_ref) = @_;

    my $root = $$tree_ref->get_root_node;

    my @descendents = $root->get_all_Descendents;
    my $L = $descendents[0];
    my $R = $descendents[1];

    my $root_L_alignment = $matrix_ref->[$root->id][$L->id];
    my $root_R_alignment = $matrix_ref->[$root->id][$R->id];
    my $L_R_alignment    = $matrix_ref->[$L->id][$R->id];

    my $root_accession_number = $loci_ref->[$root->id]->accession_number;
    my $L_accession_number    = $loci_ref->[$L->id]->accession_number;
    my $R_accession_number    = $loci_ref->[$R->id]->accession_number;

    print "$root_accession_number --- $root_R_alignment --- $R_accession_number\n";
    print "|\n";
    print "|\n";
    print "|\n";
    print "$root_L_alignment\n";
    print "|\n";
    print "|\n";
    print "|\n";
    print "$L_accession_number                           $L_R_alignment\n"

}

# Prints out the trees in a %hybridizations hash for
# each species pair. The order that the species pairs
# are printed is not guaranteed.
sub print_hybridizations {
    my ($hybridizations_ref, $loci_ref, $matrix_ref) = @_;

    while ( my ($species_pair, $tree_list_ref) = each %$hybridizations_ref ) {
        my $tree_list_size = scalar @$tree_list_ref;
        print "Hybrid trees between $species_pair: $tree_list_size\n\n";
        foreach my $tree_ref (@$tree_list_ref) {
            print_tree( $tree_ref, $loci_ref, $matrix_ref );
            print "\n";
        }
    }

    print "Summary:\n";
    while ( my ($species_pair, $tree_list_ref) = each %$hybridizations_ref ) {
        my $tree_list_size = scalar @$tree_list_ref;
        print "Hybrid trees between $species_pair: $tree_list_size\n";
    }
}


# --------------------
# Main
# --------------------

my ( $loci_ref, $species_index_lists_ref ) = get_data('./smallerlepus.gb');
my $matrix_ref = compute_alignments( $loci_ref );
my $hybridizations_ref = find_hybrids( $loci_ref, $matrix_ref, $species_index_lists_ref );

print "\n";
print_hybridizations( $hybridizations_ref, $loci_ref, $matrix_ref );


