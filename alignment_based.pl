use strict;
use warnings;

use Bio::Tools::Run::Alignment::Clustalw;

use Bio::SeqIO;

use Bio::TreeIO;

# TODO: Come up with a better name for this file.
# TODO: Add the xls output?

$| = 1; # Turn on stdout flushing

# --------------------
# Record Processing
# --------------------

sub get_data {
    my ($file_name) = @_;

    my $file_handle;
    open ($file_handle, $file_name);

    if (!$file_handle) {
        die "The file $file_name could not be opened. Make sure it exists.";
    }

    my %loci; # hash of loci by accession number
    my %species_lists; # hash of lists of accession numbers by species name

    my $in  = Bio::SeqIO->new(-file => "smallerlepus.gb" ,
                         -format => 'genbank');

    while ( my $seq = $in->next_seq() ) {

        my $accession_number = $seq->accession_number;
        my $species = $seq->species->binomial;

        $loci{$accession_number} = $seq;

        # Build lists of indices for each species so we can properly
        # pair species when we build the three-way trees
        # (1 list of indices per species)

        print "$species\n";

        if (exists $species_lists{$species}) {
            push ( $species_lists{$species}, $accession_number );
        }
        else {
            $species_lists{$species} = [ $accession_number ];
        }

    }

    return (\%loci, \%species_lists);
}



# --------------------
# Sequence Processing
# --------------------


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
    my ($accessions_ref, $loci_ref, $factory_ref) = @_;

    # Indices:
    my $A1_accession = $accessions_ref->[0];
    my $A2_accession = $accessions_ref->[1];
    my $B_accession = $accessions_ref->[2];

    # Sequences:
    my $A1 = $$loci_ref{ $A1_accession };
    my $A2 = $$loci_ref{ $A2_accession };
    my $B  = $$loci_ref{ $B_accession  };

    # Build alignent object:
    my $alignment = $$factory_ref->align([$A1, $A2, $B]);
    my $tree      = $$factory_ref->tree($alignment);

    return \$tree;
}

# Returns 0 if the provided tree is not an instance of hybridization,
# and 1 if the provided tree is an instance of hybridization.
# Assumes that two of the three descendents of the root are of the same
# species, and that the third is some different species.
sub is_hybridization {
    my ($tree_ref, $loci_ref) = @_;

    my $root = $$tree_ref->get_root_node;
    my @descendents = $root->get_all_Descendents;

    my ($A_accession) = $descendents[0]->id =~ /(.*)\//;
    my ($B_accession) = $descendents[1]->id =~ /(.*)\//;
    my ($C_accession) = $descendents[2]->id =~ /(.*)\//;

    my $A_species = $$loci_ref{$A_accession}->species->binomial;
    my $B_species = $$loci_ref{$B_accession}->species->binomial;
    my $C_species = $$loci_ref{$C_accession}->species->binomial;

    my $A_dist = $descendents[0]->branch_length;
    my $B_dist = $descendents[1]->branch_length;
    my $C_dist = $descendents[2]->branch_length;

    if ($A_species ne $B_species && $A_species ne $C_species) {
        # A is the odd one out
        my $intra_dist = abs($B_dist - $C_dist);
        my $inter_dist_1 = abs($A_dist - $B_dist);
        my $inter_dist_2 = abs($A_dist - $C_dist);
        if ($inter_dist_1 < $intra_dist || $inter_dist_2 < $intra_dist) {
            return 1;
        }

    }
    elsif ($B_species ne $A_species && $B_species ne $C_species) {
        # B is the odd one out
        my $intra_dist = abs($A_dist - $C_dist);
        my $inter_dist_1 = abs($B_dist - $A_dist);
        my $inter_dist_2 = abs($B_dist - $C_dist);
        if ($inter_dist_1 < $intra_dist || $inter_dist_2 < $intra_dist) {
            return 1;
        }
    }
    else {
        # C is the odd one out
        my $intra_dist = abs($A_dist - $B_dist);
        my $inter_dist_1 = abs($C_dist - $A_dist);
        my $inter_dist_2 = abs($C_dist - $B_dist);
        if ($inter_dist_1 < $intra_dist || $inter_dist_2 < $intra_dist) {
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
    my ( $loci_ref, $species_lists_ref ) = @_;

    my %hybridizations;

    # Sort keys alphabetically so they are in the same order every
    # time we run this program, and so the hash can be accessed
    # in a consistent way.
    my @species_list = sort keys %$species_lists_ref;

    my $fh;
    open($fh, '>output');
    my $out = Bio::TreeIO->new(-fh => $fh, -format => 'newick');

    # Factory to build alignments using clustalw
    # ktuple = word size to be used in the alignment
    # TODO: add params for how trees get built
    my @params = (   'ktuple' => 2
                    ,'matrix' => 'BLOSUM'
                    ,'quiet' => 1
                    );
    my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

    # Compare sequences for each pair of species. Since we sorted
    # the species names, species_A always comes before species_B
    # in alphabetical order (j > i in all cases).
    for (my $i = 0; $i < scalar @species_list; ++$i) {

        my $species_A = $species_list[$i];
        my $species_A_accessions_ref = $$species_lists_ref{ $species_A };

        for (my $j = $i + 1; $j < scalar @species_list; ++$j) {

            my $species_B = $species_list[$j];
            my $species_B_accessions_ref = $$species_lists_ref{ $species_B };

            my @AB_hybridizations;

            # TODO: Handle edge case where there aren't enough sequences
            #       for a three-way comparison


            # First compare triples with 2 sequences from Species A
            # and 1 sequence from Species B:
            for (my $i = 0; $i < scalar @$species_A_accessions_ref; ++$i) {
                my $A1 = $species_A_accessions_ref->[$i];

                for (my $j = $i + 1; $j < scalar @$species_A_accessions_ref; ++$j) {
                    my $A2 = $species_A_accessions_ref->[$j];

                    for (my $k = 0; $k < scalar @$species_B_accessions_ref; ++$k) {
                        my $B = $species_B_accessions_ref->[$k];

                        # Build tree and check if it is
                        # an instance of hybridization.

                        # TODO: Change to using clustalw methods

                        my $tree_ref = build_tree([$A1, $A2, $B], $loci_ref, \$factory);

                        

                        if (is_hybridization( $tree_ref, $loci_ref )) {
                            print "hybrid!\n";
                            push (@AB_hybridizations, $tree_ref);

                            print $fh "Accession #s:\n";
                            my $acc_num = $$loci_ref{$A1}->accession_number;
                            print $fh "A1: $acc_num\n";
                            $acc_num = $$loci_ref{$A2}->accession_number;
                            print $fh "A2: $acc_num\n";
                            $acc_num = $$loci_ref{$B}->accession_number;
                            print $fh "B: $acc_num\n";
                            $out->write_tree($$tree_ref);
                            print $fh "\n";
                        }
                        else {
                            print "not hybrid!\n";
                        }
                    }
                }
            }

            # Then 2 From Species B and 1 from Species A
            for (my $i = 0; $i < scalar @$species_B_accessions_ref; ++$i) {
                my $B1 = $species_B_accessions_ref->[$i];

                for (my $j = $i + 1; $j < scalar @$species_B_accessions_ref; ++$j) {
                    my $B2 = $species_B_accessions_ref->[$j];

                    for (my $k = 0; $k < scalar @$species_A_accessions_ref; ++$k) {
                        my $A = $species_A_accessions_ref->[$k];

                        # Build tree and check if it is
                        # an instance of hybridization.

                        # TODO: Change to using clustalw methods

                        my $tree_ref = build_tree([$B1, $B2, $A], $loci_ref, \$factory);

                        if (is_hybridization( $tree_ref, $loci_ref )) {
                            print "hybrid!\n";
                            push (@AB_hybridizations, $tree_ref);

                            print $fh "Accession #s:\n";
                            my $acc_num = $$loci_ref{$B1}->accession_number;
                            print $fh "B1: $acc_num\n";
                            $acc_num = $$loci_ref{$B2}->accession_number;
                            print $fh "B2: $acc_num\n";
                            $acc_num = $$loci_ref{$A}->accession_number;
                            print $fh "A: $acc_num\n";
                            $out->write_tree($$tree_ref);
                            print $fh "\n";
                        }
                        else {
                            print "not hybrid!\n";
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
    my ($hybridizations_ref, $loci_ref) = @_;

    while ( my ($species_pair, $tree_list_ref) = each %$hybridizations_ref ) {
        my $tree_list_size = scalar @$tree_list_ref;
        print "Hybrid trees between $species_pair: $tree_list_size\n\n";
        foreach my $tree_ref (@$tree_list_ref) {
            print_tree( $tree_ref, $loci_ref );
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

my ( $loci_ref, $species_lists_ref ) = get_data('./smallerlepus.gb');
my $hybridizations_ref = find_hybrids( $loci_ref, $species_lists_ref );

print "\n";
#print_hybridizations( $hybridizations_ref, $loci_ref, $matrix_ref );


