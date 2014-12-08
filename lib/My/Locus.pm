# Locus object that holds the data for each sequence
package My::Locus
{
    #use namespace::autoclean;
    use Moose;

    has 'accession_number' => (is => 'ro', isa => 'Str', required => 1);
    has 'species' => (is => 'ro', isa => 'Str', required => 1);
    has 'dna' => (is => 'ro', isa => 'Str', required => 1);

    use overload q("") => sub {
        my $self = shift;
        my $accession_number = $self->accession_number;
        my $species = $self->species;
        my $dna = $self->dna;
        return "accession_number: $accession_number";
        #return "Locus:\n\taccession_number: $accession_number\n\tspecies: $species\n\tdna: $dna";
    };

    __PACKAGE__->meta->make_immutable; # allows Moose to optimize construction
}