=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

package DBH;

use DBI;


our $AUTOLOAD;

sub connect {
  my $class = shift;
  my $dsn  = shift;
  my $user = shift;
  my $pass = shift;
  my $opts = shift;

  my $self = bless {}, $class;

  $self->dbh(DBI->connect($dsn, $user, $pass,$opts));

  $dsn =~ s/^\w+:\w+://;

  my @ar = split(';', $dsn);

  foreach my $str (@ar) {
    my ($type, $value) = split(/\s*=\s*/, $str);
    $self->$type($value);
  }

  $self->user($user);
  $self->pass($pass);

  return $self;
}


sub dbh {
  my $self = shift;
  return $self->{'dbh'} = shift if(@_);
  return $self->{'dbh'};
}

sub host {
  my $self = shift;
  return $self->{'hostname'} = shift if(@_);
  return $self->{'hostname'};
}


sub user {
  my $self = shift;
  return $self->{'user'} = shift if(@_);
  return $self->{'user'};
}


sub port {
  my $self = shift;
  return $self->{'port'} = shift if(@_);
  return $self->{'port'};
}


sub pass {
  my $self = shift;
  return $self->{'password'} = shift if(@_);
  return $self->{'password'};
}


sub dbname {
  my $self = shift;
  return $self->{'dbname'} = shift if(@_);
  return $self->{'dbname'};
}


sub AUTOLOAD {
  my $self = shift;
  $AUTOLOAD =~ /(\w+)$/;
  my $method = $1;
  $self->dbh()->$method(@_);
}


sub DESTROY {

}


1;
