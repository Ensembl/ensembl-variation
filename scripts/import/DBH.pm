use strict;
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
