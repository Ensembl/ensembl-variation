#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<helpdesk.org>.

=cut


use DBI;
use Getopt::Long;
use ImportUtils qw(load);

my $config = {};

my @exclude = qw(transcript_variation_id hgvs_genomic hgvs_protein hgvs_transcript somatic codon_allele_string);

GetOptions(
    $config,
    'tmpdir=s',
    'tmpfile=s',
    'host|h=s',
    'user|u=s',
    'port|P=i',
    'password|p=s',
    'db|d=s',
    'table|t=s',
    'version|v=i',
    'pattern=s',
    ) or die "ERROR: Could not parse command line options\n";

# check options
for(qw(host user port password tmpdir tmpfile)) {
  die("ERROR: $_ not defined, use --$_\n") unless defined $config->{$_};
}

my @db_list;

if(!defined($config->{db})) {
  @db_list = @{get_species_list($config)};
}
else {
  push @db_list, $config->{db};
}

die "ERROR: no suitable databases found on host ".$config->{host}."\n" unless scalar @db_list;

if (lc($config->{table}) eq 'all') {
  $config->{table} = 'transcript_variation motif_feature_variation regulatory_feature_variation';
} else {
  $config->{table} ||= 'transcript_variation';
}
my @exclude;

my $table_input = $config->{table};

for my $source_table (split(/\s+/, $table_input)) {
  my $filtered_db_list = filter_db_list($config, \@db_list, $source_table);    

  if ($source_table eq 'transcript_variation') {
    @exclude = qw(transcript_variation_id hgvs_genomic hgvs_protein hgvs_transcript somatic codon_allele_string);
  }
  if ($source_table eq 'regulatory_feature_variation') {
    @exclude = qw(regulatory_feature_variation_id feature_type somatic);
  }
  if ($source_table eq 'motif_feature_variation') {
    @exclude = qw(motif_feature_variation_id motif_feature_id somatic motif_end);
  }

  $table = 'MTMP_'. $source_table;

  my $TMP_DIR = $config->{tmpdir};
  my $TMP_FILE = $config->{tmpfile};

  $ImportUtils::TMP_DIR = $TMP_DIR;
  $ImportUtils::TMP_FILE = $TMP_FILE;

  foreach my $db (@$filtered_db_list) {
    print "\nProcessing database $db\n";
    my $dbc = DBI->connect(
        sprintf(
          "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
          $config->{host},
          $config->{port},
          $db,
          ), $config->{user}, $config->{password}
        );

    print "Getting column definition\n";

    # get column definition from transcript_variation
    my $sth = $dbc->prepare(qq{
        SHOW CREATE TABLE $source_table
        });
    $sth->execute();

    my $create_sth = $sth->fetchall_arrayref->[0]->[1];
    $sth->finish;

    # convert set to enum
    $create_sth =~ s/^set/enum/;

    # rename table
    $create_sth =~ s/TABLE \`$source_table\`/TABLE \`$table\`/;

    # filter out some columns
    $create_sth =~ s/\`?$_.+?,// for @exclude;

    # filter out some indices
    $create_sth =~ s/AUTO_INCREMENT=\d+//;
    $create_sth =~ s/somatic_feature_idx/feature_idx/;
    $create_sth =~ s/$_.+,// for ('PRIMARY KEY', 'KEY `somatic', 'KEY `cons');
    $create_sth =~ s/,\`somatic\`//;

    # remove final comma
    $create_sth =~ s/,(\s+\))/$1/;

    print "Creating table $table\n";

    # create a new table
    $dbc->do($create_sth);

    # get columns of new table
    $sth = $dbc->prepare(qq{
        DESCRIBE $table
        });
    $sth->execute();
    my @cols = map {$_->[0]} @{$sth->fetchall_arrayref};
    $sth->finish;

    print "Dumping data from $source_table to $TMP_DIR\/$TMP_FILE\n";

    $sth = $dbc->prepare(qq{ SELECT count(*) FROM $source_table });
    $sth->execute();
    my $row_count = $sth->fetchall_arrayref->[0]->[0];
    $sth->finish;

    # populate it
    $sth = $dbc->prepare(qq{ SELECT * FROM $source_table
    }, {mysql_use_result => 1});

    open OUT, ">$TMP_DIR/$TMP_FILE";

    $sth->execute();

    while(my $row = $sth->fetchrow_hashref()) {
      my $cons = $row->{consequence_types};

      foreach my $con(split /\,/, $cons) {
        my %vals = %{$row};
        $vals{consequence_types} = $con;
        delete $vals{$_} for @exclude;
        my @values = map {defined $vals{$_} ? $vals{$_} : '\N'} @cols;

        print OUT join("\t", @values);
        print OUT "\n";

      }
    }

    $sth->finish;

    close OUT;
    print "Loading table $table from dumped data\n";

    load($dbc, $table);

    print "Done\n";
  }
}

print "All done!\n";

sub get_species_list {
  my $config = shift;

  # connect to DB
  my $dbc = DBI->connect(
      sprintf(
        "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=mysql",
        $config->{host},
        $config->{port}
        ), $config->{user}, $config->{password}
      );

  my $version = $config->{version};

  my $sth = $dbc->prepare(qq{ SHOW DATABASES LIKE '%\_variation\_$version%' });
  $sth->execute();

  my $db;
  $sth->bind_columns(\$db);

  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = $config->{pattern};
  @dbs = grep {$_ =~ /$pattern/} @dbs if defined($pattern);

  # remove version, build
  #$_ =~ s/^([a-z]+\_[a-z]+)(.+)/$1/ for @dbs;

  return \@dbs;
}

sub filter_db_list {
  my $config = shift;
  my $db_list = shift;
  my $table = shift;
  my @filtered_db_list = ();

  foreach my $db (@$db_list) {
    my $dbc = DBI->connect(
        sprintf("DBI:mysql(RaiseError=>1):host=%s;port=%s;db=%s",
          $config->{host},
          $config->{port},
          $db
          ), $config->{user}, $config->{password});
    my $sth = $dbc->prepare(qq{SELECT * FROM $table LIMIT 1;});
    $sth->execute();
    while (my @row = $sth->fetchrow_array) {
      my $count = $row[0];
      if ($count > 0) {
        push @filtered_db_list, $db;    
      }
    }
    $sth->finish();
  }
  return \@filtered_db_list;
}

