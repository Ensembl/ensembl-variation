=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DumpVEPRefSeq;

use strict;
use warnings;
use File::Path qw(make_path);
use File::Copy;
use Storable qw(nstore_fd);

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

use Bio::EnsEMBL::Variation::Utils::VEP qw(@REG_FEAT_TYPES);

sub param_defaults {
  return {
    'species_refseq' => 0,
    'species_flags'  => {},
    'overwrite'      => 0,
    'sift'           => 0,
    'polyphen'       => 0,
    'regulatory'     => 0,
    'eg'             => 0,
  };
}


sub run {
  my $self = shift;

  my $params = $self->get_vep_params();
  
  # construct command
  my $cmd = sprintf(
    '%s %s/variant_effect_predictor.pl %s --build_part t --host %s --port %i --user %s %s --species %s --assembly %s --db_version %s --dir %s %s --cache_version %s --is_multispecies %s',
    $params->{perl},
    $params->{vep_dir},
    $params->{vep},
    
    $params->{host},
    $params->{port},
    $params->{user},
    $params->{pass},
    
    $params->{species},
    $params->{assembly},
    $params->{version},
    $params->{dir},
    $params->{species_flags_cmd},

    $params->{eg_version} || $params->{version},
    $params->{is_multispecies}
  );
 
  my $finished = 0;
 
  if($params->{debug}) {
    print STDERR "$cmd\n";
  }
  else {
    open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
    my @buffer;
    while(<CMD>) {
      $finished = 1 if /Finished/;
      push @buffer, $_;
      shift @buffer if scalar @buffer > 20;
    }
    close CMD;
  
    die "ERROR: Encountered an error running VEP\n".join("", @buffer)."\n" unless $finished;

    # now we need to copy in the var and reg data from the core run
    my $core_dir = sprintf('%s/%s/%s_%s/', map {$params->{$_}} qw(dir species version assembly));
    my $refseq_dir = sprintf('%s/%s_refseq/%s_%s/', map {$params->{$_}} qw(dir species version assembly));

    opendir CORE, $core_dir or die("ERROR: Could not read from directory $core_dir\n");
    my @core_chr_dirs = grep {!/^\./ && -d $core_dir.'/'.$_} readdir CORE;
    closedir CORE;

    my $sereal = $cmd =~ /--sereal/ ? 1 : 0;
    if($sereal) {
      eval q{use Sereal::Encoder;};
      die("ERROR: Could not use Sereal module\n") if $@;
    }

    foreach my $chr(@core_chr_dirs) {
      opendir CORE, $core_dir.'/'.$chr or die("ERROR: Could not read from directory $core_dir\/$chr\n");

      my $tr_cache = {$chr => []};

      # dir exists in core but not in refseq
      make_path($refseq_dir.'/'.$chr) unless -d $refseq_dir.'/'.$chr;

      foreach my $file(grep {/var|reg/} readdir CORE) {
        copy($core_dir.'/'.$chr.'/'.$file, $refseq_dir.'/'.$chr.'/'.$file)
          or die("ERROR: Failed to copy ".$core_dir.'/'.$chr.'/'.$file.' to '.$refseq_dir.'/'.$chr.'/'.$file."\n");

        # we need to create a dummy transcript file now too
        unless(-d $refseq_dir.'/'.$chr) {
          $file =~ s/\_.+$//;
          $self->write_serialised_dump_file("$refseq_dir/$chr/$file\.gz", $tr_cache, $sereal);
        }
      }

      closedir CORE;
    }

    # there may be some dirs in refseq not in core
    # the healthcheck script will complain we have missing var/reg files unless we create some dummy empty files
    if($self->param('variation') or $self->param('regulatory')) {
      my $var = $self->param('variation');
      my $reg = $self->param('regulatory');

      opendir REF, $refseq_dir or die("ERROR: Could not read from directory $refseq_dir\n");
      foreach my $chr(grep {!/^\./ && -d $refseq_dir.'/'.$_} readdir REF) {
        next if grep {$chr eq $_} @core_chr_dirs;

        my $rf_cache;
        $rf_cache->{$chr}->{$_} = [] for @REG_FEAT_TYPES;

        # we need to create a var/reg file for each tr file
        opendir CHR, $refseq_dir.'/'.$chr;
        foreach my $file(grep {/gz/} readdir CHR) {
          my $stem = $file;
          $stem =~ s/\.gz$//;

          # var is empty file
          if($var && !-e "$refseq_dir/$chr/$stem\_var.gz") {
            open VAR, "| gzip -c -9 > $refseq_dir/$chr/$stem\_var.gz";
            close VAR;
          }
          # reg has basic hash structure
          if($reg) {
            $self->write_serialised_dump_file("$refseq_dir/$chr/$stem\_reg.gz", $rf_cache, $sereal);
          }
        }
        close CHR;
      }
      closedir REF;
    }

    # now we want to update the refseq cache's info.txt
    # need to copy any new or updated keys from the core cache's info.txt
    open INFO, "$core_dir/info.txt" or die("ERROR: failed to read from $core_dir/info.txt\n");
    my %core_info;
    while(<INFO>) {
      next if /^\#/;
      my @split = split /\t/, $_;
      my $key = shift @split;
      $core_info{$key} = join("\t", @split);
    }
    close INFO;

    open INFO, "$refseq_dir/info.txt" or die("ERROR: failed to read from $refseq_dir/info.txt\n");
    open NEW, ">$refseq_dir/info.txt_new";
    while(<INFO>) {
      if(/^\#/) {
        print NEW;
      }
      else {
        my @split = split /\t/, $_;
        my $key = shift @split;
        my $value = join("\t", @split);

        if(defined($core_info{$key})) {
          print NEW $core_info{$key} eq $value ? $_ : $core_info{$key};
          delete $core_info{$key};
        }
        else {
          print NEW;
        }
      }
    }
    close INFO;

    # add remaining keys that weren't found in refseq
    printf NEW $_."\t".$core_info{$_} for keys %core_info;
    close NEW;

    # move new one in place of old
    unlink("$refseq_dir/info.txt");
    move("$refseq_dir/info.txt_new", "$refseq_dir/info.txt");

    # run healthcheck
    $self->healthcheck_cache($params);
  
    $self->tar($self->param('species_refseq') ? 'refseq' : '');
  }
  
  return;
}

sub write_serialised_dump_file {
  my ($self, $file, $obj, $sereal) = @_;

  # we don't want to nuke any existing files!
  return if -e $file;

  my $encoder;
  if($sereal) {
    $encoder = $self->{_encoder} ||= Sereal::Encoder->new({compress => 1});
  }

  if($sereal) {
    open OUT, ">".$file or die("ERROR: Could not write to dump file $file\n");
    print OUT $encoder->encode($obj);
    close OUT;
  }
  else {
    open my $fh, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file\n";
    nstore_fd($obj, $fh);
    close $fh;
  }
}


1;
