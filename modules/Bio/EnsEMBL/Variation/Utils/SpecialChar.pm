=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::SpecialChar

=head1 DESCRIPTION

Replace special characters in strings, example phenotype descriptions, submitter names, author names.

=cut

package Bio::EnsEMBL::Variation::Utils::SpecialChar;

use strict;
use warnings;
use base qw(Exporter);


our @EXPORT_OK = qw(replace_char replace_hex);

my %special_characters = (
  'Å' => 'A',
  'Ä' => 'A',
  'Ã' => 'A',
  'Ö' => 'O',
  'Ü' => 'U',
  'ö' => 'o',
  'ô' => 'o',
  'ó' => 'o',
  'ü' => 'u',
  'ò' => 'o',
  'ü' => 'u',
  'ú' => 'u',
  'ù' => 'u',
  'ä' => 'a',
  'ã' => 'a',
  'á' => 'a',
  'í' => 'i',
  'ï' => 'i',
  'ı' => 'i',
  'é' => 'e',
  'è' => 'e',
  'ë' => 'e',
  'ê' => 'e',
  'ç' => 'c',
  'ş' => 's',
  'ğ' => 'g',
  '<' => 'less than',
  '>' => 'more than',
  '&' => 'and',
  '%' => 'percent',
  'Đ' => '-',
  '^' => '',
);


=head2 replace_char


  Example     : my $new_string = replace_char($old_string);
  Description : returns the old_string with special characters replaced in a standard format
  ReturnType  : String
  Exceptions  : None
  Caller      : General


=cut

sub replace_char {
  my $input_string = shift;

  foreach my $char (keys(%special_characters)) {
    my $new_char = $special_characters{$char};
    $input_string =~ s/$char/$new_char/g;
  }

  return $input_string;
}

=head2 replace_hex


  Example     : my $new_string = replace_hex($old_string);
  Description : returns the old_string with hexadecimal characters replaced in a standard format
  ReturnType  : String
  Exceptions  : None
  Caller      : General

=cut

sub replace_hex {
  my $text = shift;

  $text  =~ s/&#xe9;/e/g;
  $text  =~ s/&#xe8;/e/g;
  $text  =~ s/&#xc9;/E/g;
  $text  =~ s/&#xe4;/a/g;
  $text  =~ s/&#xe0;/a/g;
  $text  =~ s/&#xe1;/a/g;
  $text  =~ s/&#xe2;/a/g;
  $text  =~ s/&#xe3;/a/g;
  $text  =~ s/&#x03b1;/A/g;
  $text  =~ s/&#x03b2;/B/g;
  $text  =~ s/&#xdf;/B/g;
  $text  =~ s/&#x2212;/-/g;
  $text  =~ s/&#x2010;/-/g;
  $text  =~ s/&#x2013;/-/g;
  $text  =~ s/&#x2014;/-/g;
  $text  =~ s/&#x2018;/'/g;
  $text  =~ s/&#x2019;/'/g;
  $text  =~ s/&#x201c;/"/g;
  $text  =~ s/&#x201d;/"/g;
  $text  =~ s/&#xf3;/o/g;
  $text  =~ s/&#xf4;/o/g;
  $text  =~ s/&#xf6;/o/g;
  $text  =~ s/&#xf4;/o/g;
  $text  =~ s/&#xed;/i/g;
  $text  =~ s/&#xee;/i/g;
  $text  =~ s/&#xef;/i/g;
  $text  =~ s/&#xfa;/u/g;
  $text  =~ s/&#xfb;/u/g;
  $text  =~ s/&#xfc;/u/g;

  return $text;
}

1;
