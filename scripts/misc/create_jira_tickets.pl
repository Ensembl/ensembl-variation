# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

# Script to generate an HTML page containing the variant sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


# Documentation:
# - https://docs.atlassian.com/jira/REST/latest/#d2e86
# - https://developer.atlassian.com/jiradev/api-reference/jira-rest-apis/jira-rest-api-tutorials/jira-rest-api-example-create-issue

use strict;
use Getopt::Long;

my ($user, $password, $input_file, $version, $help);
GetOptions(
  'user|u=s'    => \$user,
  'pass|p=s'    => \$password,
  'infile|i=s'  => \$input_file,
  'version|v=i' => \$version,
  'help!'       => \$help
);

usage("An input file (--infile) is required") unless (defined($input_file));
usage("An user name (--user) is required") unless (defined($user));
usage("An password (--pass) is required") unless (defined($password));

usage("The input file '$input_file' hasn't been found") unless(-e $input_file);

usage() if ($help);

my $jira_url_root = 'https://www.ebi.ac.uk/panda/jira/rest/api/2/issue';


#-----------#
# Main task #
#-----------#
my $data = `cat $input_file`;
   $data = cleanup_json($data);

if ($version) {
  $data =~ s/##release##/$version/g;
}
else {
  $data =~ s/##release##//g;
}
$data =~ s/##user##/$user/g;


my @issues = $data =~ /"fields":/g;

my $extra = (scalar @issues > 1) ? '/bulk' : '';
my $url = "$jira_url_root$extra";

print "## DATA - ticket ##\n$data\n\n";
my $response = run_rest_call($user, $password, $data, $url);
print "## DATA - response ##\n$response\n\n";


#-----------------------#
# Sub task(s), if exist #
#-----------------------#
if ($input_file =~ /_main_task/) {
  my $input_sub_file = $input_file;
     $input_sub_file =~ s/main/sub/i;
  die ("The input file '$input_sub_file' for the sub-tasks hasn't been found") unless(-e $input_sub_file);
  
  if ($response =~ /\"key\"\:\"(.+)"\,\n?/i) {
    my $parent_key = $1;
    
    my $sub_data = `cat $input_sub_file`;
       $sub_data = cleanup_json($sub_data);
    $sub_data =~ s/##parent-key##/$parent_key/g;
    $sub_data =~ s/##user##/$user/g;
    if ($version) {
      $sub_data =~ s/##release##/$version/g;
    }
    else {
      $sub_data =~ s/##release##//g;
    }
  
    my @sub_issues = $sub_data =~ /"fields":/g;

    my $sub_extra = (scalar @sub_issues > 1) ? '/bulk' : '';
    my $sub_url = "$jira_url_root$sub_extra";
    
    print "## SUB DATA => Sub-task(s) ##\n$sub_data\n\n";
    my $sub_response = run_rest_call($user, $password, $sub_data, $sub_url);
    print "## SUB DATA - sub-response ##\n$sub_response\n\n";
  }
}


#---------#
# Methods #
#---------#

sub cleanup_json {
  my $json = shift;
  $json =~ s/\n//g;
  $json =~ s/\s*{\s*/{/g;
  $json =~ s/\s*}\s*/}/g;
  $json =~ s/\s*\[\s*/\[/g;
  $json =~ s/\s*\]\s*/\]/g;
  $json =~ s/,\s*/, /g;
  return $json;
}

sub run_rest_call {
  my $login    = shift;
  my $pwd      = shift;
  my $json     = shift;
  my $jira_url = shift;
  
  my $cmd = sprintf( qq{curl -D- --user \%s:\%s -X POST --data '\%s' -H "Content-Type: application/json" \%s},
                     $login, $pwd, $json, $jira_url
                   );
  return `$cmd`;
}

sub usage {
  my $msg = shift;
  
  print "---------\n$msg\n---------\n" if ($msg);
  print qq{
  Usage: perl create_jira_tickets.pl [OPTION]
  
  Create JIRA tickets and sub-tasks, from a JSON file, using the JIRA REST API.
  https://developer.atlassian.com/jiradev/jira-apis/jira-rest-apis/jira-rest-api-tutorials
  
  Options:

    -help           Print this message
      
    -infile  | -i  The JSON input file name (Required)      
    -user    | -u  The JIRA user login (Required)
    -pass    | -p  The JIRA user password (Required)
    -version | -v  The Ensembl release, e.g. 88 (optional)
  exit(0);
}
