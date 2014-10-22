# Script to generate an HTML page containing the variation sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($e_version,$html_file,$hlist,$user,$port,$help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'       => \$e_version,
     'o=s'       => \$html_file,
     'help!'     => \$help,
     'hlist=s'   => \$hlist,
     'user=s'    => \$user,
     'port=i'    => \$port,
     'site=s'    => \$site,
     'etype=s'   => \$etype
);

## Missing arguments ##
if (!$e_version) {
  print "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$html_file) {
  print "> Error! Please give an output file using the option '-o'\n";
  usage();
}
if (!$hlist) {
  print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}
if (!$user) {
  print "> Error! Please give user name using the option '-user'\n";
  usage();
}
usage() if ($help);


## Settings ##
my %pops = ('1000 Genomes'                   => { 'order'      => 1,
                                                  'species'    => 'Homo sapiens',
                                                  'term'       => '1000GENOMES:phase_%',
                                                  'constraint' => 'size is not null',
                                                  'url'        => 'http://www.1000genomes.org',
                                                  'evidence'   => '1000Genomes'
                                                },
            'HapMap'                         => { 'order'    => 2,
                                                  'species'  => 'Homo sapiens',
                                                  'term'     => 'CSHL-HAPMAP:%',
                                                  'url'      => 'http://hapmap.ncbi.nlm.nih.gov/index.html.en',
                                                  'evidence' => 'HapMap'
                                                },
            'Exome Sequencing Project (ESP)' => { 'order'    => 3,
                                                  'species'  => 'Homo sapiens',
                                                  'term'     => 'ESP6500:%',
                                                  'url'      => 'http://evs.gs.washington.edu/EVS/',
                                                  'evidence' => 'ESP'
                                                }
           );
           
my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my @hostnames = split /,/, $hlist;
my $database = "";
my $pswd = "";
my $db_type = 'variation';
my $default_port = 3306;
my $margin_bottom_max = '35px';

my $evidence_icon_prefix = '/i/val/evidence_';
my $evidence_icon_suffix = '.png';
my $evidence_doc_url  = '#evidence_status';

$port ||= $default_port;
$server_name = $site if ($site) ;

my $sql  = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};
my $sql2 = qq{SELECT p.name,p.description,p.size,p.freqs_from_gts,d.display_name,d.display_priority FROM population p, display_group d WHERE p.display_group_id=d.display_group_id ORDER by d.display_priority};
my $sql3 = qq{ SELECT p2.name FROM population p1, population p2, population_structure ps 
               WHERE p1.population_id=ps.super_population_id AND p2.population_id=ps.sub_population_id AND p1.population_id=?};

my $bg = '';
my $first_species;
my %species_host;


## Headers ##
my $pop_table_header = qq{
  <tr>
    <th style="width:200px">Name</th>
    <th>Size</th>
    <th>Description</th>
  </tr>
};

my $html_pop_gen_table_header = qq{
<h2 id="populations_gen">Populations genotypes</h2>
<p>In the table displayed below, we listed the populations of the main genotyping projects having genotypes at the population level (e.g. allele and genotype frequencies)</p>
<table class="ss">\n  <tr><th>Species</th><th>Project</th><th>Population name</th><th style="max-width:700px">Population description</th><th>Size</th></tr>};


##########
## Main ##
##########

## Population genotypes ##
my $html_pop_gen = '';
foreach my $hostname (@hostnames) {
  
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    next if ($dbname =~ /sample$/);
    
    print $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    
    if ($etype) { # EG site - need to filter out species
      my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
      #  print "- checking for $img_thumb ... ";
      if (! -e $img_thumb) {
        print "\t... skipping \n";
        next;
      } 
    }
    print "\n";
    
    my $label_name = ucfirst($s_name);
       $label_name =~ s/_/ /g;
    
    $species_host{$label_name} = {'host' => $hostname, 'dbname' => $dbname};
    
    my $projects_list = store_data($dbname, $hostname);
    $html_pop_gen .= table_rows($label_name,$projects_list,$first_species);
    $first_species = $label_name if (!$first_species && scalar(keys(%$projects_list)));
  }
  $sth->finish;
}
$html_pop_gen .= qq{</table>\n<div style="margin-bottom:$margin_bottom_max"></p>};


## Populations ##
my $html_pop = '';

foreach my $project (sort{ $pops{$a}{'order'} <=> $pops{$b}{'order'} } keys(%pops)) {

  my $term = $pops{$project}{'term'};
  my $spe  = $pops{$project}{'species'};
  my $constraint  = ($pops{$project}{'constraint'}) ? $pops{$project}{'constraint'}.' AND ' : '';
  my $project_word = ($project =~ /project/i) ? '' : ' Project';
  my $url = $pops{$project}{'url'};
  my $evidence = $pops{$project}{'evidence'};
  
  $html_pop .= qq{<h4>Populations from the <a href="$url" target="_blank" style="text-decoration:none">$project$project_word</a> ($spe)</h4>\n};

  my $project_id = $project;
  $project_id =~ s/ /_/g;
  $project_id = lc $project_id;
  $html_pop .= qq{<table id="$project_id" class="ss" style="margin-bottom:4px">\n  $pop_table_header\n};

  my $dbname = $species_host{$spe}{'dbname'};
  my $host   = $species_host{$spe}{'host'};
  my $stmt = qq{ SELECT population_id, name, size, description FROM population WHERE $constraint name like ? ORDER BY name};
  my $sth  = get_connection_and_query($dbname, $host, $stmt, [$term]);
  
  $bg = '';
  my %pop_data;
  my %pop_tree;
  my %sub_pops;
  while(my @data = $sth->fetchrow_array) {
  
    my ($pop_name,$pop_suffix) = split(":", $data[1]);
    $pop_name .=  ":<b>$pop_suffix</b>" if ($pop_suffix);
    $data[2] = '-' if (!$data[2]);
    my $desc = parse_desc($data[3]);
    my $size = ($data[2] && $data[2] ne '-' ) ? $data[2] : get_size($data[0], $dbname, $host);
    
    $pop_data{$data[1]} = {'id'    => $data[0],
                           'label' => $pop_name,
                           'desc'  => $desc,
                           'size'  => $size
                          };
    # Super/sub populations                      
    my $sth2 = get_connection_and_query($dbname, $host, $sql3, [$data[0]]);
    while(my ($sub_pop) = $sth2->fetchrow_array) {
      $sub_pops{$sub_pop} = 1;
      $pop_tree{$data[1]}{$sub_pop} = 1; 
    }
    $sth2->finish;
  }
  $sth->finish;
  my $pop_list = get_population_structure(\%pop_data, \%pop_tree, \%sub_pops);

  foreach my $pop (@$pop_list) {
    my $p_name = $pop_data{$pop}{'label'};
       $p_name = qq{<ul style="margin:0px"><li style="margin:0px">$p_name</li></ul>} if ($sub_pops{$pop});
    my $desc   = $pop_data{$pop}{'desc'};
    my $size   = $pop_data{$pop}{'size'};

    $html_pop .= qq{  <tr$bg>\n    <td>$p_name</td>\n    <td style="text-align:right">$size</td>\n    <td>$desc</td>\n  </tr>\n};

    $bg = set_bg($bg);
  }

  $html_pop .= "</table>\n";

  # Evidence status
  if ($pops{$project}{'evidence'}) {
    my $evidence = $pops{$project}{'evidence'};
    my $evidence_img = "$evidence_icon_prefix$evidence$evidence_icon_suffix";
    my $margin = ($pops{$project}{'order'} == scalar(keys(%pops))-1) ? '30px' : $margin_bottom_max;
    $html_pop .= qq{
<p style="margin-bottom:$margin">
Variants which have been discovered in this project have the "evidence status" <a href="$evidence_doc_url"><b>$evidence</b></a>.
On the website this corresponds to the icon <a href="$evidence_doc_url"><img class="_ht" src="$evidence_img" title="$evidence" style="vertical-align:bottom"/></a>.
</p>
    };
  }
}


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_pop;
print HTML $html_pop_gen_table_header;
print HTML $html_pop_gen;
close(HTML);


#############
## Methods ##
#############

sub store_data {
  my $dbname   = shift;
  my $hostname = shift;
  
  my %projects_list;
  
  my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
  
  while( my ($name,$desc,$size,$freq_from_gts,$display_name,$display_order) = $sth2->fetchrow_array ) {
    $desc = parse_desc($desc);
    $projects_list{$display_name}{'order'} = $display_order;
    $projects_list{$display_name}{'pop'}{$name}{'desc'} = $desc;
    $projects_list{$display_name}{'pop'}{$name}{'size'} = $size;
  }
  $sth2->finish;
  return \%projects_list;
}

sub get_population_structure {
  my $pops     = shift;
  my $pop_tree = shift;
  my $sub_pops = shift;

  my $pop_list;
  if ($pop_tree) {

    foreach my $pop_id (sort { ($a !~ /ALL/ cmp $b !~ /ALL/) || $a cmp $b } keys %$pops) {
      next if (grep { $_ eq $pop_id} @$pop_list);

      if (!$pop_tree->{$pop_id} && !$sub_pops->{$pop_id}) {
        push (@$pop_list, $pop_id);
      }
      elsif ($pop_tree->{$pop_id}) {
        push (@$pop_list, $pop_id);
        $pop_list = add_population_to_list($pop_tree,$pop_id,$pop_list);
      }
    }
  }
  else {
    my @list = sort { ($a =~ /ALL/ cmp $b =~ /ALL/) || $a cmp $b } keys %$pops;
    $pop_list = \@list;
  }
  return $pop_list;
}


sub add_population_to_list {
  my $pop_tree = shift;
  my $pop_id   = shift;
  my $pop_list = shift;

  foreach my $sub_pop_id (keys(%{$pop_tree->{$pop_id}})) {
    push (@$pop_list, $sub_pop_id);
    if ($pop_tree->{$sub_pop_id}) {
      $pop_list = add_population_to_list($pop_tree,$sub_pop_id,$pop_list);
    }
  }
  return $pop_list;
}

sub table_rows {
  my $species = shift;
  my $projects_list = shift;
  my $first_species = shift;
  
  my $html = '';
  
  my $count_projects;
  foreach my $project (sort {$projects_list->{$a}{'order'} <=> $projects_list->{$b}{'order'} } keys(%$projects_list)) {
    $count_projects += scalar(keys(%{$projects_list->{$project}{'pop'}}));
  }

  foreach my $project (sort {$projects_list->{$a}{'order'} <=> $projects_list->{$b}{'order'} } keys(%$projects_list)) {
    my $count_pops = scalar(keys(%{$projects_list->{$project}{'pop'}}));
    foreach my $pop (sort { ($a !~ /ALL/ cmp $b !~ /ALL/) || $a cmp $b } keys(%{$projects_list->{$project}{'pop'}})) {
      my $desc = $projects_list->{$project}{'pop'}{$pop}{'desc'};
      my $size = $projects_list->{$project}{'pop'}{$pop}{'size'};
      $html .= get_row($species,$project,$pop,$desc,$size,$first_species,$count_projects,$count_pops);
      $count_projects = undef;
      $count_pops     = undef;
    }
  }
  return $html;
}

sub get_row {
  my $species = shift;
  my $project = shift;
  my $name    = shift;
  my $desc    = shift;
  my $size    = shift;
  my $first_spe = shift;
  my $count_projects = shift;
  my $count_pops     = shift;
  
  my ($pop_name,$pop_suffix) = split(":", $name);
  $pop_name .=  ":<b>$pop_suffix</b>" if ($pop_suffix);
  
  my $border_top = $count_projects ? ($first_spe ? 'border-top:4px solid #AAA' : '') : ($count_pops ? 'border-top:1px solid #CCC' : '');
  
  my $species_col = ($count_projects) ? qq{\n    <td rowspan="$count_projects" style="border-right:1px solid #CCC;background-color:#FFF;font-style:italic;$border_top">$species</td>\n} : '';   
  my $project_col = ($count_pops)     ? qq{\n    <td rowspan="$count_pops" style="border-right:1px solid #CCC;background-color:#FFF;$border_top">$project</td>} : '';
  
  my $html = qq{    
  <tr$bg>$species_col$project_col
    <td style="$border_top">$pop_name</td>
    <td style="max-width:700px;$border_top">$desc</td>
    <td style="text-align:right;$border_top">$size</td>
  </tr>};
  $bg = set_bg($bg);
  
  return $html;
}

# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  my $params = shift;
  
  my ($host, $port) = split /\:/, $hname;
  
  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(@$params);
  }
  else {
    $sth->execute;
  }
  return $sth;
}

sub get_size {
  my $pop_id   = shift;
  my $dbname   = shift;
  my $hostname = shift;

  my $stmt = qq{ SELECT count(*) FROM individual_population WHERE population_id=?};
  my $sth = get_connection_and_query($dbname, $hostname, $stmt, [$pop_id]);
  my $size = ($sth->fetchrow_array)[0];
  $sth->finish;

  return ($size == 0) ? '-' : $size;
}

sub set_bg {
  my $bg_class = shift;
  return ($bg_class eq '') ? ' class="bg2"' : '';
}

sub parse_desc {
  my $content = shift;
  #$content =~ s/<[^>]+>//g; # Remove URLs
  $content = '-' if (!$content);
  my @desc = split(/\.,/, $content);
  $content = "$desc[0]. $desc[1]." if scalar(@desc > 1);
  return $content;
}

sub usage {
  
  print qq{
  Usage: perl sources2html.pl [OPTION]
  
  Put all variation sources, for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -user           MySQL user name (Required)
    -pass           MySQL password. 3306 by default (optional)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}




