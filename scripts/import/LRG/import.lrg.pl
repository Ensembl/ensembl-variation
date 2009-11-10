#!/software/bin/perl

use strict;
#use lib '/homes/wm2/Ensembl/CheckedOut/ensembl-variation/scripts/import/LRG';
#use lib '/homes/wm2/src/lib/perl5/site_perl/5.8.8/';
#use lib '/nfs/users/nfs_i/ianl/ensembl-live/ensembl/modules/';
use Getopt::Long;
use LRG;
use Bio::Seq;
use Bio::EnsEMBL::Registry;

# set default options
my $registry_file = "ensembl.registry";
my ($help);

usage() unless scalar @ARGV;

# get options from command line
GetOptions(
	   'registry_file=s' => \$registry_file,
	   'help' => \$help,
);

# connect to core database
Bio::EnsEMBL::Registry->load_all( $registry_file );
#Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
our $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
die "ERROR: Could not connect to core database\n$usage" unless defined $dbCore;

# check for an input file
my $input_file = shift @ARGV;
die "ERROR: No input file specified\n$usage" unless -e $input_file;

# create an LRG object from it
my $lrg = LRG::LRG::newFromFile($input_file, "\.$$\.temp");

# find the LRG ID
my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();

die "ERROR: Problem with LRG identifier '$lrg_name'\n" unless $lrg_name =~ /^LRG\_[0-9]+$/;

# find mapping
my $mapping_node = $lrg->findNode("updatable_annotation/annotation_set/mapping", {'assembly' => 'GRCh37'});

# if we can't get GRCh37, get any mapping node
$mapping_node = $lrg->findNode("updatable_annotation/annotation_set/mapping") unless defined $mapping_node;

die "ERROR: Could not find mapping node in XML\n" unless defined $mapping_node;

# get a slice adaptor
my $sa = $dbCore->get_SliceAdaptor();
my $slice = $sa->fetch_by_region("chromosome", $mapping_node->data->{'chr_name'});

# first we need to check if the coordinate system for this LRG exists
my $cs_id;
my $cs_name_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from coord_system WHERE name="$lrg_name"});

# if it doesn't, insert into the coord_system table
if(! $cs_name_ref->[0][0]) {
  
  # get the max rank
  my $max_rank = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT MAX(rank) FROM coord_system;})->[0][0];
  $max_rank++;
  
  $dbCore->dbc->do(qq{INSERT INTO coord_system(species_id,name,rank,attrib) values(1,"$lrg_name",$max_rank,"default_version");});
}

# now get the ID
my $cs_id_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT coord_system_id from coord_system WHERE name="$lrg_name";});
$cs_id = $cs_id_ref->[0][0];




# now we need to do the same for seq_region
my $q_seq_region_id;
my $q_seq_length = ($mapping_node->data->{'chr_end'} - $mapping_node->data->{'chr_start'}) + 1;

my $t_seq_region_id = $sa->get_seq_region_id($slice);
my $lrg_name_ref =$dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT name from seq_region WHERE name="$lrg_name"});
if (! $lrg_name_ref->[0][0]) {
  $dbCore->dbc->do(qq{INSERT INTO seq_region(name,coord_system_id,length)values("$lrg_name",$cs_id,$q_seq_length)});
}

# now get the ID
my $q_seqid_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT seq_region_id from seq_region WHERE name="$lrg_name"});
$q_seq_region_id = $q_seqid_ref->[0][0];




# now we need to add entries to meta if needed
my $meta_ref = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT * FROM meta WHERE meta_key = "assembly.mapping" AND meta_value like "$lrg_name\%";});

if(! $meta_ref->[0][0]) {
  # get the assembly name
  my $ass_name = $dbCore->dbc->db_handle->selectall_arrayref(qq{SELECT meta_value FROM meta WHERE meta_key = "assembly.name";})->[0][0];
  
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name");});
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name\#contig");});
  $dbCore->dbc->do(qq{INSERT INTO meta(species_id,meta_key,meta_value) VALUES(1,"assembly.mapping","$lrg_name\#chromosome:$ass_name\#supercontig");});
}



# now we need to add entries to the assembly table
# and to the seq_region_attrib table for mismatches
foreach  my $span(@{$mapping_node->{'nodes'}}) {
  
  next unless $span->name eq 'mapping_span';
  
  # switched coords
  $dbCore->dbc->do(qq{INSERT IGNORE INTO assembly(asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)values($q_seq_region_id,$t_seq_region_id,}.$span->data->{'lrg_start'}.",".$span->data->{'lrg_end'}.",".$span->data->{'start'}.",".$span->data->{'end'}.",".$span->data->{'strand'}."\)");
  
  # find mismatches etc
  foreach my $diff(@{$span->{'nodes'}}) {
	next unless $diff->name eq 'diff';
	
  	$dbCore->dbc->do(qq{INSERT IGNORE INTO seq_region_attrib(seq_region_id,attrib_type_id,value) VALUES($q_seq_region_id, 145, "}.$diff->data->{'lrg_start'}." ".$diff->data->{'lrg_end'}." ".$diff->data->{'lrg_sequence'}."\"\)");
  }
}





sub usage() {
	
	print
	"Usage:
	
	perl import.lrg.pl [options] LRG.xml
	
	-r || --registry_file	registry file
	-h || --help		print this message\n";
	
	die;
}