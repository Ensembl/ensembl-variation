package Bio::EnsEMBL::Variation::Pipeline::BackupDB::InitBackup;

use strict;
use warnings;
use File::Path qw(make_path);

use Bio::EnsEMBL::Hive::Utils ('go_figure_dbc');

use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;
}

sub write_output {
  my $self = shift;
    
  my $species            = $self->param('species'),
  my $backup_dir         = $self->param('backup_dir');
  my $src_db_conn        = $self->param('src_db_conn');
  my $output_file_prefix = $self->param('output_file_prefix');
  my $output_file_suffix = $self->param('output_file_suffix');
  my $extension          = $self->param('output_extension');
  my $tables             = $self->param('tables');
  my $exclude_tables     = $self->param('exclude_tables');
  my $db_name            = $self->param('db_name');
  my $info_updates       = $self->param('info_updates');
  

  die("A species (-species) needs to be specified!") unless (defined($species));
  die("A backup directory (-backup_dir) needs to be specified!") unless (defined($backup_dir));
  die("Database parameters (-src_db_conn) needs to be specified!") unless (defined($src_db_conn));
  die("An output file prefix (-output_file_prefix) needs to be specified!") unless (defined($output_file_prefix));
  die("An output file suffix (-output_file_suffix) needs to be specified!") unless (defined($output_file_suffix));

  if (!-d $backup_dir) {
    make_path "$backup_dir" or die "Failed to create directory: '$backup_dir";
  }

  my $tmp_info_updates = "$backup_dir/tmp_$info_updates";
  

  # Get directory SQL files
  my $dh;
  opendir($dh,$backup_dir);
  die("Could not process directory $backup_dir") unless (defined($dh));
  my @dir_files = readdir($dh);
  my @sql_files = grep {$_ =~ m/\.$extension$/} @dir_files;
  # Close the dir handle
  closedir($dh);


  my $dbc = go_figure_dbc($src_db_conn);
  my %new_update_list = map { $_->[0] => $_->[1]} @{$self->get_updated_tables_list($dbc,$db_name)};

  # Get list of updated tables
  if ($tables eq 'update' && -e "$backup_dir/$info_updates") {

    my $dumped_tables = $self->get_tables_list($dbc,$exclude_tables);

    my %old_update_list = %{$self->parse_info_updates("$backup_dir/$info_updates")};

    # Get the list of tables updated since the last backup
    my @updated_tables_list;
    foreach my $table (@$dumped_tables) {
      if ($old_update_list{$table} && $new_update_list{$table}) {
        if ($old_update_list{$table} ne $new_update_list{$table}) {
          push(@updated_tables_list, $table);
        }
      }
      elsif ($new_update_list{$table}) {
        push(@updated_tables_list, $table);
      }
    }
    
    if (scalar(@updated_tables_list) == 0) {
      die("No tables have been updated in the database since the last backup - End of the pipeline");
    }
    else {
      $tables = join(' ',@updated_tables_list);
    }
  }
  elsif ($tables eq 'all' || $tables eq 'update') {
    $tables = $self->get_tables_list($dbc,$exclude_tables);
  }
  else {
    my @tables_list = ref($tables) eq 'ARRAY' ? @$tables : split(' ', $tables);
    $tables = \@tables_list;
  }

  # Create a new temporary info updates file
  $self->create_tmp_info_update_file($tmp_info_updates,\%new_update_list);
  

  my @jobs = ();
  my $id = 0;
  foreach my $table (@$tables) {
  
    my $file_version = $self->get_new_dump_file_version(\@sql_files,$species,$table,$output_file_suffix,$extension);
   
    my $file_name = "$backup_dir/$output_file_prefix$table$output_file_suffix$file_version.$extension";
    die("The file '$file_name' already exist!") if (-e $file_name);
    
    $id++;
    push @jobs, {
        'id'            => $id,
        'output_file'   => $file_name,
        'src_db_conn'   => $src_db_conn,
        'exclude_ehive' => 1,
        'table_list'    => $table,
    };
  }
    
  $self->dataflow_output_id(\@jobs, 2);
  
  $self->dataflow_output_id({ 'tables' => join(' ', @$tables) }, 1);  
  
  return;
}

sub parse_info_updates {
  my $self      = shift;
  my $info_file = shift;

  my %data;

  open INFO, "< $info_file" or die $!;
  while (<INFO>) {
    chomp $_;
    my @line = split("\t",$_);
    $data{$line[0]} = $line[1];
  }
  close(INFO);

  return \%data;
}


sub get_tables_list {
  my $self            = shift;
  my $dbc             = shift;
  my $exclude_tables = shift; 
  
  my @tables_list = map{ $_->[0] } @{$dbc->db_handle->selectall_arrayref(qq{SHOW TABLES})};

  my @dumped_tables;
  my @t_excluded;
  if ($exclude_tables) {
    $exclude_tables =~ s/%//g;
    @t_excluded = ref($exclude_tables) eq 'ARRAY' ? @$exclude_tables : split(' ', $exclude_tables);
    
    # Get a shorter list
    foreach my $table (@tables_list) {
      if (!grep{ $table =~ /$_/} @t_excluded) {
        push(@dumped_tables, $table);
      }
    }
  }
  else {
    @dumped_tables = @tables_list;
  }
  return \@dumped_tables
}


sub get_updated_tables_list {
  my $self    = shift;
  my $dbc     = shift;
  my $db_name = shift;
  
  my $stmt = qq{ SELECT TABLE_NAME,UPDATE_TIME FROM information_schema.tables WHERE TABLE_SCHEMA='$db_name'};
  my $list = $dbc->db_handle->selectall_arrayref($stmt);
  
  return $list;
}

sub get_new_dump_file_version {
  my $self      = shift;
  my $files     = shift;
  my $species   = shift;
  my $table     = shift;
  my $suffix    = shift;
  my $extension = shift;
  
  my $version = 0;
  
  my @dump_files = grep { $_ =~ m/$species\_$table$suffix/ } @$files;
  
  foreach my $file (@dump_files) {
    if ($file =~ m/$species\_$table$suffix(\d+)\.$extension$/) {
      $version = $1 if ($1 > $version);
    }
  }
  $version++;
  return $version;
}

sub create_tmp_info_update_file {
  my $self          = shift;
  my $tmp_info_file = shift;
  my $info_data     = shift;
  
  open FILE, "> $tmp_info_file" or die $!;

  while (my ($table, $time) = each %{$info_data}) {
    print FILE "$table\t$time\n";
  }
  close FILE;
}

1;
