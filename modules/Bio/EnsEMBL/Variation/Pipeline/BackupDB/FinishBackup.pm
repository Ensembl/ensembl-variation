package Bio::EnsEMBL::Variation::Pipeline::BackupDB::FinishBackup;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Utils ('go_figure_dbc');
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;
  my $species            = $self->param('species'),
  my $backup_dir         = $self->param('backup_dir');
  my $src_db_conn        = $self->param('src_db_conn');
  my $db_name            = $self->param('db_name');
  my $output_file_prefix = $self->param('output_file_prefix');
  my $output_file_suffix = $self->param('output_file_suffix');
  my $extension          = $self->param('output_extension');
  my $tables             = $self->param('tables');
  my $exclude_tables     = $self->param('exclude_tables');
  my $info_updates       = $self->param('info_updates');
  

  die("A species (-species) needs to be specified!") unless (defined($species));
  die("A backup directory (-backup_dir) needs to be specified!") unless (defined($backup_dir));
  die("Database parameters (-src_db_conn) needs to be specified!") unless (defined($src_db_conn));
  die("An output file prefix (-output_file_prefix) needs to be specified!") unless (defined($output_file_prefix));
  die("An output file suffix (-output_file_suffix) needs to be specified!") unless (defined($output_file_suffix));
  die("Database name (-db_name) needs to be specified!") unless (defined($db_name));

  die("Backup directory '$backup_dir' doesn't exist!") unless (-d $backup_dir);

  my $tmp_info_updates = "$backup_dir/tmp_$info_updates";
  die("Temporary file '$tmp_info_updates' doesn't exist!") unless (-e $tmp_info_updates);


  my $grep_start_table = 'CREATE TABLE';
  my $grep_end_table_1 = 'ALTER TABLE';
  my $grep_end_table_2 = 'ENABLE KEYS';

  # Checks that the backup succeeded.

  my @t_excluded;
  
  if ($exclude_tables) {
    $exclude_tables =~ s/%//g;
    @t_excluded = ref($exclude_tables) eq 'ARRAY' ? @$exclude_tables : split(' ', $exclude_tables);
  }

  # Retrieve list of tables
  my @tables_list;
  if ($tables eq 'all') {
    my $dbc = go_figure_dbc($src_db_conn);
    @tables_list = map{ $_->[0] } @{$dbc->db_handle->selectall_arrayref(qq{SHOW TABLES})};
  }
  else {
    @tables_list = ref($tables) eq 'ARRAY' ? @$tables : split(' ', $tables);
  }

  # Get directory SQL files
  my $dh;
  opendir($dh,$backup_dir);
  die("Could not process directory $backup_dir") unless (defined($dh));
  my @dir_files = readdir($dh);
  my @sql_files = grep {$_ =~ m/\.$extension$/} @dir_files;
  # Close the dir handle
  closedir($dh);
  

  my $report_file = $output_file_prefix.'report.txt';

  open REPORT, "> $backup_dir/$report_file" || die "Failed to open $backup_dir/$report_file : $!\n";

  # Check the number of tables dumped into the SQL file
  my @dumped_tables;
  my $total_count = @tables_list;
  my $count_create = 0;
  my $count_loaded = 0;
  my $error = 0;
  foreach my $table (@tables_list) {
    if (!grep{ $table =~ /$_/} @t_excluded || !@t_excluded) {
    
      my $file_version = $self->get_dump_file_version(\@sql_files,$species,$table,$output_file_suffix,$extension);
    
      my $dump_file = "$backup_dir/$output_file_prefix$table$output_file_suffix$file_version.$extension";
      print REPORT "Backup file '$dump_file' for the table '$table' doesn't exist!\n" unless (-e $dump_file);
      print REPORT "Backup file '$dump_file' for the table '$table' is empty\n!" if (-e $dump_file && -z $dump_file);
      
      my $res = `grep '$grep_start_table \`$table\`' $dump_file`;
      if ($res =~ /$grep_start_table `$table`/) {
        $count_create++;
       
        # Check that the table has been fully dumped (i.e. the keys have been enabled)
        my $res2 = `grep '$grep_end_table_1 \`$table\` $grep_end_table_2' $dump_file`;
        if (!$res2 || $res2 !~ /$grep_end_table_1 \`$table\` $grep_end_table_2/) {
          print REPORT "The table '$table' doesn't seem to have been dumped correctly\n";
          $error = 1;
        }
      }
      else {
        print REPORT "The table '$table' is missing from the dump $dump_file\n";
        $error = 1;
      }
      push(@dumped_tables, $table);
    }
    else {
      $total_count--;
      next;
    }
  }

  if ($total_count != $count_create) {
    print REPORT "There is a difference between the expected number of tables dumped ($total_count) and the actual number of dumped tables ($count_create)\n";
    $error = 1;
  }

  # Check if the excluded tables have been dumped - Normally, they shouldn't!
  foreach my $e_table (@t_excluded) {
    my $res = `grep '$grep_start_table \`$e_table' $backup_dir/$output_file_prefix$e_table*`;
    if ($res && $res =~ /$grep_start_table `$e_table/) {
      print REPORT "The excluded table '$e_table' has been dumped. It shouldn't have been dumped!\n";
      $error = 1;
    }
  }

  if ($error == 0) {
    print REPORT "The following tables of the MySQL database '$src_db_conn' have been successfully dumped (".scalar(@dumped_tables)." tables):\n- ".join("\n- ",@dumped_tables)."\n";
  }
  close(REPORT);

  if (-e $tmp_info_updates) {
    `mv $tmp_info_updates $backup_dir/$info_updates`;
  }

}

sub get_dump_file_version {
  my $self      = shift;
  my $files     = shift;
  my $species   = shift;
  my $table     = shift;
  my $suffix    = shift;
  my $extension = shift;
  
  my $version = 1;
  
  my @dump_files = grep { $_ =~ m/$species\_$table$suffix/ } @$files;
  
  foreach my $file (@dump_files) {
    if ($file =~ m/$species\_$table$suffix(\d+)\.$extension$/) {
      $version = $1 if ($1 > $version);
    }
  }
  return $version;
}

1;
