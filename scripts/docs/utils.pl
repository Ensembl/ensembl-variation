# Determine what type data contains in the vcf file
sub get_vcf_content_types {
  my ($project) = @_;
  my @types;
  
  # this ignores the false positive sigpipe error from tabix command 
  $SIG{PIPE} = 'DEFAULT';
  
  # add if the vcf collection mentions annotation type
  push @types, $project->{annotation_type} if $project->{annotation_type};

  # if use_as_source is set then it is the main source for tracks
  push @types, "source" if $project->{use_as_source};

  # check for populations and sample genotypes
  my $file = get_random_file($project);
  my $file_full_path = $file;
  if ($project->{type} eq "local"){
    $file_full_path = $data_dir . $file_full_path;
  }
  my $chr = `tabix -D $file_full_path -l | head -n 1`;
  chop $chr;
  my $line = `tabix -D $file_full_path $chr | head -n 1`;
  my ($info_field, $format_field) = (split /\t/, $line)[7, 8];

  # populations: can be from file if INFO has AF fields
  push @types, "populations" if ($project->{populations} && ($info_field =~ /AF=/ || ($info_field =~ /AC=/ && $info_field =~ /AN=/)));
  # a hard-coded check for NCBI-ALPHA and TOPMED as they have very special field for frequency
  push @types, "populations" if (($info_field =~ /AN_SAMN/) || ($info_field =~ /TOPMED=/));

  # sample genotype: check FORMAT field
  my $genotypes = `tabix -D $file_full_path -H | grep '##FORMAT' | grep 'ID=GT'`;
  push @types, "genotype" if ($genotypes && $format_field);
  
  return @types;
}

# Check if species have any population in database that would be calculated from genotypes
sub is_freq_from_gts {
  my ($species, $project, @hostnames) = @_;

  foreach my $hostname (@hostnames) {
    my $sql = qq{SHOW DATABASES LIKE '%$species\%variation\_$e_version%'};
    my $sth = get_connection_and_query("", $hostname, $sql);

    while (my ($var_dbname) = $sth->fetchrow_array) {
      my $sql2 = qq{SELECT COUNT(freqs_from_gts) FROM population WHERE freqs_from_gts = 1};
      my $sth2 = get_connection_and_query($var_dbname, $hostname, $sql2);

      return 1 if $sth2->fetchrow_array;
    }
  }

  return 0;
}

# Check if samples from a vcf files exist in either vcf config or database
sub genotype_samples_exists {
  my ($species, $project, @hostnames) = @_;

  # Get samples from vcf file
  my @samples;
  foreach my $file (get_all_files($project)){
    push @samples, (split / /, `bcftools query -l $file | xargs | tr -d '\n'`);
  }

  # Get samples from vcf config 
  my @samples_in_config = keys %{ $project->{sample_populations} };

  # Check if any of the sample matches in config file
  foreach ( @samples_in_config ){
    if (grep /^$_$/, @samples){
      return 1;
    }
  }

  # Check if any of the sample matches in config file
  my $samples_str = join ",", (map { "'$_'" } @samples);
  foreach my $hostname (@hostnames) {
    my $sql = qq{SHOW DATABASES LIKE '%$species\%variation\_$e_version%'};
    my $sth = get_connection_and_query("", $hostname, $sql);

    while (my ($var_dbname) = $sth->fetchrow_array) {
      my $sql2 = qq{SELECT name FROM sample WHERE name IN ($samples_str) LIMIT 1};
      my $sth2 = get_connection_and_query($var_dbname, $hostname, $sql2);

      return 1 if $sth2->fetchrow_array;

      # sample can also have prefix
      my $first_sample = $samples[0];
      my $sql3 = qq{SELECT name FROM sample WHERE name LIKE "%$first_sample" LIMIT 1};
      my $sth3 = get_connection_and_query($var_dbname, $hostname, $sql3);

      my $probable_pop_name = (split(/:/, $sth3->fetchrow_array))[0] || "";
      
      my $sql4 = qq{SELECT name FROM population WHERE name = "$probable_pop_name"};
      my $sth4 = get_connection_and_query($var_dbname, $hostname, $sql4);

      return 1 if ($sth4->fetchrow_array || ($probable_pop_name . ":" eq $project->{sample_prefix}));
    }
  }

  return 0;
}


# Get number of variant from a vcf file
sub get_variant_count {
  my ($project) = @_;
  my $count;

  foreach my $file (get_all_files($project)){
    $count += `bcftools stats $file | grep -ve '^#' | grep -e 'number of records' | cut -d\$'\t' -f 4`;

    unless($count) {
      # If file is remote try downloading it and count the line number
      if ($file =~ /^http/ || $file =~ /^ftp/) {
        `wget -O vcf.gz $file`;
        $count = `bgzip -d -c vcf.gz | grep -v '^#' | wc -l`;
        `rm vcf.gz`;
      }
      # If file is local no need for downloading; just count line number 
      else{
        $count = `bgzip -d -c $file | grep -v '^#' | wc -l`;
      }
    }
  }

  return $count;
}

# Get a random file from filename template in vcf collection
sub get_random_file {
  my ($project) = @_;
  my $file;

  my $filename_template = $project->{filename_template};

  if ($filename_template =~ /###CHR###/){
    my $chromosomes = $project->{chromosomes};

    return undef unless $chromosomes;

    my $chr = @{ $chromosomes }[0];
    
    $file = $filename_template =~ s/###CHR###/$chr/gr;
  }
  else{
    $file = $filename_template
  }

  return $file;
}

# Get all files from filename template in vcf collection
sub get_all_files {
  my ($project) = @_;
  my @files;

  my $filename_template = $project->{filename_template};

  if ($filename_template =~ /###CHR###/){
    my $chromosomes = $project->{chromosomes};

    return undef unless $chromosomes;

    foreach my $chr (@{ $chromosomes }){
      my $file = $filename_template =~ s/###CHR###/$chr/gr;
      push @files, $file;
    }
  }
  else{
    push @files, $filename_template;
  }

  return @files;
}

1;