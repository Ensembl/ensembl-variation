#!/usr/bin/env nextflow

process drop_translation_mapping {
  output: stdout

  cache false

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} <<'EOF'

    DROP TABLE IF EXISTS translation_mapping;

  EOF
  """

}

process store_translation_mapping {
  input:
    path translation_mapping
    val wait
  output: stdout

  cache false

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} --local-infile=1 <<'EOF'

    CREATE TABLE IF NOT EXISTS translation_mapping (
      stable_id   VARCHAR(255),
      md5         CHAR(32),
      PRIMARY KEY (stable_id),
      KEY md5_idx (md5)
    );

    LOAD DATA LOCAL INFILE '${translation_mapping}'
    IGNORE INTO TABLE translation_mapping;

  EOF
  """
}

process clear_assemblies {
  output: stdout
  cache false

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} <<'EOF'

    DELETE FROM meta WHERE meta_key = 'assembly';

  EOF
  """
}

process store_assemblies {
  input:
    path files
    val wait
  output: stdout
  cache false

  """
  ls * | grep -Eo "GCA_[0-9]+(\\.[0-9]+)?|GRCh3[7-8]" | awk '{print "assembly\t"\$1}' > assemblies.txt

  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} --local-infile=1 <<'EOF'

    LOAD DATA LOCAL INFILE 'assemblies.txt'
    IGNORE INTO TABLE meta (meta_key, meta_value);

  EOF
  """
}

process delete_prediction_data {
  input: val value
  output: stdout

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} <<'EOF'

    DELETE pfp.*
    FROM protein_function_predictions pfp, attrib a
    WHERE pfp.analysis_attrib_id = a.attrib_id AND a.value LIKE '${value}';

    DELETE pfpa.*
    FROM protein_function_predictions_attrib pfpa, attrib a
    WHERE pfpa.analysis_attrib_id = a.attrib_id AND a.value LIKE '${value}';

  EOF
  """
}

process update_meta {
  input:
    val key
    val value
  output: stdout

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} <<'EOF'

    REPLACE into meta ( meta_key, meta_value ) VALUES ( '${key}', '${value}' );

  EOF                                                                       
  """                                                                           
}

process get_current_MD5_translations {
  input: val analysis
  output: stdout

  cache false

  """
  # --batch prints output in nontabular format
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} --batch --skip-column-names <<'EOF'   

    SELECT DISTINCT t.translation_md5
    FROM   translation_md5 t, attrib a, protein_function_predictions p
    WHERE  t.translation_md5_id = p.translation_md5_id
    AND    a.attrib_id = p.analysis_attrib_id
    AND    a.value LIKE '${analysis}';

  EOF
  """
}

process init_sqlite_db {
  output: stdout

  cache false

  """
  #!/usr/bin/perl

  use DBI;
  
  my \$dbh = DBI->connect("dbi:SQLite:dbname=${params.sqlite_db}","","");
  \$dbh->do("DROP TABLE IF EXISTS predictions");
  \$dbh->do("CREATE TABLE predictions(md5, analysis, matrix)");
  """
}

process postprocess_sqlite_db {
  input:
    val sift_run
    val polyphen_run

  output: stdout

  cache false

  """
  #!/usr/bin/perl
  
  use DBI;

  my \$dbh = DBI->connect("dbi:SQLite:dbname=${params.sqlite_db}","","");
  \$dbh->do("CREATE INDEX md5_idx ON predictions(md5)");
  """
}