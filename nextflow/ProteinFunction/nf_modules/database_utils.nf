#!/usr/bin/env nextflow

process store_translation_mapping {
  input: path translation_mapping
  output: stdout

  cache false

  """
  mysql --host=${params.host} --port=${params.port} \
        --user=${params.user} --password=${params.pass} \
        ${params.database} --local-infile=1 <<'EOF'

    DROP TABLE IF EXISTS translation_mapping;

    CREATE TABLE translation_mapping (
      stable_id   VARCHAR(255),
      md5         CHAR(32),
      PRIMARY KEY (stable_id),
      KEY md5_idx (md5)
    );

    LOAD DATA LOCAL INFILE '${translation_mapping}'
    INTO TABLE translation_mapping;

  EOF
  """
}

process delete_prediction_data {
  input: val value
  output: stdout

  cache false

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

  cache false

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
