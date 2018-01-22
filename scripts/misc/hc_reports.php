<?php
  /* 
    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2018] EMBL-European Bioinformatics Institute
 
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
 
      http://www.apache.org/licenses/LICENSE-2.0
 
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
  */
  
  /*
    This PHP script will display all the HTML reports of the HealthChecks.
    For that you need 2 things:
      1- The HTML files should be name like this <specis_or_group>_<ensembl_release>.html or <specis_or_group>_<ensembl_release_<assembly>.html
         e.g.: cow_90.html, primates_90.html or human_90_38.html
      2- This PHP file and the HTML files should be stored in the same directory
    The script converting the TXT HealthChecks reports into HTML is available here:
    https://github.com/Ensembl/ensj-healthcheck/blob/master/healthcheck_txt2html.pl
  */
  
  $hc_id = '';
  $title = 'HC reports';
  
  if ($_GET['hc']) {
    $hc_id = $_GET['hc'];
    $hc_label = ucfirst($hc_id);
    $title = "HealthChecks of ".str_replace('_', ' ', $hc_label);
  }
  
  echo <<<EOF
  <html>
    <head>
      <title>$title</title>
      <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
      <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css">
    </head>
    <body style="padding: 5px 8px">
EOF;
  
  $select_list  = array();
  $version_list = array(); 
  if ($handle = opendir('./')) {
    while (false !== ($file = readdir($handle))) {
      if (preg_match('/^(\w+)\.html$/',$file,$matches)) {
        $hc = $matches[1];
        $label_array  = split('_',$hc);
        $species = ucfirst($label_array[0]);
        $version_list[$label_array[1]] = 1;
        
        $label_array[0] = $label_array[1].' -';
        $label_array[1] = $species;
        $report_label = implode(' ',$label_array);
        
        $select_list[$hc] = $report_label;
      }
    }
    natsort($select_list);
    closedir($handle);
  }
  
  echo <<<EOF
  <div style="padding:2px 6px;background-color:#000">
    <div class="glyphicon glyphicon glyphicon-list" style="float:left;color:#FFF;padding:6px 8px 6px 2px;font-size:18px"></div>
    <div style="float:left;padding:6px 0px;color:#FFF">List of available HC reports:</div>
    <div style="float:left;padding:6px 0px;margin-left:10px">
      <form style="margin-bottom:0px">
EOF;
  $space = '        ';
  echo "$space<select name=\"hc\" onchange='this.form.submit()'>";
  echo "<option value=\"\">-</option>";
  foreach ($select_list as $report => $report_label) {
    $selected = '';
    if ($report == $hc_id) {
      $selected = ' selected';
    }
    echo "$space  <option value=\"$report\"$selected>$report_label</option>";
  }
  echo "$space</select>";
  echo <<<EOF
      </form>
    </div>
    <div style="clear:both"></div>
  </div>
EOF;
  
  
  if ($hc_id != '') {
    include("./$hc_id.html");
  }
  else {
    krsort($version_list);
    $array_length = count($select_list);
    echo "<h3>List of available HC reports ($array_length)</h3>";
    foreach ($version_list as $version => $flag) {
      echo "<h4><span class=\"glyphicon glyphicon-chevron-right\"></span> Release <b>$version</b></h4>";
      echo "<ul>";
      foreach ($select_list as $report => $report_label) {
        if (preg_match("/^$version\s-\s(\w+\s*\w+)$/",$report_label,$matches2)) {
          echo "<li><a href=\"?hc=$report\">$matches2[1]</a></li>";
        }
      }
      echo "</ul>";
    }
  }

  echo <<<EOF
    </body>
  </html>
EOF;
?>



