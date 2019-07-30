#!/bin/bash

echo 'Getting jksrc'
if [ ! -f v335_base.tar.gz ]; then
  wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
  tar xzf v335_base.tar.gz && rm -rf v335_base.tar.gz kent-335_base/java kent-335_base/python
fi

echo 'Getting Bio::DB::HTS'
if [ ! -d Bio-HTS ]; then
  git clone --branch 2.11 --depth 1 https://github.com/Ensembl/Bio-HTS.git
fi

