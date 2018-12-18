#!/bin/bash

echo 'Getting jksrc'
if [ ! -f v335_base.tar.gz ]; then
  wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
  tar xzf v335_base.tar.gz && rm -rf v335_base.tar.gz kent-335_base/java kent-335_base/python
fi
