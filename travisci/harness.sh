#!/bin/bash

export PERL5LIB=$PWD/bioperl-live-bioperl-release-1-6-1:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/ensembl-hive/modules:$PWD/modules:$PWD/scripts/import/:$PWD/ensembl-io/modules:$PWD/ensembl-funcgen/modules

export PATH=$PATH:$PWD/C_code:$PWD/htslib

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,ensembl-hive,+ignore,ensembl-io,+ignore,ensembl-funcgen' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/modules/t $SKIP_TESTS
else
  perl $PWD/ensembl-test/scripts/runtests.pl $PWD/modules/t $SKIP_TESTS
fi

rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
