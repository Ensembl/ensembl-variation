#!/bin/bash

# kent src
export MACHTYPE=$(uname -m)
export MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'`
export MYSQLLIBS=`mysql_config --libs`

# Build kent src
cd kent-335_base/src/lib
echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk
make
cd ../jkOwnLib
make
