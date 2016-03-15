#!/bin/sh
set -ex
wget http://ftp.gnu.org/pub/gnu/gsl/gsl-2.1.tar.gz
tar -xzvf gsl-2.1.tar.gz
cd gsl-2.1 && ./configure --prefix=/usr --disable-static 
make > gsl-log-file 2>&1
sudo make install
