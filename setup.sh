#!/bin/bash

mkdir -p libs
cd libs

# download and build FFTW if not already present
if [ ! -d "fftw-3.3.10" ]; then
    echo "Downloading FFTW..."
    curl -L http://www.fftw.org/fftw-3.3.10.tar.gz | tar xz
    cd fftw-3.3.10
    mkdir -p compiled
    ./configure --prefix=$(pwd)/compiled
    make && make install
    cd ..
fi

echo "setup complete"