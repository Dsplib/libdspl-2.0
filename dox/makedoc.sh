#!/bin/bash

cd ../
make
cd test/bin
for file in *.exe
do
    "./$file"
done
cd ../../dox

doxygen doxyfile_ru

cd ../
make clean
cd dox

pkill -x gnuplot

