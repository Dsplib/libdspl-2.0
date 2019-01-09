#!/bin/bash


find ../../ru -name "*.c"   -exec cp -rf  {} ../test/src \;
find ../../ru -name "*.dox" -exec cp -rf  {} ../test/dox/ru \;
find ../../ru -name "*.plt" -exec cp -rf  {} ../test/bin/gnuplot \;

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

