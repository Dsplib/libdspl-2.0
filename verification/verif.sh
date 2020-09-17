mingw32-make

rm -f bin/verification.log

cd bin
for file in *.exe
do
    "./$file"
done
cd ../
cp  bin/verification.log verification.log
rm -f bin/verification.log
pause 5

