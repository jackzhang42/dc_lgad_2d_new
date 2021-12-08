#!/bin/sh

rootdir=~/test-for-garfieldpp/root/dc_lgad/time
builddir=~/test-for-garfieldpp/dc_lgad_2d/build

cd $builddir
make clean
make

cd $rootdir
rm time_lgad*.dat

i=0
while [ "$i" -lt 10 ]
do
echo "this is No.$i loop"
cd $builddir
./dc_lgad
cp ./time_lgad.dat $rootdir/time_lgad_$i.dat
cd $rootdir
cat time_lgad_$i.dat >> time_lgad.dat
i=`expr $i + 1`
done
