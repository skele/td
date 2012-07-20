#!/bin/bash

max=`cat out1.dat | awk 'BEGIN {max = 0.0} {if ($3 > max) {max=$3}} END {print max;}'`
echo $max

for i in `seq -w 1 300`
#
do
    pngname="${i}.png"
    cat >plotscript <<__EOF__
unset key
set terminal png
set xlabel 'Wavelength (nm)'
set ylabel 'Flux (something like flux)'
set xrange [655.3:657.3]
set yrange [0:${max}]
__EOF__

echo "set output '${pngname}'" >> plotscript
#days=`echo "${i}" | bc`
echo "set title 'After ${i} days'" >> plotscript
echo "plot 'alltimes.dat' u 2:(\$1==${i}?\$3:1/0) w lp" >> plotscript
echo "exit" >> plotscript
gnuplot plotscript
done
