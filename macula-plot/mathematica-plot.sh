#!/bin/bash
#

# Build and execute macula-plot
make veryclean
make macula-plot
./macula-plot

# Plot the data using Mathematica
math='/Applications/Mathematica.app/Contents/MacOS/MathematicaScript'
file='mathematica-plot.txt'
script='script.txt'
currentfolder=`pwd`
rm -f $script
nlines=`cat $file | wc -l`
head -2 $file > $script
linea='thisfolder="'
lineb='/";'
line=`echo $linea$currentfolder$lineb`
echo $line >> $script
ntails=$(($nlines - 3))
tail -$ntails $file >> $script
$math -script $script
rm -f $script
