#!/usr/bin/bash

mkdir 01-min
mkdir 02-product
mkdir 03-mmgbsa

cp ~/bin/MMGBSA-in/min.in     01-min/.
cp ~/bin/MMGBSA-in/heat.in    01-min/.
cp ~/bin/MMGBSA-in/density.in 01-min/.
cp ~/bin/MMGBSA-in/equil.in   01-min/.

cp ~/bin/MMGBSA-in/prod.in    02-product/.

cp ~/bin/MMGBSA-in/mmgbsa.in  03-mmgbsa/.

