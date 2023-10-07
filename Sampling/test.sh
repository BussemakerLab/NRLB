#!/bin/bash

cd /Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/dbsampling/src/

make clean
make all

cd ../

./deepbind D00410.003 20 TCGTA TCGTA -20 20 .1 100 1e-4 100 .9
