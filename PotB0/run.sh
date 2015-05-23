#!/bin/bash
#-----------------------------------------------
# run MonteCarlo simulation with different input
#-----------------------------------------------
inputDIR="./input"
resultDIR="./result"
exeDIR="../"
exe=DoubleWell

for i in 0.1 1.0 10
do
echo "running: ../DoubleWell <PotB0-RefT$i.in >PotB0-RefT$i.out"

$exeDIR$exe< $inputDIR/PotB0-RefT$i.in >$resultDIR/PotB0-RefT$i.out

echo "finish successfully"

done
