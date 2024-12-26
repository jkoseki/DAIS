#!/usr/bin/bash

#01============================================================
./01P-CA_location2.bash 

sleep 5
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  echo "Now, Step 1 - CA location -"
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done


#02============================================================
./02P-Position-Detector.bash

sleep 5
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  echo "Now, Step 2 - Position Detector -"
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done


#03============================================================
Rscript 03-Case-Contl_Common.R

sleep 5
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  echo "Now, Step 3 - Compared with Case and Contlor -"
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done


#04============================================================
./04P-PU-Coord.bash

sleep 5
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  echo "Now, Step 4 - Pickup Coordinates -"
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done


#05-08========================================================
echo "Now, Step 5 to 8"
./05-08.bash

sleep 1
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done


#09============================================================
./09P-H-bonding_detecter.bash

sleep 5
test=1
PES_period=1
while [[ $test -ge 1 ]];
do
  echo "Now, Step 9 - H-bonding detecter -"
  sleep $PES_period
  test=$(ps | grep -i " R" | wc -l )
  if [ $PES_period -lt 10 ]; then
    PES_period=$((PES_period + 1))
  elif [ $PES_period -ge 10 ] && [ $PES_period -lt 60 ]; then
    PES_period=$((PES_period + 5))
  fi
done











