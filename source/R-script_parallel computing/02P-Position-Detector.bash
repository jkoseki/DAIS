#!/usr/bin/bash

numcase=1
sumtgt=$(( $numcase + 1 ))

for ((i = 1 ; i <= sumtgt ; i++))
do
	nohup Rscript 02P-Position-Detector.R $i &
done

