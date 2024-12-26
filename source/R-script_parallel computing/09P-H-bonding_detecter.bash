#!/usr/bin/bash

numcase=1

for ((i = 1 ; i <= numcase ; i++))
do
	nohup Rscript 09P-H-bonding_detecter.R $i &
done

