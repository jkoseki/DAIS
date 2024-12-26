#!/usr/bin/bash

numcase=1
numtpdf=500

nproc=25

sumtpdf=$(( ( $numcase + 1) * $numtpdf ))
numblok=$(( $sumtpdf / $nproc ))


for ((i = 1 ; i <= nproc ; i++))
do
	vin=$(( (i - 1) * $numblok + 1))
	vot=$(( i * $numblok ))
	nohup Rscript 01P-CA_location2.R $vin $vot &
done

