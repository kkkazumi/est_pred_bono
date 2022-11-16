#!/bin/bash -eu

for username in B9001 B9002 B9003
do
	for condition in CN ML
	do
		#B9001_ML B9002_CN B9002_ML B9003_CN B9003_ML
		filename="./bono_datamaker/data/sudata/"$username"_"$condition"_factor.csv"
		LEN=`wc -l $filename`
		ARR=(${LEN//,/ })
		echo $ARR
		sed -i -e "s/define N 30/define N $ARR/" estimate.cpp
		make estimate
		./estimate ${username} ${condition}
		sed -i -e "s/define N $ARR/define N 30/" estimate.cpp
	done
done
