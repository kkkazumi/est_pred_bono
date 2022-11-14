#!/bin/bash -eu

#echo N_num
#num=$[$N_num-1]
#sed -i -e "s/define N $num/define N $N_num/" estimate.cpp
for username in B9001_CN B9001_ML B9002_CN B9002_ML B9003_CN B9003_ML
do
    filename="./bono_datamaker/data/sudata/"$username"_factor.csv"
    LEN=`wc -l $filename`
	ARR=(${LEN//,/ })
	echo $ARR
	sed -i -e "s/define N 30/define N $ARR/" estimate.cpp
	make estimate
	./estimate
	sed -i -e "s/define N $ARR/define N 30/" estimate.cpp
done
