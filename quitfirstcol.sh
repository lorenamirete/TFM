#!/bin/bash

cut -d "," -f2- $i  > $i

for f in ./*.csv
do
	name=`basename $f`
	echo "processing file : $name"
    #kepp all column excep the first one of each csv file 
	cut -d"," -f2- $f > new/$name
    #files using the same names are stored in directory new/  

done
