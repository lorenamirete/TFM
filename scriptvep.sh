#!/bin/bash


for paciente in ./*.vcf ;

do

	./vep --cache --tab -i $paciente -e --force_overwrite  -o ./txt_IDPs_CVID/$paciente.txt ;
	sed -i '/^##/ d' ./txt_IDPs_CVID/$paciente.txt

done ;


cd ./txt_IDPs_CVID/

rename 's/.vcf.txt/.txt/' *.txt








