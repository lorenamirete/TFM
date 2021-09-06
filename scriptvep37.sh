#!/bin/bash


for paciente in ./*.vcf ;

do

	./vep --dir_cache /home/lmirete/ensembl-vep/ --cache --assembly GRCh37 --offline --tab -i $paciente --sift b --polyphen b --symbol --canonical --biotype --regulatory --force_overwrite -o ./txt_IDPs_CONTROL_FILTER/$paciente.txt ;
	sed -i '/^##/ d' ./txt_IDPs_CONTROL_FILTER/$paciente.txt

done ;


cd ./txt_IDPs_CONTROL_FILTER/

rename 's/_filtered.vcf.txt/.txt/' *.txt








