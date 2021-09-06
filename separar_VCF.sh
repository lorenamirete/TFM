#!/bin/bash


for i in $(cat IDlist.txt);

do i="${i%\\n}";

    vcf-subset -c $i GCAT_WGS_colobran_SHARE.vcf > ${i}.vcf

done

