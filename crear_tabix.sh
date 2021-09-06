#!/bin/bash

for individual in ./*.vcf ;

do

bgzip -c $individual > $individual.gz ;
tabix -p vcf $individual.gz

done

